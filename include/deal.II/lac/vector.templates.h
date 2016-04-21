// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2015 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef dealii__vector_templates_h
#define dealii__vector_templates_h


#include <deal.II/base/template_constraints.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/parallel.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/vectorization.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/block_vector.h>

#ifdef DEAL_II_WITH_PETSC
#  include <deal.II/lac/petsc_vector.h>
#  include <deal.II/lac/petsc_parallel_vector.h>
#endif

#ifdef DEAL_II_WITH_TRILINOS
#  include <deal.II/lac/trilinos_vector.h>
#endif


#include <cmath>
#include <cstring>
#include <algorithm>
#include <iostream>
#include <iomanip>

DEAL_II_NAMESPACE_OPEN


namespace internal
{
  typedef types::global_dof_index size_type;

  template <typename T>
  bool is_non_negative (const T &t)
  {
    return t >= 0;
  }


  template <typename T>
  bool is_non_negative (const std::complex<T> &)
  {
    Assert (false,
            ExcMessage ("Complex numbers do not have an ordering."));

    return false;
  }


  template <typename T>
  void print (const T    &t,
              const char *format)
  {
    if (format != 0)
      std::printf (format, t);
    else
      std::printf (" %5.2f", double(t));
  }



  template <typename T>
  void print (const std::complex<T> &t,
              const char            *format)
  {
    if (format != 0)
      std::printf (format, t.real(), t.imag());
    else
      std::printf (" %5.2f+%5.2fi",
                   double(t.real()), double(t.imag()));
  }

  // call std::copy, except for in
  // the case where we want to copy
  // from std::complex to a
  // non-complex type
  template <typename T, typename U>
  void copy (const T *begin,
             const T *end,
             U       *dest)
  {
    std::copy (begin, end, dest);
  }

  template <typename T, typename U>
  void copy (const std::complex<T> *begin,
             const std::complex<T> *end,
             std::complex<U>       *dest)
  {
    std::copy (begin, end, dest);
  }

  template <typename T, typename U>
  void copy (const std::complex<T> *,
             const std::complex<T> *,
             U *)
  {
    Assert (false, ExcMessage ("Can't convert a vector of complex numbers "
                               "into a vector of reals/doubles"));
  }



#ifdef DEAL_II_WITH_THREADS
  /**
   * This struct takes the loop range from the tbb parallel for loop and
   * translates it to the actual ranges of the for loop within the vector. It
   * encodes the grain size but might choose larger values of chunks than the
   * minimum grain size. The minimum grain size given to tbb is then simple
   * 1. For affinity reasons, the layout in this loop must be kept in sync
   * with the respective class for reductions further down.
   */
  template <typename Functor>
  struct TBBForFunctor
  {
    TBBForFunctor(Functor &functor,
                  const size_type vec_size)
      :
      functor(functor),
      vec_size(vec_size)
    {
      // set chunk size for sub-tasks
      const unsigned int gs = internal::Vector::minimum_parallel_grain_size;
      n_chunks = std::min(static_cast<size_type>(4*MultithreadInfo::n_threads()),
                          vec_size / gs);
      chunk_size = vec_size / n_chunks;

      // round to next multiple of 512 (or minimum grain size if that happens
      // to be smaller). this is advantageous because our accumulation
      // algorithms favor lengths of a power of 2 due to pairwise summation ->
      // at most one 'oddly' sized chunk
      if (chunk_size > 512)
        chunk_size = ((chunk_size + 511)/512)*512;
      n_chunks = (vec_size + chunk_size - 1) / chunk_size;
      AssertIndexRange((n_chunks-1)*chunk_size, vec_size);
      AssertIndexRange(vec_size, n_chunks*chunk_size+1);
    };

    void operator() (const tbb::blocked_range<size_type> &range) const
    {
      const size_type begin = range.begin()*chunk_size;
      const size_type end = std::min(range.end()*chunk_size, vec_size);
      functor(begin, end);
    }

    Functor &functor;
    const size_type vec_size;
    unsigned int n_chunks;
    size_type chunk_size;
  };
#endif

  template <typename Functor>
  void parallel_for(Functor &functor,
                    size_type vec_size,
                    std_cxx11::shared_ptr<parallel::internal::TBBPartitioner> &partitioner)
  {
#ifdef DEAL_II_WITH_THREADS
    // only go to the parallel function in case there are at least 4 parallel
    // items, otherwise the overhead is too large
    if (vec_size >= 4*internal::Vector::minimum_parallel_grain_size &&
        MultithreadInfo::n_threads() > 1)
      {
        Assert(partitioner.get() != NULL,
               ExcInternalError("Unexpected initialization of Vector that does "
                                "not set the TBB partitioner to a usable state."));
        std_cxx11::shared_ptr<tbb::affinity_partitioner> tbb_partitioner =
          partitioner->acquire_one_partitioner();

        TBBForFunctor<Functor> generic_functor(functor, vec_size);
        tbb::parallel_for (tbb::blocked_range<size_type> (0,
                                                          generic_functor.n_chunks,
                                                          1),
                           generic_functor,
                           *tbb_partitioner);
        partitioner->release_one_partitioner(tbb_partitioner);
      }
    else if (vec_size > 0)
      functor(0,vec_size);
#else
    functor(0,vec_size);
#endif
  }


  // Define the functors necessary to use SIMD with TBB. we also include the
  // simple copy and set operations

  template <typename Number>
  struct Vector_set
  {
    Number *dst;
    Number value;

    void operator() (const size_type begin, const size_type end) const
    {
      if (value == Number())
        std::memset (dst+begin,0,(end-begin)*sizeof(Number));
      else
        std::fill (dst+begin, dst+end, value);
    }
  };

  template <typename Number, typename OtherNumber>
  struct Vector_copy
  {
    const OtherNumber *src;
    Number *dst;

    void operator() (const size_type begin, const size_type end) const
    {
      if (types_are_equal<Number,OtherNumber>::value)
        std::memcpy(dst+begin, src+begin, (end-begin)*sizeof(Number));
      else
        {
          DEAL_II_OPENMP_SIMD_PRAGMA
          for (typename dealii::Vector<Number>::size_type i=begin; i<end; ++i)
            dst[i] = src[i];
        }
    }
  };

  template <typename Number>
  struct Vectorization_multiply_factor
  {
    Number *val;
    Number factor;

    void operator() (const size_type begin, const size_type end) const
    {
      if (parallel::internal::EnableOpenMPSimdFor<Number>::value)
        {
          DEAL_II_OPENMP_SIMD_PRAGMA
          for (size_type i=begin; i<end; ++i)
            val[i] *= factor;
        }
      else
        {
          for (size_type i=begin; i<end; ++i)
            val[i] *= factor;
        }
    }
  };

  template <typename Number>
  struct Vectorization_add_av
  {
    Number *val;
    Number *v_val;
    Number factor;

    void operator() (const size_type begin, const size_type end) const
    {
      if (parallel::internal::EnableOpenMPSimdFor<Number>::value)
        {
          DEAL_II_OPENMP_SIMD_PRAGMA
          for (size_type i=begin; i<end; ++i)
            val[i] += factor*v_val[i];
        }
      else
        {
          for (size_type i=begin; i<end; ++i)
            val[i] += factor*v_val[i];
        }
    }
  };

  template <typename Number>
  struct Vectorization_sadd_xav
  {
    Number *val;
    Number *v_val;
    Number a;
    Number x;

    void operator() (const size_type begin, const size_type end) const
    {
      if (parallel::internal::EnableOpenMPSimdFor<Number>::value)
        {
          DEAL_II_OPENMP_SIMD_PRAGMA
          for (size_type i=begin; i<end; ++i)
            val[i] = x*val[i] + a*v_val[i];
        }
      else
        {
          for (size_type i=begin; i<end; ++i)
            val[i] = x*val[i] + a*v_val[i];
        }
    }
  };

  template <typename Number>
  struct Vectorization_subtract_v
  {
    Number *val;
    Number *v_val;

    void operator() (const size_type begin, const size_type end) const
    {
      if (parallel::internal::EnableOpenMPSimdFor<Number>::value)
        {
          DEAL_II_OPENMP_SIMD_PRAGMA
          for (size_type i=begin; i<end; ++i)
            val[i] -= v_val[i];
        }
      else
        {
          for (size_type i=begin; i<end; ++i)
            val[i] -= v_val[i];
        }
    }
  };

  template <typename Number>
  struct Vectorization_add_factor
  {
    Number *val;
    Number factor;

    void operator() (const size_type begin, const size_type end) const
    {
      if (parallel::internal::EnableOpenMPSimdFor<Number>::value)
        {
          DEAL_II_OPENMP_SIMD_PRAGMA
          for (size_type i=begin; i<end; ++i)
            val[i] += factor;
        }
      else
        {
          for (size_type i=begin; i<end; ++i)
            val[i] += factor;
        }
    }
  };

  template <typename Number>
  struct Vectorization_add_v
  {
    Number *val;
    Number *v_val;

    void operator() (const size_type begin, const size_type end) const
    {
      if (parallel::internal::EnableOpenMPSimdFor<Number>::value)
        {
          DEAL_II_OPENMP_SIMD_PRAGMA
          for (size_type i=begin; i<end; ++i)
            val[i] += v_val[i];
        }
      else
        {
          for (size_type i=begin; i<end; ++i)
            val[i] += v_val[i];
        }
    }
  };

  template <typename Number>
  struct Vectorization_add_avpbw
  {
    Number *val;
    Number *v_val;
    Number *w_val;
    Number a;
    Number b;

    void operator() (const size_type begin, const size_type end) const
    {
      if (parallel::internal::EnableOpenMPSimdFor<Number>::value)
        {
          DEAL_II_OPENMP_SIMD_PRAGMA
          for (size_type i=begin; i<end; ++i)
            val[i] = val[i] + a*v_val[i] + b*w_val[i];
        }
      else
        {
          for (size_type i=begin; i<end; ++i)
            val[i] = val[i] + a*v_val[i] + b*w_val[i];
        }
    }
  };

  template <typename Number>
  struct Vectorization_sadd_xv
  {
    Number *val;
    Number *v_val;
    Number x;

    void operator() (const size_type begin, const size_type end) const
    {
      if (parallel::internal::EnableOpenMPSimdFor<Number>::value)
        {
          DEAL_II_OPENMP_SIMD_PRAGMA
          for (size_type i=begin; i<end; ++i)
            val[i] = x*val[i] + v_val[i];
        }
      else
        {
          for (size_type i=begin; i<end; ++i)
            val[i] = x*val[i] + v_val[i];
        }
    }
  };

  template <typename Number>
  struct Vectorization_sadd_xavbw
  {
    Number *val;
    Number *v_val;
    Number *w_val;
    Number x;
    Number a;
    Number b;

    void operator() (const size_type begin, const size_type end) const
    {
      if (parallel::internal::EnableOpenMPSimdFor<Number>::value)
        {
          DEAL_II_OPENMP_SIMD_PRAGMA
          for (size_type i=begin; i<end; ++i)
            val[i] = x*val[i] + a*v_val[i] + b*w_val[i];
        }
      else
        {
          for (size_type i=begin; i<end; ++i)
            val[i] = x*val[i] + a*v_val[i] + b*w_val[i];
        }
    }
  };

  template <typename Number>
  struct Vectorization_scale
  {
    Number *val;
    Number *v_val;

    void operator() (const size_type begin, const size_type end) const
    {
      if (parallel::internal::EnableOpenMPSimdFor<Number>::value)
        {
          DEAL_II_OPENMP_SIMD_PRAGMA
          for (size_type i=begin; i<end; ++i)
            val[i] *= v_val[i];
        }
      else
        {
          for (size_type i=begin; i<end; ++i)
            val[i] *= v_val[i];
        }
    }
  };

  template <typename Number>
  struct Vectorization_equ_au
  {
    Number *val;
    Number *u_val;
    Number a;

    void operator() (const size_type begin, const size_type end) const
    {
      if (parallel::internal::EnableOpenMPSimdFor<Number>::value)
        {
          DEAL_II_OPENMP_SIMD_PRAGMA
          for (size_type i=begin; i<end; ++i)
            val[i] = a*u_val[i];
        }
      else
        {
          for (size_type i=begin; i<end; ++i)
            val[i] = a*u_val[i];
        }
    }
  };

  template <typename Number>
  struct Vectorization_equ_aubv
  {
    Number *val;
    Number *u_val;
    Number *v_val;
    Number a;
    Number b;

    void operator() (const size_type begin, const size_type end) const
    {
      if (parallel::internal::EnableOpenMPSimdFor<Number>::value)
        {
          DEAL_II_OPENMP_SIMD_PRAGMA
          for (size_type i=begin; i<end; ++i)
            val[i] = a*u_val[i] + b*v_val[i];
        }
      else
        {
          for (size_type i=begin; i<end; ++i)
            val[i] = a*u_val[i] + b*v_val[i];
        }
    }
  };

  template <typename Number>
  struct Vectorization_equ_aubvcw
  {
    Number *val;
    Number *u_val;
    Number *v_val;
    Number *w_val;
    Number a;
    Number b;
    Number c;

    void operator() (const size_type begin, const size_type end) const
    {
      if (parallel::internal::EnableOpenMPSimdFor<Number>::value)
        {
          DEAL_II_OPENMP_SIMD_PRAGMA
          for (size_type i=begin; i<end; ++i)
            val[i] = a*u_val[i] + b*v_val[i] + c*w_val[i];
        }
      else
        {
          for (size_type i=begin; i<end; ++i)
            val[i] = a*u_val[i] + b*v_val[i] + c*w_val[i];
        }
    }
  };

  template <typename Number>
  struct Vectorization_ratio
  {
    Number *val;
    Number *a_val;
    Number *b_val;

    void operator() (const size_type begin, const size_type end) const
    {
      if (parallel::internal::EnableOpenMPSimdFor<Number>::value)
        {
          DEAL_II_OPENMP_SIMD_PRAGMA
          for (size_type i=begin; i<end; ++i)
            val[i] = a_val[i]/b_val[i];
        }
      else
        {
          for (size_type i=begin; i<end; ++i)
            val[i] = a_val[i]/b_val[i];
        }
    }
  };



  // All sums over all the vector entries (l2-norm, inner product, etc.) are
  // performed with the same code, using a templated operation defined
  // here. There are always two versions defined, a standard one that covers
  // most cases and a vectorized one which is only for equal types and float
  // and double.
  template <typename Number, typename Number2>
  struct Dot
  {
    static const bool vectorizes = types_are_equal<Number,Number2>::value &&
                                   (VectorizedArray<Number>::n_array_elements > 1);

    Number
    operator() (const size_type i) const
    {
      return X[i] * Number(numbers::NumberTraits<Number2>::conjugate(Y[i]));
    }

    VectorizedArray<Number>
    do_vectorized(const size_type i) const
    {
      VectorizedArray<Number> x, y;
      x.load(X+i);
      y.load(Y+i);
      return x * y;
    }

    const Number  *X;
    const Number2 *Y;
  };

  template <typename Number, typename RealType>
  struct Norm2
  {
    static const bool vectorizes = VectorizedArray<Number>::n_array_elements > 1;

    RealType
    operator() (const size_type i) const
    {
      return numbers::NumberTraits<Number>::abs_square(X[i]);
    }

    VectorizedArray<Number>
    do_vectorized(const size_type i) const
    {
      VectorizedArray<Number> x;
      x.load(X+i);
      return x * x;
    }

    const Number *X;
  };

  template <typename Number, typename RealType>
  struct Norm1
  {
    static const bool vectorizes = VectorizedArray<Number>::n_array_elements > 1;

    RealType
    operator() (const size_type i) const
    {
      return numbers::NumberTraits<Number>::abs(X[i]);
    }

    VectorizedArray<Number>
    do_vectorized(const size_type i) const
    {
      VectorizedArray<Number> x;
      x.load(X+i);
      return std::abs(x);
    }

    const Number *X;
  };

  template <typename Number, typename RealType>
  struct NormP
  {
    static const bool vectorizes = VectorizedArray<Number>::n_array_elements > 1;

    RealType
    operator() (const size_type i) const
    {
      return std::pow(numbers::NumberTraits<Number>::abs(X[i]), p);
    }

    VectorizedArray<Number>
    do_vectorized(const size_type i) const
    {
      VectorizedArray<Number> x;
      x.load(X+i);
      return std::pow(std::abs(x),p);
    }

    const Number *X;
    RealType p;
  };

  template <typename Number>
  struct MeanValue
  {
    static const bool vectorizes = VectorizedArray<Number>::n_array_elements > 1;

    Number
    operator() (const size_type i) const
    {
      return X[i];
    }

    VectorizedArray<Number>
    do_vectorized(const size_type i) const
    {
      VectorizedArray<Number> x;
      x.load(X+i);
      return x;
    }

    const Number *X;
  };

  template <typename Number>
  struct AddAndDot
  {
    static const bool vectorizes = VectorizedArray<Number>::n_array_elements > 1;

    Number
    operator() (const size_type i) const
    {
      X[i] += a * V[i];
      return X[i] * Number(numbers::NumberTraits<Number>::conjugate(W[i]));
    }

    VectorizedArray<Number>
    do_vectorized(const size_type i) const
    {
      VectorizedArray<Number> x, w, v;
      x.load(X+i);
      v.load(V+i);
      x += a * v;
      x.store(X+i);
      // may only load from W after storing in X because the pointers might
      // point to the same memory
      w.load(W+i);
      return x * w;
    }

    Number *X;
    const Number *V, *W;
    Number a;
  };



  // this is the main working loop for all vector sums using the templated
  // operation above. it accumulates the sums using a block-wise summation
  // algorithm with post-update. this blocked algorithm has been proposed in
  // a similar form by Castaldo, Whaley and Chronopoulos (SIAM
  // J. Sci. Comput. 31, 1156-1174, 2008) and we use the smallest possible
  // block size, 2. Sometimes it is referred to as pairwise summation. The
  // worst case error made by this algorithm is on the order O(eps *
  // log2(vec_size)), whereas a naive summation is O(eps * vec_size). Even
  // though the Kahan summation is even more accurate with an error O(eps)
  // by carrying along remainders not captured by the main sum, that involves
  // additional costs which are not worthwhile. See the Wikipedia article on
  // the Kahan summation algorithm.

  // The algorithm implemented here has the additional benefit that it is
  // easily parallelized without changing the order of how the elements are
  // added (floating point addition is not associative). For the same vector
  // size and minimum_parallel_grainsize, the blocks are always the
  // same and added pairwise.

  // The depth of recursion is controlled by the 'magic' parameter
  // vector_accumulation_recursion_threshold: If the length is below
  // vector_accumulation_recursion_threshold * 32 (32 is the part of code we
  // unroll), a straight loop instead of recursion will be used.  At the
  // innermost level, eight values are added consecutively in order to better
  // balance multiplications and additions.

  // The code returns the result as the last argument in order to make
  // spawning tasks simpler and use automatic template deduction.

  const unsigned int vector_accumulation_recursion_threshold = 128;

  template <typename Operation, typename ResultType>
  void accumulate_recursive (const Operation   &op,
                             const size_type    first,
                             const size_type    last,
                             ResultType        &result)
  {
    const size_type vec_size = last - first;
    if (vec_size <= vector_accumulation_recursion_threshold * 32)
      {
        // the vector is short enough so we perform the summation. first
        // work on the regular part. The innermost 32 values are expanded in
        // order to obtain known loop bounds for most of the work.
        size_type index = first;
        ResultType outer_results [vector_accumulation_recursion_threshold];
        size_type n_chunks = vec_size / 32;
        const size_type remainder = vec_size % 32;
        Assert (remainder == 0 || n_chunks < vector_accumulation_recursion_threshold,
                ExcInternalError());

        // Select between the regular version and vectorized version based
        // on the number types we are given. To choose the vectorized
        // version often enough, we need to have all tasks but the last one
        // to be divisible by the vectorization length
        accumulate_regular(op, n_chunks, index, outer_results,
                           internal::bool2type<Operation::vectorizes>());

        // now work on the remainder, i.e., the last up to 32 values. Use
        // switch statement with fall-through to work on these values.
        if (remainder > 0)
          {
            AssertIndexRange(n_chunks, vector_accumulation_recursion_threshold+1);
            const size_type inner_chunks = remainder / 8;
            Assert (inner_chunks <= 3, ExcInternalError());
            const size_type remainder_inner = remainder % 8;
            ResultType r0 = ResultType(), r1 = ResultType(),
                       r2 = ResultType();
            switch (inner_chunks)
              {
              case 3:
                r2 = op(index++);
                for (size_type j=1; j<8; ++j)
                  r2 += op(index++);
              // no break
              case 2:
                r1 = op(index++);
                for (size_type j=1; j<8; ++j)
                  r1 += op(index++);
                r1 += r2;
              // no break
              case 1:
                r2 = op(index++);
                for (size_type j=1; j<8; ++j)
                  r2 += op(index++);
              // no break
              default:
                for (size_type j=0; j<remainder_inner; ++j)
                  r0 += op(index++);
                r0 += r2;
                r0 += r1;
                if (n_chunks == vector_accumulation_recursion_threshold)
                  outer_results[vector_accumulation_recursion_threshold-1] += r0;
                else
                  {
                    outer_results[n_chunks] = r0;
                    n_chunks++;
                  }
                break;
              }
          }
        AssertDimension(index, last);

        // now sum the results from the chunks
        // recursively
        while (n_chunks > 1)
          {
            if (n_chunks % 2 == 1)
              outer_results[n_chunks++] = ResultType();
            for (size_type i=0; i<n_chunks; i+=2)
              outer_results[i/2] = outer_results[i] + outer_results[i+1];
            n_chunks /= 2;
          }
        result = outer_results[0];
      }
    else
      {
        // split vector into four pieces and work on the pieces
        // recursively. Make pieces (except last) divisible by one fourth the
        // recursion threshold.
        const size_type new_size =
          (vec_size / (vector_accumulation_recursion_threshold * 32)) *
          vector_accumulation_recursion_threshold * 8;
        ResultType r0, r1, r2, r3;
        accumulate_recursive (op, first, first+new_size, r0);
        accumulate_recursive (op, first+new_size, first+2*new_size, r1);
        accumulate_recursive (op, first+2*new_size, first+3*new_size, r2);
        accumulate_recursive (op, first+3*new_size, last, r3);
        r0 += r1;
        r2 += r3;
        result = r0 + r2;
      }
  }


  // this is the inner working routine for the accumulation loops
  // below. This is the standard case where the loop bounds are known. We
  // pulled this function out of the regular accumulate routine because we
  // might do this thing vectorized (see specialized function below)
  template <typename Operation, typename ResultType>
  void
  accumulate_regular(const Operation &op,
                     size_type       &n_chunks,
                     size_type       &index,
                     ResultType (&outer_results)[vector_accumulation_recursion_threshold],
                     internal::bool2type<false>)
  {
    for (size_type i=0; i<n_chunks; ++i)
      {
        ResultType r0 = op(index);
        ResultType r1 = op(index+1);
        ResultType r2 = op(index+2);
        ResultType r3 = op(index+3);
        index += 4;
        for (size_type j=1; j<8; ++j, index += 4)
          {
            r0 += op(index);
            r1 += op(index+1);
            r2 += op(index+2);
            r3 += op(index+3);
          }
        r0 += r1;
        r2 += r3;
        outer_results[i] = r0 + r2;
      }
  }



  // this is the inner working routine for the accumulation loops
  // below. This is the specialized case where the loop bounds are known and
  // where we can vectorize. In that case, we request the 'do_vectorized'
  // routine of the operation instead of the regular one which does several
  // operations at once.
  template <typename Operation, typename Number>
  void
  accumulate_regular(const Operation &op,
                     size_type       &n_chunks,
                     size_type       &index,
                     Number (&outer_results)[vector_accumulation_recursion_threshold],
                     internal::bool2type<true>)
  {
    const unsigned int nvecs = VectorizedArray<Number>::n_array_elements;
    const size_type regular_chunks = n_chunks/nvecs;
    for (size_type i=0; i<regular_chunks; ++i)
      {
        VectorizedArray<Number> r0 = op.do_vectorized(index);
        VectorizedArray<Number> r1 = op.do_vectorized(index+nvecs);
        VectorizedArray<Number> r2 = op.do_vectorized(index+2*nvecs);
        VectorizedArray<Number> r3 = op.do_vectorized(index+3*nvecs);
        index += nvecs*4;
        for (size_type j=1; j<8; ++j, index += nvecs*4)
          {
            r0 += op.do_vectorized(index);
            r1 += op.do_vectorized(index+nvecs);
            r2 += op.do_vectorized(index+2*nvecs);
            r3 += op.do_vectorized(index+3*nvecs);
          }
        r0 += r1;
        r2 += r3;
        r0 += r2;
        r0.store(&outer_results[i*VectorizedArray<Number>::n_array_elements]);
      }

    // If we are treating a case where the vector length is not divisible by
    // the vectorization length, need a cleanup loop
    AssertIndexRange(VectorizedArray<Number>::n_array_elements,
                     17);
    if (n_chunks % VectorizedArray<Number>::n_array_elements != 0)
      {
        VectorizedArray<Number> r0 = VectorizedArray<Number>(),
                                r1 = VectorizedArray<Number>();
        const size_type start_irreg = regular_chunks * nvecs;
        for (size_type c=start_irreg; c<n_chunks; ++c)
          for (size_type j=0; j<32; j+=2*nvecs, index+=2*nvecs)
            {
              r0 += op.do_vectorized(index);
              r1 += op.do_vectorized(index+nvecs);
            }
        r0 += r1;
        r0.store(&outer_results[start_irreg]);
        n_chunks = start_irreg + VectorizedArray<Number>::n_array_elements;
      }
  }



#ifdef DEAL_II_WITH_THREADS
  /**
   * This struct takes the loop range from the tbb parallel for loop and
   * translates it to the actual ranges of the reduction loop inside the
   * vector. It encodes the grain size but might choose larger values of
   * chunks than the minimum grain size. The minimum grain size given to tbb
   * is 1. For affinity reasons, the layout in this loop must be kept in sync
   * with the respective class for plain for loops further up.
   *
   * Due to this construction, TBB usually only sees a loop of length
   * 4*num_threads with grain size 1. The actual ranges inside the vector are
   * computed outside of TBB because otherwise TBB would split the ranges in
   * some unpredictable position which destroys exact bitwise
   * reproducibility. An important part of this is that inside
   * TBBReduceFunctor::operator() the recursive calls to accumulate are done
   * sequentially on one item a time (even though we could directly run it on
   * the whole range given through the tbb::blocked_range times the chunk size
   * - but that would be unpredictable). Thus, the values we cannot control
   * are the positions in the array that gets filled - but up to that point
   * the algorithm TBB sees is just a parallel for and nothing unpredictable
   * can happen.
   *
   * To sum up: Once the number of threads and the vector size are fixed, we
   * have an exact layout of how the calls into the recursive function will
   * happen. Inside the recursive function, we again only depend on the
   * length. Finally, the concurrent threads write into different positions in
   * a result vector in a thread-safe way and the addition in the short array
   * is again serial.
   */
  template <typename Operation, typename ResultType>
  struct TBBReduceFunctor
  {
    static const unsigned int threshold_array_allocate = 512;

    TBBReduceFunctor(const Operation   &op,
                     const size_type    vec_size)
      :
      op(op),
      vec_size(vec_size)
    {
      // set chunk size for sub-tasks
      const unsigned int gs = internal::Vector::minimum_parallel_grain_size;
      n_chunks = std::min(static_cast<size_type>(4*MultithreadInfo::n_threads()),
                          vec_size / gs);
      chunk_size = vec_size / n_chunks;

      // round to next multiple of 512 (or leave it at the minimum grain size
      // if that happens to be smaller). this is advantageous because our
      // algorithm favors lengths of a power of 2 due to pairwise summation ->
      // at most one 'oddly' sized chunk
      if (chunk_size > 512)
        chunk_size = ((chunk_size + 511)/512)*512;
      n_chunks = (vec_size + chunk_size - 1) / chunk_size;
      AssertIndexRange((n_chunks-1)*chunk_size, vec_size);
      AssertIndexRange(vec_size, n_chunks*chunk_size+1);

      if (n_chunks > threshold_array_allocate)
        {
          large_array.resize(n_chunks);
          array_ptr = &large_array[0];
        }
      else
        array_ptr = &small_array[0];
    };

    void operator() (const tbb::blocked_range<size_type> &range) const
    {
      for (size_type i = range.begin(); i < range.end(); ++i)
        accumulate_recursive(op, i*chunk_size, std::min((i+1)*chunk_size, vec_size),
                             array_ptr[i]);
    }

    ResultType do_sum() const
    {
      while (n_chunks > 1)
        {
          if (n_chunks % 2 == 1)
            array_ptr[n_chunks++] = ResultType();
          for (size_type i=0; i<n_chunks; i+=2)
            array_ptr[i/2] = array_ptr[i] + array_ptr[i+1];
          n_chunks /= 2;
        }
      return array_ptr[0];
    }

    const Operation &op;
    const size_type vec_size;

    mutable unsigned int n_chunks;
    unsigned int chunk_size;
    ResultType small_array [threshold_array_allocate];
    std::vector<ResultType> large_array;
    // this variable either points to small_array or large_array depending on
    // the number of threads we want to feed
    mutable ResultType *array_ptr;
  };
#endif



  /**
   * This is the general caller for parallel reduction operations that work in
   * parallel.
   */
  template <typename Operation, typename ResultType>
  void parallel_reduce (const Operation   &op,
                        const size_type    vec_size,
                        ResultType        &result,
                        std_cxx11::shared_ptr<parallel::internal::TBBPartitioner> &partitioner)
  {
#ifdef DEAL_II_WITH_THREADS
    // only go to the parallel function in case there are at least 4 parallel
    // items, otherwise the overhead is too large
    if (vec_size >= 4*internal::Vector::minimum_parallel_grain_size &&
        MultithreadInfo::n_threads() > 1)
      {
        Assert(partitioner.get() != NULL,
               ExcInternalError("Unexpected initialization of Vector that does "
                                "not set the TBB partitioner to a usable state."));
        std_cxx11::shared_ptr<tbb::affinity_partitioner> tbb_partitioner =
          partitioner->acquire_one_partitioner();

        TBBReduceFunctor<Operation,ResultType> generic_functor(op, vec_size);
        tbb::parallel_for (tbb::blocked_range<size_type> (0,
                                                          generic_functor.n_chunks,
                                                          1),
                           generic_functor,
                           *tbb_partitioner);
        partitioner->release_one_partitioner(tbb_partitioner);
        result = generic_functor.do_sum();
      }
    else if (vec_size > 0)
      accumulate_recursive(op,0,vec_size,result);
#else
    accumulate_recursive(op,0,vec_size,result);
#endif
  }
}



template <typename Number>
Vector<Number>::Vector (const Vector<Number> &v)
  :
  Subscriptor(),
  vec_size(v.size()),
  max_vec_size(v.size()),
  val(0)
{
  if (vec_size != 0)
    {
      allocate();
      *this = v;
    }
}



#ifdef DEAL_II_WITH_CXX11
template <typename Number>
Vector<Number>::Vector (Vector<Number> &&v)
  :
  Subscriptor(std::move(v)),
  vec_size(v.vec_size),
  max_vec_size(v.max_vec_size),
  val(v.val),
  thread_loop_partitioner(std::move(v.thread_loop_partitioner))
{
  v.vec_size = 0;
  v.max_vec_size = 0;
  v.val = nullptr;
}
#endif



#ifndef DEAL_II_EXPLICIT_CONSTRUCTOR_BUG
template <typename Number>
template <typename OtherNumber>
Vector<Number>::Vector (const Vector<OtherNumber> &v)
  :
  Subscriptor(),
  vec_size(v.size()),
  max_vec_size(v.size()),
  val(0)
{
  if (vec_size != 0)
    {
      allocate();
      *this = v;
    }
}
#endif


#ifdef DEAL_II_WITH_PETSC

template <typename Number>
Vector<Number>::Vector (const PETScWrappers::Vector &v)
  :
  Subscriptor(),
  vec_size(v.size()),
  max_vec_size(v.size()),
  val(0)
{
  if (vec_size != 0)
    {
      allocate();

      // get a representation of the vector
      // and copy it
      PetscScalar *start_ptr;
      int ierr = VecGetArray (static_cast<const Vec &>(v), &start_ptr);
      AssertThrow (ierr == 0, ExcPETScError(ierr));

      internal::copy (start_ptr, start_ptr+vec_size, begin());

      // restore the representation of the
      // vector
      ierr = VecRestoreArray (static_cast<const Vec &>(v), &start_ptr);
      AssertThrow (ierr == 0, ExcPETScError(ierr));
    }
}



template <typename Number>
Vector<Number>::Vector (const PETScWrappers::MPI::Vector &v)
  :
  Subscriptor(),
  vec_size(0),
  max_vec_size(0),
  val(0)
{
  if (v.size() != 0)
    {
      // do this in a two-stage process:
      // first convert to a sequential petsc
      // vector, then copy that
      PETScWrappers::Vector seq (v);
      *this = seq;
    }
}

#endif


#ifdef DEAL_II_WITH_TRILINOS

template <typename Number>
Vector<Number>::Vector (const TrilinosWrappers::MPI::Vector &v)
  :
  Subscriptor(),
  vec_size(v.size()),
  max_vec_size(v.size()),
  val(0)
{
  if (vec_size != 0)
    {
      allocate();

      // Copy the distributed vector to
      // a local one at all
      // processors. TODO: There could
      // be a better solution than
      // this, but it has not yet been
      // found.
      TrilinosWrappers::Vector localized_vector (v);

      // get a representation of the vector
      // and copy it
      TrilinosScalar **start_ptr;

      int ierr = localized_vector.trilinos_vector().ExtractView (&start_ptr);
      AssertThrow (ierr == 0, ExcTrilinosError(ierr));

      std::copy (start_ptr[0], start_ptr[0]+vec_size, begin());
    }
}



template <typename Number>
Vector<Number>::Vector (const TrilinosWrappers::Vector &v)
  :
  Subscriptor(),
  vec_size(v.size()),
  max_vec_size(v.size()),
  val(0)
{
  if (vec_size != 0)
    {
      allocate();

      // get a representation of the vector
      // and copy it
      TrilinosScalar **start_ptr;

      int ierr = v.trilinos_vector().ExtractView (&start_ptr);
      AssertThrow (ierr == 0, ExcTrilinosError(ierr));

      std::copy (start_ptr[0], start_ptr[0]+vec_size, begin());
    }
}

#endif


template <typename Number>
inline
Vector<Number> &
Vector<Number>::operator= (const Vector<Number> &v)
{
  if (PointerComparison::equal(this, &v))
    return *this;

  thread_loop_partitioner = v.thread_loop_partitioner;
  if (vec_size != v.vec_size)
    reinit (v, true);

  dealii::internal::Vector_copy<Number,Number> copier;
  copier.dst = val;
  copier.src = v.val;
  internal::parallel_for(copier,vec_size,thread_loop_partitioner);

  return *this;
}



#ifdef DEAL_II_WITH_CXX11
template <typename Number>
inline
Vector<Number> &
Vector<Number>::operator= (Vector<Number> &&v)
{
  Subscriptor::operator=(std::move(v));

  if (val) deallocate();

  vec_size = v.vec_size;
  max_vec_size = v.max_vec_size;
  val = v.val;
  thread_loop_partitioner = std::move(v.thread_loop_partitioner);

  v.vec_size = 0;
  v.max_vec_size = 0;
  v.val = nullptr;

  return *this;
}
#endif



template <typename Number>
template <typename Number2>
inline
Vector<Number> &
Vector<Number>::operator= (const Vector<Number2> &v)
{
  thread_loop_partitioner = v.thread_loop_partitioner;
  if (vec_size != v.vec_size)
    reinit (v, true);

  dealii::internal::Vector_copy<Number,Number2> copier;
  copier.dst = val;
  copier.src = v.val;
  internal::parallel_for(copier,vec_size,thread_loop_partitioner);

  return *this;
}



template <typename Number>
inline
void Vector<Number>::reinit (const size_type n,
                             const bool omit_zeroing_entries)
{
  if (n==0)
    {
      if (val) deallocate();
      val = 0;
      max_vec_size = vec_size = 0;
      thread_loop_partitioner.reset(new parallel::internal::TBBPartitioner());
      return;
    };

  if (n>max_vec_size)
    {
      if (val) deallocate();
      max_vec_size = n;
      allocate();
    };

  if (vec_size != n)
    {
      vec_size = n;

      // only reset the partitioner if we actually expect a significant vector
      // size
      if (vec_size >= 4*internal::Vector::minimum_parallel_grain_size)
        thread_loop_partitioner.reset(new parallel::internal::TBBPartitioner());
    }

  if (omit_zeroing_entries == false)
    *this = static_cast<Number>(0);
}



template <typename Number>
template <typename Number2>
void Vector<Number>::reinit (const Vector<Number2> &v,
                             const bool omit_zeroing_entries)
{
  thread_loop_partitioner = v.thread_loop_partitioner;

  if (v.vec_size==0)
    {
      if (val) deallocate();
      val = 0;
      max_vec_size = vec_size = 0;
      return;
    };

  if (v.vec_size>max_vec_size)
    {
      if (val) deallocate();
      max_vec_size = v.vec_size;
      allocate();
    };
  vec_size = v.vec_size;
  if (omit_zeroing_entries == false)
    *this = static_cast<Number>(0);
}



template <typename Number>
bool
Vector<Number>::all_zero () const
{
  Assert (vec_size!=0, ExcEmptyObject());

  for (size_type i=0; i<vec_size; ++i)
    if (val[i] != Number(0))
      return false;
  return true;
}



template <typename Number>
bool
Vector<Number>::is_non_negative () const
{
  Assert (vec_size!=0, ExcEmptyObject());

  for (size_type i=0; i<vec_size; ++i)
    if ( ! internal::is_non_negative (val[i]))
      return false;

  return true;
}



template <typename Number>
Vector<Number> &
Vector<Number>::operator= (const Number s)
{
  AssertIsFinite(s);
  if (s != Number())
    Assert (vec_size!=0, ExcEmptyObject());

  internal::Vector_set<Number> setter;
  setter.dst = val;
  setter.value = s;

  internal::parallel_for(setter,vec_size,thread_loop_partitioner);

  return *this;
}



#ifdef DEAL_II_BOOST_BIND_COMPILER_BUG
template <>
Vector<std::complex<float> > &
Vector<std::complex<float> >::operator= (const std::complex<float> s)
{
  AssertIsFinite(s);
  if (s != std::complex<float>())
    Assert (vec_size!=0, ExcEmptyObject());
  if (vec_size!=0)
    std::fill (begin(), end(), s);

  return *this;
}
#endif



template <typename Number>
Vector<Number> &Vector<Number>::operator *= (const Number factor)
{
  AssertIsFinite(factor);

  Assert (vec_size!=0, ExcEmptyObject());

  internal::Vectorization_multiply_factor<Number> vector_multiply;
  vector_multiply.val = val;
  vector_multiply.factor = factor;

  internal::parallel_for(vector_multiply,vec_size,thread_loop_partitioner);

  return *this;
}



template <typename Number>
void
Vector<Number>::add (const Number a,
                     const Vector<Number> &v)
{
  AssertIsFinite(a);

  Assert (vec_size!=0, ExcEmptyObject());
  Assert (vec_size == v.vec_size, ExcDimensionMismatch(vec_size, v.vec_size));

  internal::Vectorization_add_av<Number> vector_add_av;
  vector_add_av.val = val;
  vector_add_av.v_val = v.val;
  vector_add_av.factor = a;
  internal::parallel_for(vector_add_av,vec_size,thread_loop_partitioner);
}



template <typename Number>
void
Vector<Number>::sadd (const Number x,
                      const Number a,
                      const Vector<Number> &v)
{
  AssertIsFinite(x);
  AssertIsFinite(a);

  Assert (vec_size!=0, ExcEmptyObject());
  Assert (vec_size == v.vec_size, ExcDimensionMismatch(vec_size, v.vec_size));

  internal::Vectorization_sadd_xav<Number> vector_sadd_xav;
  vector_sadd_xav.val = val;
  vector_sadd_xav.v_val = v.val;
  vector_sadd_xav.a = a;
  vector_sadd_xav.x = x;
  internal::parallel_for(vector_sadd_xav,vec_size,thread_loop_partitioner);
}



namespace internal
{
  namespace Vector
  {
  }
}



template <typename Number>
template <typename Number2>
Number Vector<Number>::operator * (const Vector<Number2> &v) const
{
  Assert (vec_size!=0, ExcEmptyObject());

  if (PointerComparison::equal (this, &v))
    return norm_sqr();

  Assert (vec_size == v.size(),
          ExcDimensionMismatch(vec_size, v.size()));

  Number sum;
  internal::Dot<Number,Number2> dot;
  dot.X = val;
  dot.Y = v.val;
  internal::parallel_reduce (dot, vec_size, sum, thread_loop_partitioner);
  AssertIsFinite(sum);

  return sum;
}



template <typename Number>
typename Vector<Number>::real_type
Vector<Number>::norm_sqr () const
{
  Assert (vec_size!=0, ExcEmptyObject());

  real_type sum;
  internal::Norm2<Number,real_type> norm2;
  norm2.X = val;
  internal::parallel_reduce (norm2, vec_size, sum, thread_loop_partitioner);

  AssertIsFinite(sum);

  return sum;
}



template <typename Number>
Number Vector<Number>::mean_value () const
{
  Assert (vec_size!=0, ExcEmptyObject());

  Number sum;
  internal::MeanValue<Number> mean;
  mean.X = val;
  internal::parallel_reduce (mean, vec_size, sum, thread_loop_partitioner);

  return sum / real_type(size());
}



template <typename Number>
typename Vector<Number>::real_type
Vector<Number>::l1_norm () const
{
  Assert (vec_size!=0, ExcEmptyObject());

  real_type sum;
  internal::Norm1<Number, real_type> norm1;
  norm1.X = val;
  internal::parallel_reduce (norm1, vec_size, sum, thread_loop_partitioner);

  return sum;
}



template <typename Number>
typename Vector<Number>::real_type
Vector<Number>::l2_norm () const
{
  // if l2_norm()^2 is finite and non-zero, the answer is computed as
  // std::sqrt(norm_sqr()). If norm_sqr() is infinite or zero, the l2 norm
  // might still be finite. In that case, recompute it (this is a rare case,
  // so working on the vector twice is uncritical and paid off by the extended
  // precision) using the BLAS approach with a weight, see e.g. dnrm2.f.
  Assert (vec_size!=0, ExcEmptyObject());

  real_type norm_square;
  internal::Norm2<Number, real_type> norm2;
  norm2.X = val;
  internal::parallel_reduce (norm2, vec_size, norm_square,
                             thread_loop_partitioner);
  if (numbers::is_finite(norm_square) &&
      norm_square >= std::numeric_limits<real_type>::min())
    return std::sqrt(norm_square);
  else
    {
      real_type scale = 0.;
      real_type sum = 1.;
      for (size_type i=0; i<vec_size; ++i)
        {
          if (val[i] != Number())
            {
              const real_type abs_x =
                numbers::NumberTraits<Number>::abs(val[i]);
              if (scale < abs_x)
                {
                  sum = 1. + sum * (scale/abs_x) * (scale/abs_x);
                  scale = abs_x;
                }
              else
                sum += (abs_x/scale) * (abs_x/scale);
            }
        }
      AssertIsFinite(scale*std::sqrt(sum));
      return scale * std::sqrt(sum);
    }
}



template <typename Number>
typename Vector<Number>::real_type
Vector<Number>::lp_norm (const real_type p) const
{
  Assert (vec_size!=0, ExcEmptyObject());

  if (p == 1.)
    return l1_norm();
  else if (p == 2.)
    return l2_norm();

  real_type sum;
  internal::NormP<Number, real_type> normp;
  normp.X = val;
  normp.p = p;
  internal::parallel_reduce (normp, vec_size, sum, thread_loop_partitioner);

  if (numbers::is_finite(sum) && sum >= std::numeric_limits<real_type>::min())
    return std::pow(sum, static_cast<real_type>(1./p));
  else
    {
      real_type scale = 0.;
      real_type sum = 1.;
      for (size_type i=0; i<vec_size; ++i)
        {
          if (val[i] != Number())
            {
              const real_type abs_x =
                numbers::NumberTraits<Number>::abs(val[i]);
              if (scale < abs_x)
                {
                  sum = 1. + sum * std::pow(scale/abs_x, p);
                  scale = abs_x;
                }
              else
                sum += std::pow(abs_x/scale, p);
            }
        }
      return scale * std::pow(sum, static_cast<real_type>(1./p));
    }
}



template <>
Vector<int>::real_type
Vector<int>::lp_norm (const real_type) const
{
  Assert(false, ExcMessage("No lp norm for integer vectors"));
  return -1;
}



template <typename Number>
typename Vector<Number>::real_type
Vector<Number>::linfty_norm () const
{
  Assert (vec_size!=0, ExcEmptyObject());

  real_type max = 0.;

  for (size_type i=0; i<vec_size; ++i)
    max = std::max (numbers::NumberTraits<Number>::abs(val[i]), max);

  return max;
}



template <typename Number>
Number
Vector<Number>::add_and_dot (const Number          a,
                             const Vector<Number> &V,
                             const Vector<Number> &W)
{
  Assert (vec_size!=0, ExcEmptyObject());
  AssertDimension (vec_size, V.size());
  AssertDimension (vec_size, W.size());

  Number sum;
  internal::AddAndDot<Number> adder;
  adder.X = val;
  adder.a = a;
  adder.V = V.val;
  adder.W = W.val;
  internal::parallel_reduce (adder, vec_size, sum, thread_loop_partitioner);
  AssertIsFinite(sum);

  return sum;
}



template <typename Number>
Vector<Number> &Vector<Number>::operator += (const Vector<Number> &v)
{
  Assert (vec_size!=0, ExcEmptyObject());

  add (v);
  return *this;
}



template <typename Number>
Vector<Number> &Vector<Number>::operator -= (const Vector<Number> &v)
{
  Assert (vec_size!=0, ExcEmptyObject());
  Assert (vec_size == v.vec_size, ExcDimensionMismatch(vec_size, v.vec_size));

  internal::Vectorization_subtract_v<Number> vector_subtract;
  vector_subtract.val = val;
  vector_subtract.v_val = v.val;
  internal::parallel_for(vector_subtract,vec_size,thread_loop_partitioner);

  return *this;
}



template <typename Number>
void Vector<Number>::add (const Number v)
{
  Assert (vec_size!=0, ExcEmptyObject());

  internal::Vectorization_add_factor<Number> vector_add;
  vector_add.val = val;
  vector_add.factor = v;
  internal::parallel_for(vector_add,vec_size,thread_loop_partitioner);
}



template <typename Number>
void Vector<Number>::add (const Vector<Number> &v)
{
  Assert (vec_size!=0, ExcEmptyObject());
  Assert (vec_size == v.vec_size, ExcDimensionMismatch(vec_size, v.vec_size));

  internal::Vectorization_add_v<Number> vector_add;
  vector_add.val = val;
  vector_add.v_val = v.val;
  internal::parallel_for(vector_add,vec_size,thread_loop_partitioner);
}



template <typename Number>
void Vector<Number>::add (const Number a, const Vector<Number> &v,
                          const Number b, const Vector<Number> &w)
{
  AssertIsFinite(a);
  AssertIsFinite(b);

  Assert (vec_size!=0, ExcEmptyObject());
  Assert (vec_size == v.vec_size, ExcDimensionMismatch(vec_size, v.vec_size));
  Assert (vec_size == w.vec_size, ExcDimensionMismatch(vec_size, w.vec_size));

  internal::Vectorization_add_avpbw<Number> vector_add;
  vector_add.val = val;
  vector_add.v_val = v.val;
  vector_add.w_val = w.val;
  vector_add.a = a;
  vector_add.b = b;
  internal::parallel_for(vector_add,vec_size,thread_loop_partitioner);
}



template <typename Number>
void Vector<Number>::sadd (const Number x,
                           const Vector<Number> &v)
{
  AssertIsFinite(x);

  Assert (vec_size!=0, ExcEmptyObject());
  Assert (vec_size == v.vec_size, ExcDimensionMismatch(vec_size, v.vec_size));

  internal::Vectorization_sadd_xv<Number> vector_sadd;
  vector_sadd.val = val;
  vector_sadd.v_val = v.val;
  vector_sadd.x = x;
  internal::parallel_for(vector_sadd,vec_size,thread_loop_partitioner);
}



template <typename Number>
void Vector<Number>::sadd (const Number x, const Number a,
                           const Vector<Number> &v, const Number b,
                           const Vector<Number> &w)
{
  AssertIsFinite(x);
  AssertIsFinite(a);
  AssertIsFinite(b);

  Assert (vec_size!=0, ExcEmptyObject());
  Assert (vec_size == v.vec_size, ExcDimensionMismatch(vec_size, v.vec_size));
  Assert (vec_size == w.vec_size, ExcDimensionMismatch(vec_size, w.vec_size));

  internal::Vectorization_sadd_xavbw<Number> vector_sadd;
  vector_sadd.val = val;
  vector_sadd.v_val = v.val;
  vector_sadd.w_val = w.val;
  vector_sadd.x = x;
  vector_sadd.a = a;
  vector_sadd.b = b;
  internal::parallel_for(vector_sadd,vec_size,thread_loop_partitioner);
}


template <typename Number>
void Vector<Number>::sadd (const Number x, const Number a,
                           const Vector<Number> &v, const Number b,
                           const Vector<Number> &w, const Number c,
                           const Vector<Number> &y)
{
  sadd (x, a, v, b, w);
  add (c, y);
}



template <typename Number>
void Vector<Number>::scale (const Vector<Number> &s)
{
  Assert (vec_size!=0, ExcEmptyObject());
  Assert (vec_size == s.vec_size, ExcDimensionMismatch(vec_size, s.vec_size));

  internal::Vectorization_scale<Number> vector_scale;
  vector_scale.val = val;
  vector_scale.v_val = s.val;
  internal::parallel_for(vector_scale,vec_size,thread_loop_partitioner);
}



template <typename Number>
template <typename Number2>
void Vector<Number>::scale (const Vector<Number2> &s)
{
  Assert (vec_size!=0, ExcEmptyObject());
  Assert (vec_size == s.vec_size, ExcDimensionMismatch(vec_size, s.vec_size));

  for (size_type i=0; i<vec_size; ++i)
    val[i] *= Number(s.val[i]);
}



template <typename Number>
void Vector<Number>::equ (const Number a,
                          const Vector<Number> &u)
{
  AssertIsFinite(a);

  Assert (vec_size!=0, ExcEmptyObject());
  Assert (vec_size == u.vec_size, ExcDimensionMismatch(vec_size, u.vec_size));

  internal::Vectorization_equ_au<Number> vector_equ;
  vector_equ.val = val;
  vector_equ.u_val = u.val;
  vector_equ.a = a;
  internal::parallel_for(vector_equ,vec_size,thread_loop_partitioner);
}



template <typename Number>
template <typename Number2>
void Vector<Number>::equ (const Number a,
                          const Vector<Number2> &u)
{
  AssertIsFinite(a);

  Assert (vec_size!=0, ExcEmptyObject());
  Assert (vec_size == u.vec_size, ExcDimensionMismatch(vec_size, u.vec_size));

  // set the result vector to a*u. we have to
  // convert the elements of u to the type of
  // the result vector. this is necessary
  // because
  // operator*(complex<float>,complex<double>)
  // is not defined by default
  for (size_type i=0; i<vec_size; ++i)
    val[i] = a * Number(u.val[i]);
}



template <typename Number>
void Vector<Number>::equ (const Number a, const Vector<Number> &u,
                          const Number b, const Vector<Number> &v)
{
  AssertIsFinite(a);
  AssertIsFinite(b);

  Assert (vec_size!=0, ExcEmptyObject());
  Assert (vec_size == u.vec_size, ExcDimensionMismatch(vec_size, u.vec_size));
  Assert (vec_size == v.vec_size, ExcDimensionMismatch(vec_size, v.vec_size));

  internal::Vectorization_equ_aubv<Number> vector_equ;
  vector_equ.val = val;
  vector_equ.u_val = u.val;
  vector_equ.v_val = v.val;
  vector_equ.a = a;
  vector_equ.b = b;
  internal::parallel_for(vector_equ,vec_size,thread_loop_partitioner);
}


template <typename Number>
void Vector<Number>::equ (const Number a, const Vector<Number> &u,
                          const Number b, const Vector<Number> &v,
                          const Number c, const Vector<Number> &w)
{
  Assert (vec_size!=0, ExcEmptyObject());
  Assert (vec_size == u.vec_size, ExcDimensionMismatch(vec_size, u.vec_size));
  Assert (vec_size == v.vec_size, ExcDimensionMismatch(vec_size, v.vec_size));
  Assert (vec_size == w.vec_size, ExcDimensionMismatch(vec_size, w.vec_size));

  internal::Vectorization_equ_aubvcw<Number> vector_equ;
  vector_equ.val = val;
  vector_equ.u_val = u.val;
  vector_equ.v_val = v.val;
  vector_equ.w_val = w.val;
  vector_equ.a = a;
  vector_equ.b = b;
  vector_equ.c = c;
  internal::parallel_for(vector_equ,vec_size,thread_loop_partitioner);
}


template <typename Number>
void Vector<Number>::ratio (const Vector<Number> &a,
                            const Vector<Number> &b)
{
  Assert (vec_size!=0, ExcEmptyObject());
  Assert (a.vec_size == b.vec_size,
          ExcDimensionMismatch (a.vec_size, b.vec_size));

  // no need to reinit with zeros, since
  // we overwrite them anyway
  reinit (a.size(), true);

  internal::Vectorization_ratio<Number> vector_ratio;
  vector_ratio.val = val;
  vector_ratio.a_val = a.val;
  vector_ratio.b_val = b.val;
  internal::parallel_for(vector_ratio,vec_size,thread_loop_partitioner);
}



template <typename Number>
Vector<Number> &
Vector<Number>::operator= (const BlockVector<Number> &v)
{
  if (v.size() != vec_size)
    reinit (v.size(), true);

  size_type this_index = 0;
  for (size_type b=0; b<v.n_blocks(); ++b)
    for (size_type i=0; i<v.block(b).size(); ++i, ++this_index)
      val[this_index] = v.block(b)(i);

  return *this;
}



#ifdef DEAL_II_WITH_PETSC

template <typename Number>
Vector<Number> &
Vector<Number>::operator= (const PETScWrappers::Vector &v)
{
  if (v.size() != vec_size)
    reinit (v.size(), true);
  if (vec_size != 0)
    {
      // get a representation of the vector
      // and copy it
      PetscScalar *start_ptr;
      int ierr = VecGetArray (static_cast<const Vec &>(v), &start_ptr);
      AssertThrow (ierr == 0, ExcPETScError(ierr));

      internal::copy (start_ptr, start_ptr+vec_size, begin());

      // restore the representation of the
      // vector
      ierr = VecRestoreArray (static_cast<const Vec &>(v), &start_ptr);
      AssertThrow (ierr == 0, ExcPETScError(ierr));
    }

  return *this;
}



template <typename Number>
Vector<Number> &
Vector<Number>::operator= (const PETScWrappers::MPI::Vector &v)
{
  // do this in a two-stage process:
  // first convert to a sequential petsc
  // vector, then copy that
  PETScWrappers::Vector seq (v);
  *this = seq;

  return *this;
}

#endif


#ifdef DEAL_II_WITH_TRILINOS

template <typename Number>
Vector<Number> &
Vector<Number>::operator= (const TrilinosWrappers::MPI::Vector &v)
{
  // Generate a localized version
  // of the Trilinos vectors and
  // then call the other =
  // operator.
  TrilinosWrappers::Vector localized_vector (v);
  *this = localized_vector;
  return *this;
}



template <typename Number>
Vector<Number> &
Vector<Number>::operator= (const TrilinosWrappers::Vector &v)
{
  if (v.size() != vec_size)
    reinit (v.size(), true);
  if (vec_size != 0)
    {
      // get a representation of the vector
      // and copy it
      TrilinosScalar **start_ptr;
      int ierr = v.trilinos_vector().ExtractView (&start_ptr);
      AssertThrow (ierr == 0, ExcTrilinosError(ierr));

      std::copy (start_ptr[0], start_ptr[0]+vec_size, begin());
    }

  return *this;
}

#endif

template <typename Number>
template <typename Number2>
bool
Vector<Number>::operator== (const Vector<Number2> &v) const
{
  Assert (vec_size!=0, ExcEmptyObject());
  Assert (vec_size == v.size(), ExcDimensionMismatch(vec_size, v.size()));

  // compare the two vector. we have to
  // convert the elements of v to the type of
  // the result vector. this is necessary
  // because
  // operator==(complex<float>,complex<double>)
  // is not defined by default
  for (size_type i=0; i<vec_size; ++i)
    if (val[i] != Number(v.val[i]))
      return false;

  return true;
}



template <typename Number>
void Vector<Number>::print (const char *format) const
{
  Assert (vec_size!=0, ExcEmptyObject());

  for (size_type j=0; j<size(); ++j)
    internal::print (val[j], format);
  std::printf ("\n");
}



template <typename Number>
void Vector<Number>::print (std::ostream      &out,
                            const unsigned int precision,
                            const bool         scientific,
                            const bool         across) const
{
  Assert (vec_size!=0, ExcEmptyObject());
  AssertThrow (out, ExcIO());

  std::ios::fmtflags old_flags = out.flags();
  unsigned int old_precision = out.precision (precision);

  out.precision (precision);
  if (scientific)
    out.setf (std::ios::scientific, std::ios::floatfield);
  else
    out.setf (std::ios::fixed, std::ios::floatfield);

  if (across)
    for (size_type i=0; i<size(); ++i)
      out << val[i] << ' ';
  else
    for (size_type i=0; i<size(); ++i)
      out << val[i] << std::endl;
  out << std::endl;

  AssertThrow (out, ExcIO());
  // reset output format
  out.flags (old_flags);
  out.precision(old_precision);
}



template <typename Number>
void
Vector<Number>::print (LogStream &out, const unsigned int width, const bool across) const
{
  Assert (vec_size!=0, ExcEmptyObject());

  if (across)
    for (size_type i=0; i<size(); ++i)
      out << std::setw(width) << val[i] << ' ';
  else
    for (size_type i=0; i<size(); ++i)
      out << val[i] << std::endl;
  out << std::endl;
}


template <typename Number>
void Vector<Number>::block_write (std::ostream &out) const
{
  AssertThrow (out, ExcIO());

  // other version of the following
  //  out << size() << std::endl << '[';
  // reason: operator<< seems to use
  // some resources that lead to
  // problems in a multithreaded
  // environment
  const size_type sz = size();
  char buf[16];

#ifdef DEAL_II_WITH_64BIT_INDICES
  std::sprintf(buf, "%llu", sz);
#else
  std::sprintf(buf, "%u", sz);
#endif
  std::strcat(buf, "\n[");

  out.write(buf, std::strlen(buf));
  out.write (reinterpret_cast<const char *>(begin()),
             reinterpret_cast<const char *>(end())
             - reinterpret_cast<const char *>(begin()));

  // out << ']';
  const char outro = ']';
  out.write (&outro, 1);

  AssertThrow (out, ExcIO());
}



template <typename Number>
void Vector<Number>::block_read (std::istream &in)
{
  AssertThrow (in, ExcIO());

  size_type sz;

  char buf[16];


  in.getline(buf,16,'\n');
  sz=std::atoi(buf);

  // fast initialization, since the
  // data elements are overwritten anyway
  reinit (sz, true);

  char c;
  //  in >> c;
  in.read (&c, 1);
  AssertThrow (c=='[', ExcIO());

  in.read (reinterpret_cast<char *>(begin()),
           reinterpret_cast<const char *>(end())
           - reinterpret_cast<const char *>(begin()));

  //  in >> c;
  in.read (&c, 1);
  AssertThrow (c==']', ExcIO());
}



template <typename Number>
IndexSet
Vector<Number>::locally_owned_elements() const
{
  return complete_index_set(size());
}



template <typename Number>
std::size_t
Vector<Number>::memory_consumption () const
{
  return sizeof(*this) + (max_vec_size * sizeof(Number));
}



template <typename Number>
void
Vector<Number>::allocate()
{
  // make sure that we don't create a memory leak
  Assert (val == 0, ExcInternalError());

  // then allocate memory with the proper alignment requirements of 64 bytes
  Utilities::System::posix_memalign ((void **)&val, 64, sizeof(Number)*max_vec_size);
}



template <typename Number>
void
Vector<Number>::deallocate()
{
  free(val);
  val = 0;
}

DEAL_II_NAMESPACE_CLOSE

#endif
