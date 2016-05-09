// ---------------------------------------------------------------------
//
// Copyright (C) 2016 by the deal.II authors
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


#ifndef dealii__vector_operations_internal_h
#define dealii__vector_operations_internal_h

#include <deal.II/base/multithread_info.h>
#include <deal.II/base/parallel.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/base/vectorization.h>

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
    (void)partitioner;
#endif
  }


  // Define the functors necessary to use SIMD with TBB. we also include the
  // simple copy and set operations

  template <typename Number>
  struct Vector_set
  {
    Vector_set(Number value, Number *dst)
      :
      value(value),
      dst(dst)
    {}

    void operator() (const size_type begin, const size_type end) const
    {
      if (value == Number())
        std::memset (dst+begin,0,(end-begin)*sizeof(Number));
      else
        std::fill (dst+begin, dst+end, value);
    }

    Number value;
    Number *dst;
  };

  template <typename Number, typename OtherNumber>
  struct Vector_copy
  {
    Vector_copy(const OtherNumber *src, Number *dst)
      :
      src(src),
      dst(dst)
    {}

    void operator() (const size_type begin, const size_type end) const
    {
      if (types_are_equal<Number,OtherNumber>::value)
        std::memcpy(dst+begin, src+begin, (end-begin)*sizeof(Number));
      else
        {
          DEAL_II_OPENMP_SIMD_PRAGMA
          for (size_type i=begin; i<end; ++i)
            dst[i] = src[i];
        }
    }

    const OtherNumber *src;
    Number *dst;
  };

  template <typename Number>
  struct Vectorization_multiply_factor
  {
    Vectorization_multiply_factor(Number *val, Number factor)
      :
      val(val),
      factor(factor)
    {}

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

    Number *val;
    Number factor;
  };

  template <typename Number>
  struct Vectorization_add_av
  {
    Vectorization_add_av(Number *val, Number *v_val, Number factor)
      :
      val(val),
      v_val(v_val),
      factor(factor)
    {}

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

    Number *val;
    Number *v_val;
    Number factor;
  };

  template <typename Number>
  struct Vectorization_sadd_xav
  {
    Vectorization_sadd_xav(Number *val, Number *v_val, Number a, Number x)
      :
      val(val),
      v_val(v_val),
      a(a),
      x(x)
    {}

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

    Number *val;
    Number *v_val;
    Number a;
    Number x;
  };

  template <typename Number>
  struct Vectorization_subtract_v
  {
    Vectorization_subtract_v(Number *val, Number *v_val)
      :
      val(val),
      v_val(v_val)
    {}

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

    Number *val;
    Number *v_val;
  };

  template <typename Number>
  struct Vectorization_add_factor
  {
    Vectorization_add_factor(Number *val, Number factor)
      :
      val(val),
      factor(factor)
    {}

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

    Number *val;
    Number factor;
  };

  template <typename Number>
  struct Vectorization_add_v
  {
    Vectorization_add_v(Number *val, Number *v_val)
      :
      val(val),
      v_val(v_val)
    {}

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

    Number *val;
    Number *v_val;
  };

  template <typename Number>
  struct Vectorization_add_avpbw
  {
    Vectorization_add_avpbw(Number *val, Number *v_val, Number *w_val, Number a, Number b)
      :
      val(val),
      v_val(v_val),
      w_val(w_val),
      a(a),
      b(b)
    {}

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

    Number *val;
    Number *v_val;
    Number *w_val;
    Number a;
    Number b;
  };

  template <typename Number>
  struct Vectorization_sadd_xv
  {
    Vectorization_sadd_xv(Number *val, Number *v_val, Number x)
      :
      val(val),
      v_val(v_val),
      x(x)
    {}

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

    Number *val;
    Number *v_val;
    Number x;
  };

  template <typename Number>
  struct Vectorization_sadd_xavbw
  {
    Vectorization_sadd_xavbw(Number *val, Number *v_val, Number *w_val,
                             Number x, Number a, Number b)
      :
      val(val),
      v_val(v_val),
      w_val(w_val),
      x(x),
      a(a),
      b(b)
    {}

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

    Number *val;
    Number *v_val;
    Number *w_val;
    Number x;
    Number a;
    Number b;
  };

  template <typename Number>
  struct Vectorization_scale
  {
    Vectorization_scale(Number *val, Number *v_val)
      :
      val(val),
      v_val(v_val)
    {}

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

    Number *val;
    Number *v_val;
  };

  template <typename Number>
  struct Vectorization_equ_au
  {
    Vectorization_equ_au(Number *val, Number *u_val, Number a)
      :
      val(val),
      u_val(u_val),
      a(a)
    {}

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

    Number *val;
    Number *u_val;
    Number a;
  };

  template <typename Number>
  struct Vectorization_equ_aubv
  {
    Vectorization_equ_aubv(Number *val, Number *u_val, Number *v_val,
                           Number a, Number b)
      :
      val(val),
      u_val(u_val),
      v_val(v_val),
      a(a),
      b(b)
    {}

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

    Number *val;
    Number *u_val;
    Number *v_val;
    Number a;
    Number b;
  };

  template <typename Number>
  struct Vectorization_equ_aubvcw
  {
    Vectorization_equ_aubvcw(Number *val, Number *u_val, Number *v_val,
                             Number *w_val, Number a, Number b, Number c)
      :
      val(val),
      u_val(u_val),
      v_val(v_val),
      w_val(w_val),
      a(a),
      b(b),
      c(c)
    {}

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

    Number *val;
    Number *u_val;
    Number *v_val;
    Number *w_val;
    Number a;
    Number b;
    Number c;
  };

  template <typename Number>
  struct Vectorization_ratio
  {
    Vectorization_ratio(Number *val, Number *a_val, Number *b_val)
      :
      val(val),
      a_val(a_val),
      b_val(b_val)
    {}

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

    Number *val;
    Number *a_val;
    Number *b_val;
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

    Dot(const Number *X, const Number2 *Y)
      :
      X(X),
      Y(Y)
    {}

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

    Norm2(const Number *X)
      :
      X(X)
    {}

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

    Norm1(const Number *X)
      :
      X(X)
    {}

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

    NormP(const Number *X, RealType p)
      :
      X(X),
      p(p)
    {}

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

    MeanValue(const Number *X)
      :
      X(X)
    {}

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

    AddAndDot(Number *X, const Number *V, const Number *W, Number a)
      :
      X(X),
      V(V),
      W(W),
      a(a)
    {}

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
    (void)partitioner;
#endif
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif
