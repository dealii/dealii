// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


#ifndef dealii_vector_operations_internal_h
#define dealii_vector_operations_internal_h

#include <deal.II/base/config.h>

#include <deal.II/base/memory_space.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/parallel.h>
#include <deal.II/base/types.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/lac/cuda_kernels.h>
#include <deal.II/lac/cuda_kernels.templates.h>
#include <deal.II/lac/vector_operation.h>

#include <cstdio>
#include <cstring>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace VectorOperations
  {
    using size_type = types::global_dof_index;

    template <typename T>
    bool
    is_non_negative(const T &t)
    {
      return t >= 0;
    }


    template <typename T>
    bool
    is_non_negative(const std::complex<T> &)
    {
      Assert(false, ExcMessage("Complex numbers do not have an ordering."));

      return false;
    }


    // call std::copy, except for in
    // the case where we want to copy
    // from std::complex to a
    // non-complex type
    template <typename T, typename U>
    void
    copy(const T *begin, const T *end, U *dest)
    {
      std::copy(begin, end, dest);
    }

    template <typename T, typename U>
    void
    copy(const std::complex<T> *begin,
         const std::complex<T> *end,
         std::complex<U> *      dest)
    {
      std::copy(begin, end, dest);
    }

    template <typename T, typename U>
    void
    copy(const std::complex<T> *, const std::complex<T> *, U *)
    {
      Assert(false,
             ExcMessage("Can't convert a vector of complex numbers "
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
      TBBForFunctor(Functor &       functor,
                    const size_type start,
                    const size_type end)
        : functor(functor)
        , start(start)
        , end(end)
      {
        const size_type vec_size = end - start;
        // set chunk size for sub-tasks
        const unsigned int gs =
          internal::VectorImplementation::minimum_parallel_grain_size;
        n_chunks =
          std::min(static_cast<size_type>(4 * MultithreadInfo::n_threads()),
                   vec_size / gs);
        chunk_size = vec_size / n_chunks;

        // round to next multiple of 512 (or minimum grain size if that happens
        // to be smaller). this is advantageous because our accumulation
        // algorithms favor lengths of a power of 2 due to pairwise summation ->
        // at most one 'oddly' sized chunk
        if (chunk_size > 512)
          chunk_size = ((chunk_size + 511) / 512) * 512;
        n_chunks = (vec_size + chunk_size - 1) / chunk_size;
        AssertIndexRange((n_chunks - 1) * chunk_size, vec_size);
        AssertIndexRange(vec_size, n_chunks * chunk_size + 1);
      }

      void
      operator()(const tbb::blocked_range<size_type> &range) const
      {
        const size_type r_begin = start + range.begin() * chunk_size;
        const size_type r_end = std::min(start + range.end() * chunk_size, end);
        functor(r_begin, r_end);
      }

      Functor &       functor;
      const size_type start;
      const size_type end;
      unsigned int    n_chunks;
      size_type       chunk_size;
    };
#endif

    template <typename Functor>
    void
    parallel_for(
      Functor &       functor,
      const size_type start,
      const size_type end,
      const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
        &partitioner)
    {
#ifdef DEAL_II_WITH_THREADS
      const size_type vec_size = end - start;
      // only go to the parallel function in case there are at least 4 parallel
      // items, otherwise the overhead is too large
      if (vec_size >=
            4 * internal::VectorImplementation::minimum_parallel_grain_size &&
          MultithreadInfo::n_threads() > 1)
        {
          Assert(partitioner.get() != nullptr,
                 ExcInternalError(
                   "Unexpected initialization of Vector that does "
                   "not set the TBB partitioner to a usable state."));
          std::shared_ptr<tbb::affinity_partitioner> tbb_partitioner =
            partitioner->acquire_one_partitioner();

          TBBForFunctor<Functor> generic_functor(functor, start, end);
          // We use a minimum grain size of 1 here since the grains at this
          // stage of dividing the work refer to the number of vector chunks
          // that are processed by (possibly different) threads in the
          // parallelized for loop (i.e., they do not refer to individual
          // vector entries). The number of chunks here is calculated inside
          // TBBForFunctor. See also GitHub issue #2496 for further discussion
          // of this strategy.
          ::dealii::parallel::internal::parallel_for(
            static_cast<size_type>(0),
            static_cast<size_type>(generic_functor.n_chunks),
            generic_functor,
            1,
            tbb_partitioner);
          partitioner->release_one_partitioner(tbb_partitioner);
        }
      else if (vec_size > 0)
        functor(start, end);
#else
      functor(start, end);
      (void)partitioner;
#endif
    }


    // Define the functors necessary to use SIMD with TBB. we also include the
    // simple copy and set operations

    template <typename Number>
    struct Vector_set
    {
      Vector_set(const Number value, Number *const dst)
        : value(value)
        , dst(dst)
      {
        Assert(dst != nullptr, ExcInternalError());
      }

      void
      operator()(const size_type begin, const size_type end) const
      {
        Assert(end >= begin, ExcInternalError());

        if (value == Number())
          {
#ifdef DEAL_II_HAVE_CXX17
            if constexpr (std::is_trivial<Number>::value)
#else
            if (std::is_trivial<Number>::value)
#endif
              {
                std::memset(dst + begin, 0, sizeof(Number) * (end - begin));
                return;
              }
          }
        std::fill(dst + begin, dst + end, value);
      }

      const Number  value;
      Number *const dst;
    };

    template <typename Number, typename OtherNumber>
    struct Vector_copy
    {
      Vector_copy(const OtherNumber *const src, Number *const dst)
        : src(src)
        , dst(dst)
      {
        Assert(src != nullptr, ExcInternalError());
        Assert(dst != nullptr, ExcInternalError());
      }

      void
      operator()(const size_type begin, const size_type end) const
      {
        Assert(end >= begin, ExcInternalError());

#if __GNUG__ && __GNUC__ < 5
        if (__has_trivial_copy(Number) &&
            std::is_same<Number, OtherNumber>::value)
#else
#  ifdef DEAL_II_HAVE_CXX17
        if constexpr (std::is_trivially_copyable<Number>() &&
                      std::is_same<Number, OtherNumber>::value)
#  else
        if (std::is_trivially_copyable<Number>() &&
            std::is_same<Number, OtherNumber>::value)
#  endif
#endif
          std::memcpy(dst + begin, src + begin, (end - begin) * sizeof(Number));
        else
          {
            DEAL_II_OPENMP_SIMD_PRAGMA
            for (size_type i = begin; i < end; ++i)
              dst[i] = src[i];
          }
      }

      const OtherNumber *const src;
      Number *const            dst;
    };

    template <typename Number>
    struct Vectorization_multiply_factor
    {
      Vectorization_multiply_factor(Number *const val, const Number factor)
        : val(val)
        , factor(factor)
      {}

      void
      operator()(const size_type begin, const size_type end) const
      {
        if (::dealii::parallel::internal::EnableOpenMPSimdFor<Number>::value)
          {
            DEAL_II_OPENMP_SIMD_PRAGMA
            for (size_type i = begin; i < end; ++i)
              val[i] *= factor;
          }
        else
          {
            for (size_type i = begin; i < end; ++i)
              val[i] *= factor;
          }
      }

      Number *const val;
      const Number  factor;
    };

    template <typename Number>
    struct Vectorization_add_av
    {
      Vectorization_add_av(Number *const       val,
                           const Number *const v_val,
                           const Number        factor)
        : val(val)
        , v_val(v_val)
        , factor(factor)
      {}

      void
      operator()(const size_type begin, const size_type end) const
      {
        if (::dealii::parallel::internal::EnableOpenMPSimdFor<Number>::value)
          {
            DEAL_II_OPENMP_SIMD_PRAGMA
            for (size_type i = begin; i < end; ++i)
              val[i] += factor * v_val[i];
          }
        else
          {
            for (size_type i = begin; i < end; ++i)
              val[i] += factor * v_val[i];
          }
      }

      Number *const       val;
      const Number *const v_val;
      const Number        factor;
    };

    template <typename Number>
    struct Vectorization_sadd_xav
    {
      Vectorization_sadd_xav(Number *            val,
                             const Number *const v_val,
                             const Number        a,
                             const Number        x)
        : val(val)
        , v_val(v_val)
        , a(a)
        , x(x)
      {}

      void
      operator()(const size_type begin, const size_type end) const
      {
        if (::dealii::parallel::internal::EnableOpenMPSimdFor<Number>::value)
          {
            DEAL_II_OPENMP_SIMD_PRAGMA
            for (size_type i = begin; i < end; ++i)
              val[i] = x * val[i] + a * v_val[i];
          }
        else
          {
            for (size_type i = begin; i < end; ++i)
              val[i] = x * val[i] + a * v_val[i];
          }
      }

      Number *const       val;
      const Number *const v_val;
      const Number        a;
      const Number        x;
    };

    template <typename Number>
    struct Vectorization_subtract_v
    {
      Vectorization_subtract_v(Number *val, const Number *const v_val)
        : val(val)
        , v_val(v_val)
      {}

      void
      operator()(const size_type begin, const size_type end) const
      {
        if (::dealii::parallel::internal::EnableOpenMPSimdFor<Number>::value)
          {
            DEAL_II_OPENMP_SIMD_PRAGMA
            for (size_type i = begin; i < end; ++i)
              val[i] -= v_val[i];
          }
        else
          {
            for (size_type i = begin; i < end; ++i)
              val[i] -= v_val[i];
          }
      }

      Number *const       val;
      const Number *const v_val;
    };

    template <typename Number>
    struct Vectorization_add_factor
    {
      Vectorization_add_factor(Number *const val, const Number factor)
        : val(val)
        , factor(factor)
      {}

      void
      operator()(const size_type begin, const size_type end) const
      {
        if (::dealii::parallel::internal::EnableOpenMPSimdFor<Number>::value)
          {
            DEAL_II_OPENMP_SIMD_PRAGMA
            for (size_type i = begin; i < end; ++i)
              val[i] += factor;
          }
        else
          {
            for (size_type i = begin; i < end; ++i)
              val[i] += factor;
          }
      }

      Number *const val;
      const Number  factor;
    };

    template <typename Number>
    struct Vectorization_add_v
    {
      Vectorization_add_v(Number *const val, const Number *const v_val)
        : val(val)
        , v_val(v_val)
      {}

      void
      operator()(const size_type begin, const size_type end) const
      {
        if (::dealii::parallel::internal::EnableOpenMPSimdFor<Number>::value)
          {
            DEAL_II_OPENMP_SIMD_PRAGMA
            for (size_type i = begin; i < end; ++i)
              val[i] += v_val[i];
          }
        else
          {
            for (size_type i = begin; i < end; ++i)
              val[i] += v_val[i];
          }
      }

      Number *const       val;
      const Number *const v_val;
    };

    template <typename Number>
    struct Vectorization_add_avpbw
    {
      Vectorization_add_avpbw(Number *const       val,
                              const Number *const v_val,
                              const Number *const w_val,
                              const Number        a,
                              const Number        b)
        : val(val)
        , v_val(v_val)
        , w_val(w_val)
        , a(a)
        , b(b)
      {}

      void
      operator()(const size_type begin, const size_type end) const
      {
        if (::dealii::parallel::internal::EnableOpenMPSimdFor<Number>::value)
          {
            DEAL_II_OPENMP_SIMD_PRAGMA
            for (size_type i = begin; i < end; ++i)
              val[i] = val[i] + a * v_val[i] + b * w_val[i];
          }
        else
          {
            for (size_type i = begin; i < end; ++i)
              val[i] = val[i] + a * v_val[i] + b * w_val[i];
          }
      }

      Number *const       val;
      const Number *const v_val;
      const Number *const w_val;
      const Number        a;
      const Number        b;
    };

    template <typename Number>
    struct Vectorization_sadd_xv
    {
      Vectorization_sadd_xv(Number *const       val,
                            const Number *const v_val,
                            const Number        x)
        : val(val)
        , v_val(v_val)
        , x(x)
      {}

      void
      operator()(const size_type begin, const size_type end) const
      {
        if (::dealii::parallel::internal::EnableOpenMPSimdFor<Number>::value)
          {
            DEAL_II_OPENMP_SIMD_PRAGMA
            for (size_type i = begin; i < end; ++i)
              val[i] = x * val[i] + v_val[i];
          }
        else
          {
            for (size_type i = begin; i < end; ++i)
              val[i] = x * val[i] + v_val[i];
          }
      }

      Number *const       val;
      const Number *const v_val;
      const Number        x;
    };

    template <typename Number>
    struct Vectorization_sadd_xavbw
    {
      Vectorization_sadd_xavbw(Number *      val,
                               const Number *v_val,
                               const Number *w_val,
                               Number        x,
                               Number        a,
                               Number        b)
        : val(val)
        , v_val(v_val)
        , w_val(w_val)
        , x(x)
        , a(a)
        , b(b)
      {}

      void
      operator()(const size_type begin, const size_type end) const
      {
        if (::dealii::parallel::internal::EnableOpenMPSimdFor<Number>::value)
          {
            DEAL_II_OPENMP_SIMD_PRAGMA
            for (size_type i = begin; i < end; ++i)
              val[i] = x * val[i] + a * v_val[i] + b * w_val[i];
          }
        else
          {
            for (size_type i = begin; i < end; ++i)
              val[i] = x * val[i] + a * v_val[i] + b * w_val[i];
          }
      }

      Number *const       val;
      const Number *const v_val;
      const Number *const w_val;
      const Number        x;
      const Number        a;
      const Number        b;
    };

    template <typename Number>
    struct Vectorization_scale
    {
      Vectorization_scale(Number *const val, const Number *const v_val)
        : val(val)
        , v_val(v_val)
      {}

      void
      operator()(const size_type begin, const size_type end) const
      {
        if (::dealii::parallel::internal::EnableOpenMPSimdFor<Number>::value)
          {
            DEAL_II_OPENMP_SIMD_PRAGMA
            for (size_type i = begin; i < end; ++i)
              val[i] *= v_val[i];
          }
        else
          {
            for (size_type i = begin; i < end; ++i)
              val[i] *= v_val[i];
          }
      }

      Number *const       val;
      const Number *const v_val;
    };

    template <typename Number>
    struct Vectorization_equ_au
    {
      Vectorization_equ_au(Number *const       val,
                           const Number *const u_val,
                           const Number        a)
        : val(val)
        , u_val(u_val)
        , a(a)
      {}

      void
      operator()(const size_type begin, const size_type end) const
      {
        if (::dealii::parallel::internal::EnableOpenMPSimdFor<Number>::value)
          {
            DEAL_II_OPENMP_SIMD_PRAGMA
            for (size_type i = begin; i < end; ++i)
              val[i] = a * u_val[i];
          }
        else
          {
            for (size_type i = begin; i < end; ++i)
              val[i] = a * u_val[i];
          }
      }

      Number *const       val;
      const Number *const u_val;
      const Number        a;
    };

    template <typename Number>
    struct Vectorization_equ_aubv
    {
      Vectorization_equ_aubv(Number *const       val,
                             const Number *const u_val,
                             const Number *const v_val,
                             const Number        a,
                             const Number        b)
        : val(val)
        , u_val(u_val)
        , v_val(v_val)
        , a(a)
        , b(b)
      {}

      void
      operator()(const size_type begin, const size_type end) const
      {
        if (::dealii::parallel::internal::EnableOpenMPSimdFor<Number>::value)
          {
            DEAL_II_OPENMP_SIMD_PRAGMA
            for (size_type i = begin; i < end; ++i)
              val[i] = a * u_val[i] + b * v_val[i];
          }
        else
          {
            for (size_type i = begin; i < end; ++i)
              val[i] = a * u_val[i] + b * v_val[i];
          }
      }

      Number *const       val;
      const Number *const u_val;
      const Number *const v_val;
      const Number        a;
      const Number        b;
    };

    template <typename Number>
    struct Vectorization_equ_aubvcw
    {
      Vectorization_equ_aubvcw(Number *      val,
                               const Number *u_val,
                               const Number *v_val,
                               const Number *w_val,
                               const Number  a,
                               const Number  b,
                               const Number  c)
        : val(val)
        , u_val(u_val)
        , v_val(v_val)
        , w_val(w_val)
        , a(a)
        , b(b)
        , c(c)
      {}

      void
      operator()(const size_type begin, const size_type end) const
      {
        if (::dealii::parallel::internal::EnableOpenMPSimdFor<Number>::value)
          {
            DEAL_II_OPENMP_SIMD_PRAGMA
            for (size_type i = begin; i < end; ++i)
              val[i] = a * u_val[i] + b * v_val[i] + c * w_val[i];
          }
        else
          {
            for (size_type i = begin; i < end; ++i)
              val[i] = a * u_val[i] + b * v_val[i] + c * w_val[i];
          }
      }

      Number *const       val;
      const Number *const u_val;
      const Number *const v_val;
      const Number *const w_val;
      const Number        a;
      const Number        b;
      const Number        c;
    };

    template <typename Number>
    struct Vectorization_ratio
    {
      Vectorization_ratio(Number *val, const Number *a_val, const Number *b_val)
        : val(val)
        , a_val(a_val)
        , b_val(b_val)
      {}

      void
      operator()(const size_type begin, const size_type end) const
      {
        if (::dealii::parallel::internal::EnableOpenMPSimdFor<Number>::value)
          {
            DEAL_II_OPENMP_SIMD_PRAGMA
            for (size_type i = begin; i < end; ++i)
              val[i] = a_val[i] / b_val[i];
          }
        else
          {
            for (size_type i = begin; i < end; ++i)
              val[i] = a_val[i] / b_val[i];
          }
      }

      Number *const       val;
      const Number *const a_val;
      const Number *const b_val;
    };



    // All sums over all the vector entries (l2-norm, inner product, etc.) are
    // performed with the same code, using a templated operation defined
    // here. There are always two versions defined, a standard one that covers
    // most cases and a vectorized one which is only for equal types and float
    // and double.
    template <typename Number, typename Number2>
    struct Dot
    {
      static constexpr bool vectorizes = std::is_same<Number, Number2>::value &&
                                         (VectorizedArray<Number>::size() > 1);

      Dot(const Number *const X, const Number2 *const Y)
        : X(X)
        , Y(Y)
      {}

      Number
      operator()(const size_type i) const
      {
        return X[i] * Number(numbers::NumberTraits<Number2>::conjugate(Y[i]));
      }

      VectorizedArray<Number>
      do_vectorized(const size_type i) const
      {
        VectorizedArray<Number> x, y;
        x.load(X + i);
        y.load(Y + i);

        // the following operation in VectorizedArray does an element-wise
        // scalar product without taking into account complex values and
        // the need to take the complex-conjugate of one argument. this
        // may be a bug, but because all VectorizedArray classes only
        // work on real scalars, it doesn't really matter very much.
        // in any case, assert that we really don't get here for
        // complex-valued objects
        static_assert(numbers::NumberTraits<Number>::is_complex == false,
                      "This operation is not correctly implemented for "
                      "complex-valued objects.");
        return x * y;
      }

      const Number *const  X;
      const Number2 *const Y;
    };

    template <typename Number, typename RealType>
    struct Norm2
    {
      static const bool vectorizes = VectorizedArray<Number>::size() > 1;

      Norm2(const Number *const X)
        : X(X)
      {}

      RealType
      operator()(const size_type i) const
      {
        return numbers::NumberTraits<Number>::abs_square(X[i]);
      }

      VectorizedArray<Number>
      do_vectorized(const size_type i) const
      {
        VectorizedArray<Number> x;
        x.load(X + i);
        return x * x;
      }

      const Number *const X;
    };

    template <typename Number, typename RealType>
    struct Norm1
    {
      static const bool vectorizes = VectorizedArray<Number>::size() > 1;

      Norm1(const Number *X)
        : X(X)
      {}

      RealType
      operator()(const size_type i) const
      {
        return numbers::NumberTraits<Number>::abs(X[i]);
      }

      VectorizedArray<Number>
      do_vectorized(const size_type i) const
      {
        VectorizedArray<Number> x;
        x.load(X + i);
        return std::abs(x);
      }

      const Number *X;
    };

    template <typename Number, typename RealType>
    struct NormP
    {
      static const bool vectorizes = VectorizedArray<Number>::size() > 1;

      NormP(const Number *X, RealType p)
        : X(X)
        , p(p)
      {}

      RealType
      operator()(const size_type i) const
      {
        return std::pow(numbers::NumberTraits<Number>::abs(X[i]), p);
      }

      VectorizedArray<Number>
      do_vectorized(const size_type i) const
      {
        VectorizedArray<Number> x;
        x.load(X + i);
        return std::pow(std::abs(x), p);
      }

      const Number * X;
      const RealType p;
    };

    template <typename Number>
    struct MeanValue
    {
      static const bool vectorizes = VectorizedArray<Number>::size() > 1;

      MeanValue(const Number *X)
        : X(X)
      {}

      Number
      operator()(const size_type i) const
      {
        return X[i];
      }

      VectorizedArray<Number>
      do_vectorized(const size_type i) const
      {
        VectorizedArray<Number> x;
        x.load(X + i);
        return x;
      }

      const Number *X;
    };

    template <typename Number>
    struct AddAndDot
    {
      static const bool vectorizes = VectorizedArray<Number>::size() > 1;

      AddAndDot(Number *const       X,
                const Number *const V,
                const Number *const W,
                const Number        a)
        : X(X)
        , V(V)
        , W(W)
        , a(a)
      {}

      Number
      operator()(const size_type i) const
      {
        X[i] += a * V[i];
        return X[i] * Number(numbers::NumberTraits<Number>::conjugate(W[i]));
      }

      VectorizedArray<Number>
      do_vectorized(const size_type i) const
      {
        VectorizedArray<Number> x, w, v;
        x.load(X + i);
        v.load(V + i);
        x += a * v;
        x.store(X + i);
        // may only load from W after storing in X because the pointers might
        // point to the same memory
        w.load(W + i);

        // the following operation in VectorizedArray does an element-wise
        // scalar product without taking into account complex values and
        // the need to take the complex-conjugate of one argument. this
        // may be a bug, but because all VectorizedArray classes only
        // work on real scalars, it doesn't really matter very much.
        // in any case, assert that we really don't get here for
        // complex-valued objects
        static_assert(numbers::NumberTraits<Number>::is_complex == false,
                      "This operation is not correctly implemented for "
                      "complex-valued objects.");
        return x * w;
      }

      Number *const       X;
      const Number *const V;
      const Number *const W;
      const Number        a;
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

    // Loops are unrolled as follows: the range [first,last) is broken into
    // @p n_chunks each of size 32 plus the @p remainder.
    // accumulate_regular() does the work on 32*n_chunks elements employing SIMD
    // if possible and stores the result of the operation for each chunk in @p outer_results.

    // The code returns the result as the last argument in order to make
    // spawning tasks simpler and use automatic template deduction.


    /**
     * The minimum number of chunks (each of size 32) to divide the range
     * [first,last) into two (second part of the if branch in
     * accumulate_recursive).
     */
    const unsigned int vector_accumulation_recursion_threshold = 128;

    template <typename Operation, typename ResultType>
    void
    accumulate_recursive(const Operation &op,
                         const size_type  first,
                         const size_type  last,
                         ResultType &     result)
    {
      const size_type vec_size = last - first;
      if (vec_size <= vector_accumulation_recursion_threshold * 32)
        {
          // the vector is short enough so we perform the summation. first
          // work on the regular part. The innermost 32 values are expanded in
          // order to obtain known loop bounds for most of the work.
          size_type  index = first;
          ResultType outer_results[vector_accumulation_recursion_threshold];

          // set the zeroth element to zero to correctly handle the case where
          // vec_size == 0
          outer_results[0] = ResultType();

          // the variable serves two purposes: (i)  number of chunks (each 32
          // indices) for the given size; all results are stored in
          // outer_results[0,n_chunks) (ii) in the SIMD case n_chunks is also a
          // next free index in outer_results[] to which we can write after
          // accumulate_regular() is executed.
          size_type       n_chunks  = vec_size / 32;
          const size_type remainder = vec_size % 32;
          Assert(remainder == 0 ||
                   n_chunks < vector_accumulation_recursion_threshold,
                 ExcInternalError());

          // Select between the regular version and vectorized version based
          // on the number types we are given. To choose the vectorized
          // version often enough, we need to have all tasks but the last one
          // to be divisible by the vectorization length
          accumulate_regular(
            op,
            n_chunks,
            index,
            outer_results,
            std::integral_constant<bool, Operation::vectorizes>());

          // now work on the remainder, i.e., the last up to 32 values. Use
          // switch statement with fall-through to work on these values.
          if (remainder > 0)
            {
              // if we got here, it means that (vec_size <=
              // vector_accumulation_recursion_threshold * 32), which is to say
              // that the domain can be split into n_chunks <=
              // vector_accumulation_recursion_threshold:
              AssertIndexRange(n_chunks,
                               vector_accumulation_recursion_threshold + 1);
              // split the remainder into chunks of 8, there could be up to 3
              // such chunks since remainder < 32.
              // Work on those chunks without any SIMD, that is we call
              // op(index).
              const size_type inner_chunks = remainder / 8;
              Assert(inner_chunks <= 3, ExcInternalError());
              const size_type remainder_inner = remainder % 8;
              ResultType      r0 = ResultType(), r1 = ResultType(),
                         r2 = ResultType();
              switch (inner_chunks)
                {
                  case 3:
                    r2 = op(index++);
                    for (size_type j = 1; j < 8; ++j)
                      r2 += op(index++);
                    DEAL_II_FALLTHROUGH;
                  case 2:
                    r1 = op(index++);
                    for (size_type j = 1; j < 8; ++j)
                      r1 += op(index++);
                    r1 += r2;
                    DEAL_II_FALLTHROUGH;
                  case 1:
                    r2 = op(index++);
                    for (size_type j = 1; j < 8; ++j)
                      r2 += op(index++);
                    DEAL_II_FALLTHROUGH;
                  default:
                    for (size_type j = 0; j < remainder_inner; ++j)
                      r0 += op(index++);
                    r0 += r2;
                    r0 += r1;
                    if (n_chunks == vector_accumulation_recursion_threshold)
                      outer_results[vector_accumulation_recursion_threshold -
                                    1] += r0;
                    else
                      {
                        outer_results[n_chunks] = r0;
                        n_chunks++;
                      }
                    break;
                }
            }
          // make sure we worked through all indices
          AssertDimension(index, last);

          // now sum the results from the chunks stored in
          // outer_results[0,n_chunks) recursively
          while (n_chunks > 1)
            {
              if (n_chunks % 2 == 1)
                outer_results[n_chunks++] = ResultType();
              for (size_type i = 0; i < n_chunks; i += 2)
                outer_results[i / 2] = outer_results[i] + outer_results[i + 1];
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
          Assert(first + 3 * new_size < last, ExcInternalError());
          ResultType r0, r1, r2, r3;
          accumulate_recursive(op, first, first + new_size, r0);
          accumulate_recursive(op, first + new_size, first + 2 * new_size, r1);
          accumulate_recursive(op,
                               first + 2 * new_size,
                               first + 3 * new_size,
                               r2);
          accumulate_recursive(op, first + 3 * new_size, last, r3);
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
    accumulate_regular(
      const Operation &op,
      const size_type &n_chunks,
      size_type &      index,
      ResultType (&outer_results)[vector_accumulation_recursion_threshold],
      std::integral_constant<bool, false>)
    {
      // note that each chunk is chosen to have a width of 32, thereby the index
      // is incremented by 4*8 for each @p i.
      for (size_type i = 0; i < n_chunks; ++i)
        {
          ResultType r0 = op(index);
          ResultType r1 = op(index + 1);
          ResultType r2 = op(index + 2);
          ResultType r3 = op(index + 3);
          index += 4;
          for (size_type j = 1; j < 8; ++j, index += 4)
            {
              r0 += op(index);
              r1 += op(index + 1);
              r2 += op(index + 2);
              r3 += op(index + 3);
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
    accumulate_regular(
      const Operation &op,
      size_type &      n_chunks,
      size_type &      index,
      Number (&outer_results)[vector_accumulation_recursion_threshold],
      std::integral_constant<bool, true>)
    {
      // we start from @p index and workout @p n_chunks each of size 32.
      // in order employ SIMD and work on @p nvecs at a time, we split this
      // loop yet again:
      // First we work on (n_chunks/nvecs) chunks, where each chunk processes
      // nvecs*(4*8) elements.

      constexpr unsigned int nvecs          = VectorizedArray<Number>::size();
      const size_type        regular_chunks = n_chunks / nvecs;
      for (size_type i = 0; i < regular_chunks; ++i)
        {
          VectorizedArray<Number> r0 = op.do_vectorized(index);
          VectorizedArray<Number> r1 = op.do_vectorized(index + nvecs);
          VectorizedArray<Number> r2 = op.do_vectorized(index + 2 * nvecs);
          VectorizedArray<Number> r3 = op.do_vectorized(index + 3 * nvecs);
          index += nvecs * 4;
          for (size_type j = 1; j < 8; ++j, index += nvecs * 4)
            {
              r0 += op.do_vectorized(index);
              r1 += op.do_vectorized(index + nvecs);
              r2 += op.do_vectorized(index + 2 * nvecs);
              r3 += op.do_vectorized(index + 3 * nvecs);
            }
          r0 += r1;
          r2 += r3;
          r0 += r2;
          r0.store(&outer_results[i * nvecs]);
        }

      // If we are treating a case where the vector length is not divisible by
      // the vectorization length, need a cleanup loop
      // The remaining chunks are processed one by one starting from
      // regular_chunks * nvecs; We do as much as possible with 2 SIMD
      // operations within each chunk. Here we assume that nvecs < 32/2 = 16 as
      // well as 16%nvecs==0.
      static_assert(
        VectorizedArray<Number>::size() <= 16 &&
          16 % VectorizedArray<Number>::size() == 0,
        "VectorizedArray::size() must be a power of 2 and not more than 16");
      Assert(16 % nvecs == 0, ExcInternalError());
      if (n_chunks % nvecs != 0)
        {
          VectorizedArray<Number> r0  = VectorizedArray<Number>(),
                                  r1  = VectorizedArray<Number>();
          const size_type start_irreg = regular_chunks * nvecs;
          for (size_type c = start_irreg; c < n_chunks; ++c)
            for (size_type j = 0; j < 32; j += 2 * nvecs, index += 2 * nvecs)
              {
                r0 += op.do_vectorized(index);
                r1 += op.do_vectorized(index + nvecs);
              }
          r0 += r1;
          r0.store(&outer_results[start_irreg]);
          // update n_chunks to denote unused element in outer_results[] from
          // which we can keep writing.
          n_chunks = start_irreg + VectorizedArray<Number>::size();
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

      TBBReduceFunctor(const Operation &op,
                       const size_type  start,
                       const size_type  end)
        : op(op)
        , start(start)
        , end(end)
      {
        const size_type vec_size = end - start;
        // set chunk size for sub-tasks
        const unsigned int gs =
          internal::VectorImplementation::minimum_parallel_grain_size;
        n_chunks =
          std::min(static_cast<size_type>(4 * MultithreadInfo::n_threads()),
                   vec_size / gs);
        chunk_size = vec_size / n_chunks;

        // round to next multiple of 512 (or leave it at the minimum grain size
        // if that happens to be smaller). this is advantageous because our
        // algorithm favors lengths of a power of 2 due to pairwise summation ->
        // at most one 'oddly' sized chunk
        if (chunk_size > 512)
          chunk_size = ((chunk_size + 511) / 512) * 512;
        n_chunks = (vec_size + chunk_size - 1) / chunk_size;
        AssertIndexRange((n_chunks - 1) * chunk_size, vec_size);
        AssertIndexRange(vec_size, n_chunks * chunk_size + 1);

        if (n_chunks > threshold_array_allocate)
          {
            // make sure we allocate an even number of elements,
            // access to the new last element is needed in do_sum()
            large_array.resize(2 * ((n_chunks + 1) / 2));
            array_ptr = large_array.data();
          }
        else
          array_ptr = &small_array[0];
      }

      /**
       * An operator used by TBB to work on a given @p range of chunks
       * [range.begin(), range.end()).
       */
      void
      operator()(const tbb::blocked_range<size_type> &range) const
      {
        for (size_type i = range.begin(); i < range.end(); ++i)
          accumulate_recursive(op,
                               start + i * chunk_size,
                               std::min(start + (i + 1) * chunk_size, end),
                               array_ptr[i]);
      }

      ResultType
      do_sum() const
      {
        while (n_chunks > 1)
          {
            if (n_chunks % 2 == 1)
              array_ptr[n_chunks++] = ResultType();
            for (size_type i = 0; i < n_chunks; i += 2)
              array_ptr[i / 2] = array_ptr[i] + array_ptr[i + 1];
            n_chunks /= 2;
          }
        return array_ptr[0];
      }

      const Operation &op;
      const size_type  start;
      const size_type  end;

      mutable unsigned int    n_chunks;
      unsigned int            chunk_size;
      ResultType              small_array[threshold_array_allocate];
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
    void
    parallel_reduce(
      const Operation &op,
      const size_type  start,
      const size_type  end,
      ResultType &     result,
      const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
        &partitioner)
    {
#ifdef DEAL_II_WITH_THREADS
      const size_type vec_size = end - start;
      // only go to the parallel function in case there are at least 4 parallel
      // items, otherwise the overhead is too large
      if (vec_size >=
            4 * internal::VectorImplementation::minimum_parallel_grain_size &&
          MultithreadInfo::n_threads() > 1)
        {
          Assert(partitioner.get() != nullptr,
                 ExcInternalError(
                   "Unexpected initialization of Vector that does "
                   "not set the TBB partitioner to a usable state."));
          std::shared_ptr<tbb::affinity_partitioner> tbb_partitioner =
            partitioner->acquire_one_partitioner();

          TBBReduceFunctor<Operation, ResultType> generic_functor(op,
                                                                  start,
                                                                  end);
          // We use a minimum grain size of 1 here since the grains at this
          // stage of dividing the work refer to the number of vector chunks
          // that are processed by (possibly different) threads in the
          // parallelized for loop (i.e., they do not refer to individual
          // vector entries). The number of chunks here is calculated inside
          // TBBForFunctor. See also GitHub issue #2496 for further discussion
          // of this strategy.
          ::dealii::parallel::internal::parallel_for(
            static_cast<size_type>(0),
            static_cast<size_type>(generic_functor.n_chunks),
            generic_functor,
            1,
            tbb_partitioner);
          partitioner->release_one_partitioner(tbb_partitioner);
          result = generic_functor.do_sum();
        }
      else
        accumulate_recursive(op, start, end, result);
#else
      accumulate_recursive(op, start, end, result);
      (void)partitioner;
#endif
    }


    template <typename Number, typename Number2, typename MemorySpace>
    struct functions
    {
      static void
      copy(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &
        /*thread_loop_partitioner*/,
        const size_type /*size*/,
        const ::dealii::MemorySpace::MemorySpaceData<Number2, MemorySpace>
          & /*v_data*/,
        ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpace> & /*data*/)
      {
        static_assert(
          std::is_same<MemorySpace, ::dealii::MemorySpace::CUDA>::value &&
            std::is_same<Number, Number2>::value,
          "For the CUDA MemorySpace Number and Number2 should be the same type");
      }

      static void
      set(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &
        /*thread_loop_partitioner*/,
        const size_type /*size*/,
        const Number /*s*/,
        ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpace> & /*data*/)
      {}

      static void
      add_vector(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &
        /*thread_loop_partitioner*/,
        const size_type /*size*/,
        const ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpace>
          & /*v_data*/,
        ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpace> & /*data*/)
      {}

      static void
      subtract_vector(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &
        /*thread_loop_partitioner*/,
        const size_type /*size*/,
        const ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpace>
          & /*v_data*/,
        ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpace> & /*data*/)
      {}

      static void
      add_factor(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &
        /*thread_loop_partitioner*/,
        const size_type /*size*/,
        Number /*a*/,
        ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpace> & /*data*/)
      {}

      static void
      add_av(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &
        /*thread_loop_partitioner*/,
        const size_type /*size*/,
        const Number /*a*/,
        const ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpace>
          & /*v_data*/,
        ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpace> & /*data*/)
      {}

      static void
      add_avpbw(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &
        /*thread_loop_partitioner*/,
        const size_type /*size*/,
        const Number /*a*/,
        const Number /*b*/,
        const ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpace>
          & /*v_data*/,
        const ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpace>
          & /*w_data*/,
        ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpace> & /*data*/)
      {}

      static void
      sadd_xv(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &
        /*thread_loop_partitioner*/,
        const size_type /*size*/,
        const Number /*x*/,
        const ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpace>
          & /*v_data*/,
        ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpace> & /*data*/)
      {}

      static void
      sadd_xav(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &
        /*thread_loop_partitioner*/,
        const size_type /*size*/,
        const Number /*x*/,
        const Number /*a*/,
        const ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpace>
          & /*v_data*/,
        ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpace> & /*data*/)
      {}

      static void
      sadd_xavbw(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &
        /*thread_loop_partitioner*/,
        const size_type /*size*/,
        const Number /*x*/,
        const Number /*a*/,
        const Number /*b*/,
        const ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpace>
          & /*v_data*/,
        const ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpace>
          & /*w_data*/,
        ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpace> & /*data*/)
      {}

      static void
      multiply_factor(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &
        /*thread_loop_partitioner*/,
        const size_type /*size*/,
        const Number /*factor*/,
        ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpace> & /*data*/)
      {}

      static void
      scale(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &
        /*thread_loop_partitioner*/,
        const size_type /*size*/,
        const ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpace>
          & /*v_data*/,
        ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpace> & /*data*/)
      {}

      static void
      equ_au(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &
        /*thread_loop_partitioner*/,
        const size_type /*size*/,
        const Number /*a*/,
        const ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpace>
          & /*v_data*/,
        ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpace> & /*data*/)
      {}

      static void
      equ_aubv(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &
        /*thread_loop_partitioner*/,
        const size_type /*size*/,
        const Number /*a*/,
        const Number /*b*/,
        const ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpace>
          & /*v_data*/,
        const ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpace>
          & /*w_data*/,
        ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpace> & /*data*/)
      {}

      static Number
      dot(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &
        /*thread_loop_partitioner*/,
        const size_type /*size*/,
        const ::dealii::MemorySpace::MemorySpaceData<Number2, MemorySpace>
          & /*v_data*/,
        ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpace> & /*data*/)
      {
        return Number();
      }

      template <typename real_type>
      static void
      norm_2(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &
        /*thread_loop_partitioner*/,
        const size_type /*size*/,
        real_type & /*sum*/,
        const ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpace>
          & /*v_data*/,
        ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpace> & /*data*/)
      {}

      static Number
      mean_value(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &
        /*thread_loop_partitioner*/,
        const size_type /*size*/,
        const ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpace>
          & /*data*/)
      {
        return Number();
      }

      template <typename real_type>
      static void
      norm_1(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &
        /*thread_loop_partitioner*/,
        const size_type /*size*/,
        real_type & /*sum*/,
        Number * /*values*/,
        Number * /*values_dev*/)
      {}

      template <typename real_type>
      static void
      norm_p(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &
        /*thread_loop_partitioner*/,
        const size_type /*size*/,
        real_type & /*sum*/,
        real_type /*p*/,
        ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpace> & /*data*/)
      {}

      static Number
      add_and_dot(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &
        /*thread_loop_partitioner*/,
        const size_type /*size*/,
        const Number /*a*/,
        const ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpace>
          & /*v_data*/,
        const ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpace>
          & /*w_data*/,
        ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpace> & /*data*/)
      {
        return Number();
      }

      template <typename MemorySpace2>
      static void
      import(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &
        /*thread_loop_partitioner*/,
        const size_type /*size*/,
        VectorOperation::values /*operation*/,
        const ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpace2>
          & /*v_data*/,
        ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpace> & /*data*/)
      {}
    };



    template <typename Number, typename Number2>
    struct functions<Number, Number2, ::dealii::MemorySpace::Host>
    {
      static void
      copy(const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
             &             thread_loop_partitioner,
           const size_type size,
           const ::dealii::MemorySpace::
             MemorySpaceData<Number2, ::dealii::MemorySpace::Host> &v_data,
           ::dealii::MemorySpace::MemorySpaceData<Number,
                                                  ::dealii::MemorySpace::Host>
             &data)
      {
        Vector_copy<Number, Number2> copier(v_data.values.get(),
                                            data.values.get());
        parallel_for(copier, 0, size, thread_loop_partitioner);
      }

      static void
      set(const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
            &             thread_loop_partitioner,
          const size_type size,
          const Number    s,
          ::dealii::MemorySpace::MemorySpaceData<Number,
                                                 ::dealii::MemorySpace::Host>
            &data)
      {
        Vector_set<Number> setter(s, data.values.get());
        parallel_for(setter, 0, size, thread_loop_partitioner);
      }

      static void
      add_vector(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
          &             thread_loop_partitioner,
        const size_type size,
        const ::dealii::MemorySpace::
          MemorySpaceData<Number, ::dealii::MemorySpace::Host> &v_data,
        ::dealii::MemorySpace::MemorySpaceData<Number,
                                               ::dealii::MemorySpace::Host>
          &data)
      {
        Vectorization_add_v<Number> vector_add(data.values.get(),
                                               v_data.values.get());
        parallel_for(vector_add, 0, size, thread_loop_partitioner);
      }

      static void
      subtract_vector(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
          &             thread_loop_partitioner,
        const size_type size,
        const ::dealii::MemorySpace::
          MemorySpaceData<Number, ::dealii::MemorySpace::Host> &v_data,
        ::dealii::MemorySpace::MemorySpaceData<Number,
                                               ::dealii::MemorySpace::Host>
          &data)
      {
        Vectorization_subtract_v<Number> vector_subtract(data.values.get(),
                                                         v_data.values.get());
        parallel_for(vector_subtract, 0, size, thread_loop_partitioner);
      }

      static void
      add_factor(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
          &             thread_loop_partitioner,
        const size_type size,
        Number          a,
        ::dealii::MemorySpace::MemorySpaceData<Number,
                                               ::dealii::MemorySpace::Host>
          &data)
      {
        Vectorization_add_factor<Number> vector_add(data.values.get(), a);
        parallel_for(vector_add, 0, size, thread_loop_partitioner);
      }

      static void
      add_av(const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
               &             thread_loop_partitioner,
             const size_type size,
             const Number    a,
             const ::dealii::MemorySpace::
               MemorySpaceData<Number, ::dealii::MemorySpace::Host> &v_data,
             ::dealii::MemorySpace::MemorySpaceData<Number,
                                                    ::dealii::MemorySpace::Host>
               &data)
      {
        Vectorization_add_av<Number> vector_add(data.values.get(),
                                                v_data.values.get(),
                                                a);
        parallel_for(vector_add, 0, size, thread_loop_partitioner);
      }

      static void
      add_avpbw(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
          &             thread_loop_partitioner,
        const size_type size,
        const Number    a,
        const Number    b,
        const ::dealii::MemorySpace::
          MemorySpaceData<Number, ::dealii::MemorySpace::Host> &v_data,
        const ::dealii::MemorySpace::
          MemorySpaceData<Number, ::dealii::MemorySpace::Host> &w_data,
        ::dealii::MemorySpace::MemorySpaceData<Number,
                                               ::dealii::MemorySpace::Host>
          &data)
      {
        Vectorization_add_avpbw<Number> vector_add(
          data.values.get(), v_data.values.get(), w_data.values.get(), a, b);
        parallel_for(vector_add, 0, size, thread_loop_partitioner);
      }

      static void
      sadd_xv(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
          &             thread_loop_partitioner,
        const size_type size,
        const Number    x,
        const ::dealii::MemorySpace::
          MemorySpaceData<Number, ::dealii::MemorySpace::Host> &v_data,
        ::dealii::MemorySpace::MemorySpaceData<Number,
                                               ::dealii::MemorySpace::Host>
          &data)
      {
        Vectorization_sadd_xv<Number> vector_sadd(data.values.get(),
                                                  v_data.values.get(),
                                                  x);
        parallel_for(vector_sadd, 0, size, thread_loop_partitioner);
      }

      static void
      sadd_xav(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
          &             thread_loop_partitioner,
        const size_type size,
        const Number    x,
        const Number    a,
        const ::dealii::MemorySpace::
          MemorySpaceData<Number, ::dealii::MemorySpace::Host> &v_data,
        ::dealii::MemorySpace::MemorySpaceData<Number,
                                               ::dealii::MemorySpace::Host>
          &data)
      {
        Vectorization_sadd_xav<Number> vector_sadd(data.values.get(),
                                                   v_data.values.get(),
                                                   a,
                                                   x);
        parallel_for(vector_sadd, 0, size, thread_loop_partitioner);
      }

      static void
      sadd_xavbw(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
          &             thread_loop_partitioner,
        const size_type size,
        const Number    x,
        const Number    a,
        const Number    b,
        const ::dealii::MemorySpace::
          MemorySpaceData<Number, ::dealii::MemorySpace::Host> &v_data,
        const ::dealii::MemorySpace::
          MemorySpaceData<Number, ::dealii::MemorySpace::Host> &w_data,
        ::dealii::MemorySpace::MemorySpaceData<Number,
                                               ::dealii::MemorySpace::Host>
          &data)
      {
        Vectorization_sadd_xavbw<Number> vector_sadd(
          data.values.get(), v_data.values.get(), w_data.values.get(), x, a, b);
        parallel_for(vector_sadd, 0, size, thread_loop_partitioner);
      }

      static void
      multiply_factor(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
          &             thread_loop_partitioner,
        const size_type size,
        const Number    factor,
        ::dealii::MemorySpace::MemorySpaceData<Number,
                                               ::dealii::MemorySpace::Host>
          &data)
      {
        Vectorization_multiply_factor<Number> vector_multiply(data.values.get(),
                                                              factor);
        parallel_for(vector_multiply, 0, size, thread_loop_partitioner);
      }

      static void
      scale(const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
              &             thread_loop_partitioner,
            const size_type size,
            const ::dealii::MemorySpace::
              MemorySpaceData<Number, ::dealii::MemorySpace::Host> &v_data,
            ::dealii::MemorySpace::MemorySpaceData<Number,
                                                   ::dealii::MemorySpace::Host>
              &data)
      {
        Vectorization_scale<Number> vector_scale(data.values.get(),
                                                 v_data.values.get());
        parallel_for(vector_scale, 0, size, thread_loop_partitioner);
      }

      static void
      equ_au(const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
               &             thread_loop_partitioner,
             const size_type size,
             const Number    a,
             const ::dealii::MemorySpace::
               MemorySpaceData<Number, ::dealii::MemorySpace::Host> &v_data,
             ::dealii::MemorySpace::MemorySpaceData<Number,
                                                    ::dealii::MemorySpace::Host>
               &data)
      {
        Vectorization_equ_au<Number> vector_equ(data.values.get(),
                                                v_data.values.get(),
                                                a);
        parallel_for(vector_equ, 0, size, thread_loop_partitioner);
      }

      static void
      equ_aubv(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
          &             thread_loop_partitioner,
        const size_type size,
        const Number    a,
        const Number    b,
        const ::dealii::MemorySpace::
          MemorySpaceData<Number, ::dealii::MemorySpace::Host> &v_data,
        const ::dealii::MemorySpace::
          MemorySpaceData<Number, ::dealii::MemorySpace::Host> &w_data,
        ::dealii::MemorySpace::MemorySpaceData<Number,
                                               ::dealii::MemorySpace::Host>
          &data)
      {
        Vectorization_equ_aubv<Number> vector_equ(
          data.values.get(), v_data.values.get(), w_data.values.get(), a, b);
        parallel_for(vector_equ, 0, size, thread_loop_partitioner);
      }

      static Number
      dot(const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
            &             thread_loop_partitioner,
          const size_type size,
          const ::dealii::MemorySpace::
            MemorySpaceData<Number2, ::dealii::MemorySpace::Host> &v_data,
          ::dealii::MemorySpace::MemorySpaceData<Number,
                                                 ::dealii::MemorySpace::Host>
            &data)
      {
        Number                                                   sum;
        dealii::internal::VectorOperations::Dot<Number, Number2> dot(
          data.values.get(), v_data.values.get());
        dealii::internal::VectorOperations::parallel_reduce(
          dot, 0, size, sum, thread_loop_partitioner);
        AssertIsFinite(sum);

        return sum;
      }

      template <typename real_type>
      static void
      norm_2(const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
               &             thread_loop_partitioner,
             const size_type size,
             real_type &     sum,
             ::dealii::MemorySpace::MemorySpaceData<Number,
                                                    ::dealii::MemorySpace::Host>
               &data)
      {
        Norm2<Number, real_type> norm2(data.values.get());
        parallel_reduce(norm2, 0, size, sum, thread_loop_partitioner);
      }

      static Number
      mean_value(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
          &             thread_loop_partitioner,
        const size_type size,
        const ::dealii::MemorySpace::
          MemorySpaceData<Number, ::dealii::MemorySpace::Host> &data)
      {
        Number            sum;
        MeanValue<Number> mean(data.values.get());
        parallel_reduce(mean, 0, size, sum, thread_loop_partitioner);

        return sum;
      }

      template <typename real_type>
      static void
      norm_1(const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
               &             thread_loop_partitioner,
             const size_type size,
             real_type &     sum,
             ::dealii::MemorySpace::MemorySpaceData<Number,
                                                    ::dealii::MemorySpace::Host>
               &data)
      {
        Norm1<Number, real_type> norm1(data.values.get());
        parallel_reduce(norm1, 0, size, sum, thread_loop_partitioner);
      }

      template <typename real_type>
      static void
      norm_p(const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
               &             thread_loop_partitioner,
             const size_type size,
             real_type &     sum,
             const real_type p,
             ::dealii::MemorySpace::MemorySpaceData<Number,
                                                    ::dealii::MemorySpace::Host>
               &data)
      {
        NormP<Number, real_type> normp(data.values.get(), p);
        parallel_reduce(normp, 0, size, sum, thread_loop_partitioner);
      }

      static Number
      add_and_dot(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
          &             thread_loop_partitioner,
        const size_type size,
        const Number    a,
        const ::dealii::MemorySpace::
          MemorySpaceData<Number, ::dealii::MemorySpace::Host> &v_data,
        const ::dealii::MemorySpace::
          MemorySpaceData<Number, ::dealii::MemorySpace::Host> &w_data,
        ::dealii::MemorySpace::MemorySpaceData<Number,
                                               ::dealii::MemorySpace::Host>
          &data)
      {
        Number            sum;
        AddAndDot<Number> adder(data.values.get(),
                                v_data.values.get(),
                                w_data.values.get(),
                                a);
        parallel_reduce(adder, 0, size, sum, thread_loop_partitioner);

        return sum;
      }

      template <typename MemorySpace2>
      static void
      import(const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
               &                     thread_loop_partitioner,
             const size_type         size,
             VectorOperation::values operation,
             const ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpace2>
               &v_data,
             ::dealii::MemorySpace::MemorySpaceData<Number,
                                                    ::dealii::MemorySpace::Host>
               &data,
             typename std::enable_if<
               std::is_same<MemorySpace2, dealii::MemorySpace::Host>::value,
               int>::type = 0)
      {
        if (operation == VectorOperation::insert)
          {
            copy(thread_loop_partitioner, size, v_data, data);
          }
        else if (operation == VectorOperation::add)
          {
            add_vector(thread_loop_partitioner, size, v_data, data);
          }
        else
          {
            AssertThrow(false, ExcNotImplemented());
          }
      }

#ifdef DEAL_II_COMPILER_CUDA_AWARE
      template <typename MemorySpace2>
      static void
      import(const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
               & /*thread_loop_partitioner*/,
             const size_type         size,
             VectorOperation::values operation,
             const ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpace2>
               &v_data,
             ::dealii::MemorySpace::MemorySpaceData<Number,
                                                    ::dealii::MemorySpace::Host>
               &data,
             typename std::enable_if<
               std::is_same<MemorySpace2, ::dealii::MemorySpace::CUDA>::value,
               int>::type = 0)
      {
        if (operation == VectorOperation::insert)
          {
            cudaError_t cuda_error_code = cudaMemcpy(data.values.get(),
                                                     v_data.values_dev.get(),
                                                     size * sizeof(Number),
                                                     cudaMemcpyDeviceToHost);
            AssertCuda(cuda_error_code);
          }
        else
          {
            AssertThrow(false, ExcNotImplemented());
          }
      }
#endif
    };



#ifdef DEAL_II_COMPILER_CUDA_AWARE
    template <typename Number>
    struct functions<Number, Number, ::dealii::MemorySpace::CUDA>
    {
      static const int block_size =
        ::dealii::LinearAlgebra::CUDAWrappers::kernel::block_size;
      static const int chunk_size =
        ::dealii::LinearAlgebra::CUDAWrappers::kernel::chunk_size;

      static void
      copy(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &,
        const size_type size,
        const ::dealii::MemorySpace::
          MemorySpaceData<Number, ::dealii::MemorySpace::CUDA> &v_data,
        ::dealii::MemorySpace::MemorySpaceData<Number,
                                               ::dealii::MemorySpace::CUDA>
          &data)
      {
        cudaError_t cuda_error_code = cudaMemcpy(data.values_dev.get(),
                                                 v_data.values_dev.get(),
                                                 size * sizeof(Number),
                                                 cudaMemcpyDeviceToDevice);
        AssertCuda(cuda_error_code);
      }

      static void
      set(const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &,
          const size_type size,
          const Number    s,
          ::dealii::MemorySpace::MemorySpaceData<Number,
                                                 ::dealii::MemorySpace::CUDA>
            &data)
      {
        const int n_blocks = 1 + size / (chunk_size * block_size);
        ::dealii::LinearAlgebra::CUDAWrappers::kernel::set<Number>
          <<<n_blocks, block_size>>>(data.values_dev.get(), s, size);
        AssertCudaKernel();
      }

      static void
      add_vector(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &,
        const size_type size,
        const ::dealii::MemorySpace::
          MemorySpaceData<Number, ::dealii::MemorySpace::CUDA> &v_data,
        ::dealii::MemorySpace::MemorySpaceData<Number,
                                               ::dealii::MemorySpace::CUDA>
          &data)
      {
        const int n_blocks = 1 + size / (chunk_size * block_size);
        ::dealii::LinearAlgebra::CUDAWrappers::kernel::add_aV<Number>
          <<<n_blocks, block_size>>>(data.values_dev.get(),
                                     1.,
                                     v_data.values_dev.get(),
                                     size);
        AssertCudaKernel();
      }

      static void
      subtract_vector(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &,
        const size_type size,
        const ::dealii::MemorySpace::
          MemorySpaceData<Number, ::dealii::MemorySpace::CUDA> &v_data,
        ::dealii::MemorySpace::MemorySpaceData<Number,
                                               ::dealii::MemorySpace::CUDA>
          &data)
      {
        const int n_blocks = 1 + size / (chunk_size * block_size);
        ::dealii::LinearAlgebra::CUDAWrappers::kernel::add_aV<Number>
          <<<n_blocks, block_size>>>(data.values_dev.get(),
                                     -1.,
                                     v_data.values_dev.get(),
                                     size);
        AssertCudaKernel();
      }

      static void
      add_factor(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &,
        const size_type size,
        Number          a,
        ::dealii::MemorySpace::MemorySpaceData<Number,
                                               ::dealii::MemorySpace::CUDA>
          &data)
      {
        const int n_blocks = 1 + size / (chunk_size * block_size);
        ::dealii::LinearAlgebra::CUDAWrappers::kernel::vec_add<Number>
          <<<n_blocks, block_size>>>(data.values_dev.get(), a, size);
        AssertCudaKernel();
      }

      static void
      add_av(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &,
        const size_type size,
        const Number    a,
        const ::dealii::MemorySpace::
          MemorySpaceData<Number, ::dealii::MemorySpace::CUDA> &v_data,
        ::dealii::MemorySpace::MemorySpaceData<Number,
                                               ::dealii::MemorySpace::CUDA>
          &data)
      {
        const int n_blocks = 1 + size / (chunk_size * block_size);
        ::dealii::LinearAlgebra::CUDAWrappers::kernel::add_aV<Number>
          <<<n_blocks, block_size>>>(data.values_dev.get(),
                                     a,
                                     v_data.values_dev.get(),
                                     size);
        AssertCudaKernel();
      }

      static void
      add_avpbw(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &,
        const size_type size,
        const Number    a,
        const Number    b,
        const ::dealii::MemorySpace::
          MemorySpaceData<Number, ::dealii::MemorySpace::CUDA> &v_data,
        const ::dealii::MemorySpace::
          MemorySpaceData<Number, ::dealii::MemorySpace::CUDA> &w_data,
        ::dealii::MemorySpace::MemorySpaceData<Number,
                                               ::dealii::MemorySpace::CUDA>
          &data)
      {
        const int n_blocks = 1 + size / (chunk_size * block_size);
        ::dealii::LinearAlgebra::CUDAWrappers::kernel::add_aVbW<Number>
          <<<dim3(n_blocks, 1), dim3(block_size)>>>(data.values_dev.get(),
                                                    a,
                                                    v_data.values_dev.get(),
                                                    b,
                                                    w_data.values_dev.get(),
                                                    size);
        AssertCudaKernel();
      }

      static void
      sadd_xv(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &,
        const size_type size,
        const Number    x,
        const ::dealii::MemorySpace::
          MemorySpaceData<Number, ::dealii::MemorySpace::CUDA> &v_data,
        ::dealii::MemorySpace::MemorySpaceData<Number,
                                               ::dealii::MemorySpace::CUDA>
          &data)
      {
        const int n_blocks = 1 + size / (chunk_size * block_size);
        ::dealii::LinearAlgebra::CUDAWrappers::kernel::sadd<Number>
          <<<dim3(n_blocks, 1), dim3(block_size)>>>(
            x, data.values_dev.get(), 1., v_data.values_dev.get(), size);
        AssertCudaKernel();
      }

      static void
      sadd_xav(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &,
        const size_type size,
        const Number    x,
        const Number    a,
        const ::dealii::MemorySpace::
          MemorySpaceData<Number, ::dealii::MemorySpace::CUDA> &v_data,
        ::dealii::MemorySpace::MemorySpaceData<Number,
                                               ::dealii::MemorySpace::CUDA>
          &data)
      {
        const int n_blocks = 1 + size / (chunk_size * block_size);
        ::dealii::LinearAlgebra::CUDAWrappers::kernel::sadd<Number>
          <<<dim3(n_blocks, 1), dim3(block_size)>>>(
            x, data.values_dev.get(), a, v_data.values_dev.get(), size);
        AssertCudaKernel();
      }

      static void
      sadd_xavbw(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &,
        const size_type size,
        const Number    x,
        const Number    a,
        const Number    b,
        const ::dealii::MemorySpace::
          MemorySpaceData<Number, ::dealii::MemorySpace::CUDA> &v_data,
        const ::dealii::MemorySpace::
          MemorySpaceData<Number, ::dealii::MemorySpace::CUDA> &w_data,
        ::dealii::MemorySpace::MemorySpaceData<Number,
                                               ::dealii::MemorySpace::CUDA>
          &data)
      {
        const int n_blocks = 1 + size / (chunk_size * block_size);
        ::dealii::LinearAlgebra::CUDAWrappers::kernel::sadd<Number>
          <<<dim3(n_blocks, 1), dim3(block_size)>>>(x,
                                                    data.values_dev.get(),
                                                    a,
                                                    v_data.values_dev.get(),
                                                    b,
                                                    w_data.values_dev.get(),
                                                    size);
        AssertCudaKernel();
      }

      static void
      multiply_factor(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &,
        const size_type size,
        const Number    factor,
        ::dealii::MemorySpace::MemorySpaceData<Number,
                                               ::dealii::MemorySpace::CUDA>
          &data)
      {
        const int n_blocks = 1 + size / (chunk_size * block_size);
        ::dealii::LinearAlgebra::CUDAWrappers::kernel::vec_scale<Number>
          <<<n_blocks, block_size>>>(data.values_dev.get(), factor, size);
        AssertCudaKernel();
      }

      static void
      scale(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &,
        const size_type size,
        const ::dealii::MemorySpace::
          MemorySpaceData<Number, ::dealii::MemorySpace::CUDA> &v_data,
        ::dealii::MemorySpace::MemorySpaceData<Number,
                                               ::dealii::MemorySpace::CUDA>
          &data)
      {
        const int n_blocks = 1 + size / (chunk_size * block_size);
        ::dealii::LinearAlgebra::CUDAWrappers::kernel::scale<Number>
          <<<dim3(n_blocks, 1), dim3(block_size)>>>(data.values_dev.get(),
                                                    v_data.values_dev.get(),
                                                    size);
        AssertCudaKernel();
      }

      static void
      equ_au(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &,
        const size_type size,
        const Number    a,
        const ::dealii::MemorySpace::
          MemorySpaceData<Number, ::dealii::MemorySpace::CUDA> &v_data,
        ::dealii::MemorySpace::MemorySpaceData<Number,
                                               ::dealii::MemorySpace::CUDA>
          &data)
      {
        const int n_blocks = 1 + size / (chunk_size * block_size);
        ::dealii::LinearAlgebra::CUDAWrappers::kernel::equ<Number>
          <<<dim3(n_blocks, 1), dim3(block_size)>>>(data.values_dev.get(),
                                                    a,
                                                    v_data.values_dev.get(),
                                                    size);
        AssertCudaKernel();
      }

      static void
      equ_aubv(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &,
        const size_type size,
        const Number    a,
        const Number    b,
        const ::dealii::MemorySpace::
          MemorySpaceData<Number, ::dealii::MemorySpace::CUDA> &v_data,
        const ::dealii::MemorySpace::
          MemorySpaceData<Number, ::dealii::MemorySpace::CUDA> &w_data,
        ::dealii::MemorySpace::MemorySpaceData<Number,
                                               ::dealii::MemorySpace::CUDA>
          &data)
      {
        const int n_blocks = 1 + size / (chunk_size * block_size);
        ::dealii::LinearAlgebra::CUDAWrappers::kernel::equ<Number>
          <<<dim3(n_blocks, 1), dim3(block_size)>>>(data.values_dev.get(),
                                                    a,
                                                    v_data.values_dev.get(),
                                                    b,
                                                    w_data.values_dev.get(),
                                                    size);
        AssertCudaKernel();
      }

      static Number
      dot(const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &,
          const size_type size,
          const ::dealii::MemorySpace::
            MemorySpaceData<Number, ::dealii::MemorySpace::CUDA> &v_data,
          ::dealii::MemorySpace::MemorySpaceData<Number,
                                                 ::dealii::MemorySpace::CUDA>
            &data)
      {
        Number *    result_device;
        cudaError_t error_code = cudaMalloc(&result_device, sizeof(Number));
        AssertCuda(error_code);
        error_code = cudaMemset(result_device, 0, sizeof(Number));
        AssertCuda(error_code);

        const int n_blocks = 1 + size / (chunk_size * block_size);
        ::dealii::LinearAlgebra::CUDAWrappers::kernel::double_vector_reduction<
          Number,
          ::dealii::LinearAlgebra::CUDAWrappers::kernel::DotProduct<Number>>
          <<<dim3(n_blocks, 1), dim3(block_size)>>>(result_device,
                                                    data.values_dev.get(),
                                                    v_data.values_dev.get(),
                                                    static_cast<unsigned int>(
                                                      size));
        AssertCudaKernel();

        // Copy the result back to the host
        Number result;
        error_code = cudaMemcpy(&result,
                                result_device,
                                sizeof(Number),
                                cudaMemcpyDeviceToHost);
        AssertCuda(error_code);
        // Free the memory on the device
        error_code = cudaFree(result_device);
        AssertCuda(error_code);

        AssertIsFinite(result);

        return result;
      }

      template <typename real_type>
      static void
      norm_2(const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
               &             thread_loop_partitioner,
             const size_type size,
             real_type &     sum,
             ::dealii::MemorySpace::MemorySpaceData<Number,
                                                    ::dealii::MemorySpace::CUDA>
               &data)
      {
        sum = dot(thread_loop_partitioner, size, data, data);
      }

      static Number
      mean_value(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &,
        const size_type size,
        const ::dealii::MemorySpace::
          MemorySpaceData<Number, ::dealii::MemorySpace::CUDA> &data)
      {
        Number *    result_device;
        cudaError_t error_code = cudaMalloc(&result_device, sizeof(Number));
        AssertCuda(error_code);
        error_code = cudaMemset(result_device, 0, sizeof(Number));

        const int n_blocks = 1 + size / (chunk_size * block_size);
        ::dealii::LinearAlgebra::CUDAWrappers::kernel::reduction<
          Number,
          ::dealii::LinearAlgebra::CUDAWrappers::kernel::ElemSum<Number>>
          <<<dim3(n_blocks, 1), dim3(block_size)>>>(result_device,
                                                    data.values_dev.get(),
                                                    size);

        // Copy the result back to the host
        Number result;
        error_code = cudaMemcpy(&result,
                                result_device,
                                sizeof(Number),
                                cudaMemcpyDeviceToHost);
        AssertCuda(error_code);
        // Free the memory on the device
        error_code = cudaFree(result_device);
        AssertCuda(error_code);

        return result;
      }

      template <typename real_type>
      static void
      norm_1(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &,
        const size_type size,
        real_type &     sum,
        ::dealii::MemorySpace::MemorySpaceData<Number,
                                               ::dealii::MemorySpace::CUDA>
          &data)
      {
        Number *    result_device;
        cudaError_t error_code = cudaMalloc(&result_device, sizeof(Number));
        AssertCuda(error_code);
        error_code = cudaMemset(result_device, 0, sizeof(Number));

        const int n_blocks = 1 + size / (chunk_size * block_size);
        ::dealii::LinearAlgebra::CUDAWrappers::kernel::reduction<
          Number,
          ::dealii::LinearAlgebra::CUDAWrappers::kernel::L1Norm<Number>>
          <<<dim3(n_blocks, 1), dim3(block_size)>>>(result_device,
                                                    data.values_dev.get(),
                                                    size);

        // Copy the result back to the host
        error_code = cudaMemcpy(&sum,
                                result_device,
                                sizeof(Number),
                                cudaMemcpyDeviceToHost);
        AssertCuda(error_code);
        // Free the memory on the device
        error_code = cudaFree(result_device);
        AssertCuda(error_code);
      }

      template <typename real_type>
      static void
      norm_p(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &,
        const size_type,
        real_type &,
        real_type,
        ::dealii::MemorySpace::MemorySpaceData<Number,
                                               ::dealii::MemorySpace::CUDA> &)
      {
        Assert(false, ExcNotImplemented());
      }

      static Number
      add_and_dot(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &,
        const size_type size,
        const Number    a,
        const ::dealii::MemorySpace::
          MemorySpaceData<Number, ::dealii::MemorySpace::CUDA> &v_data,
        const ::dealii::MemorySpace::
          MemorySpaceData<Number, ::dealii::MemorySpace::CUDA> &w_data,
        ::dealii::MemorySpace::MemorySpaceData<Number,
                                               ::dealii::MemorySpace::CUDA>
          &data)
      {
        Number *    res_d;
        cudaError_t error_code = cudaMalloc(&res_d, sizeof(Number));
        AssertCuda(error_code);
        error_code = cudaMemset(res_d, 0, sizeof(Number));
        AssertCuda(error_code);

        const int n_blocks = 1 + size / (chunk_size * block_size);
        ::dealii::LinearAlgebra::CUDAWrappers::kernel::add_and_dot<Number>
          <<<dim3(n_blocks, 1), dim3(block_size)>>>(res_d,
                                                    data.values_dev.get(),
                                                    v_data.values_dev.get(),
                                                    w_data.values_dev.get(),
                                                    a,
                                                    size);

        Number res;
        error_code =
          cudaMemcpy(&res, res_d, sizeof(Number), cudaMemcpyDeviceToHost);
        AssertCuda(error_code);
        error_code = cudaFree(res_d);

        return res;
      }

      template <typename MemorySpace2>
      static void
      import(const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
               &                     thread_loop_partitioner,
             const size_type         size,
             VectorOperation::values operation,
             const ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpace2>
               &v_data,
             ::dealii::MemorySpace::MemorySpaceData<Number,
                                                    ::dealii::MemorySpace::CUDA>
               &data,
             typename std::enable_if<
               std::is_same<MemorySpace2, ::dealii::MemorySpace::CUDA>::value,
               int>::type = 0)
      {
        if (operation == VectorOperation::insert)
          {
            copy(thread_loop_partitioner, size, v_data, data);
          }
        else if (operation == VectorOperation::add)
          {
            add_vector(thread_loop_partitioner, size, v_data, data);
          }
        else
          {
            AssertThrow(false, ExcNotImplemented());
          }
      }

      template <typename MemorySpace2>
      static void
      import(const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
               & /*thread_loop_partitioner*/,
             const size_type         size,
             VectorOperation::values operation,
             const ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpace2>
               &v_data,
             ::dealii::MemorySpace::MemorySpaceData<Number,
                                                    ::dealii::MemorySpace::CUDA>
               &data,
             typename std::enable_if<
               std::is_same<MemorySpace2, ::dealii::MemorySpace::Host>::value,
               int>::type = 0)
      {
        if (operation == VectorOperation::insert)
          {
            cudaError_t cuda_error_code = cudaMemcpy(data.values_dev.get(),
                                                     v_data.values.get(),
                                                     size * sizeof(Number),
                                                     cudaMemcpyHostToDevice);
            AssertCuda(cuda_error_code);
          }
        else
          {
            AssertThrow(false, ExcNotImplemented());
          }
      }
    };
#endif
  } // namespace VectorOperations
} // namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif
