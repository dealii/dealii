// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#ifndef dealii_vector_operations_internal_h
#define dealii_vector_operations_internal_h

#include <deal.II/base/config.h>

#include <deal.II/base/memory_space.h>
#include <deal.II/base/memory_space_data.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/parallel.h>
#include <deal.II/base/types.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/lac/vector_operation.h>

#include <Kokkos_Core.hpp>

#include <cstdio>
#include <cstring>

#ifdef DEAL_II_WITH_TBB
#  include <tbb/blocked_range.h>
#  include <tbb/partitioner.h>
#endif


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
         std::complex<U>       *dest)
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



#ifdef DEAL_II_WITH_TBB
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
      TBBForFunctor(Functor        &functor,
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

      Functor        &functor;
      const size_type start;
      const size_type end;
      unsigned int    n_chunks;
      size_type       chunk_size;
    };
#endif

    template <typename Functor>
    void
    parallel_for(
      Functor        &functor,
      const size_type start,
      const size_type end,
      const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
        &partitioner)
    {
#ifdef DEAL_II_WITH_TBB
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
          std::fill(dst + begin, dst + end, Number());
        else
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

        if constexpr (std::is_trivially_copyable<Number>() &&
                      std::is_same_v<Number, OtherNumber>)
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
        , stored_factor(factor)
      {}

      void
      operator()(const size_type begin, const size_type end) const
      {
        // create a local copy of the variable to help the compiler with the
        // aliasing analysis
        const Number factor = stored_factor;

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
      const Number  stored_factor;
    };

    template <typename Number>
    struct Vectorization_add_av
    {
      Vectorization_add_av(Number *const       val,
                           const Number *const v_val,
                           const Number        factor)
        : val(val)
        , v_val(v_val)
        , stored_factor(factor)
      {}

      void
      operator()(const size_type begin, const size_type end) const
      {
        // create a local copy of the variable to help the compiler with the
        // aliasing analysis
        const Number factor = stored_factor;
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
      const Number        stored_factor;
    };

    template <typename Number>
    struct Vectorization_sadd_xav
    {
      Vectorization_sadd_xav(Number             *val,
                             const Number *const v_val,
                             const Number        a,
                             const Number        x)
        : val(val)
        , v_val(v_val)
        , stored_a(a)
        , stored_x(x)
      {}

      void
      operator()(const size_type begin, const size_type end) const
      {
        // create a local copy of the variable to help the compiler with the
        // aliasing analysis
        const Number x = stored_x, a = stored_a;

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
      const Number        stored_a;
      const Number        stored_x;
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
        , stored_factor(factor)
      {}

      void
      operator()(const size_type begin, const size_type end) const
      {
        const Number factor = stored_factor;

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
      const Number  stored_factor;
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
        , stored_a(a)
        , stored_b(b)
      {}

      void
      operator()(const size_type begin, const size_type end) const
      {
        const Number a = stored_a, b = stored_b;

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
      const Number        stored_a;
      const Number        stored_b;
    };

    template <typename Number>
    struct Vectorization_sadd_xv
    {
      Vectorization_sadd_xv(Number *const       val,
                            const Number *const v_val,
                            const Number        x)
        : val(val)
        , v_val(v_val)
        , stored_x(x)
      {}

      void
      operator()(const size_type begin, const size_type end) const
      {
        const Number x = stored_x;

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
      const Number        stored_x;
    };

    template <typename Number>
    struct Vectorization_sadd_xavbw
    {
      Vectorization_sadd_xavbw(Number       *val,
                               const Number *v_val,
                               const Number *w_val,
                               Number        x,
                               Number        a,
                               Number        b)
        : val(val)
        , v_val(v_val)
        , w_val(w_val)
        , stored_x(x)
        , stored_a(a)
        , stored_b(b)
      {}

      void
      operator()(const size_type begin, const size_type end) const
      {
        const Number x = stored_x, a = stored_a, b = stored_b;

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
      const Number        stored_x;
      const Number        stored_a;
      const Number        stored_b;
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
        , stored_a(a)
      {}

      void
      operator()(const size_type begin, const size_type end) const
      {
        const Number a = stored_a;

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
      const Number        stored_a;
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
        , stored_a(a)
        , stored_b(b)
      {}

      void
      operator()(const size_type begin, const size_type end) const
      {
        const Number a = stored_a, b = stored_b;

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
      const Number        stored_a;
      const Number        stored_b;
    };

    template <typename Number>
    struct Vectorization_equ_aubvcw
    {
      Vectorization_equ_aubvcw(Number       *val,
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
        , stored_a(a)
        , stored_b(b)
        , stored_c(c)
      {}

      void
      operator()(const size_type begin, const size_type end) const
      {
        const Number a = stored_a, b = stored_b, c = stored_c;

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
      const Number        stored_a;
      const Number        stored_b;
      const Number        stored_c;
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
      static constexpr bool vectorizes = std::is_same_v<Number, Number2> &&
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

      const Number  *X;
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
                         ResultType      &result)
    {
      if (first == last)
        {
          result = ResultType();
          return;
        }

      const size_type vec_size = last - first;
      if (vec_size <= vector_accumulation_recursion_threshold * 32)
        {
          // The vector is short enough so we perform the summation.  We store
          // the number of chunks (each 32 indices) for the given vector
          // length; all results are stored in outer_results[0,n_chunks). We
          // keep twice the number around to be able to do the pairwise
          // summation with a single for loop (see the loop over j below)
          ResultType outer_results[vector_accumulation_recursion_threshold * 2];

          // Select between the regular version and vectorized version based
          // on the number types we are given. To choose the vectorized
          // version often enough, we need to have all tasks but the last one
          // to be divisible by the vectorization length
          size_type n_chunks =
            do_accumulate(op,
                          vec_size,
                          first,
                          outer_results,
                          std::bool_constant<Operation::vectorizes>());

          AssertIndexRange(n_chunks,
                           vector_accumulation_recursion_threshold + 1);

          // now sum the results from the chunks stored in
          // outer_results[0,n_chunks) recursively
          unsigned int           j       = 0;
          constexpr unsigned int n_lanes = VectorizedArray<ResultType>::size();
          for (; j + 2 * n_lanes - 1 < n_chunks;
               j += 2 * n_lanes, n_chunks += n_lanes)
            {
              VectorizedArray<ResultType> a, b;
              a.load(outer_results + j);
              b.load(outer_results + j + n_lanes);
              a += b;
              a.store(outer_results + n_chunks);
            }

          // In the vectorized case, we know the loop bounds and can do things
          // more efficiently
          if (Operation::vectorizes)
            {
              AssertDimension(j + n_lanes, n_chunks);
              AssertIndexRange(n_chunks,
                               2 * vector_accumulation_recursion_threshold + 1);
              ResultType *result_ptr = outer_results + j;
              if (n_lanes >= 16)
                for (unsigned int i = 0; i < 8; ++i)
                  result_ptr[i] = result_ptr[i] + result_ptr[i + 8];
              if (n_lanes >= 8)
                for (unsigned int i = 0; i < 4; ++i)
                  result_ptr[i] = result_ptr[i] + result_ptr[i + 4];
              if (n_lanes >= 4)
                for (unsigned int i = 0; i < 2; ++i)
                  result_ptr[i] = result_ptr[i] + result_ptr[i + 2];
              result = result_ptr[0] + result_ptr[1];
            }
          else
            {
              // Without vectorization, we do not know the exact bounds, so we
              // need to continue the variable-length pairwise summation loop
              // from above
              for (; j + 1 < n_chunks; j += 2, ++n_chunks)
                outer_results[n_chunks] =
                  outer_results[j] + outer_results[j + 1];

              AssertIndexRange(n_chunks,
                               2 * vector_accumulation_recursion_threshold + 1);
              Assert(n_chunks > 0, ExcInternalError());
              result = outer_results[n_chunks - 1];
            }
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
          result = (r0 + r1) + (r2 + r3);
        }
    }


    // this is the inner working routine for the accumulation loops below. We
    // pulled this part out of the regular accumulate routine because we might
    // do this thing vectorized (see specialized function below; this is the
    // un-vectorized version). As opposed to the vector add functions above,
    // we here pass the functor 'op' by value, because we cannot create a copy
    // of the scalar inline, and instead make sure that the numbers get local
    // (and thus definitely not aliased) for the compiler
    template <typename Operation, typename ResultType>
    size_type
    do_accumulate(const Operation op,
                  const size_type vec_size,
                  const size_type start_index,
                  ResultType     *outer_results,
                  std::bool_constant<false>)
    {
      // Create local copy to indicate no aliasing to the compiler
      size_type index = start_index;

      // choose each chunk to have a width of 32, thereby the index
      // is incremented by 4*8 for each @p i.
      size_type n_chunks = vec_size / 32;
      for (size_type i = 0; i < n_chunks; ++i)
        {
          ResultType r = {};
          for (unsigned int k = 0; k < 2; ++k)
            {
              ResultType r0 = op(index);
              ResultType r1 = op(index + 1);
              ResultType r2 = op(index + 2);
              ResultType r3 = op(index + 3);
              index += 4;
              for (size_type j = 1; j < 4; ++j, index += 4)
                {
                  r0 += op(index);
                  r1 += op(index + 1);
                  r2 += op(index + 2);
                  r3 += op(index + 3);
                }
              r += (r0 + r1) + (r2 + r3);
            }
          outer_results[i] = r;
        }

      if (n_chunks * 32 < vec_size)
        {
          const size_type remainder       = vec_size - n_chunks * 32;
          const size_type inner_chunks    = remainder / 8;
          const size_type remainder_inner = remainder % 8;
          ResultType r0 = ResultType(), r1 = ResultType(), r2 = ResultType();
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
                outer_results[n_chunks++] = (r0 + r2) + r1;
                break;
            }
        }

      // make sure we worked through all indices
      AssertDimension(index, start_index + vec_size);

      return n_chunks;
    }



    // this is the inner working routine for the accumulation loops
    // below. This is the specialized case where we can vectorize. We request
    // the 'do_vectorized' routine of the operation instead of the regular one
    // which does several operations at once. As above, pass in the functor by
    // value to create a local copy of the scalar factors in the function (if
    // there are any).
    template <typename Operation, typename Number>
    size_type
    do_accumulate(const Operation op,
                  const size_type vec_size,
                  const size_type start_index,
                  Number         *outer_results,
                  std::bool_constant<true>)
    {
      // Create local copy to indicate no aliasing to the compiler
      size_type index = start_index;

      // we start from @p index and workout @p n_chunks each of size 32.
      // in order employ SIMD and work on @p nvecs at a time, we split this
      // loop yet again:
      // First we work on (n_chunks/nvecs) chunks, where each chunk processes
      // nvecs*(4*8) elements.

      constexpr size_type n_lanes        = VectorizedArray<Number>::size();
      const size_type     regular_chunks = vec_size / (32 * n_lanes);
      for (size_type i = 0; i < regular_chunks; ++i)
        {
          VectorizedArray<Number> r = {};
          for (unsigned int k = 0; k < 2; ++k)
            {
              VectorizedArray<Number> r0 = op.do_vectorized(index);
              VectorizedArray<Number> r1 = op.do_vectorized(index + n_lanes);
              VectorizedArray<Number> r2 =
                op.do_vectorized(index + 2 * n_lanes);
              VectorizedArray<Number> r3 =
                op.do_vectorized(index + 3 * n_lanes);
              index += n_lanes * 4;
              for (size_type j = 1; j < 4; ++j, index += n_lanes * 4)
                {
                  r0 += op.do_vectorized(index);
                  r1 += op.do_vectorized(index + n_lanes);
                  r2 += op.do_vectorized(index + 2 * n_lanes);
                  r3 += op.do_vectorized(index + 3 * n_lanes);
                }
              r += (r0 + r1) + (r2 + r3);
            }
          r.store(&outer_results[i * n_lanes]);
        }

      // If we are treating a case where the vector length is not divisible by
      // the vectorization length, need a cleanup loop
      // The remaining chunks are processed one by one starting from
      // regular_chunks * n_lanes; We do as much as possible with 2 SIMD
      // operations within each chunk. Here we assume that n_lanes < 32/2 = 16
      // as well as 16 % n_lanes == 0.
      static_assert(n_lanes <= 16 && 16 % n_lanes == 0,
                    "VectorizedArray::size() must be 1, 2, 4, 8, or 16");
      size_type       n_chunks        = regular_chunks * n_lanes;
      const size_type start_irregular = regular_chunks * n_lanes * 32;
      if (start_irregular < vec_size)
        {
          VectorizedArray<Number> r0  = VectorizedArray<Number>(),
                                  r1  = VectorizedArray<Number>();
          const size_type remainder   = vec_size - start_irregular;
          const size_type loop_length = remainder / (2 * n_lanes);
          for (size_type j = 0; j < loop_length; ++j, index += 2 * n_lanes)
            {
              r0 += op.do_vectorized(index);
              r1 += op.do_vectorized(index + n_lanes);
            }
          Number    scalar_part = Number();
          size_type last        = remainder % (2 * n_lanes);
          if (last > 0)
            {
              if (last >= n_lanes)
                {
                  r0 += op.do_vectorized(index);
                  index += n_lanes;
                  last -= n_lanes;
                }
              for (unsigned int i = 0; i < last; ++i)
                scalar_part += op(index++);
            }

          r0 += r1;
          r0.store(&outer_results[n_chunks]);
          outer_results[n_chunks] += scalar_part;

          // update n_chunks to denote range of entries to sum up in
          // outer_results[].
          n_chunks += n_lanes;
        }

      // make sure we worked through all indices
      AssertDimension(index, start_index + vec_size);

      return n_chunks;
    }



#ifdef DEAL_II_WITH_TBB
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
#ifndef DEBUG
    DEAL_II_ALWAYS_INLINE
#endif
      inline void
      parallel_reduce(
        const Operation &op,
        const size_type  start,
        const size_type  end,
        ResultType      &result,
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
          &partitioner)
    {
#ifdef DEAL_II_WITH_TBB
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
          std::is_same_v<MemorySpace, ::dealii::MemorySpace::Default> &&
            std::is_same_v<Number, Number2>,
          "For the Default MemorySpace Number and Number2 should be the same type");
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
        Number * /*values*/)
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
      import_elements(
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
                          &thread_loop_partitioner,
           const size_type size,
           const ::dealii::MemorySpace::
             MemorySpaceData<Number2, ::dealii::MemorySpace::Host> &v_data,
           ::dealii::MemorySpace::MemorySpaceData<Number,
                                                  ::dealii::MemorySpace::Host>
             &data)
      {
        Vector_copy<Number, Number2> copier(v_data.values.data(),
                                            data.values.data());
        parallel_for(copier, 0, size, thread_loop_partitioner);
      }

      static void
      set(const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
                         &thread_loop_partitioner,
          const size_type size,
          const Number    s,
          ::dealii::MemorySpace::MemorySpaceData<Number,
                                                 ::dealii::MemorySpace::Host>
            &data)
      {
        Vector_set<Number> setter(s, data.values.data());
        parallel_for(setter, 0, size, thread_loop_partitioner);
      }

      static void
      add_vector(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
                       &thread_loop_partitioner,
        const size_type size,
        const ::dealii::MemorySpace::
          MemorySpaceData<Number, ::dealii::MemorySpace::Host> &v_data,
        ::dealii::MemorySpace::MemorySpaceData<Number,
                                               ::dealii::MemorySpace::Host>
          &data)
      {
        Vectorization_add_v<Number> vector_add(data.values.data(),
                                               v_data.values.data());
        parallel_for(vector_add, 0, size, thread_loop_partitioner);
      }

      static void
      subtract_vector(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
                       &thread_loop_partitioner,
        const size_type size,
        const ::dealii::MemorySpace::
          MemorySpaceData<Number, ::dealii::MemorySpace::Host> &v_data,
        ::dealii::MemorySpace::MemorySpaceData<Number,
                                               ::dealii::MemorySpace::Host>
          &data)
      {
        Vectorization_subtract_v<Number> vector_subtract(data.values.data(),
                                                         v_data.values.data());
        parallel_for(vector_subtract, 0, size, thread_loop_partitioner);
      }

      static void
      add_factor(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
                       &thread_loop_partitioner,
        const size_type size,
        Number          a,
        ::dealii::MemorySpace::MemorySpaceData<Number,
                                               ::dealii::MemorySpace::Host>
          &data)
      {
        Vectorization_add_factor<Number> vector_add(data.values.data(), a);
        parallel_for(vector_add, 0, size, thread_loop_partitioner);
      }

      static void
      add_av(const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
                            &thread_loop_partitioner,
             const size_type size,
             const Number    a,
             const ::dealii::MemorySpace::
               MemorySpaceData<Number, ::dealii::MemorySpace::Host> &v_data,
             ::dealii::MemorySpace::MemorySpaceData<Number,
                                                    ::dealii::MemorySpace::Host>
               &data)
      {
        Vectorization_add_av<Number> vector_add(data.values.data(),
                                                v_data.values.data(),
                                                a);
        parallel_for(vector_add, 0, size, thread_loop_partitioner);
      }

      static void
      add_avpbw(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
                       &thread_loop_partitioner,
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
          data.values.data(), v_data.values.data(), w_data.values.data(), a, b);
        parallel_for(vector_add, 0, size, thread_loop_partitioner);
      }

      static void
      sadd_xv(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
                       &thread_loop_partitioner,
        const size_type size,
        const Number    x,
        const ::dealii::MemorySpace::
          MemorySpaceData<Number, ::dealii::MemorySpace::Host> &v_data,
        ::dealii::MemorySpace::MemorySpaceData<Number,
                                               ::dealii::MemorySpace::Host>
          &data)
      {
        Vectorization_sadd_xv<Number> vector_sadd(data.values.data(),
                                                  v_data.values.data(),
                                                  x);
        parallel_for(vector_sadd, 0, size, thread_loop_partitioner);
      }

      static void
      sadd_xav(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
                       &thread_loop_partitioner,
        const size_type size,
        const Number    x,
        const Number    a,
        const ::dealii::MemorySpace::
          MemorySpaceData<Number, ::dealii::MemorySpace::Host> &v_data,
        ::dealii::MemorySpace::MemorySpaceData<Number,
                                               ::dealii::MemorySpace::Host>
          &data)
      {
        Vectorization_sadd_xav<Number> vector_sadd(data.values.data(),
                                                   v_data.values.data(),
                                                   a,
                                                   x);
        parallel_for(vector_sadd, 0, size, thread_loop_partitioner);
      }

      static void
      sadd_xavbw(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
                       &thread_loop_partitioner,
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
        Vectorization_sadd_xavbw<Number> vector_sadd(data.values.data(),
                                                     v_data.values.data(),
                                                     w_data.values.data(),
                                                     x,
                                                     a,
                                                     b);
        parallel_for(vector_sadd, 0, size, thread_loop_partitioner);
      }

      static void
      multiply_factor(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
                       &thread_loop_partitioner,
        const size_type size,
        const Number    factor,
        ::dealii::MemorySpace::MemorySpaceData<Number,
                                               ::dealii::MemorySpace::Host>
          &data)
      {
        Vectorization_multiply_factor<Number> vector_multiply(
          data.values.data(), factor);
        parallel_for(vector_multiply, 0, size, thread_loop_partitioner);
      }

      static void
      scale(const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
                           &thread_loop_partitioner,
            const size_type size,
            const ::dealii::MemorySpace::
              MemorySpaceData<Number, ::dealii::MemorySpace::Host> &v_data,
            ::dealii::MemorySpace::MemorySpaceData<Number,
                                                   ::dealii::MemorySpace::Host>
              &data)
      {
        Vectorization_scale<Number> vector_scale(data.values.data(),
                                                 v_data.values.data());
        parallel_for(vector_scale, 0, size, thread_loop_partitioner);
      }

      static void
      equ_au(const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
                            &thread_loop_partitioner,
             const size_type size,
             const Number    a,
             const ::dealii::MemorySpace::
               MemorySpaceData<Number, ::dealii::MemorySpace::Host> &v_data,
             ::dealii::MemorySpace::MemorySpaceData<Number,
                                                    ::dealii::MemorySpace::Host>
               &data)
      {
        Vectorization_equ_au<Number> vector_equ(data.values.data(),
                                                v_data.values.data(),
                                                a);
        parallel_for(vector_equ, 0, size, thread_loop_partitioner);
      }

      static void
      equ_aubv(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
                       &thread_loop_partitioner,
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
          data.values.data(), v_data.values.data(), w_data.values.data(), a, b);
        parallel_for(vector_equ, 0, size, thread_loop_partitioner);
      }

      static Number
      dot(const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
                         &thread_loop_partitioner,
          const size_type size,
          const ::dealii::MemorySpace::
            MemorySpaceData<Number2, ::dealii::MemorySpace::Host> &v_data,
          ::dealii::MemorySpace::MemorySpaceData<Number,
                                                 ::dealii::MemorySpace::Host>
            &data)
      {
        Number                                                   sum;
        dealii::internal::VectorOperations::Dot<Number, Number2> dot(
          data.values.data(), v_data.values.data());
        dealii::internal::VectorOperations::parallel_reduce(
          dot, 0, size, sum, thread_loop_partitioner);
        AssertIsFinite(sum);

        return sum;
      }

      template <typename real_type>
      static void
      norm_2(const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
                            &thread_loop_partitioner,
             const size_type size,
             real_type      &sum,
             ::dealii::MemorySpace::MemorySpaceData<Number,
                                                    ::dealii::MemorySpace::Host>
               &data)
      {
        Norm2<Number, real_type> norm2(data.values.data());
        parallel_reduce(norm2, 0, size, sum, thread_loop_partitioner);
      }

      static Number
      mean_value(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
                       &thread_loop_partitioner,
        const size_type size,
        const ::dealii::MemorySpace::
          MemorySpaceData<Number, ::dealii::MemorySpace::Host> &data)
      {
        Number            sum;
        MeanValue<Number> mean(data.values.data());
        parallel_reduce(mean, 0, size, sum, thread_loop_partitioner);

        return sum;
      }

      template <typename real_type>
      static void
      norm_1(const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
                            &thread_loop_partitioner,
             const size_type size,
             real_type      &sum,
             ::dealii::MemorySpace::MemorySpaceData<Number,
                                                    ::dealii::MemorySpace::Host>
                            &data,
             const size_type optional_offset = 0)
      {
        Norm1<Number, real_type> norm1(data.values.data());
        parallel_reduce(norm1,
                        optional_offset,
                        optional_offset + size,
                        sum,
                        thread_loop_partitioner);
      }

      template <typename real_type>
      static void
      norm_p(const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
                            &thread_loop_partitioner,
             const size_type size,
             real_type      &sum,
             const real_type p,
             ::dealii::MemorySpace::MemorySpaceData<Number,
                                                    ::dealii::MemorySpace::Host>
               &data)
      {
        NormP<Number, real_type> normp(data.values.data(), p);
        parallel_reduce(normp, 0, size, sum, thread_loop_partitioner);
      }

      static Number
      add_and_dot(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
                       &thread_loop_partitioner,
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
        AddAndDot<Number> adder(data.values.data(),
                                v_data.values.data(),
                                w_data.values.data(),
                                a);
        parallel_reduce(adder, 0, size, sum, thread_loop_partitioner);

        return sum;
      }

      template <typename MemorySpace2>
      static void
      import_elements(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
                               &thread_loop_partitioner,
        const size_type         size,
        VectorOperation::values operation,
        const ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpace2>
          &v_data,
        ::dealii::MemorySpace::MemorySpaceData<Number,
                                               ::dealii::MemorySpace::Host>
          &data,
        std::enable_if_t<
          std::is_same_v<MemorySpace2, dealii::MemorySpace::Host>,
          int> = 0)
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
      import_elements(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
          & /*thread_loop_partitioner*/,
        const size_type         size,
        VectorOperation::values operation,
        const ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpace2>
          &v_data,
        ::dealii::MemorySpace::MemorySpaceData<Number,
                                               ::dealii::MemorySpace::Host>
          &data,
        std::enable_if_t<
          std::is_same_v<MemorySpace2, ::dealii::MemorySpace::Default>,
          int> = 0)
      {
        if (operation == VectorOperation::insert)
          {
            Kokkos::deep_copy(
              Kokkos::subview(data.values,
                              Kokkos::pair<size_type, size_type>(0, size)),
              Kokkos::subview(v_data.values,
                              Kokkos::pair<size_type, size_type>(0, size)));
          }
        else
          {
            AssertThrow(false, ExcNotImplemented());
          }
      }
    };



    template <typename Number>
    struct functions<Number, Number, ::dealii::MemorySpace::Default>
    {
      static void
      copy(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &,
        const size_type size,
        const ::dealii::MemorySpace::
          MemorySpaceData<Number, ::dealii::MemorySpace::Default> &v_data,
        ::dealii::MemorySpace::MemorySpaceData<Number,
                                               ::dealii::MemorySpace::Default>
          &data)
      {
        Kokkos::deep_copy(
          Kokkos::subview(data.values,
                          Kokkos::pair<size_type, size_type>(0, size)),
          Kokkos::subview(v_data.values,
                          Kokkos::pair<size_type, size_type>(0, size)));
      }

      static void
      set(const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &,
          const size_type size,
          const Number    s,
          ::dealii::MemorySpace::MemorySpaceData<Number,
                                                 ::dealii::MemorySpace::Default>
            &data)
      {
        Kokkos::deep_copy(
          Kokkos::subview(data.values,
                          Kokkos::pair<size_type, size_type>(0, size)),
          s);
      }

      static void
      add_vector(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &,
        const size_type size,
        const ::dealii::MemorySpace::
          MemorySpaceData<Number, ::dealii::MemorySpace::Default> &v_data,
        ::dealii::MemorySpace::MemorySpaceData<Number,
                                               ::dealii::MemorySpace::Default>
          &data)
      {
        auto exec = typename ::dealii::MemorySpace::Default::kokkos_space::
          execution_space{};
        Kokkos::parallel_for(
          "dealii::add_vector",
          Kokkos::RangePolicy<
            ::dealii::MemorySpace::Default::kokkos_space::execution_space>(
            exec, 0, size),
          KOKKOS_LAMBDA(int i) { data.values(i) += v_data.values(i); });
        exec.fence();
      }

      static void
      subtract_vector(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &,
        const size_type size,
        const ::dealii::MemorySpace::
          MemorySpaceData<Number, ::dealii::MemorySpace::Default> &v_data,
        ::dealii::MemorySpace::MemorySpaceData<Number,
                                               ::dealii::MemorySpace::Default>
          &data)
      {
        auto exec = typename ::dealii::MemorySpace::Default::kokkos_space::
          execution_space{};
        Kokkos::parallel_for(
          "dealii::subtract_vector",
          Kokkos::RangePolicy<
            ::dealii::MemorySpace::Default::kokkos_space::execution_space>(
            exec, 0, size),
          KOKKOS_LAMBDA(size_type i) { data.values(i) -= v_data.values(i); });
        exec.fence();
      }

      static void
      add_factor(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &,
        const size_type size,
        Number          a,
        ::dealii::MemorySpace::MemorySpaceData<Number,
                                               ::dealii::MemorySpace::Default>
          &data)
      {
        auto exec = typename ::dealii::MemorySpace::Default::kokkos_space::
          execution_space{};
        Kokkos::parallel_for(
          "dealii::add_factor",
          Kokkos::RangePolicy<
            ::dealii::MemorySpace::Default::kokkos_space::execution_space>(
            exec, 0, size),
          KOKKOS_LAMBDA(size_type i) { data.values(i) += a; });
        exec.fence();
      }

      static void
      add_av(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &,
        const size_type size,
        const Number    a,
        const ::dealii::MemorySpace::
          MemorySpaceData<Number, ::dealii::MemorySpace::Default> &v_data,
        ::dealii::MemorySpace::MemorySpaceData<Number,
                                               ::dealii::MemorySpace::Default>
          &data)
      {
        auto exec = typename ::dealii::MemorySpace::Default::kokkos_space::
          execution_space{};
        Kokkos::parallel_for(
          "dealii::add_av",
          Kokkos::RangePolicy<
            ::dealii::MemorySpace::Default::kokkos_space::execution_space>(
            exec, 0, size),
          KOKKOS_LAMBDA(size_type i) {
            data.values(i) += a * v_data.values(i);
          });
        exec.fence();
      }

      static void
      add_avpbw(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &,
        const size_type size,
        const Number    a,
        const Number    b,
        const ::dealii::MemorySpace::
          MemorySpaceData<Number, ::dealii::MemorySpace::Default> &v_data,
        const ::dealii::MemorySpace::
          MemorySpaceData<Number, ::dealii::MemorySpace::Default> &w_data,
        ::dealii::MemorySpace::MemorySpaceData<Number,
                                               ::dealii::MemorySpace::Default>
          &data)
      {
        auto exec = typename ::dealii::MemorySpace::Default::kokkos_space::
          execution_space{};
        Kokkos::parallel_for(
          "dealii::add_avpbw",
          Kokkos::RangePolicy<
            ::dealii::MemorySpace::Default::kokkos_space::execution_space>(
            exec, 0, size),
          KOKKOS_LAMBDA(size_type i) {
            data.values(i) += a * v_data.values(i) + b * w_data.values(i);
          });
        exec.fence();
      }

      static void
      sadd_xv(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &,
        const size_type size,
        const Number    x,
        const ::dealii::MemorySpace::
          MemorySpaceData<Number, ::dealii::MemorySpace::Default> &v_data,
        ::dealii::MemorySpace::MemorySpaceData<Number,
                                               ::dealii::MemorySpace::Default>
          &data)
      {
        auto exec = typename ::dealii::MemorySpace::Default::kokkos_space::
          execution_space{};
        Kokkos::parallel_for(
          "dealii::sadd_xv",
          Kokkos::RangePolicy<
            ::dealii::MemorySpace::Default::kokkos_space::execution_space>(
            exec, 0, size),
          KOKKOS_LAMBDA(size_type i) {
            data.values(i) = x * data.values(i) + v_data.values(i);
          });
        exec.fence();
      }

      static void
      sadd_xav(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &,
        const size_type size,
        const Number    x,
        const Number    a,
        const ::dealii::MemorySpace::
          MemorySpaceData<Number, ::dealii::MemorySpace::Default> &v_data,
        ::dealii::MemorySpace::MemorySpaceData<Number,
                                               ::dealii::MemorySpace::Default>
          &data)
      {
        auto exec = typename ::dealii::MemorySpace::Default::kokkos_space::
          execution_space{};
        Kokkos::parallel_for(
          "dealii::sadd_xav",
          Kokkos::RangePolicy<
            ::dealii::MemorySpace::Default::kokkos_space::execution_space>(
            exec, 0, size),
          KOKKOS_LAMBDA(size_type i) {
            data.values(i) = x * data.values(i) + a * v_data.values(i);
          });
        exec.fence();
      }

      static void
      sadd_xavbw(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &,
        const size_type size,
        const Number    x,
        const Number    a,
        const Number    b,
        const ::dealii::MemorySpace::
          MemorySpaceData<Number, ::dealii::MemorySpace::Default> &v_data,
        const ::dealii::MemorySpace::
          MemorySpaceData<Number, ::dealii::MemorySpace::Default> &w_data,
        ::dealii::MemorySpace::MemorySpaceData<Number,
                                               ::dealii::MemorySpace::Default>
          &data)
      {
        auto exec = typename ::dealii::MemorySpace::Default::kokkos_space::
          execution_space{};
        Kokkos::parallel_for(
          "dealii::sadd_xavbw",
          Kokkos::RangePolicy<
            ::dealii::MemorySpace::Default::kokkos_space::execution_space>(
            exec, 0, size),
          KOKKOS_LAMBDA(size_type i) {
            data.values(i) =
              x * data.values(i) + a * v_data.values(i) + b * w_data.values(i);
          });
        exec.fence();
      }

      static void
      multiply_factor(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &,
        const size_type size,
        const Number    factor,
        ::dealii::MemorySpace::MemorySpaceData<Number,
                                               ::dealii::MemorySpace::Default>
          &data)
      {
        auto exec = typename ::dealii::MemorySpace::Default::kokkos_space::
          execution_space{};
        Kokkos::parallel_for(
          "dealii::multiply_factor",
          Kokkos::RangePolicy<
            ::dealii::MemorySpace::Default::kokkos_space::execution_space>(
            exec, 0, size),
          KOKKOS_LAMBDA(size_type i) { data.values(i) *= factor; });
        exec.fence();
      }

      static void
      scale(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &,
        const size_type size,
        const ::dealii::MemorySpace::
          MemorySpaceData<Number, ::dealii::MemorySpace::Default> &v_data,
        ::dealii::MemorySpace::MemorySpaceData<Number,
                                               ::dealii::MemorySpace::Default>
          &data)
      {
        auto exec = typename ::dealii::MemorySpace::Default::kokkos_space::
          execution_space{};
        Kokkos::parallel_for(
          "dealii::scale",
          Kokkos::RangePolicy<
            ::dealii::MemorySpace::Default::kokkos_space::execution_space>(
            exec, 0, size),
          KOKKOS_LAMBDA(size_type i) { data.values(i) *= v_data.values(i); });
        exec.fence();
      }

      static void
      equ_au(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &,
        const size_type size,
        const Number    a,
        const ::dealii::MemorySpace::
          MemorySpaceData<Number, ::dealii::MemorySpace::Default> &v_data,
        ::dealii::MemorySpace::MemorySpaceData<Number,
                                               ::dealii::MemorySpace::Default>
          &data)
      {
        auto exec = typename ::dealii::MemorySpace::Default::kokkos_space::
          execution_space{};
        Kokkos::parallel_for(
          "dealii::equ_au",
          Kokkos::RangePolicy<
            ::dealii::MemorySpace::Default::kokkos_space::execution_space>(
            exec, 0, size),
          KOKKOS_LAMBDA(size_type i) {
            data.values(i) = a * v_data.values(i);
          });
        exec.fence();
      }

      static void
      equ_aubv(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &,
        const size_type size,
        const Number    a,
        const Number    b,
        const ::dealii::MemorySpace::
          MemorySpaceData<Number, ::dealii::MemorySpace::Default> &v_data,
        const ::dealii::MemorySpace::
          MemorySpaceData<Number, ::dealii::MemorySpace::Default> &w_data,
        ::dealii::MemorySpace::MemorySpaceData<Number,
                                               ::dealii::MemorySpace::Default>
          &data)
      {
        auto exec = typename ::dealii::MemorySpace::Default::kokkos_space::
          execution_space{};
        Kokkos::parallel_for(
          "dealii::equ_aubv",
          Kokkos::RangePolicy<
            ::dealii::MemorySpace::Default::kokkos_space::execution_space>(
            exec, 0, size),
          KOKKOS_LAMBDA(size_type i) {
            data.values(i) = a * v_data.values(i) + b * w_data.values(i);
          });
        exec.fence();
      }

      static Number
      dot(const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &,
          const size_type size,
          const ::dealii::MemorySpace::
            MemorySpaceData<Number, ::dealii::MemorySpace::Default> &v_data,
          ::dealii::MemorySpace::MemorySpaceData<Number,
                                                 ::dealii::MemorySpace::Default>
            &data)
      {
        Number result;

        auto exec = typename ::dealii::MemorySpace::Default::kokkos_space::
          execution_space{};
        Kokkos::parallel_reduce(
          "dealii::dot",
          Kokkos::RangePolicy<
            ::dealii::MemorySpace::Default::kokkos_space::execution_space>(
            exec, 0, size),
          KOKKOS_LAMBDA(size_type i, Number & update) {
            update += data.values(i) * v_data.values(i);
          },
          result);

        AssertIsFinite(result);
        return result;
      }

      template <typename real_type>
      static void
      norm_2(const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
                            &thread_loop_partitioner,
             const size_type size,
             real_type      &sum,
             ::dealii::MemorySpace::
               MemorySpaceData<Number, ::dealii::MemorySpace::Default> &data)
      {
        sum = dot(thread_loop_partitioner, size, data, data);
      }

      static Number
      mean_value(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &,
        const size_type size,
        const ::dealii::MemorySpace::
          MemorySpaceData<Number, ::dealii::MemorySpace::Default> &data)
      {
        Number result;

        auto exec = typename ::dealii::MemorySpace::Default::kokkos_space::
          execution_space{};
        Kokkos::parallel_reduce(
          "dealii::mean_value",
          Kokkos::RangePolicy<
            ::dealii::MemorySpace::Default::kokkos_space::execution_space>(
            exec, 0, size),
          KOKKOS_LAMBDA(size_type i, Number & update) {
            update += data.values(i);
          },
          result);

        AssertIsFinite(result);
        return result;
      }

      template <typename real_type>
      static void
      norm_1(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &,
        const size_type size,
        real_type      &sum,
        ::dealii::MemorySpace::MemorySpaceData<Number,
                                               ::dealii::MemorySpace::Default>
                       &data,
        const size_type optional_offset = 0)
      {
        auto exec = typename ::dealii::MemorySpace::Default::kokkos_space::
          execution_space{};
        Kokkos::parallel_reduce(
          "dealii::norm_1",
          Kokkos::RangePolicy<
            ::dealii::MemorySpace::Default::kokkos_space::execution_space>(
            exec, optional_offset, optional_offset + size),
          KOKKOS_LAMBDA(size_type i, Number & update) {
#if DEAL_II_KOKKOS_VERSION_GTE(3, 7, 0)
            update += Kokkos::abs(data.values(i));
#else
            update += Kokkos::Experimental::fabs(data.values(i));
#endif
          },
          sum);
      }

      template <typename real_type>
      static void
      norm_p(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &,
        const size_type size,
        real_type      &sum,
        real_type       exp,
        ::dealii::MemorySpace::MemorySpaceData<Number,
                                               ::dealii::MemorySpace::Default>
          &data)
      {
        auto exec = typename ::dealii::MemorySpace::Default::kokkos_space::
          execution_space{};
        Kokkos::parallel_reduce(
          "dealii::norm_p",
          Kokkos::RangePolicy<
            ::dealii::MemorySpace::Default::kokkos_space::execution_space>(
            exec, 0, size),
          KOKKOS_LAMBDA(size_type i, Number & update) {
#if DEAL_II_KOKKOS_VERSION_GTE(3, 7, 0)
            update += Kokkos::pow(Kokkos::abs(data.values(i)), exp);
#else
            update += Kokkos::Experimental::pow(
              Kokkos::Experimental::fabs(data.values(i)), exp);
#endif
          },
          sum);
      }

      static Number
      add_and_dot(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner> &,
        const size_type size,
        const Number    a,
        const ::dealii::MemorySpace::
          MemorySpaceData<Number, ::dealii::MemorySpace::Default> &v_data,
        const ::dealii::MemorySpace::
          MemorySpaceData<Number, ::dealii::MemorySpace::Default> &w_data,
        ::dealii::MemorySpace::MemorySpaceData<Number,
                                               ::dealii::MemorySpace::Default>
          &data)
      {
        Number res;

        auto exec = typename ::dealii::MemorySpace::Default::kokkos_space::
          execution_space{};
        Kokkos::parallel_reduce(
          "dealii::add_and_dot",
          Kokkos::RangePolicy<
            ::dealii::MemorySpace::Default::kokkos_space::execution_space>(
            exec, 0, size),
          KOKKOS_LAMBDA(size_type i, Number & update) {
            data.values(i) += a * v_data.values(i);
            update +=
              data.values(i) * Number(numbers::NumberTraits<Number>::conjugate(
                                 w_data.values(i)));
          },
          res);

        return res;
      }

      template <typename MemorySpace2>
      static void
      import_elements(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
                               &thread_loop_partitioner,
        const size_type         size,
        VectorOperation::values operation,
        const ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpace2>
          &v_data,
        ::dealii::MemorySpace::MemorySpaceData<Number,
                                               ::dealii::MemorySpace::Default>
          &data,
        std::enable_if_t<
          std::is_same_v<MemorySpace2, ::dealii::MemorySpace::Default>,
          int> = 0)
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
      import_elements(
        const std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
          & /*thread_loop_partitioner*/,
        const size_type         size,
        VectorOperation::values operation,
        const ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpace2>
          &v_data,
        ::dealii::MemorySpace::MemorySpaceData<Number,
                                               ::dealii::MemorySpace::Default>
          &data,
        std::enable_if_t<
          std::is_same_v<MemorySpace2, ::dealii::MemorySpace::Host>,
          int> = 0)
      {
        if (operation == VectorOperation::insert)
          {
            Kokkos::deep_copy(
              Kokkos::subview(data.values,
                              Kokkos::pair<size_type, size_type>(0, size)),
              Kokkos::subview(v_data.values,
                              Kokkos::pair<size_type, size_type>(0, size)));
          }
        else
          {
            AssertThrow(false, ExcNotImplemented());
          }
      }
    };
  } // namespace VectorOperations
} // namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif
