// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


#include <deal.II/base/memory_space.h>

#include <deal.II/lac/solver_cg.h>

#include <Kokkos_Core.hpp>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace SolverCG
  {
    /**
     * Value type for the reduction of the seven inner products in
     * do_reductions_device().
     */
    template <typename Number>
    struct InnerProducts
    {
      Number values[7];

      KOKKOS_INLINE_FUNCTION
      InnerProducts()
      {
        for (Number &a : values)
          a = Number();
      }

      KOKKOS_INLINE_FUNCTION InnerProducts &
      operator+=(const InnerProducts &other)
      {
        for (unsigned int i = 0; i < 7; ++i)
          values[i] += other.values[i];
        return *this;
      }
    };



    /**
     * Custom Kokkos reducer for InnerProducts. We use a reducer class rather
     * than a specialization of Kokkos::reduction_identity because the latter
     * would need to be placed in namespace Kokkos, outside of namespace
     * dealii, which is not compatible with the conversion of this file into a
     * C++20 module partition.
     */
    template <typename Number>
    struct InnerProductsReducer
    {
      using reducer    = InnerProductsReducer;
      using value_type = InnerProducts<Number>;
      using result_view_type =
        Kokkos::View<value_type,
                     Kokkos::HostSpace,
                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

      KOKKOS_INLINE_FUNCTION
      InnerProductsReducer(value_type &result)
        : value(&result)
      {}

      KOKKOS_INLINE_FUNCTION void
      join(value_type &dst, const value_type &src) const
      {
        dst += src;
      }

      KOKKOS_INLINE_FUNCTION void
      init(value_type &val) const
      {
        val = value_type();
      }

      KOKKOS_INLINE_FUNCTION value_type &
      reference() const
      {
        return *value;
      }

      KOKKOS_INLINE_FUNCTION result_view_type
      view() const
      {
        return result_view_type(value);
      }

      KOKKOS_INLINE_FUNCTION bool
      references_scalar() const
      {
        return true;
      }

    private:
      value_type *value;
    };



    using ExecutionSpace = MemorySpace::Default::kokkos_space::execution_space;
    using Policy =
      Kokkos::RangePolicy<ExecutionSpace, Kokkos::IndexType<std::size_t>>;



    /**
     * Apply the diagonal preconditioner to the value @p value of entry
     * @p j, or return the value unchanged for the identity preconditioner.
     */
    template <bool has_diagonal, typename Number>
    KOKKOS_FORCEINLINE_FUNCTION Number
    apply_preconditioner(const Number     *diagonal,
                         const std::size_t j,
                         const Number      value)
    {
      if constexpr (has_diagonal)
        return diagonal[j] * value;
      else
        {
          (void)diagonal;
          (void)j;
          return value;
        }
    }



    template <bool has_diagonal, typename Number>
    void
    do_vector_updates_device(const unsigned int iteration_index,
                             const std::size_t  size,
                             const Number       alpha,
                             const Number       beta,
                             const Number       previous_alpha,
                             const Number       previous_beta,
                             const Number      *diagonal,
                             Number            *x,
                             Number            *r,
                             Number            *p,
                             const Number      *v)
    {
      ExecutionSpace exec;
      if (iteration_index == 1)
        {
          Kokkos::parallel_for(
            "dealii::cg_update_vectors_first",
            Policy(exec, 0, size),
            KOKKOS_LAMBDA(const std::size_t j) {
              p[j] = apply_preconditioner<has_diagonal>(diagonal, j, r[j]);
            });
        }
      else if (iteration_index % 2 == 0 && beta != Number())
        {
          Kokkos::parallel_for(
            "dealii::cg_update_vectors_even",
            Policy(exec, 0, size),
            KOKKOS_LAMBDA(const std::size_t j) {
              const Number rj = r[j] - alpha * v[j];
              r[j]            = rj;
              p[j] = beta * p[j] +
                     apply_preconditioner<has_diagonal>(diagonal, j, rj);
            });
        }
      else if (iteration_index % 2 == 0 && beta == Number())
        {
          // Case where beta is zero: we cannot reconstruct p_{j-1} in the
          // next iteration, and must hence compute x here. This can happen
          // before termination
          Kokkos::parallel_for(
            "dealii::cg_update_vectors_even_zero_beta",
            Policy(exec, 0, size),
            KOKKOS_LAMBDA(const std::size_t j) {
              const Number rj = r[j] - alpha * v[j];
              r[j]            = rj;
              x[j] += alpha * p[j];
              p[j] = apply_preconditioner<has_diagonal>(diagonal, j, rj);
            });
        }
      else
        {
          const Number alpha_plus_previous_alpha_over_beta =
            alpha + previous_alpha / previous_beta;
          const Number previous_alpha_over_beta =
            previous_alpha / previous_beta;
          Kokkos::parallel_for(
            "dealii::cg_update_vectors_odd",
            Policy(exec, 0, size),
            KOKKOS_LAMBDA(const std::size_t j) {
              const Number pj     = p[j];
              const Number rj     = r[j];
              const Number vj     = v[j];
              const Number new_rj = rj - alpha * vj;
              const Number xj_update =
                alpha_plus_previous_alpha_over_beta * pj -
                previous_alpha_over_beta *
                  apply_preconditioner<has_diagonal>(diagonal, j, rj);
              p[j] = beta * pj +
                     apply_preconditioner<has_diagonal>(diagonal, j, new_rj);
              r[j] = new_rj;
              x[j] += xj_update;
            });
        }
    }



    template <typename Number>
    void
    do_vector_updates_device(const unsigned int iteration_index,
                             const std::size_t  size,
                             const Number       alpha,
                             const Number       beta,
                             const Number       previous_alpha,
                             const Number       previous_beta,
                             const Number      *diagonal,
                             Number            *x,
                             Number            *r,
                             Number            *p,
                             const Number      *v)
    {
      if (size == 0)
        return;
      else if (diagonal != nullptr)
        do_vector_updates_device<true>(iteration_index,
                                       size,
                                       alpha,
                                       beta,
                                       previous_alpha,
                                       previous_beta,
                                       diagonal,
                                       x,
                                       r,
                                       p,
                                       v);
      else
        do_vector_updates_device<false>(iteration_index,
                                        size,
                                        alpha,
                                        beta,
                                        previous_alpha,
                                        previous_beta,
                                        diagonal,
                                        x,
                                        r,
                                        p,
                                        v);
    }



    // Reductions with identity preconditioner -> specialization of the
    // diagonal code further down
    template <typename Number>
    std::array<Number, 7>
    do_reductions_device(const std::size_t size,
                         const Number     *p,
                         const Number     *v,
                         const Number     *r)
    {
      ExecutionSpace        exec;
      InnerProducts<Number> sums;
      Kokkos::parallel_reduce(
        "dealii::cg_inner_products",
        Policy(exec, 0, size),
        KOKKOS_LAMBDA(const std::size_t j, InnerProducts<Number> &s) {
          const Number pj = p[j];
          const Number vj = v[j];
          const Number rj = r[j];
          s.values[0] += pj * vj;
          s.values[1] += vj * vj;
          s.values[2] += rj * vj;
          s.values[3] += rj * rj;
        },
        InnerProductsReducer<Number>(sums));

      // Finalize the reduction by setting the values we did not compute for
      // identity preconditioner but are expected in the code
      std::array<Number, 7> result;
      for (unsigned int i = 0; i < 4; ++i)
        result[i] = sums.values[i];
      result[4] = sums.values[2];
      result[5] = sums.values[1];
      result[6] = sums.values[3];
      return result;
    }



    template <typename Number>
    std::array<Number, 7>
    do_reductions_device(const std::size_t size,
                         const Number     *diagonal,
                         const Number     *p,
                         const Number     *v,
                         const Number     *r)
    {
      if (size == 0)
        return std::array<Number, 7>{};
      else if (diagonal == nullptr)
        return do_reductions_device(size, p, v, r);
      else
        {
          ExecutionSpace        exec;
          InnerProducts<Number> sums;
          Kokkos::parallel_reduce(
            "dealii::cg_inner_products",
            Policy(exec, 0, size),
            KOKKOS_LAMBDA(const std::size_t j, InnerProducts<Number> &s) {
              const Number pj      = p[j];
              const Number vj      = v[j];
              const Number rj      = r[j];
              const Number prec_vj = diagonal[j] * vj;
              const Number prec_rj = diagonal[j] * rj;
              s.values[0] += pj * vj;
              s.values[1] += vj * vj;
              s.values[2] += rj * vj;
              s.values[3] += rj * rj;
              s.values[4] += rj * prec_vj;
              s.values[5] += vj * prec_vj;
              s.values[6] += rj * prec_rj;
            },
            InnerProductsReducer<Number>(sums));

          // The reduction into a scalar on the host is synchronous, so the
          // result is available here
          std::array<Number, 7> result;
          for (unsigned int i = 0; i < 7; ++i)
            result[i] = sums.values[i];
          return result;
        }
    }



    template <typename Number>
    void
    do_finalize_device(const std::size_t size,
                       const Number      alpha_plus_previous_alpha_over_beta,
                       const Number      previous_alpha_over_beta,
                       const Number     *diagonal,
                       Number           *x,
                       const Number     *r,
                       const Number     *p)
    {
      if (size == 0)
        return;

      ExecutionSpace exec;
      if (diagonal != nullptr)
        Kokkos::parallel_for(
          "dealii::cg_finalize",
          Policy(exec, 0, size),
          KOKKOS_LAMBDA(const std::size_t j) {
            x[j] += alpha_plus_previous_alpha_over_beta * p[j] -
                    previous_alpha_over_beta * diagonal[j] * r[j];
          });
      else
        Kokkos::parallel_for(
          "dealii::cg_finalize",
          Policy(exec, 0, size),
          KOKKOS_LAMBDA(const std::size_t j) {
            x[j] += alpha_plus_previous_alpha_over_beta * p[j] -
                    previous_alpha_over_beta * r[j];
          });
    }
  } // namespace SolverCG
} // namespace internal

#include "lac/solver_cg.inst"

DEAL_II_NAMESPACE_CLOSE
