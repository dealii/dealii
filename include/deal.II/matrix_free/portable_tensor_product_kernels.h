// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#ifndef dealii__tensor_product_kernels_h
#define dealii__tensor_product_kernels_h

#include <deal.II/base/config.h>

#include <deal.II/base/memory_space.h>
#include <deal.II/base/utilities.h>

#include <Kokkos_Core.hpp>


DEAL_II_NAMESPACE_OPEN


namespace Portable
{
  namespace internal
  {
    /**
     * In this namespace, the evaluator routines that evaluate the tensor
     * products are implemented.
     */
    // TODO: for now only the general variant is implemented
    enum EvaluatorVariant
    {
      evaluate_general,
      evaluate_symmetric,
      evaluate_evenodd
    };



    /**
     * Helper function that copies or adds the first N entries of src to
     * dst, depending on the template argument "add".
     */
    template <bool add, typename ViewTypeIn, typename ViewTypeOut>
    DEAL_II_HOST_DEVICE void
    populate_view(
      const Kokkos::TeamPolicy<
        MemorySpace::Default::kokkos_space::execution_space>::member_type
                      &team_member,
      ViewTypeOut      dst,
      const ViewTypeIn src,
      const int        N)
    {
      Assert(dst.size() >= static_cast<unsigned int>(N), ExcInternalError());
      Assert(src.size() >= static_cast<unsigned int>(N), ExcInternalError());
      Kokkos::parallel_for(Kokkos::TeamVectorRange(team_member, N),
                           [&](const int i) {
                             if constexpr (add)
                               Kokkos::atomic_add(&dst(i), src(i));
                             else
                               dst(i) = src(i);
                           });

      team_member.team_barrier();
    }



#if DEAL_II_KOKKOS_VERSION_GTE(4, 0, 0)
    /**
     * Helper function for apply() in 1D
     */
    template <int n_rows,
              int n_columns,
              int direction,
              typename Number,
              bool contract_over_rows,
              bool add,
              typename ViewTypeIn,
              typename ViewTypeOut>
    DEAL_II_HOST_DEVICE void
    apply_1d(const Kokkos::TeamPolicy<
               MemorySpace::Default::kokkos_space::execution_space>::member_type
               &team_member,
             const Kokkos::View<Number *, MemorySpace::Default::kokkos_space>
                              shape_data,
             const ViewTypeIn in,
             ViewTypeOut      out)
    {
      constexpr int Nk = (contract_over_rows ? n_rows : n_columns),
                    Nq = (contract_over_rows ? n_columns : n_rows);

      Assert(shape_data.size() == n_rows * n_columns, ExcInternalError());
      Assert(in.size() >= Nk, ExcInternalError());
      Assert(out.size() >= Nq, ExcInternalError());

      Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, Nq),
                           [&](const int q) {
                             Number sum = 0;
                             for (int k = 0; k < Nk; ++k)
                               {
                                 const int shape_idx =
                                   (contract_over_rows ? q + k * Nq :
                                                         k + q * Nk);
                                 sum += shape_data(shape_idx) * in(k);
                               }

                             if constexpr (add)
                               Kokkos::atomic_add(&out(q), sum);
                             else
                               out(q) = sum;
                           });

      team_member.team_barrier();
    }



    /**
     * Helper function for apply() in 2D
     */
    template <int n_rows,
              int n_columns,
              int direction,
              typename Number,
              bool contract_over_rows,
              bool add,
              typename ViewTypeIn,
              typename ViewTypeOut>
    DEAL_II_HOST_DEVICE void
    apply_2d(const Kokkos::TeamPolicy<
               MemorySpace::Default::kokkos_space::execution_space>::member_type
               &team_member,
             const Kokkos::View<Number *, MemorySpace::Default::kokkos_space>
                              shape_data,
             const ViewTypeIn in,
             ViewTypeOut      out)
    {
      using TeamType = Kokkos::TeamPolicy<
        MemorySpace::Default::kokkos_space::execution_space>::member_type;

      // Sizes of the input and output vectors:
      // -----------------------------------------------------------
      //   direction  |  contract_over_rows  |  !contract_over_rows
      // -----------------------------------------------------------
      //       0      |    m x m -> n x m    |     n x m -> m x m
      // -----------------------------------------------------------
      //       1      |    n x m -> n x n    |     n x n -> n x m
      // -----------------------------------------------------------
      //
      // Directions of the cycle indices:
      // -----------------------------
      //   direction  |   j   |  q/k
      // -----------------------------
      //       0      |   1   |   0
      // -----------------------------
      //       1      |   0   |   1
      // -----------------------------
      constexpr int Nj = (direction < 1 ? n_rows : n_columns),
                    Nk = (contract_over_rows ? n_rows : n_columns),
                    Nq = (contract_over_rows ? n_columns : n_rows);

      Assert(shape_data.size() == n_rows * n_columns, ExcInternalError());
      Assert(in.size() >= Nj * Nk, ExcInternalError());
      Assert(out.size() >= Nj * Nq, ExcInternalError());

      auto thread_policy =
        Kokkos::TeamThreadMDRange<Kokkos::Rank<2>, TeamType>(team_member,
                                                             Nj,
                                                             Nq);
      Kokkos::parallel_for(thread_policy, [&](const int j, const int q) {
        const int base_shape   = contract_over_rows ? q : q * n_columns;
        const int stride_shape = contract_over_rows ? n_columns : 1;

        const int base_in   = (direction == 0 ? j * Nk : j);
        const int stride_in = Utilities::pow(n_columns, direction);

        Number sum = shape_data(base_shape) * in(base_in);
        for (int k = 1; k < Nk; ++k)
          sum += shape_data(base_shape + k * stride_shape) *
                 in(base_in + k * stride_in);

        const int index_out = (direction == 0 ? j * Nq + q : j + q * Nj);

        if constexpr (add)
          Kokkos::atomic_add(&out(index_out), sum);
        else
          out(index_out) = sum;
      });

      team_member.team_barrier();
    }



    /**
     * Helper function for apply() in 3D
     */
    template <int n_rows,
              int n_columns,
              int direction,
              typename Number,
              bool contract_over_rows,
              bool add,
              typename ViewTypeIn,
              typename ViewTypeOut>
    DEAL_II_HOST_DEVICE void
    apply_3d(const Kokkos::TeamPolicy<
               MemorySpace::Default::kokkos_space::execution_space>::member_type
               &team_member,
             const Kokkos::View<Number *, MemorySpace::Default::kokkos_space>
                              shape_data,
             const ViewTypeIn in,
             ViewTypeOut      out)
    {
      using TeamType = Kokkos::TeamPolicy<
        MemorySpace::Default::kokkos_space::execution_space>::member_type;

      // Sizes of the input and output vectors:
      // ------------------------------------------------------------------
      //   direction  |    contract_over_rows    |   !contract_over_rows
      // ------------------------------------------------------------------
      //       0      |  m x m x m -> n x m x m  |  n x m x m -> m x m x m
      // ------------------------------------------------------------------
      //       1      |  n x m x m -> n x n x m  |  n x n x m -> n x m x m
      // ------------------------------------------------------------------
      //       2      |  n x n x m -> n x n x n  |  n x n x n -> n x n x m
      // ------------------------------------------------------------------
      //
      // Directions of the cycle indices:
      // -------------------------------------
      //   direction  |   i   |   j   |  q/k
      // -------------------------------------
      //       0      |   2   |   1   |   0
      // -------------------------------------
      //       1      |   0   |   2   |   1
      // -------------------------------------
      //       2      |   1   |   0   |   2
      // -------------------------------------
      constexpr int Ni = (direction < 1 ? n_rows : n_columns),
                    Nj = (direction < 2 ? n_rows : n_columns),
                    Nk = (contract_over_rows ? n_rows : n_columns),
                    Nq = (contract_over_rows ? n_columns : n_rows);

      Assert(shape_data.size() == n_rows * n_columns, ExcInternalError());
      Assert(in.size() >= Ni * Nj * Nk, ExcInternalError());
      Assert(out.size() >= Ni * Nj * Nq, ExcInternalError());

      auto thread_policy = Kokkos::TeamThreadMDRange<Kokkos::Rank<3>, TeamType>(
        team_member, Ni, Nj, Nq);
      Kokkos::parallel_for(
        thread_policy, [&](const int i, const int j, const int q) {
          const int base_shape   = contract_over_rows ? q : q * n_columns;
          const int stride_shape = contract_over_rows ? n_columns : 1;

          const int base_in =
            (direction == 0 ? (i * Nj + j) * Nk :
                              (direction == 1 ? i + j * Ni * Nk : i * Nj + j));
          const int stride_in = Utilities::pow(n_columns, direction);

          Number sum = shape_data(base_shape) * in(base_in);
          for (int k = 1; k < Nk; ++k)
            sum += shape_data(base_shape + k * stride_shape) *
                   in(base_in + k * stride_in);

          const int index_out =
            (direction == 0 ? (i * Nj + j) * Nq + q :
                              (direction == 1 ? i + (j * Nq + q) * Ni :
                                                (i + q * Ni) * Nj + j));

          if constexpr (add)
            Kokkos::atomic_add(&out(index_out), sum);
          else
            out(index_out) = sum;
        });

      team_member.team_barrier();
    }
#endif



    template <int dim,
              int n_rows,
              int n_columns,
              typename Number,
              int  direction,
              bool contract_over_rows,
              bool add,
              typename ViewTypeIn,
              typename ViewTypeOut>
    DEAL_II_HOST_DEVICE void
    apply(const Kokkos::TeamPolicy<
            MemorySpace::Default::kokkos_space::execution_space>::member_type
            &team_member,
          const Kokkos::View<Number *, MemorySpace::Default::kokkos_space>
                           shape_data,
          const ViewTypeIn in,
          ViewTypeOut      out)
    {
#if DEAL_II_KOKKOS_VERSION_GTE(4, 0, 0)
      if constexpr (dim == 1)
        apply_1d<n_rows, n_columns, direction, Number, contract_over_rows, add>(
          team_member, shape_data, in, out);
      if constexpr (dim == 2)
        apply_2d<n_rows, n_columns, direction, Number, contract_over_rows, add>(
          team_member, shape_data, in, out);
      if constexpr (dim == 3)
        apply_3d<n_rows, n_columns, direction, Number, contract_over_rows, add>(
          team_member, shape_data, in, out);
#else
      // I: [0, m^{dim - direction - 1})
      // J: [0, n^direction)
      constexpr int NI = Utilities::pow(n_rows, dim - direction - 1);
      constexpr int NJ = Utilities::pow(n_columns, direction);

      constexpr int Nk = contract_over_rows ? n_rows : n_columns;
      constexpr int Nq = contract_over_rows ? n_columns : n_rows;

      Assert(shape_data.size() == n_rows * n_columns, ExcInternalError());
      Assert(in.size() >= NI * NJ * Nk, ExcInternalError());
      Assert(out.size() >= NI * NJ * Nq, ExcInternalError());

      constexpr int N      = NI * NJ * Nq;
      constexpr int stride = Utilities::pow(n_columns, direction);

      Kokkos::parallel_for(
        Kokkos::TeamThreadRange(team_member, N), [&](const int index_out) {
          // index_in  = (I Nk + k) n^direction + J
          // index_out = (I Nq + q) n^direction + J
          const int q = (index_out / stride) % Nq;
          const int I = (index_out / stride) / Nq;
          const int J = index_out % stride;

          const int base_shape   = contract_over_rows ? q : q * n_columns;
          const int stride_shape = contract_over_rows ? n_columns : 1;
          const int base_in      = I * Nk * stride + J;

          Number sum = shape_data(base_shape) * in(base_in);
          for (int k = 1; k < Nk; ++k)
            {
              const int index_in = (I * Nk + k) * stride + J;
              sum += shape_data(base_shape + k * stride_shape) * in(index_in);
            }

          if constexpr (add)
            Kokkos::atomic_add(&out(index_out), sum);
          else
            out(index_out) = sum;
        });

      team_member.team_barrier();
#endif
    }



    /**
     * Generic evaluator framework.
     */
    template <EvaluatorVariant variant,
              int              dim,
              int              n_rows,
              int              n_columns,
              typename Number>
    struct EvaluatorTensorProduct
    {};



    /**
     * Internal evaluator for 1d-3d shape function using the tensor product form
     * of the basis functions.
     */
    template <int dim, int n_rows, int n_columns, typename Number>
    struct EvaluatorTensorProduct<evaluate_general,
                                  dim,
                                  n_rows,
                                  n_columns,
                                  Number>
    {
    public:
      using TeamHandle = Kokkos::TeamPolicy<
        MemorySpace::Default::kokkos_space::execution_space>::member_type;

      using SharedView = Kokkos::View<Number *,
                                      MemorySpace::Default::kokkos_space::
                                        execution_space::scratch_memory_space,
                                      Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

      DEAL_II_HOST_DEVICE
      EvaluatorTensorProduct(
        const TeamHandle                                          &team_member,
        Kokkos::View<Number *, MemorySpace::Default::kokkos_space> shape_values,
        Kokkos::View<Number *, MemorySpace::Default::kokkos_space>
          shape_gradients,
        Kokkos::View<Number *, MemorySpace::Default::kokkos_space>
                   co_shape_gradients,
        SharedView temp);

      /**
       * Evaluate/integrate the values of a finite element function at the
       * quadrature points for a given @p direction.
       */
      template <int  direction,
                bool dof_to_quad,
                bool add,
                bool in_place,
                typename ViewTypeIn,
                typename ViewTypeOut>
      DEAL_II_HOST_DEVICE void
      values(const ViewTypeIn in, ViewTypeOut out) const;

      /**
       * Evaluate/integrate the gradient of a finite element function at the
       * quadrature points for a given @p direction.
       */
      template <int  direction,
                bool dof_to_quad,
                bool add,
                bool in_place,
                typename ViewTypeIn,
                typename ViewTypeOut>
      DEAL_II_HOST_DEVICE void
      gradients(const ViewTypeIn in, ViewTypeOut out) const;

      /**
       * Evaluate the gradient of a finite element function at the quadrature
       * points for a given @p direction for collocation methods.
       */
      template <int  direction,
                bool dof_to_quad,
                bool add,
                bool in_place,
                typename ViewTypeIn,
                typename ViewTypeOut>
      DEAL_II_HOST_DEVICE void
      co_gradients(const ViewTypeIn in, ViewTypeOut out) const;

      /**
       * TeamPolicy handle.
       */
      const TeamHandle &team_member;

      /**
       * Values of the shape functions.
       */
      Kokkos::View<Number *, MemorySpace::Default::kokkos_space> shape_values;

      /**
       * Values of the shape function gradients.
       */
      Kokkos::View<Number *, MemorySpace::Default::kokkos_space>
        shape_gradients;

      /**
       * Values of the shape function gradients for collocation methods.
       */
      Kokkos::View<Number *, MemorySpace::Default::kokkos_space>
        co_shape_gradients;

      /**
       * Temporary storage for in-place evaluations.
       */
      SharedView temp;
    };



    template <int dim, int n_rows, int n_columns, typename Number>
    DEAL_II_HOST_DEVICE
    EvaluatorTensorProduct<evaluate_general, dim, n_rows, n_columns, Number>::
      EvaluatorTensorProduct(
        const TeamHandle                                          &team_member,
        Kokkos::View<Number *, MemorySpace::Default::kokkos_space> shape_values,
        Kokkos::View<Number *, MemorySpace::Default::kokkos_space>
          shape_gradients,
        Kokkos::View<Number *, MemorySpace::Default::kokkos_space>
                   co_shape_gradients,
        SharedView temp)
      : team_member(team_member)
      , shape_values(shape_values)
      , shape_gradients(shape_gradients)
      , co_shape_gradients(co_shape_gradients)
      , temp(temp)
    {}



    template <int dim, int n_rows, int n_columns, typename Number>
    template <int  direction,
              bool dof_to_quad,
              bool add,
              bool in_place,
              typename ViewTypeIn,
              typename ViewTypeOut>
    DEAL_II_HOST_DEVICE void
    EvaluatorTensorProduct<evaluate_general, dim, n_rows, n_columns, Number>::
      values(const ViewTypeIn in, ViewTypeOut out) const
    {
      if constexpr (in_place)
        {
          apply<dim, n_rows, n_columns, Number, direction, dof_to_quad, false>(
            team_member, shape_values, in, temp);

          populate_view<add>(team_member, out, temp, out.extent(0));
        }
      else
        apply<dim, n_rows, n_columns, Number, direction, dof_to_quad, add>(
          team_member, shape_values, in, out);
    }



    template <int dim, int n_rows, int n_columns, typename Number>
    template <int  direction,
              bool dof_to_quad,
              bool add,
              bool in_place,
              typename ViewTypeIn,
              typename ViewTypeOut>
    DEAL_II_HOST_DEVICE void
    EvaluatorTensorProduct<evaluate_general, dim, n_rows, n_columns, Number>::
      gradients(const ViewTypeIn in, ViewTypeOut out) const
    {
      if constexpr (in_place)
        {
          apply<dim, n_rows, n_columns, Number, direction, dof_to_quad, false>(
            team_member, shape_gradients, in, temp);

          populate_view<add>(team_member, out, temp, out.extent(0));
        }
      else
        apply<dim, n_rows, n_columns, Number, direction, dof_to_quad, add>(
          team_member, shape_gradients, in, out);
    }



    template <int dim, int n_rows, int n_columns, typename Number>
    template <int  direction,
              bool dof_to_quad,
              bool add,
              bool in_place,
              typename ViewTypeIn,
              typename ViewTypeOut>
    DEAL_II_HOST_DEVICE void
    EvaluatorTensorProduct<evaluate_general, dim, n_rows, n_columns, Number>::
      co_gradients(const ViewTypeIn in, ViewTypeOut out) const
    {
      if constexpr (in_place)
        {
          apply<dim,
                n_columns,
                n_columns,
                Number,
                direction,
                dof_to_quad,
                false>(team_member, co_shape_gradients, in, temp);

          populate_view<add>(team_member, out, temp, out.extent(0));
        }
      else
        apply<dim, n_columns, n_columns, Number, direction, dof_to_quad, add>(
          team_member, co_shape_gradients, in, out);
    }
  } // namespace internal
} // namespace Portable

DEAL_II_NAMESPACE_CLOSE

#endif
