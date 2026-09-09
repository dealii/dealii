// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2017 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


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
              typename ShapeDataMemorySpace,
              bool contract_over_rows,
              bool add,
              typename ViewTypeIn,
              typename ViewTypeOut>
    DEAL_II_HOST_DEVICE void
    apply_1d(const Kokkos::TeamPolicy<
               MemorySpace::Default::kokkos_space::execution_space>::member_type
                                                               &team_member,
             const Kokkos::View<Number *, ShapeDataMemorySpace> shape_data,
             const ViewTypeIn                                   in,
             ViewTypeOut                                        out)
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
              typename ShapeDataMemorySpace,
              bool contract_over_rows,
              bool add,
              typename ViewTypeIn,
              typename ViewTypeOut>
    DEAL_II_HOST_DEVICE void
    apply_2d(const Kokkos::TeamPolicy<
               MemorySpace::Default::kokkos_space::execution_space>::member_type
                                                               &team_member,
             const Kokkos::View<Number *, ShapeDataMemorySpace> shape_data,
             const ViewTypeIn                                   in,
             ViewTypeOut                                        out)
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
              typename ShapeDataMemorySpace,
              bool contract_over_rows,
              bool add,
              typename ViewTypeIn,
              typename ViewTypeOut>
    DEAL_II_HOST_DEVICE void
    apply_3d(const Kokkos::TeamPolicy<
               MemorySpace::Default::kokkos_space::execution_space>::member_type
                                                               &team_member,
             const Kokkos::View<Number *, ShapeDataMemorySpace> shape_data,
             const ViewTypeIn                                   in,
             ViewTypeOut                                        out)
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
              typename ShapeDataMemorySpace,
              int  direction,
              bool contract_over_rows,
              bool add,
              typename ViewTypeIn,
              typename ViewTypeOut>
    DEAL_II_HOST_DEVICE void
    apply(const Kokkos::TeamPolicy<
            MemorySpace::Default::kokkos_space::execution_space>::member_type
                                                            &team_member,
          const Kokkos::View<Number *, ShapeDataMemorySpace> shape_data,
          const ViewTypeIn                                   in,
          ViewTypeOut                                        out)
    {
      // We have two implementations for this apply() function. The first
      // requires Kokkos version 4.0 or later, uses the modern
      // TeamThreadMDRange, and is simpler. The second option (below)
      // performs the i,j,k loop manually.
      //
      // The second implementation turns out to be slightly faster,
      // at least until Kokkos 5.2. For now, always use the second
      // implementation. The step-104 throughput for vmult() is 33 instead of 24
      // MDoFs/s for an AMD W7800 and 135 instead of 96 MDoFs/s for an NVIDIA
      // RTX 6000.
#if 0
      if constexpr (dim == 1)
        apply_1d<n_rows,
                 n_columns,
                 direction,
                 Number,
                 ShapeDataMemorySpace,
                 contract_over_rows,
                 add>(team_member, shape_data, in, out);
      if constexpr (dim == 2)
        apply_2d<n_rows,
                 n_columns,
                 direction,
                 Number,
                 ShapeDataMemorySpace,
                 contract_over_rows,
                 add>(team_member, shape_data, in, out);
      if constexpr (dim == 3)
        apply_3d<n_rows,
                 n_columns,
                 direction,
                 Number,
                 ShapeDataMemorySpace,
                 contract_over_rows,
                 add>(team_member, shape_data, in, out);
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
              typename Number,
              typename ShapeDataMemorySpace =
                MemorySpace::Default::kokkos_space>
    struct EvaluatorTensorProduct
    {};



    /**
     * Internal evaluator for 1d-3d shape function using the tensor product form
     * of the basis functions.
     */
    template <int dim,
              int n_rows,
              int n_columns,
              typename Number,
              typename ShapeDataMemorySpace>
    struct EvaluatorTensorProduct<evaluate_general,
                                  dim,
                                  n_rows,
                                  n_columns,
                                  Number,
                                  ShapeDataMemorySpace>
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
        const TeamHandle                            &team_member,
        Kokkos::View<Number *, ShapeDataMemorySpace> shape_values,
        Kokkos::View<Number *, ShapeDataMemorySpace> shape_gradients,
        Kokkos::View<Number *, ShapeDataMemorySpace> co_shape_gradients,
        SharedView                                   temp);

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
      Kokkos::View<Number *, ShapeDataMemorySpace> shape_values;

      /**
       * Values of the shape function gradients.
       */
      Kokkos::View<Number *, ShapeDataMemorySpace> shape_gradients;

      /**
       * Values of the shape function gradients for collocation methods.
       */
      Kokkos::View<Number *, ShapeDataMemorySpace> co_shape_gradients;

      /**
       * Temporary storage for in-place evaluations.
       */
      SharedView temp;
    };



    template <int dim,
              int n_rows,
              int n_columns,
              typename Number,
              typename ShapeDataMemorySpace>
    DEAL_II_HOST_DEVICE
    EvaluatorTensorProduct<evaluate_general,
                           dim,
                           n_rows,
                           n_columns,
                           Number,
                           ShapeDataMemorySpace>::
      EvaluatorTensorProduct(
        const TeamHandle                            &team_member,
        Kokkos::View<Number *, ShapeDataMemorySpace> shape_values,
        Kokkos::View<Number *, ShapeDataMemorySpace> shape_gradients,
        Kokkos::View<Number *, ShapeDataMemorySpace> co_shape_gradients,
        SharedView                                   temp)
      : team_member(team_member)
      , shape_values(shape_values)
      , shape_gradients(shape_gradients)
      , co_shape_gradients(co_shape_gradients)
      , temp(temp)
    {}



    template <int dim,
              int n_rows,
              int n_columns,
              typename Number,
              typename ShapeDataMemorySpace>
    template <int  direction,
              bool dof_to_quad,
              bool add,
              bool in_place,
              typename ViewTypeIn,
              typename ViewTypeOut>
    DEAL_II_HOST_DEVICE void
    EvaluatorTensorProduct<evaluate_general,
                           dim,
                           n_rows,
                           n_columns,
                           Number,
                           ShapeDataMemorySpace>::values(const ViewTypeIn in,
                                                         ViewTypeOut out) const
    {
      if constexpr (in_place)
        {
          apply<dim,
                n_rows,
                n_columns,
                Number,
                ShapeDataMemorySpace,
                direction,
                dof_to_quad,
                false>(team_member, shape_values, in, temp);

          populate_view<add>(team_member, out, temp, out.extent(0));
        }
      else
        apply<dim,
              n_rows,
              n_columns,
              Number,
              ShapeDataMemorySpace,
              direction,
              dof_to_quad,
              add>(team_member, shape_values, in, out);
    }



    template <int dim,
              int n_rows,
              int n_columns,
              typename Number,
              typename ShapeDataMemorySpace>
    template <int  direction,
              bool dof_to_quad,
              bool add,
              bool in_place,
              typename ViewTypeIn,
              typename ViewTypeOut>
    DEAL_II_HOST_DEVICE void
    EvaluatorTensorProduct<evaluate_general,
                           dim,
                           n_rows,
                           n_columns,
                           Number,
                           ShapeDataMemorySpace>::gradients(const ViewTypeIn in,
                                                            ViewTypeOut out)
      const
    {
      if constexpr (in_place)
        {
          apply<dim,
                n_rows,
                n_columns,
                Number,
                ShapeDataMemorySpace,
                direction,
                dof_to_quad,
                false>(team_member, shape_gradients, in, temp);

          populate_view<add>(team_member, out, temp, out.extent(0));
        }
      else
        apply<dim,
              n_rows,
              n_columns,
              Number,
              ShapeDataMemorySpace,
              direction,
              dof_to_quad,
              add>(team_member, shape_gradients, in, out);
    }



    template <int dim,
              int n_rows,
              int n_columns,
              typename Number,
              typename ShapeDataMemorySpace>
    template <int  direction,
              bool dof_to_quad,
              bool add,
              bool in_place,
              typename ViewTypeIn,
              typename ViewTypeOut>
    DEAL_II_HOST_DEVICE void
    EvaluatorTensorProduct<
      evaluate_general,
      dim,
      n_rows,
      n_columns,
      Number,
      ShapeDataMemorySpace>::co_gradients(const ViewTypeIn in,
                                          ViewTypeOut      out) const
    {
      if constexpr (in_place)
        {
          apply<dim,
                n_columns,
                n_columns,
                Number,
                ShapeDataMemorySpace,
                direction,
                dof_to_quad,
                false>(team_member, co_shape_gradients, in, temp);

          populate_view<add>(team_member, out, temp, out.extent(0));
        }
      else
        apply<dim,
              n_columns,
              n_columns,
              Number,
              ShapeDataMemorySpace,
              direction,
              dof_to_quad,
              add>(team_member, co_shape_gradients, in, out);
    }

    namespace batched
    {
      /**
       * One-dimensional kernel for use by the generic tensor product
       * interpolation as provided by the class batched::EvaluatorTensorProduct,
       * implementing a matrix-vector product along this dimension, controlled
       * by the number of rows and columns and the stride in the input and
       * output arrays, which are embedded into some lexicographic ordering of
       * unknowns in a tensor-product arrangement.
       */
      template <int  n_rows,
                int  n_columns,
                bool contract_over_rows,
                bool add,
                int  stride_in,
                int  stride_out,
                typename Number>
      DEAL_II_HOST_DEVICE inline void
      apply_matrix_vector_product(const Number *matrix,
                                  const Number *in,
                                  Number       *out)
      {
        constexpr int mm = contract_over_rows ? n_rows : n_columns;
        constexpr int nn = contract_over_rows ? n_columns : n_rows;

        // Cache the input in registers once, rather than re-reading
        // shared memory for every output index below.
        Number r_in[mm];
        for (int k = 0; k < mm; ++k)
          r_in[k] = in[k * stride_in];

        for (int q = 0; q < nn; ++q)
          {
            Number sum = 0;
            for (int k = 0; k < mm; ++k)
              {
                const int row = contract_over_rows ? k : q;
                const int col = contract_over_rows ? q : k;
                sum += matrix[row * n_columns + col] * r_in[k];
              }

            if constexpr (add)
              out[q * stride_out] += sum;
            else
              out[q * stride_out] = sum;
          }
      }

      /**
       * Specialized version of apply_matrix_vector_product() that takes the
       * strides as arguments, rather than as template parameters.
       */
      template <int  n_rows,
                int  n_columns,
                bool contract_over_rows,
                bool add,
                typename Number>
      DEAL_II_HOST_DEVICE inline void
      apply_matrix_vector_product(const Number *matrix,
                                  const Number *in,
                                  Number       *out,
                                  const int     stride_in,
                                  const int     stride_out)
      {
        constexpr int mm = contract_over_rows ? n_rows : n_columns;
        constexpr int nn = contract_over_rows ? n_columns : n_rows;

        // Cache the input in registers once, rather than re-reading
        // shared memory for every output index below.
        Number r_in[mm];
        for (int k = 0; k < mm; ++k)
          r_in[k] = in[k * stride_in];

        for (int q = 0; q < nn; ++q)
          {
            Number sum = 0;
            for (int k = 0; k < mm; ++k)
              {
                const int row = contract_over_rows ? k : q;
                const int col = contract_over_rows ? q : k;
                sum += matrix[row * n_columns + col] * r_in[k];
              }

            if constexpr (add)
              out[q * stride_out] += sum;
            else
              out[q * stride_out] = sum;
          }
      }


      /**
       * Kokkos::View-based overload of apply_matrix_vector_product()
       * above with `matrix`/`in`/`out` passed as Kokkos::Views.
       */
      template <
        int  n_rows,
        int  n_columns,
        bool contract_over_rows,
        bool add,
        int  stride_in,
        int  stride_out,
        typename ViewTypeMatrix,
        typename ViewTypeIn,
        typename ViewTypeOut,
        typename = std::enable_if_t<Kokkos::is_view<ViewTypeOut>::value>>
      DEAL_II_HOST_DEVICE inline void
      apply_matrix_vector_product(const ViewTypeMatrix matrix,
                                  const ViewTypeIn     in,
                                  ViewTypeOut          out)
      {
        using Number = typename ViewTypeOut::non_const_value_type;

        constexpr int mm = contract_over_rows ? n_rows : n_columns;
        constexpr int nn = contract_over_rows ? n_columns : n_rows;

        // Cache the input in registers once, rather than re-reading
        // shared memory for every output index below.
        Number r_in[mm];
        for (int k = 0; k < mm; ++k)
          r_in[k] = in(k * stride_in);

        for (int q = 0; q < nn; ++q)
          {
            Number sum = 0;
            for (int k = 0; k < mm; ++k)
              {
                const int row = contract_over_rows ? k : q;
                const int col = contract_over_rows ? q : k;
                sum += matrix(row * n_columns + col) * r_in[k];
              }

            if constexpr (add)
              out(q * stride_out) += sum;
            else
              out(q * stride_out) = sum;
          }
      }

      /**
       * View-based overload of the runtime-stride apply_matrix_vector_product()
       * above.
       */
      template <
        int  n_rows,
        int  n_columns,
        bool contract_over_rows,
        bool add,
        typename ViewTypeMatrix,
        typename ViewTypeIn,
        typename ViewTypeOut,
        typename = std::enable_if_t<Kokkos::is_view<ViewTypeOut>::value>>
      DEAL_II_HOST_DEVICE inline void
      apply_matrix_vector_product(const ViewTypeMatrix matrix,
                                  const ViewTypeIn     in,
                                  ViewTypeOut          out,
                                  const int            stride_in,
                                  const int            stride_out)
      {
        using Number = typename ViewTypeOut::non_const_value_type;

        constexpr int mm = contract_over_rows ? n_rows : n_columns;
        constexpr int nn = contract_over_rows ? n_columns : n_rows;

        // Cache the input in registers once, rather than re-reading
        // shared memory for every output index below.
        Number r_in[mm];
        for (int k = 0; k < mm; ++k)
          r_in[k] = in(k * stride_in);

        for (int q = 0; q < nn; ++q)
          {
            Number sum = 0;
            for (int k = 0; k < mm; ++k)
              {
                const int row = contract_over_rows ? k : q;
                const int col = contract_over_rows ? q : k;
                sum += matrix(row * n_columns + col) * r_in[k];
              }

            if constexpr (add)
              out(q * stride_out) += sum;
            else
              out(q * stride_out) = sum;
          }
      }

      /**
       * Helper function that applies sum factorization in a specified direction
       * using batched kernel and apply_matrix_vector_product().
       *
       * Sizes of the input and output vectors in 1D:
       * -----------------------------------------------------------
       *   direction  |  contract_over_rows  |  !contract_over_rows
       * -----------------------------------------------------------
       *       0      |        m -> n        |       n   -> m
       * -----------------------------------------------------------
       *
       * Sizes of the input and output vectors in 2D:
       * -----------------------------------------------------------
       *   direction  |  contract_over_rows  |  !contract_over_rows
       * -----------------------------------------------------------
       *       0      |    m x m -> n x m    |     n x m -> m x m
       * -----------------------------------------------------------
       *       1      |    n x m -> n x n    |     n x n -> n x m
       * -----------------------------------------------------------
       *
       * Sizes of the input and output vectors in 3D:
       * ------------------------------------------------------------------
       *   direction  |    contract_over_rows    |   !contract_over_rows
       * ------------------------------------------------------------------
       *       0      |  m x m x m -> n x m x m  |  n x m x m -> m x m x m
       * ------------------------------------------------------------------
       *       1      |  n x m x m -> n x n x m  |  n x n x m -> n x m x m
       * ------------------------------------------------------------------
       *       2      |  n x n x m -> n x n x n  |  n x n x n -> n x n x m
       * ------------------------------------------------------------------
       */
      template <
        int  dim,
        int  direction,
        int  n_rows,
        int  n_columns,
        bool contract_over_rows,
        bool add,
        typename ViewTypeMatrix,
        typename ViewTypeIn,
        typename ViewTypeOut,
        typename = std::enable_if_t<Kokkos::is_view<ViewTypeOut>::value>>
      DEAL_II_HOST_DEVICE inline void
      apply(const Kokkos::TeamPolicy<
              MemorySpace::Default::kokkos_space::execution_space>::member_type
                                &team_member,
            const ViewTypeMatrix shape_data,
            const ViewTypeIn     in,
            ViewTypeOut          out,
            const int            batch_size = 1)
      {
        static_assert(direction >= 0 && direction < dim,
                      "direction must be in [0, dim)");
        Assert(shape_data.size() == n_rows * n_columns, ExcInternalError());

        constexpr int mm = contract_over_rows ? n_rows : n_columns;
        constexpr int nn = contract_over_rows ? n_columns : n_rows;

        // combined extent of the already-transformed axes
        // (role < direction, extent n_columns each)
        constexpr int n_blocks1 = Utilities::pow(n_columns, direction);

        // n_blocks2: combined extent ofthe not-yet-transformed
        // axes (role > direction, extent n_rows  each)
        constexpr int n_blocks2 = Utilities::pow(n_rows, dim - direction - 1);

        constexpr int n_in_per_elmt  = n_blocks1 * mm * n_blocks2;
        constexpr int n_out_per_elmt = n_blocks1 * nn * n_blocks2;

        // Unlike apply_1d()/apply_2d()/apply_3d() above (one cell per call),
        // this is the batched variant -- batch_size cells'
        // worth of data, laid out contiguously per element, share one in/out
        // view, so the required size scales with the batch size too.
        Assert(in.size() >=
                 static_cast<std::size_t>(batch_size) * n_in_per_elmt,
               ExcInternalError());
        Assert(out.size() >=
                 static_cast<std::size_t>(batch_size) * n_out_per_elmt,
               ExcInternalError());

        Kokkos::parallel_for(
          Kokkos::TeamVectorRange(team_member,
                                  batch_size * n_blocks1 * n_blocks2),
          [&](const int tid) {
            const int e   = tid / (n_blocks1 * n_blocks2);
            const int rem = tid % (n_blocks1 * n_blocks2);
            const int i2  = rem / n_blocks1;
            const int i1  = rem % n_blocks1;

            const int in_offset = e * n_in_per_elmt + i2 * n_blocks1 * mm + i1;
            const int out_offset =
              e * n_out_per_elmt + i2 * n_blocks1 * nn + i1;

            apply_matrix_vector_product<n_rows,
                                        n_columns,
                                        contract_over_rows,
                                        add,
                                        n_blocks1,
                                        n_blocks1>(
              shape_data,
              Kokkos::subview(
                in,
                Kokkos::make_pair(in_offset, static_cast<int>(in.extent(0)))),
              Kokkos::subview(out,
                              Kokkos::make_pair(
                                out_offset, static_cast<int>(out.extent(0)))));
          });

        team_member.team_barrier();
      }

      /**
       * Helper function that copies or adds the first N entries of src to
       * dst, depending on the template argument "add".
       */
      template <
        bool add,
        typename ViewTypeOut,
        typename ViewTypeIn,
        typename = std::enable_if_t<Kokkos::is_view<ViewTypeOut>::value>>
      DEAL_II_HOST_DEVICE inline void
      populate_view(
        const Kokkos::TeamPolicy<
          MemorySpace::Default::kokkos_space::execution_space>::member_type
                        &team_member,
        ViewTypeOut      dst,
        const ViewTypeIn src,
        const int        N)
      {
        Kokkos::parallel_for(Kokkos::TeamVectorRange(team_member, N),
                             [&](const int tid) {
                               if constexpr (add)
                                 dst(tid) += src(tid);
                               else
                                 dst(tid) = src(tid);
                             });

        team_member.team_barrier();
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
       * Internal evaluator for 1d-3d shape function using the tensor product
       * form of the basis functions.
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

        using ShapeDataType =
          Kokkos::View<Number *,
                       MemorySpace::Default::kokkos_space::execution_space::
                         scratch_memory_space,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

        using SharedView =
          Kokkos::View<Number *,
                       MemorySpace::Default::kokkos_space::execution_space::
                         scratch_memory_space,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>;


        DEAL_II_HOST_DEVICE
        EvaluatorTensorProduct(const TeamHandle &team_member,
                               ShapeDataType     shape_values,
                               ShapeDataType     shape_gradients,
                               ShapeDataType     co_shape_gradients,
                               SharedView        temp,
                               const int         batch_size);

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
         * Evaluate (transpose = false) or integrate (transpose = true) the
         * full gradient (all `dim` components at once) in the
         * collocation space. It fuses what would otherwise be `dim` separate
         * direction-by-direction co_gradients() calls into a single pass
         * with one register-cached read of `in` per output quadrature point.
         * transpose = false broadcasts one scalar field into `dim` gradient
         * components (the forward, evaluate_gradients() case);
         * transpose = true is the adjoint. i.e.,  it reduces `dim` gradient
         * components back into one scalar field (the integrate_gradients()
         * case). This is faster alternative to dim consecutive exuction of the
         * above function.
         */
        template <bool transpose,
                  bool add,
                  typename ViewTypeIn,
                  typename ViewTypeOut>
        DEAL_II_HOST_DEVICE void
        co_gradients(const ViewTypeIn in, ViewTypeOut out) const;

      private:
        const TeamHandle &team_member;
        ShapeDataType     shape_values;
        ShapeDataType     shape_gradients;
        ShapeDataType     co_shape_gradients;
        SharedView        temp;
        const int         batch_size;
      };

      template <int dim, int n_rows, int n_columns, typename Number>
      DEAL_II_HOST_DEVICE
      EvaluatorTensorProduct<evaluate_general, dim, n_rows, n_columns, Number>::
        EvaluatorTensorProduct(const TeamHandle &team_member,
                               ShapeDataType     shape_values,
                               ShapeDataType     shape_gradients,
                               ShapeDataType     co_shape_gradients,
                               SharedView        temp,
                               const int         batch_size)
        : team_member(team_member)
        , shape_values(shape_values)
        , shape_gradients(shape_gradients)
        , co_shape_gradients(co_shape_gradients)
        , temp(temp)
        , batch_size(batch_size)
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
            apply<dim, direction, n_rows, n_columns, dof_to_quad, false>(
              team_member, shape_values, in, temp, batch_size);

            constexpr int nn        = dof_to_quad ? n_columns : n_rows;
            constexpr int n_blocks1 = Utilities::pow(n_columns, direction);
            constexpr int n_blocks2 =
              Utilities::pow(n_rows, dim - direction - 1);

            populate_view<add>(team_member,
                               out,
                               temp,
                               batch_size * n_blocks1 * nn * n_blocks2);
          }
        else
          {
            apply<dim, direction, n_rows, n_columns, dof_to_quad, add>(
              team_member, shape_values, in, out, batch_size);
          }
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
            apply<dim, direction, n_rows, n_columns, dof_to_quad, false>(
              team_member, shape_gradients, in, temp, batch_size);

            constexpr int nn        = dof_to_quad ? n_columns : n_rows;
            constexpr int n_blocks1 = Utilities::pow(n_columns, direction);
            constexpr int n_blocks2 =
              Utilities::pow(n_rows, dim - direction - 1);

            populate_view<add>(team_member,
                               out,
                               temp,
                               batch_size * n_blocks1 * nn * n_blocks2);
          }
        else
          {
            apply<dim, direction, n_rows, n_columns, dof_to_quad, add>(
              team_member, shape_gradients, in, out, batch_size);
          }
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
            apply<dim, direction, n_columns, n_columns, dof_to_quad, false>(
              team_member, co_shape_gradients, in, temp, batch_size);

            constexpr int n_blocks1 = Utilities::pow(n_columns, direction);
            constexpr int n_blocks2 =
              Utilities::pow(n_columns, dim - direction - 1);

            populate_view<add>(team_member,
                               out,
                               temp,
                               batch_size * n_blocks1 * n_columns * n_blocks2);
          }
        else
          {
            apply<dim, direction, n_columns, n_columns, dof_to_quad, add>(
              team_member, co_shape_gradients, in, out, batch_size);
          }
      }

      template <int dim, int n_rows, int n_columns, typename Number>
      template <bool transpose,
                bool add,
                typename ViewTypeIn,
                typename ViewTypeOut>
      DEAL_II_HOST_DEVICE void
      EvaluatorTensorProduct<evaluate_general, dim, n_rows, n_columns, Number>::
        co_gradients(const ViewTypeIn in, ViewTypeOut out) const
      {
        static_assert(dim >= 1, "dim must be at least 1");
        static_assert(
          ViewTypeIn::rank == (transpose ? 2 : 1),
          "in must be the values (1D view) for evaluate_gradients (transpose = "
          "false), or the gradients (2D view) for integrate_gradients (transpose "
          "= true).");
        static_assert(
          ViewTypeOut::rank == (transpose ? 1 : 2),
          "out must be the gradients (2D view) for evaluate_gradients (transpose "
          "= false), or the values (1D view) for integrate_gradients (transpose "
          "= true).");

        constexpr int n_q_points        = Utilities::pow(n_columns, dim);
        constexpr int co_dimension_size = Utilities::pow(n_columns, dim - 1);

        Kokkos::parallel_for(
          Kokkos::TeamVectorRange(team_member, batch_size * co_dimension_size),
          [&](const int tid) {
            const int elmnt_idx = tid / co_dimension_size;
            const int reminder  = tid % co_dimension_size;

            Kokkos::Array<int, dim - 1> idx_d, stride_d;
            Number                      reg[dim][n_columns];

            for (int d = 0; d < dim - 1; ++d)
              {
                stride_d[d] = Utilities::pow(n_columns, d);
                idx_d[d]    = (reminder / stride_d[d]) % n_columns;
                for (int n = 0; n < n_columns; ++n)
                  {
                    if constexpr (!transpose)
                      reg[d][n] = co_shape_gradients(n * n_columns + idx_d[d]);
                    else
                      reg[d][n] = co_shape_gradients(idx_d[d] * n_columns + n);
                  }
              }

            for (int n = 0; n < n_columns; ++n)
              {
                if constexpr (!transpose)
                  reg[dim - 1][n] = in(elmnt_idx * n_q_points + reminder +
                                       n * co_dimension_size);
                else
                  reg[dim - 1][n] = in(elmnt_idx * n_q_points + reminder +
                                         n * co_dimension_size,
                                       dim - 1);
              }

            for (int last = 0; last < n_columns; ++last)
              {
                const int q_point = reminder + last * co_dimension_size;

                if constexpr (!transpose)
                  {
                    Number result[dim];
                    for (int d = 0; d < dim - 1; ++d)
                      {
                        const int q_point_base =
                          q_point - idx_d[d] * stride_d[d];
                        const int in_base =
                          elmnt_idx * n_q_points + q_point_base;

                        Number res_d = 0;
                        for (int n = 0; n < n_columns; ++n)
                          res_d += reg[d][n] * in(in_base + n * stride_d[d]);
                        result[d] = res_d;
                      }
                    {
                      Number res_d = 0;
                      for (int n = 0; n < n_columns; ++n)
                        res_d += co_shape_gradients(n * n_columns + last) *
                                 reg[dim - 1][n];
                      result[dim - 1] = res_d;
                    }

                    for (int d = 0; d < dim; ++d)
                      {
                        if constexpr (add)
                          out(elmnt_idx * n_q_points + q_point, d) += result[d];
                        else
                          out(elmnt_idx * n_q_points + q_point, d) = result[d];
                      }
                  }
                else
                  {
                    Number result = 0;
                    for (int d = 0; d < dim - 1; ++d)
                      {
                        const int point_base = q_point - idx_d[d] * stride_d[d];
                        const int grad_row =
                          elmnt_idx * n_q_points + point_base;

                        for (int n = 0; n < n_columns; ++n)
                          result +=
                            in(grad_row + n * stride_d[d], d) * reg[d][n];
                      }
                    for (int n = 0; n < n_columns; ++n)
                      result += reg[dim - 1][n] *
                                co_shape_gradients(last * n_columns + n);

                    if constexpr (add)
                      out(elmnt_idx * n_q_points + q_point) += result;
                    else
                      out(elmnt_idx * n_q_points + q_point) = result;
                  }
              }
          });

        team_member.team_barrier();
      }


    } // namespace batched
  }   // namespace internal
} // namespace Portable

DEAL_II_NAMESPACE_CLOSE

#endif
