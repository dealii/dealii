// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_portable_mg_global_coarsening_templates_h
#define dealii_portable_mg_global_coarsening_templates_h

#include <deal.II/base/config.h>

#include <deal.II/base/memory_space.h>
#include <deal.II/base/mpi_consensus_algorithms.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/matrix_free/hanging_nodes_internal.h>
#include <deal.II/matrix_free/portable_hanging_nodes_internal.h>
#include <deal.II/matrix_free/portable_tensor_product_kernels.h>

#include <deal.II/multigrid/mg_transfer_matrix_free.h>
#include <deal.II/multigrid/mg_transfer_matrix_free.templates.h>
#include <deal.II/multigrid/portable_mg_transfer_global_coarsening.h>

DEAL_II_NAMESPACE_OPEN

namespace Portable
{
  namespace internal
  {

    /**
     * Helper class to select the right templated implementation.
     *
     * @note This class is similar to internal::FEEvaluationFactory
     */
    class CellTransferFactory
    {
    public:
      static const unsigned int max_degree = 4;

      CellTransferFactory(const unsigned int degree_fine,
                          const unsigned int degree_coarse)
        : degree_fine(degree_fine)
        , degree_coarse(degree_coarse)
      {}

      template <typename KernelType>
      DEAL_II_HOST_DEVICE bool
      run(KernelType &kernel) const
      {
        // h-MG can use a faster path
        if (degree_fine == 2 && degree_coarse == 1)
          {
            kernel.template run<2, 1>();
            return true; // fast path taken
          }
        if (degree_fine == 4 && degree_coarse == 2)
          {
            kernel.template run<4, 2>();
            return true; // fast path taken
          }
        if (degree_fine == 6 && degree_coarse == 3)
          {
            kernel.template run<6, 3>();
            return true; // fast path taken
          }
        if (degree_fine == 8 && degree_coarse == 4)
          {
            kernel.template run<8, 4>();
            return true; // fast path taken
          }
        if (degree_fine == 10 && degree_coarse == 5)
          {
            kernel.template run<10, 5>();
            return true; // fast path taken
          }
        if (degree_fine == 12 && degree_coarse == 6)
          {
            kernel.template run<12, 6>();
            return true; // fast path taken
          }

        // due to limited functionality (e.g., FE_Q only) at the moment this
        // should be unreachable due to the checks in the initialization phase
        return false; // indicate that slow path has been taken
      }

    private:
      const unsigned int degree_fine;
      const unsigned int degree_coarse;
    };



    /**
     * Dispatch resolve_hanging_nodes for a run-time coarse degree.
     */
    class HangingNodeDegreeFactory
    {
    public:
      static const unsigned int max_degree = 9;

      explicit HangingNodeDegreeFactory(const unsigned int degree)
        : degree(degree)
      {}

      template <int dim, bool transpose, typename Number, typename ViewType>
      DEAL_II_HOST_DEVICE bool
      run(const Kokkos::TeamPolicy<
            MemorySpace::Default::kokkos_space::execution_space>::member_type
            &team_member,
          Kokkos::View<Number *, MemorySpace::Default::kokkos_space>
            constraint_weights,
          const dealii::internal::MatrixFreeFunctions::ConstraintKinds
                   constraint_mask,
          ViewType values) const
      {
        // resolve_hanging_nodes() only implements dim == 2 and dim == 3.
        if constexpr (dim < 2)
          {
            (void)team_member;
            (void)constraint_weights;
            (void)constraint_mask;
            (void)values;
            return false;
          }
        else
          {
#define DEAL_II_HANGING_NODE_DEGREE_CASE(degree_value)             \
  case degree_value:                                               \
    {                                                              \
      resolve_hanging_nodes<dim, degree_value, transpose, Number>( \
        team_member, constraint_weights, constraint_mask, values); \
      return true;                                                 \
    }

            switch (degree)
              {
                DEAL_II_HANGING_NODE_DEGREE_CASE(1)
                DEAL_II_HANGING_NODE_DEGREE_CASE(2)
                DEAL_II_HANGING_NODE_DEGREE_CASE(3)
                DEAL_II_HANGING_NODE_DEGREE_CASE(4)
                DEAL_II_HANGING_NODE_DEGREE_CASE(5)
                DEAL_II_HANGING_NODE_DEGREE_CASE(6)
                DEAL_II_HANGING_NODE_DEGREE_CASE(7)
                DEAL_II_HANGING_NODE_DEGREE_CASE(8)
                DEAL_II_HANGING_NODE_DEGREE_CASE(9)
                default:
                  return false;
              }
#undef DEAL_II_HANGING_NODE_DEGREE_CASE
          }
      }

    private:
      const unsigned int degree;
    };



    /**
     * Helper class containing the cell-wise prolongation operation.
     */
    template <int dim, typename VectorType>
    class CellProlongator
    {
    public:
      using Number = typename VectorType::value_type;

      using TeamHandle = Kokkos::TeamPolicy<
        MemorySpace::Default::kokkos_space::execution_space>::member_type;

      using SharedViewScratchPad =
        Kokkos::View<Number *,
                     MemorySpace::Default::kokkos_space::execution_space::
                       scratch_memory_space,
                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

      DEAL_II_HOST_DEVICE
      CellProlongator(
        const typename MGTwoLevelTransfer<dim, VectorType>::TransferCellData
          *cell_data)
        : cell_data(cell_data)
      {}

      template <int degree_fine, int degree_coarse>
      DEAL_II_HOST_DEVICE void
      run() const
      {
        const auto &team_member = cell_data->team_member;

        const auto &prolongation_matrix_scratch =
          cell_data->prolongation_matrix_device;

        auto &values_coarse = cell_data->values_coarse;
        auto &values_fine   = cell_data->values_fine;
        auto &scratch_pad   = cell_data->scratch_pad;

        // the evaluator does not need temporary storage since no in-place
        // operation takes place in this function
        auto scratch_for_eval =
          Kokkos::subview(scratch_pad, Kokkos::make_pair(0, 0));

        internal::EvaluatorTensorProduct<
          EvaluatorVariant::evaluate_general,
          dim,
          degree_coarse + 1,
          degree_fine + 1,
          Number,
          MemorySpace::Default::kokkos_space::execution_space::
            scratch_memory_space>
          eval(team_member,
               prolongation_matrix_scratch,
               SharedViewScratchPad(),
               SharedViewScratchPad(),
               scratch_for_eval);

        // apply kernel in each direction
        if constexpr (dim == 1)
          {
            eval.template values<0, true, false, false>(values_coarse,
                                                        values_fine);
          }
        else if constexpr (dim == 2)
          {
            const int temp_size = (degree_coarse + 1) * (degree_fine + 1);
            auto      temp =
              Kokkos::subview(scratch_pad, Kokkos::make_pair(0, temp_size));

            eval.template values<0, true, false, false>(values_coarse, temp);
            eval.template values<1, true, false, false>(temp, values_fine);
          }
        else if constexpr (dim == 3)
          {
            const int temp1_size =
              Utilities::pow(degree_coarse + 1, 2) * (degree_fine + 1);
            const int temp2_size =
              Utilities::pow(degree_fine + 1, 2) * (degree_coarse + 1);
            auto temp1 =
              Kokkos::subview(scratch_pad, Kokkos::make_pair(0, temp1_size));
            auto temp2 =
              Kokkos::subview(scratch_pad,
                              Kokkos::make_pair(temp1_size,
                                                temp1_size + temp2_size));

            eval.template values<0, true, false, false>(values_coarse, temp1);
            eval.template values<1, true, false, false>(temp1, temp2);
            eval.template values<2, true, false, false>(temp2, values_fine);
          }
      }

    private:
      const typename MGTwoLevelTransfer<dim, VectorType>::TransferCellData
        *cell_data;
    };

    /**
     * Helper class containing the cell-wise restriction operation.
     */
    template <int dim, typename VectorType>
    class CellRestrictor
    {
    public:
      using Number = typename VectorType::value_type;

      using TeamHandle = Kokkos::TeamPolicy<
        MemorySpace::Default::kokkos_space::execution_space>::member_type;

      using SharedViewScratchPad =
        Kokkos::View<Number *,
                     MemorySpace::Default::kokkos_space::execution_space::
                       scratch_memory_space,
                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

      DEAL_II_HOST_DEVICE
      CellRestrictor(
        const typename MGTwoLevelTransfer<dim, VectorType>::TransferCellData
          *cell_data)
        : cell_data(cell_data)
      {}

      template <int degree_fine, int degree_coarse>
      DEAL_II_HOST_DEVICE void
      run() const
      {
        const auto &team_member = cell_data->team_member;

        const auto &prolongation_matrix_scratch =
          cell_data->prolongation_matrix_device;

        auto &values_coarse = cell_data->values_coarse;
        auto &values_fine   = cell_data->values_fine;
        auto &scratch_pad   = cell_data->scratch_pad;

        // the evaluator does not need temporary storage since no in-place
        // operation takes place in this function
        auto scratch_for_eval =
          Kokkos::subview(scratch_pad, Kokkos::make_pair(0, 0));

        internal::EvaluatorTensorProduct<
          EvaluatorVariant::evaluate_general,
          dim,
          degree_coarse + 1,
          degree_fine + 1,
          Number,
          MemorySpace::Default::kokkos_space::execution_space::
            scratch_memory_space>
          eval(team_member,
               prolongation_matrix_scratch,
               SharedViewScratchPad(),
               SharedViewScratchPad(),
               scratch_for_eval);

        // apply kernel in each direction
        if constexpr (dim == 1)
          {
            eval.template values<0, false, false, false>(values_fine,
                                                         values_coarse);
          }
        else if constexpr (dim == 2)
          {
            const int temp_size = (degree_coarse + 1) * (degree_fine + 1);
            auto      temp =
              Kokkos::subview(scratch_pad, Kokkos::make_pair(0, temp_size));

            eval.template values<1, false, false, false>(values_fine, temp);
            eval.template values<0, false, false, false>(temp, values_coarse);
          }
        else if constexpr (dim == 3)
          {
            const int temp1_size =
              Utilities::pow(degree_fine + 1, 2) * (degree_coarse + 1);
            const int temp2_size =
              Utilities::pow(degree_coarse + 1, 2) * (degree_fine + 1);
            auto temp1 =
              Kokkos::subview(scratch_pad, Kokkos::make_pair(0, temp1_size));
            auto temp2 =
              Kokkos::subview(scratch_pad,
                              Kokkos::make_pair(temp1_size,
                                                temp1_size + temp2_size));

            eval.template values<2, false, false, false>(values_fine, temp1);
            eval.template values<1, false, false, false>(temp1, temp2);
            eval.template values<0, false, false, false>(temp2, values_coarse);
          }
      }

    private:
      const typename MGTwoLevelTransfer<dim, VectorType>::TransferCellData
        *cell_data;
    };

    class MGTwoLevelTransferImplementation
    {
      /**
       * Compute weights.
       */
      template <int dim, typename Number>
      static void
      setup_weights(
        const dealii::AffineConstraints<Number> &constraints_fine,
        MGTwoLevelTransfer<
          dim,
          LinearAlgebra::distributed::Vector<Number, MemorySpace::Default>>
                  &transfer,
        const bool is_feq)
      {
        if (transfer.fine_element_is_continuous == false)
          return; // nothing to do

        Assert(is_feq, ExcNotImplemented());

        // compute weights globally
        LinearAlgebra::distributed::Vector<Number> weight_vector;
        weight_vector.reinit(transfer.transfer_cpu.partitioner_fine);

        // ... compute valence of DoFs
        for (const auto i :
             transfer.transfer_cpu.constraint_info_fine.dof_indices)
          weight_vector.local_element(i) += 1.0;
        weight_vector.compress(VectorOperation::add);

        // ... invert valence
        for (unsigned int i = 0; i < weight_vector.locally_owned_size(); ++i)
          if (weight_vector.local_element(i) > 0)
            weight_vector.local_element(i) =
              1.0 / weight_vector.local_element(i);

        // ... clear constrained indices
        for (const auto &constrained_dofs : constraints_fine.get_lines())
          if (weight_vector.locally_owned_elements().is_element(
                constrained_dofs.index))
            weight_vector[constrained_dofs.index] = 0.0;

        weight_vector.update_ghost_values();

        unsigned int cell_counter = 0;
        unsigned int scheme_index = 0;
        for (auto &scheme : transfer.schemes)
          {
            scheme.weights =
              Kokkos::View<Number **, MemorySpace::Default::kokkos_space>(
                Kokkos::view_alloc("weights_h_transfer_scheme_" +
                                     std::to_string(scheme_index),
                                   Kokkos::WithoutInitializing),
                scheme.n_dofs_per_cell_fine,
                scheme.n_coarse_cells);

            auto weights_view_host = Kokkos::create_mirror_view(scheme.weights);

            const unsigned int first_cell = cell_counter;
            for (unsigned int cell = 0; cell < scheme.n_coarse_cells; ++cell)
              {
                const unsigned int  coarse_cell_index = first_cell + cell;
                const unsigned int *dof_indices_fine =
                  transfer.transfer_cpu.constraint_info_fine.dof_indices
                    .data() +
                  transfer.transfer_cpu.constraint_info_fine
                    .row_starts[coarse_cell_index]
                    .first;

                for (unsigned int i = 0; i < scheme.n_dofs_per_cell_fine;
                     ++dof_indices_fine, ++i)
                  {
                    weights_view_host(i, cell) =
                      weight_vector.local_element(*dof_indices_fine);
                  }
                ++cell_counter;
              }

            Kokkos::deep_copy(scheme.weights, weights_view_host);
            Kokkos::fence();
          }
      }

    public:
      template <int dim, typename Number>
      static void
      reinit_geometric_transfer(
        const DoFHandler<dim>                   &dof_handler_fine,
        const DoFHandler<dim>                   &dof_handler_coarse,
        const dealii::AffineConstraints<Number> &constraints_fine,
        const dealii::AffineConstraints<Number> &constraints_coarse,
        const unsigned int                       mg_level_fine,
        const unsigned int                       mg_level_coarse,
        MGTwoLevelTransfer<
          dim,
          LinearAlgebra::distributed::Vector<Number, MemorySpace::Default>>
          &transfer)
      {
        AssertDimension(constraints_fine.n_inhomogeneities(), 0);
        AssertDimension(constraints_coarse.n_inhomogeneities(), 0);

        transfer.transfer_cpu.reinit_geometric_transfer(dof_handler_fine,
                                                        dof_handler_coarse,
                                                        constraints_fine,
                                                        constraints_coarse,
                                                        mg_level_fine,
                                                        mg_level_coarse);

        const auto &transfer_cpu = transfer.transfer_cpu;

        transfer.dof_handler_fine = &dof_handler_fine;
        transfer.mg_level_fine    = mg_level_fine;

        transfer.schemes.resize(transfer_cpu.schemes.size());

        transfer.n_components = transfer_cpu.n_components;

        const auto &fe_fine = dof_handler_fine.get_fe();

        const bool is_feq = fe_fine.n_base_elements() == 1 &&
                            ((dynamic_cast<const FE_Q<dim> *>(
                                &fe_fine.base_element(0)) != nullptr));

        Assert(is_feq, ExcNotImplemented());

        transfer.fine_element_is_continuous =
          transfer_cpu.fine_element_is_continuous;

        Assert(transfer.fine_element_is_continuous, ExcNotImplemented());

        for (unsigned int scheme_index = 0;
             scheme_index < transfer.schemes.size();
             ++scheme_index)
          {
            auto       &scheme     = transfer.schemes[scheme_index];
            const auto &scheme_cpu = transfer_cpu.schemes[scheme_index];

            scheme.n_dofs_per_cell_coarse = scheme_cpu.n_dofs_per_cell_coarse;
            scheme.n_dofs_per_cell_fine   = scheme_cpu.n_dofs_per_cell_fine;

            scheme.degree_coarse = scheme_cpu.degree_coarse;
            scheme.degree_fine   = scheme_cpu.degree_fine;

            scheme.n_coarse_cells = scheme_cpu.n_coarse_cells;
          }

        {
          transfer.partitioner_coarse = transfer_cpu.partitioner_coarse;
          transfer.vec_coarse.reinit(transfer.partitioner_coarse);

          transfer.partitioner_fine = transfer_cpu.partitioner_fine;
          transfer.vec_fine.reinit(transfer.partitioner_fine);
        }

        // -------- prolongation matrix (0) -> identity matrix --------
        //
        // nothing to do since for identity prolongation matrices a short-cut
        // code path is used during prolongation/restriction

        // --------------prolongation matrix (i = 1 ... n)--------------
        {
          for (unsigned int transfer_scheme_index = 1;
               transfer_scheme_index < transfer.schemes.size();
               ++transfer_scheme_index)
            {
              const auto prolongation_matrix_host =
                transfer_cpu.schemes[transfer_scheme_index].prolongation_matrix;

              if (prolongation_matrix_host.size() > 0)
                {
                  transfer.schemes[transfer_scheme_index].prolongation_matrix =
                    Kokkos::View<Number *, MemorySpace::Default::kokkos_space>(
                      Kokkos::view_alloc(
                        "prolongation_matrix_h_transfer_scheme_" +
                          std::to_string(transfer_scheme_index),
                        Kokkos::WithoutInitializing),
                      prolongation_matrix_host.size());

                  auto prolongation_matrix_mirror_host =
                    Kokkos::create_mirror_view(
                      transfer.schemes[transfer_scheme_index]
                        .prolongation_matrix);

                  for (unsigned int i = 0; i < prolongation_matrix_host.size();
                       ++i)
                    prolongation_matrix_mirror_host[i] =
                      prolongation_matrix_host[i];

                  Kokkos::deep_copy(
                    transfer.schemes[transfer_scheme_index].prolongation_matrix,
                    prolongation_matrix_mirror_host);
                  Kokkos::fence();
                }
            }
        }

        // ------------------- hanging-node weights -------------------
        const auto &constraint_info_coarse =
          transfer_cpu.constraint_info_coarse;
        const bool has_hanging_nodes =
          constraint_info_coarse.hanging_node_constraint_masks.empty() == false;

        if (has_hanging_nodes)
          {
            Assert(constraint_info_coarse.shape_infos.empty() == false,
                   ExcInternalError());
            const auto &subface_matrix =
              constraint_info_coarse.shape_infos.front()
                .data.front()
                .subface_interpolation_matrices[0];

            transfer.hanging_node_constraint_weights =
              Kokkos::View<Number *, MemorySpace::Default::kokkos_space>(
                Kokkos::view_alloc("h_transfer_hanging_node_weights",
                                   Kokkos::WithoutInitializing),
                subface_matrix.size());

            auto weights_host = Kokkos::create_mirror_view(
              transfer.hanging_node_constraint_weights);
            for (unsigned int i = 0; i < subface_matrix.size(); ++i)
              weights_host[i] = subface_matrix[i];
            Kokkos::deep_copy(transfer.hanging_node_constraint_weights,
                              weights_host);
            Kokkos::fence();
          }

        // ---------- coarse constraint CRS + fine dof indices ----------
        {
          unsigned int scheme_counter = 0;
          unsigned int cell_counter   = 0;
          for (auto &scheme : transfer.schemes)
            {
              if (scheme.n_coarse_cells == 0)
                {
                  ++scheme_counter;
                  continue;
                }

              const unsigned int n_rows =
                scheme.n_coarse_cells * scheme.n_dofs_per_cell_coarse;

              std::vector<unsigned int> constraint_offsets(n_rows + 1);
              std::vector<unsigned int> constraint_indices;
              std::vector<Number>       constraint_weights;
              constraint_indices.reserve(n_rows);
              constraint_weights.reserve(n_rows);

              scheme.dof_indices_fine =
                Kokkos::View<unsigned int **,
                             MemorySpace::Default::kokkos_space>(
                  Kokkos::view_alloc("h_transfer_dof_indices_fine_scheme_" +
                                       std::to_string(scheme_counter),
                                     Kokkos::WithoutInitializing),
                  scheme.n_dofs_per_cell_fine,
                  scheme.n_coarse_cells);

              auto dof_indices_fine_host =
                Kokkos::create_mirror_view(scheme.dof_indices_fine);

              typename std::remove_reference_t<
                decltype(scheme.hanging_node_masks)>::host_mirror_type
                hanging_node_masks_host;
              if (has_hanging_nodes)
                {
                  scheme.hanging_node_masks = Kokkos::View<
                    dealii::internal::MatrixFreeFunctions::ConstraintKinds *,
                    MemorySpace::Default::kokkos_space>(
                    "h_transfer_hanging_node_masks_scheme_" +
                      std::to_string(scheme_counter),
                    scheme.n_coarse_cells);
                  hanging_node_masks_host =
                    Kokkos::create_mirror_view(scheme.hanging_node_masks);
                }

              const unsigned int first_cell = cell_counter;
              for (unsigned int cell = 0; cell < scheme.n_coarse_cells; ++cell)
                {
                  const unsigned int coarse_cell_index = first_cell + cell;

                  // expand ConstraintInfo into a CRS gather/scatter map
                  {
                    const unsigned int *dof_indices_coarse =
                      constraint_info_coarse.dof_indices.data() +
                      constraint_info_coarse.row_starts[coarse_cell_index]
                        .first;
                    unsigned int index_indicators =
                      constraint_info_coarse.row_starts[coarse_cell_index]
                        .second;
                    const unsigned int next_index_indicators =
                      constraint_info_coarse.row_starts[coarse_cell_index + 1]
                        .second;

                    unsigned int ind_local = 0;
                    for (; index_indicators != next_index_indicators;
                         ++index_indicators)
                      {
                        const std::pair<unsigned short, unsigned short>
                          indicator = constraint_info_coarse
                                        .constraint_indicator[index_indicators];

                        for (unsigned int j = 0; j < indicator.first; ++j)
                          {
                            const unsigned int row =
                              cell * scheme.n_dofs_per_cell_coarse + ind_local +
                              j;
                            constraint_offsets[row] = constraint_indices.size();
                            constraint_indices.push_back(dof_indices_coarse[j]);
                            constraint_weights.push_back(Number(1.));
                          }
                        ind_local += indicator.first;
                        dof_indices_coarse += indicator.first;

                        const unsigned int row =
                          cell * scheme.n_dofs_per_cell_coarse + ind_local;
                        constraint_offsets[row] = constraint_indices.size();

                        const unsigned int pool_row = indicator.second;
                        const Number      *data_val =
                          constraint_info_coarse.constraint_pool_data.empty() ?
                                 nullptr :
                                 constraint_info_coarse.constraint_pool_data.data() +
                              constraint_info_coarse
                                .constraint_pool_row_index[pool_row];
                        const Number *end_pool =
                          constraint_info_coarse.constraint_pool_data.empty() ?
                            nullptr :
                            constraint_info_coarse.constraint_pool_data.data() +
                              constraint_info_coarse
                                .constraint_pool_row_index[pool_row + 1];

                        for (; data_val != end_pool;
                             ++data_val, ++dof_indices_coarse)
                          {
                            constraint_indices.push_back(*dof_indices_coarse);
                            constraint_weights.push_back(*data_val);
                          }

                        ++ind_local;
                      }

                    for (; ind_local < scheme.n_dofs_per_cell_coarse;
                         ++dof_indices_coarse, ++ind_local)
                      {
                        const unsigned int row =
                          cell * scheme.n_dofs_per_cell_coarse + ind_local;
                        constraint_offsets[row] = constraint_indices.size();
                        constraint_indices.push_back(*dof_indices_coarse);
                        constraint_weights.push_back(Number(1.));
                      }
                  }

                  // fill fine indices
                  {
                    const unsigned int *dof_indices_fine =
                      transfer_cpu.constraint_info_fine.dof_indices.data() +
                      transfer_cpu.constraint_info_fine
                        .row_starts[coarse_cell_index]
                        .first;

                    for (unsigned int j = 0; j < scheme.n_dofs_per_cell_fine;
                         ++dof_indices_fine, ++j)
                      dof_indices_fine_host(j, cell) = *dof_indices_fine;
                  }

                  if (has_hanging_nodes)
                    hanging_node_masks_host(cell) =
                      dealii::internal::MatrixFreeFunctions::decompress(
                        constraint_info_coarse
                          .hanging_node_constraint_masks[coarse_cell_index],
                        dim);

                  ++cell_counter;
                }

              constraint_offsets[n_rows] = constraint_indices.size();

              scheme.coarse_constraint_offsets =
                Kokkos::View<unsigned int *,
                             MemorySpace::Default::kokkos_space>(
                  Kokkos::view_alloc("h_transfer_coarse_constraint_offsets_" +
                                       std::to_string(scheme_counter),
                                     Kokkos::WithoutInitializing),
                  constraint_offsets.size());
              scheme.coarse_constraint_indices =
                Kokkos::View<unsigned int *,
                             MemorySpace::Default::kokkos_space>(
                  Kokkos::view_alloc("h_transfer_coarse_constraint_indices_" +
                                       std::to_string(scheme_counter),
                                     Kokkos::WithoutInitializing),
                  constraint_indices.size());
              scheme.coarse_constraint_weights =
                Kokkos::View<Number *, MemorySpace::Default::kokkos_space>(
                  Kokkos::view_alloc("h_transfer_coarse_constraint_weights_" +
                                       std::to_string(scheme_counter),
                                     Kokkos::WithoutInitializing),
                  constraint_weights.size());

              {
                auto offsets_host =
                  Kokkos::create_mirror_view(scheme.coarse_constraint_offsets);
                auto indices_host =
                  Kokkos::create_mirror_view(scheme.coarse_constraint_indices);
                auto weights_host =
                  Kokkos::create_mirror_view(scheme.coarse_constraint_weights);

                for (unsigned int i = 0; i < constraint_offsets.size(); ++i)
                  offsets_host[i] = constraint_offsets[i];
                for (unsigned int i = 0; i < constraint_indices.size(); ++i)
                  indices_host[i] = constraint_indices[i];
                for (unsigned int i = 0; i < constraint_weights.size(); ++i)
                  weights_host[i] = constraint_weights[i];

                Kokkos::deep_copy(scheme.coarse_constraint_offsets,
                                  offsets_host);
                Kokkos::deep_copy(scheme.coarse_constraint_indices,
                                  indices_host);
                Kokkos::deep_copy(scheme.coarse_constraint_weights,
                                  weights_host);
                Kokkos::deep_copy(scheme.dof_indices_fine,
                                  dof_indices_fine_host);
                if (has_hanging_nodes)
                  Kokkos::deep_copy(scheme.hanging_node_masks,
                                    hanging_node_masks_host);
                Kokkos::fence();
              }

              ++scheme_counter;
            }
        }

        // ------------------- weights-------------------
        if (transfer.fine_element_is_continuous)
          setup_weights(constraints_fine, transfer, is_feq);
      }

      template <int dim, typename VectorType>
      static void
      prolongate_and_add_internal(
        const MGTwoLevelTransfer<dim, VectorType> &transfer,
        VectorType                                &dst,
        const VectorType                          &src)
      {
        using TeamPolicy = Kokkos::TeamPolicy<
          MemorySpace::Default::kokkos_space::execution_space>;

        using TeamHandle = Kokkos::TeamPolicy<
          MemorySpace::Default::kokkos_space::execution_space>::member_type;

        using Number = typename VectorType::value_type;

        using SharedViewValues =
          Kokkos::View<Number **,
                       Kokkos::LayoutLeft,
                       MemorySpace::Default::kokkos_space::execution_space::
                         scratch_memory_space,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

        using SharedViewScratchPad =
          Kokkos::View<Number *,
                       MemorySpace::Default::kokkos_space::execution_space::
                         scratch_memory_space,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

        MemorySpace::Default::kokkos_space::execution_space exec;

        DeviceVector<Number> dst_device(dst.get_values(),
                                        dst.locally_owned_size());
        DeviceVector<Number> src_device(src.get_values(),
                                        src.locally_owned_size());

        unsigned int scheme_index = 0;

        const unsigned int n_components = transfer.n_components;
        const auto         hanging_node_weights =
          transfer.hanging_node_constraint_weights;

        for (const auto &scheme : transfer.schemes)
          {
            if (scheme.n_coarse_cells == 0)
              continue;

            const unsigned int n_scalar_dofs_per_cell_coarse =
              scheme.n_dofs_per_cell_coarse / n_components;
            const unsigned int n_scalar_dofs_per_cell_fine =
              scheme.n_dofs_per_cell_fine / n_components;

            const bool needs_interpolation =
              scheme.prolongation_matrix.extent(0) > 0;
            const bool has_hanging_nodes =
              scheme.hanging_node_masks.extent(0) > 0;

            auto team_policy =
              TeamPolicy(exec, scheme.n_coarse_cells, Kokkos::AUTO);

            CellTransferFactory      cell_transfer(scheme.degree_fine,
                                              scheme.degree_coarse);
            HangingNodeDegreeFactory hanging_node_factory(scheme.degree_coarse);

            unsigned int scratch_pad_size         = 0;
            unsigned int prolongation_matrix_size = 0;
            if (needs_interpolation)
              {
                prolongation_matrix_size =
                  (scheme.degree_coarse + 1) * (scheme.degree_fine + 1);
                if constexpr (dim == 2)
                  scratch_pad_size =
                    (scheme.degree_coarse + 1) * (scheme.degree_fine + 1);
                else if constexpr (dim == 3)
                  scratch_pad_size =
                    (scheme.degree_coarse + 1) * (scheme.degree_fine + 1) *
                    (scheme.degree_coarse + scheme.degree_fine + 2);
              }

            const auto team_shmem_size =
              SharedViewValues::shmem_size(n_scalar_dofs_per_cell_coarse,
                                           n_components) +
              SharedViewValues::shmem_size(n_scalar_dofs_per_cell_fine,
                                           n_components) +
              SharedViewScratchPad::shmem_size(prolongation_matrix_size) +
              SharedViewScratchPad::shmem_size(scratch_pad_size);

            team_policy.set_scratch_size(0, Kokkos::PerTeam(team_shmem_size));

            const auto coarse_constraint_offsets =
              scheme.coarse_constraint_offsets;
            const auto coarse_constraint_indices =
              scheme.coarse_constraint_indices;
            const auto coarse_constraint_weights =
              scheme.coarse_constraint_weights;
            const auto         hanging_node_masks  = scheme.hanging_node_masks;
            const auto         dof_indices_fine    = scheme.dof_indices_fine;
            const auto         weights             = scheme.weights;
            const auto         prolongation_matrix = scheme.prolongation_matrix;
            const unsigned int n_dofs_coarse = scheme.n_dofs_per_cell_coarse;
            const unsigned int n_dofs_fine   = scheme.n_dofs_per_cell_fine;

            Kokkos::parallel_for(
              "prolongate_and_add_geometric_transfer_scheme_" +
                std::to_string(scheme_index),
              team_policy,
              KOKKOS_LAMBDA(const TeamHandle &team_member) {
                const int coarse_cell_index = team_member.league_rank();

                SharedViewScratchPad prolongation_matrix_device(
                  team_member.team_shmem(), prolongation_matrix_size);

                SharedViewValues values_coarse(team_member.team_shmem(),
                                               n_scalar_dofs_per_cell_coarse,
                                               n_components);

                SharedViewValues values_fine(team_member.team_shmem(),
                                             n_scalar_dofs_per_cell_fine,
                                             n_components);

                SharedViewScratchPad scratch_pad(team_member.team_shmem(),
                                                 scratch_pad_size);

                if (needs_interpolation)
                  {
                    Kokkos::parallel_for(
                      Kokkos::TeamThreadRange(team_member,
                                              prolongation_matrix_size),
                      [&](const int &i) {
                        prolongation_matrix_device(i) = prolongation_matrix(i);
                      });
                    team_member.team_barrier();
                  }

                // read coarse dof values (AffineConstraints via CRS)
                Kokkos::parallel_for(
                  Kokkos::TeamThreadRange(team_member, n_dofs_coarse),
                  [&](const int &thread_id) {
                    const int component =
                      thread_id / n_scalar_dofs_per_cell_coarse;
                    const int local_dof =
                      thread_id % n_scalar_dofs_per_cell_coarse;

                    Number             value = 0.;
                    const unsigned int row =
                      coarse_cell_index * n_dofs_coarse + thread_id;
                    for (unsigned int k = coarse_constraint_offsets(row);
                         k < coarse_constraint_offsets(row + 1);
                         ++k)
                      value += coarse_constraint_weights(k) *
                               src_device[coarse_constraint_indices(k)];

                    values_coarse(local_dof, component) = value;
                  });
                team_member.team_barrier();

                // fast hanging-node interpolation on the coarse cell
                if (has_hanging_nodes)
                  {
                    const auto mask = hanging_node_masks(coarse_cell_index);
                    if (mask != dealii::internal::MatrixFreeFunctions::
                                  ConstraintKinds::unconstrained)
                      for (unsigned int c = 0; c < n_components; ++c)
                        {
                          SharedViewScratchPad values_coarse_component =
                            Kokkos::subview(values_coarse, Kokkos::ALL, c);
                          hanging_node_factory.template run<dim, false, Number>(
                            team_member,
                            hanging_node_weights,
                            mask,
                            values_coarse_component);
                        }
                  }

                if (needs_interpolation)
                  {
                    for (unsigned int c = 0; c < n_components; ++c)
                      {
                        SharedViewScratchPad values_coarse_component =
                          Kokkos::subview(values_coarse, Kokkos::ALL, c);
                        SharedViewScratchPad values_fine_component =
                          Kokkos::subview(values_fine, Kokkos::ALL, c);

                        typename MGTwoLevelTransfer<dim, VectorType>::
                          TransferCellData cell_data{team_member,
                                                     prolongation_matrix_device,
                                                     values_coarse_component,
                                                     values_fine_component,
                                                     scratch_pad};

                        CellProlongator<dim, VectorType> cell_prolongator(
                          &cell_data);

                        cell_transfer.run(cell_prolongator);
                      }
                  }
                else
                  {
                    // identity scheme: cell exists on both levels
                    Kokkos::parallel_for(
                      Kokkos::TeamThreadRange(team_member, n_dofs_coarse),
                      [&](const int &thread_id) {
                        const int component =
                          thread_id / n_scalar_dofs_per_cell_coarse;
                        const int local_dof =
                          thread_id % n_scalar_dofs_per_cell_coarse;
                        values_fine(local_dof, component) =
                          values_coarse(local_dof, component);
                      });
                    team_member.team_barrier();
                  }

                // apply weights and add fine values into results vector
                Kokkos::parallel_for(
                  Kokkos::TeamThreadRange(team_member, n_dofs_fine),
                  [&](const int &thread_id) {
                    const int component =
                      thread_id / n_scalar_dofs_per_cell_fine;
                    const int local_dof =
                      thread_id % n_scalar_dofs_per_cell_fine;

                    const unsigned int global_dof =
                      dof_indices_fine(thread_id, coarse_cell_index);

                    const Number weight = weights(thread_id, coarse_cell_index);

                    Kokkos::atomic_add(&dst_device[global_dof],
                                       weight *
                                         values_fine(local_dof, component));
                  });
                team_member.team_barrier();
              });

            ++scheme_index;
          }
      }

      template <int dim, typename VectorType>
      static void
      restrict_and_add_internal(
        const MGTwoLevelTransfer<dim, VectorType> &transfer,
        VectorType                                &dst,
        const VectorType                          &src)
      {
        using TeamPolicy = Kokkos::TeamPolicy<
          MemorySpace::Default::kokkos_space::execution_space>;

        using TeamHandle = Kokkos::TeamPolicy<
          MemorySpace::Default::kokkos_space::execution_space>::member_type;

        using Number = typename VectorType::value_type;

        using SharedViewValues =
          Kokkos::View<Number **,
                       Kokkos::LayoutLeft,
                       MemorySpace::Default::kokkos_space::execution_space::
                         scratch_memory_space,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

        using SharedViewScratchPad =
          Kokkos::View<Number *,
                       MemorySpace::Default::kokkos_space::execution_space::
                         scratch_memory_space,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

        MemorySpace::Default::kokkos_space::execution_space exec;

        DeviceVector<Number> dst_device(dst.get_values(),
                                        dst.locally_owned_size());
        DeviceVector<Number> src_device(src.get_values(),
                                        src.locally_owned_size());

        const unsigned int n_components = transfer.n_components;
        const auto         hanging_node_weights =
          transfer.hanging_node_constraint_weights;

        unsigned int scheme_index = 0;

        for (const auto &scheme : transfer.schemes)
          {
            if (scheme.n_coarse_cells == 0)
              continue;

            const unsigned int n_scalar_dofs_per_cell_coarse =
              scheme.n_dofs_per_cell_coarse / n_components;
            const unsigned int n_scalar_dofs_per_cell_fine =
              scheme.n_dofs_per_cell_fine / n_components;

            const bool needs_interpolation =
              scheme.prolongation_matrix.extent(0) > 0;
            const bool has_hanging_nodes =
              scheme.hanging_node_masks.extent(0) > 0;

            auto team_policy =
              TeamPolicy(exec, scheme.n_coarse_cells, Kokkos::AUTO);

            CellTransferFactory      cell_transfer(scheme.degree_fine,
                                              scheme.degree_coarse);
            HangingNodeDegreeFactory hanging_node_factory(scheme.degree_coarse);

            unsigned int scratch_pad_size         = 0;
            unsigned int prolongation_matrix_size = 0;
            if (needs_interpolation)
              {
                prolongation_matrix_size =
                  (scheme.degree_coarse + 1) * (scheme.degree_fine + 1);
                if constexpr (dim == 2)
                  scratch_pad_size =
                    (scheme.degree_coarse + 1) * (scheme.degree_fine + 1);
                else if constexpr (dim == 3)
                  scratch_pad_size =
                    (scheme.degree_coarse + 1) * (scheme.degree_fine + 1) *
                    (scheme.degree_coarse + scheme.degree_fine + 2);
              }

            const auto team_shmem_size =
              SharedViewValues::shmem_size(n_scalar_dofs_per_cell_coarse,
                                           n_components) +
              SharedViewValues::shmem_size(n_scalar_dofs_per_cell_fine,
                                           n_components) +
              SharedViewScratchPad::shmem_size(prolongation_matrix_size) +
              SharedViewScratchPad::shmem_size(scratch_pad_size);

            team_policy.set_scratch_size(0, Kokkos::PerTeam(team_shmem_size));

            const auto coarse_constraint_offsets =
              scheme.coarse_constraint_offsets;
            const auto coarse_constraint_indices =
              scheme.coarse_constraint_indices;
            const auto coarse_constraint_weights =
              scheme.coarse_constraint_weights;
            const auto         hanging_node_masks  = scheme.hanging_node_masks;
            const auto         dof_indices_fine    = scheme.dof_indices_fine;
            const auto         weights             = scheme.weights;
            const auto         prolongation_matrix = scheme.prolongation_matrix;
            const unsigned int n_dofs_coarse = scheme.n_dofs_per_cell_coarse;
            const unsigned int n_dofs_fine   = scheme.n_dofs_per_cell_fine;

            Kokkos::parallel_for(
              "restrict_and_add_geometric_transfer_scheme_" +
                std::to_string(scheme_index),
              team_policy,
              KOKKOS_LAMBDA(const TeamHandle &team_member) {
                const int coarse_cell_index = team_member.league_rank();

                SharedViewScratchPad prolongation_matrix_device(
                  team_member.team_shmem(), prolongation_matrix_size);

                SharedViewValues values_coarse(team_member.team_shmem(),
                                               n_scalar_dofs_per_cell_coarse,
                                               n_components);

                SharedViewValues values_fine(team_member.team_shmem(),
                                             n_scalar_dofs_per_cell_fine,
                                             n_components);

                SharedViewScratchPad scratch_pad(team_member.team_shmem(),
                                                 scratch_pad_size);

                if (needs_interpolation)
                  {
                    Kokkos::parallel_for(
                      Kokkos::TeamThreadRange(team_member,
                                              prolongation_matrix_size),
                      [&](const int &i) {
                        prolongation_matrix_device(i) = prolongation_matrix(i);
                      });
                    team_member.team_barrier();
                  }

                // read fine dof values and apply weights
                Kokkos::parallel_for(
                  Kokkos::TeamThreadRange(team_member, n_dofs_fine),
                  [&](const int &thread_id) {
                    const int component =
                      thread_id / n_scalar_dofs_per_cell_fine;
                    const int local_dof =
                      thread_id % n_scalar_dofs_per_cell_fine;

                    const Number weight = weights(thread_id, coarse_cell_index);

                    const unsigned int global_dof =
                      dof_indices_fine(thread_id, coarse_cell_index);

                    values_fine(local_dof, component) =
                      weight * src_device[global_dof];
                  });
                team_member.team_barrier();

                if (needs_interpolation)
                  {
                    for (unsigned int c = 0; c < n_components; ++c)
                      {
                        SharedViewScratchPad values_coarse_component =
                          Kokkos::subview(values_coarse, Kokkos::ALL, c);
                        SharedViewScratchPad values_fine_component =
                          Kokkos::subview(values_fine, Kokkos::ALL, c);

                        typename MGTwoLevelTransfer<dim, VectorType>::
                          TransferCellData cell_data{team_member,
                                                     prolongation_matrix_device,
                                                     values_coarse_component,
                                                     values_fine_component,
                                                     scratch_pad};

                        CellRestrictor<dim, VectorType> cell_restrictor(
                          &cell_data);

                        cell_transfer.run(cell_restrictor);
                      }
                  }
                else
                  {
                    Kokkos::parallel_for(
                      Kokkos::TeamThreadRange(team_member, n_dofs_coarse),
                      [&](const int &thread_id) {
                        const int component =
                          thread_id / n_scalar_dofs_per_cell_coarse;
                        const int local_dof =
                          thread_id % n_scalar_dofs_per_cell_coarse;
                        values_coarse(local_dof, component) =
                          values_fine(local_dof, component);
                      });
                    team_member.team_barrier();
                  }

                // transpose hanging-node interpolation before scattering
                if (has_hanging_nodes)
                  {
                    const auto mask = hanging_node_masks(coarse_cell_index);
                    if (mask != dealii::internal::MatrixFreeFunctions::
                                  ConstraintKinds::unconstrained)
                      for (unsigned int c = 0; c < n_components; ++c)
                        {
                          SharedViewScratchPad values_coarse_component =
                            Kokkos::subview(values_coarse, Kokkos::ALL, c);
                          hanging_node_factory.template run<dim, true, Number>(
                            team_member,
                            hanging_node_weights,
                            mask,
                            values_coarse_component);
                        }
                  }

                // distribute coarse dof values (AffineConstraints via CRS)
                Kokkos::parallel_for(
                  Kokkos::TeamThreadRange(team_member, n_dofs_coarse),
                  [&](const int &thread_id) {
                    const int component =
                      thread_id / n_scalar_dofs_per_cell_coarse;
                    const int local_dof =
                      thread_id % n_scalar_dofs_per_cell_coarse;

                    const Number value = values_coarse(local_dof, component);
                    const unsigned int row =
                      coarse_cell_index * n_dofs_coarse + thread_id;
                    for (unsigned int k = coarse_constraint_offsets(row);
                         k < coarse_constraint_offsets(row + 1);
                         ++k)
                      Kokkos::atomic_add(
                        &dst_device[coarse_constraint_indices(k)],
                        coarse_constraint_weights(k) * value);
                  });
                team_member.team_barrier();
              });

            ++scheme_index;
          }
      }
    };
  } // namespace internal

  template <int dim, typename VectorType>
  void
  MGTwoLevelTransfer<dim, VectorType>::prolongate_and_add_internal(
    VectorType       &dst,
    const VectorType &src) const
  {
    internal::MGTwoLevelTransferImplementation::prolongate_and_add_internal(
      *this, dst, src);
  }

  template <int dim, typename VectorType>
  void
  MGTwoLevelTransfer<dim, VectorType>::restrict_and_add_internal(
    VectorType       &dst,
    const VectorType &src) const
  {
    internal::MGTwoLevelTransferImplementation::restrict_and_add_internal(*this,
                                                                          dst,
                                                                          src);
  }

  template <int dim, typename VectorType>
  std::pair<const DoFHandler<dim> *, unsigned int>
  MGTwoLevelTransfer<dim, VectorType>::get_dof_handler_fine() const
  {
    return {this->dof_handler_fine, this->mg_level_fine};
  }

  template <int dim, typename VectorType>
  std::size_t
  MGTwoLevelTransfer<dim, VectorType>::memory_consumption() const
  {
    Assert(false, ExcNotImplemented());

    return 0;
  }

  template <int dim, typename VectorType>
  void
  MGTwoLevelTransfer<dim, VectorType>::interpolate(VectorType       &dst,
                                                   const VectorType &src) const
  {
    (void)dst;
    (void)src;
    Assert(false, ExcNotImplemented());
  }

  template <int dim, typename VectorType>
  std::pair<bool, bool>
  MGTwoLevelTransfer<dim, VectorType>::enable_inplace_operations_if_possible(
    const std::shared_ptr<const Utilities::MPI::Partitioner>
      &external_partitioner_coarse,
    const std::shared_ptr<const Utilities::MPI::Partitioner>
      &external_partitioner_fine)
  {
    return this->internal_enable_inplace_operations_if_possible(
      external_partitioner_coarse,
      external_partitioner_fine,
      this->vec_fine_needs_ghost_update,
      transfer_cpu.constraint_info_coarse,
      transfer_cpu.constraint_info_fine.dof_indices);
  }

  template <int dim, typename VectorType>
  void
  MGTwoLevelTransfer<dim, VectorType>::reinit_geometric_transfer(
    const DoFHandler<dim>           &dof_handler_fine,
    const DoFHandler<dim>           &dof_handler_coarse,
    const AffineConstraints<Number> &constraints_fine,
    const AffineConstraints<Number> &constraints_coarse,
    const unsigned int               mg_level_fine,
    const unsigned int               mg_level_coarse)
  {
    internal::MGTwoLevelTransferImplementation::reinit_geometric_transfer(
      dof_handler_fine,
      dof_handler_coarse,
      constraints_fine,
      constraints_coarse,
      mg_level_fine,
      mg_level_coarse,
      *this);
  }
} // namespace Portable

DEAL_II_NAMESPACE_CLOSE

#endif
