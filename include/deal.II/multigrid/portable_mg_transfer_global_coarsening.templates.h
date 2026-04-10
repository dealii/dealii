// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_mg_portable_global_coarsening_templates_h
#define dealii_mg_portable_global_coarsening_templates_h

#include <deal.II/base/config.h>

#include <deal.II/base/memory_space.h>
#include <deal.II/base/mpi_consensus_algorithms.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/multigrid/mg_transfer_matrix_free.h>
#include <deal.II/multigrid/mg_transfer_matrix_free.templates.h>
#include <deal.II/multigrid/portable_mg_transfer_global_coarsening.h>

DEAL_II_NAMESPACE_OPEN

namespace Portable
{
  namespace internal
  {

    template <int dim, typename VectorType, typename Functor>
    struct ApplyCellKernel
    {
      using Number = typename VectorType::value_type;

      using TeamHandle = Kokkos::TeamPolicy<
        MemorySpace::Default::kokkos_space::execution_space>::member_type;

      using SharedViewValues =
        Kokkos::View<Number *,
                     MemorySpace::Default::kokkos_space::execution_space::
                       scratch_memory_space,
                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

      ApplyCellKernel(
        Functor func,
        const typename MGTwoLevelTransfer<dim, VectorType>::MGTransferScheme
                          transfer_scheme,
        const VectorType &src,
        const VectorType &dst)
        : func(func)
        , transfer_scheme(transfer_scheme)
        , src(src.get_values(), src.locally_owned_size())
        , dst(dst.get_values(), dst.locally_owned_size())
      {}

      Functor func;

      const typename MGTwoLevelTransfer<dim, VectorType>::MGTransferScheme
        transfer_scheme;

      const DeviceVector<Number> src;
      const DeviceVector<Number> dst;

      // Provide the shared memory capacity. This function takes the team_size
      // as an argument, which allows team_size dependent allocations.
      std::size_t
      team_shmem_size(int /*team_size*/) const
      {
        return SharedViewValues::shmem_size(
          transfer_scheme.n_dofs_per_cell_coarse + // coarse dof values
          transfer_scheme.n_dofs_per_cell_fine +   // fine dof values
          2 * transfer_scheme
                .n_dofs_per_cell_fine // at most two tmp vectors of at
                                      // most n_dofs_per_cell_fine size
          + (transfer_scheme.degree_coarse + 1) *
              (transfer_scheme.degree_fine + 1) // prolongation matrix
        );
      }

      DEAL_II_HOST_DEVICE
      void
      operator()(const TeamHandle &team_member) const
      {
        const int coarse_cell_index = team_member.league_rank();

        SharedViewValues values_coarse(team_member.team_shmem(),
                                       transfer_scheme.n_dofs_per_cell_coarse);

        SharedViewValues values_fine(team_member.team_shmem(),
                                     transfer_scheme.n_dofs_per_cell_fine);

        SharedViewValues prolongation_matrix_device(
          team_member.team_shmem(),
          (transfer_scheme.degree_coarse + 1) *
            (transfer_scheme.degree_fine + 1));

        SharedViewValues scratch_pad(team_member.team_shmem(),
                                     transfer_scheme.n_dofs_per_cell_fine * 2);

        Kokkos::parallel_for(
          Kokkos::TeamThreadRange(team_member,
                                  (transfer_scheme.degree_coarse + 1) *
                                    (transfer_scheme.degree_fine + 1)),
          [&](const int &i) {
            prolongation_matrix_device(i) =
              transfer_scheme.prolongation_matrix(i);
          });
        team_member.team_barrier();

        typename MGTwoLevelTransfer<dim, VectorType>::TransferCellData data{
          team_member,
          coarse_cell_index,
          transfer_scheme,
          prolongation_matrix_device,
          values_coarse,
          values_fine,
          scratch_pad};

        func(&data, src, dst);
      }
    };

    template <int dim, typename VectorType>
    class CellProlongationKernel : public EnableObserverPointer
    {
    public:
      using Number = typename VectorType::value_type;

      using TeamHandle = Kokkos::TeamPolicy<
        MemorySpace::Default::kokkos_space::execution_space>::member_type;

      using SharedViewValues =
        Kokkos::View<Number *,
                     MemorySpace::Default::kokkos_space::execution_space::
                       scratch_memory_space,
                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

      CellProlongationKernel() = default;

      DEAL_II_HOST_DEVICE void
      operator()(
        const typename MGTwoLevelTransfer<dim, VectorType>::TransferCellData
                                   *cell_data,
        const DeviceVector<Number> &src,
        const DeviceVector<Number> &dst) const;
    };

    template <int dim, typename VectorType>
    DEAL_II_HOST_DEVICE void
    CellProlongationKernel<dim, VectorType>::operator()(
      const typename MGTwoLevelTransfer<dim, VectorType>::TransferCellData
                                 *cell_data,
      const DeviceVector<Number> &src,
      const DeviceVector<Number> &dst) const
    {
      const int   coarse_cell_index = cell_data->coarse_cell_index;
      const auto &team_member       = cell_data->team_member;

      const auto &prolongation_matrix_scratch =
        cell_data->prolongation_matrix_device;

      const auto &transfer_data = cell_data->transfer_data;

      const auto &dof_indices_coarse = transfer_data.dof_indices_coarse;
      const auto &dof_indices_fine   = transfer_data.dof_indices_fine;

      const unsigned int degree_coarse = transfer_data.degree_coarse;
      const unsigned int degree_fine   = transfer_data.degree_fine;

      auto &values_coarse = cell_data->values_coarse;
      auto &values_fine   = cell_data->values_fine;
      auto &scratch_pad   = cell_data->scratch_pad;

      // read coarse dof values
      Kokkos::parallel_for(Kokkos::TeamThreadRange(
                             team_member, transfer_data.n_dofs_per_cell_coarse),
                           [&](const int &i) {
                             const unsigned int dof_index =
                               dof_indices_coarse(i, coarse_cell_index);
                             if (dof_index != numbers::invalid_unsigned_int)
                               values_coarse(i) = src[dof_index];
                             else
                               values_coarse(i) = 0.;
                           });

      team_member.team_barrier();

#if DEAL_II_KOKKOS_VERSION_GTE(4, 0, 0)

      // apply kernel in each direction
      if constexpr (dim == 2)
        {
          const int temp_size = (degree_coarse + 1) * (degree_fine + 1);
          auto      tmp =
            Kokkos::subview(scratch_pad, Kokkos::make_pair(0, temp_size));

          {
            const int Ni = degree_coarse + 1;
            const int Nj = degree_fine + 1;
            const int Nk = degree_coarse + 1;

            auto thread_policy =
              Kokkos::TeamThreadMDRange<Kokkos::Rank<2>, TeamHandle>(
                team_member, Ni, Nj);
            Kokkos::parallel_for(thread_policy, [&](const int i, const int j) {
              const int base_kernel   = j;
              const int stride_kernel = degree_fine + 1;

              const int base_coarse   = i * Nk;
              const int stride_coarse = 1;

              Number sum = prolongation_matrix_scratch(base_kernel) *
                           values_coarse(base_coarse);

              for (int k = 1; k < Nk; ++k)
                sum +=
                  prolongation_matrix_scratch(base_kernel + k * stride_kernel) *
                  values_coarse(base_coarse + k * stride_coarse);

              const int index_tmp = i * Nj + j;

              tmp(index_tmp) = sum;
            });
          }

          team_member.team_barrier();

          {
            const int Ni = degree_fine + 1;
            const int Nj = degree_fine + 1;
            const int Nk = degree_coarse + 1;

            auto thread_policy =
              Kokkos::TeamThreadMDRange<Kokkos::Rank<2>, TeamHandle>(
                team_member, Ni, Nj);
            Kokkos::parallel_for(thread_policy, [&](const int i, const int j) {
              const int base_kernel   = j;
              const int stride_kernel = degree_fine + 1;

              const int base_tmp   = i;
              const int stride_tmp = degree_fine + 1;

              Number sum =
                prolongation_matrix_scratch(base_kernel) * tmp(base_tmp);

              for (int k = 1; k < Nk; ++k)
                sum +=
                  prolongation_matrix_scratch(base_kernel + k * stride_kernel) *
                  tmp(base_tmp + k * stride_tmp);

              const int index_fine    = i + j * Ni;
              values_fine(index_fine) = sum;
            });
          }

          team_member.team_barrier();
        }
      else if constexpr (dim == 3)
        {
          const int tmp1_size =
            Utilities::pow(degree_coarse + 1, 2) * (degree_fine + 1);
          const int tmp2_size =
            Utilities::pow(degree_fine + 1, 2) * (degree_coarse + 1);
          auto tmp1 =
            Kokkos::subview(scratch_pad, Kokkos::make_pair(0, tmp1_size));
          auto tmp2 = Kokkos::subview(scratch_pad,
                                      Kokkos::make_pair(tmp1_size,
                                                        tmp1_size + tmp2_size));

          {
            const int Ni = degree_coarse + 1;
            const int Nj = degree_coarse + 1;
            const int Nm = degree_fine + 1;
            const int Nk = degree_coarse + 1;

            auto thread_policy =
              Kokkos::TeamThreadMDRange<Kokkos::Rank<3>, TeamHandle>(
                team_member, Ni, Nj, Nm);
            Kokkos::parallel_for(
              thread_policy, [&](const int i, const int j, const int m) {
                const int base_kernel   = m;
                const int stride_kernel = degree_fine + 1;

                const int base_coarse   = (i * Nj + j) * Nk;
                const int stride_coarse = 1;

                Number sum = prolongation_matrix_scratch(base_kernel) *
                             values_coarse(base_coarse);

                for (int k = 1; k < Nk; ++k)
                  sum += prolongation_matrix_scratch(base_kernel +
                                                     k * stride_kernel) *
                         values_coarse(base_coarse + k * stride_coarse);

                const int index_tmp1 = (i * Nj + j) * Nm + m;
                tmp1(index_tmp1)     = sum;
              });
          }

          team_member.team_barrier();

          {
            const int Ni = degree_fine + 1;
            const int Nj = degree_coarse + 1;
            const int Nm = degree_fine + 1;
            const int Nk = degree_coarse + 1;

            auto thread_policy =
              Kokkos::TeamThreadMDRange<Kokkos::Rank<3>, TeamHandle>(
                team_member, Ni, Nj, Nm);
            Kokkos::parallel_for(thread_policy,
                                 [&](const int i, const int j, const int m) {
                                   const int base_kernel   = m;
                                   const int stride_kernel = degree_fine + 1;

                                   const int base_tmp1   = i + j * Ni * Nk;
                                   const int stride_tmp1 = degree_fine + 1;

                                   Number sum =
                                     prolongation_matrix_scratch(base_kernel) *
                                     tmp1(base_tmp1);

                                   for (int k = 1; k < Nk; ++k)
                                     sum += prolongation_matrix_scratch(
                                              base_kernel + k * stride_kernel) *
                                            tmp1(base_tmp1 + k * stride_tmp1);

                                   const int index_tmp2 = i + (j * Nm + m) * Ni;
                                   tmp2(index_tmp2)     = sum;
                                 });
          }

          team_member.team_barrier();

          {
            const int Ni = degree_fine + 1;
            const int Nj = degree_fine + 1;
            const int Nm = degree_fine + 1;
            const int Nk = degree_coarse + 1;

            auto thread_policy =
              Kokkos::TeamThreadMDRange<Kokkos::Rank<3>, TeamHandle>(
                team_member, Ni, Nj, Nm);
            Kokkos::parallel_for(
              thread_policy, [&](const int i, const int j, const int m) {
                const int base_kernel   = m;
                const int stride_kernel = degree_fine + 1;

                const int base_tmp2   = i * Nj + j;
                const int stride_tmp2 = Utilities::pow(degree_fine + 1, 2);
                Number    sum =
                  prolongation_matrix_scratch(base_kernel) * tmp2(base_tmp2);

                for (int k = 1; k < Nk; ++k)
                  sum += prolongation_matrix_scratch(base_kernel +
                                                     k * stride_kernel) *
                         tmp2(base_tmp2 + k * stride_tmp2);

                const int index_fine    = (i + m * Ni) * Nj + j;
                values_fine(index_fine) = sum;
              });
          }
          team_member.team_barrier();
        }
#else
      AssertThrow(false,
                  ExcMessage(
                    "This feature requires Kokkos version 4.0.0 or newer."));
#endif
      // apply weights
      Kokkos::parallel_for(Kokkos::TeamThreadRange(
                             team_member, transfer_data.n_dofs_per_cell_fine),
                           [&](const int &i) {
                             values_fine(i) *=
                               transfer_data.weights(i, coarse_cell_index);
                           });
      team_member.team_barrier();

      Kokkos::parallel_for(
        Kokkos::TeamThreadRange(team_member,
                                transfer_data.n_dofs_per_cell_fine),
        [&](const int &i) {
          const unsigned int dof_index = dof_indices_fine(i, coarse_cell_index);
          Kokkos::atomic_add(&dst[dof_index], values_fine(i));
        });
      team_member.team_barrier();
    }

    template <int dim, typename VectorType>
    class CellRestrictionKernel : public EnableObserverPointer
    {
    public:
      using Number = typename VectorType::value_type;

      using TeamHandle = Kokkos::TeamPolicy<
        MemorySpace::Default::kokkos_space::execution_space>::member_type;

      using SharedView = Kokkos::View<Number *,
                                      MemorySpace::Default::kokkos_space::
                                        execution_space::scratch_memory_space,
                                      Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

      CellRestrictionKernel() = default;

      DEAL_II_HOST_DEVICE void
      operator()(
        const typename MGTwoLevelTransfer<dim, VectorType>::TransferCellData
                                   *cell_data,
        const DeviceVector<Number> &src,
        const DeviceVector<Number> &dst) const;
    };

    template <int dim, typename VectorType>
    DEAL_II_HOST_DEVICE void
    CellRestrictionKernel<dim, VectorType>::operator()(
      const typename MGTwoLevelTransfer<dim, VectorType>::TransferCellData
                                 *cell_data,
      const DeviceVector<Number> &src,
      const DeviceVector<Number> &dst) const
    {
      const int   coarse_cell_index = cell_data->coarse_cell_index;
      const auto &team_member       = cell_data->team_member;

      const auto &prolongation_matrix_scratch =
        cell_data->prolongation_matrix_device;

      const auto &transfer_data = cell_data->transfer_data;

      const auto &dof_indices_coarse = transfer_data.dof_indices_coarse;
      const auto &dof_indices_fine   = transfer_data.dof_indices_fine;

      const unsigned int degree_coarse = transfer_data.degree_coarse;
      const unsigned int degree_fine   = transfer_data.degree_fine;

      auto &values_coarse = cell_data->values_coarse;
      auto &values_fine   = cell_data->values_fine;
      auto &scratch_pad   = cell_data->scratch_pad;

      // read fine dof values
      Kokkos::parallel_for(Kokkos::TeamThreadRange(
                             team_member, transfer_data.n_dofs_per_cell_fine),
                           [&](const int &i) {
                             values_fine(i) =
                               src[dof_indices_fine(i, coarse_cell_index)];
                           });
      team_member.team_barrier();

      // apply weights
      Kokkos::parallel_for(Kokkos::TeamThreadRange(
                             team_member, transfer_data.n_dofs_per_cell_fine),
                           [&](const int &i) {
                             values_fine(i) *=
                               transfer_data.weights(i, coarse_cell_index);
                           });
      team_member.team_barrier();

#if DEAL_II_KOKKOS_VERSION_GTE(4, 0, 0)

      // apply kernel in each direction
      if constexpr (dim == 2)
        {
          const int tmp_size = (degree_coarse + 1) * (degree_fine + 1);

          auto tmp =
            Kokkos::subview(scratch_pad, Kokkos::make_pair(0, tmp_size));
          {
            const int Ni = degree_fine + 1;
            const int Nj = degree_coarse + 1;
            const int Nk = degree_fine + 1;

            auto thread_policy =
              Kokkos::TeamThreadMDRange<Kokkos::Rank<2>, TeamHandle>(
                team_member, Ni, Nj);
            Kokkos::parallel_for(thread_policy, [&](const int i, const int j) {
              const int base_kernel   = j * (degree_fine + 1);
              const int stride_kernel = 1;

              const int base_fine   = i;
              const int stride_fine = degree_fine + 1;

              Number sum = prolongation_matrix_scratch(base_kernel) *
                           values_fine(base_fine);

              for (int k = 1; k < Nk; ++k)
                sum +=
                  prolongation_matrix_scratch(base_kernel + k * stride_kernel) *
                  values_fine(base_fine + k * stride_fine);

              const int index_tmp = i + j * Ni;

              tmp(index_tmp) = sum;
            });
          }

          team_member.team_barrier();

          {
            const int Ni = degree_coarse + 1;
            const int Nj = degree_coarse + 1;
            const int Nk = degree_fine + 1;

            auto thread_policy =
              Kokkos::TeamThreadMDRange<Kokkos::Rank<2>, TeamHandle>(
                team_member, Ni, Nj);
            Kokkos::parallel_for(thread_policy, [&](const int i, const int j) {
              const int base_kernel   = j * (degree_fine + 1);
              const int stride_kernel = 1;

              const int base_tmp   = i * Nk;
              const int stride_tmp = 1;

              Number sum =
                prolongation_matrix_scratch(base_kernel) * tmp(base_tmp);

              for (int k = 1; k < Nk; ++k)
                sum +=
                  prolongation_matrix_scratch(base_kernel + k * stride_kernel) *
                  tmp(base_tmp + k * stride_tmp);

              const int index_coarse = i * Nj + j;

              values_coarse(index_coarse) = sum;
            });
          }

          team_member.team_barrier();
        }
      else if constexpr (dim == 3)
        {
          const int tmp1_size =
            Utilities::pow(degree_fine + 1, 2) * (degree_coarse + 1);
          const int tmp2_size =
            Utilities::pow(degree_coarse + 1, 2) * (degree_fine + 1);

          auto tmp1 =
            Kokkos::subview(scratch_pad, Kokkos::make_pair(0, tmp1_size));
          auto tmp2 = Kokkos::subview(scratch_pad,
                                      Kokkos::make_pair(tmp1_size,
                                                        tmp1_size + tmp2_size));
          {
            const int Ni = degree_fine + 1;
            const int Nj = degree_fine + 1;
            const int Nm = degree_coarse + 1;
            const int Nk = degree_fine + 1;

            auto thread_policy =
              Kokkos::TeamThreadMDRange<Kokkos::Rank<3>, TeamHandle>(
                team_member, Ni, Nj, Nm);
            Kokkos::parallel_for(
              thread_policy, [&](const int i, const int j, const int m) {
                const int base_kernel   = m * (degree_fine + 1);
                const int stride_kernel = 1;

                const int base_fine   = i * Nj + j;
                const int stride_fine = Utilities::pow(degree_fine + 1, 2);

                Number sum = prolongation_matrix_scratch(base_kernel) *
                             values_fine(base_fine);

                for (int k = 1; k < Nk; ++k)
                  sum += prolongation_matrix_scratch(base_kernel +
                                                     k * stride_kernel) *
                         values_fine(base_fine + k * stride_fine);

                const int index_tmp1 = (i + m * Ni) * Nj + j;
                tmp1(index_tmp1)     = sum;
              });
          }

          team_member.team_barrier();

          {
            const int Ni = degree_fine + 1;
            const int Nj = degree_coarse + 1;
            const int Nm = degree_coarse + 1;
            const int Nk = degree_fine + 1;

            auto thread_policy =
              Kokkos::TeamThreadMDRange<Kokkos::Rank<3>, TeamHandle>(
                team_member, Ni, Nj, Nm);
            Kokkos::parallel_for(
              thread_policy, [&](const int i, const int j, const int m) {
                const int base_kernel   = m * (degree_fine + 1);
                const int stride_kernel = 1;

                const int base_tmp1   = i + j * Ni * Nk;
                const int stride_tmp1 = degree_fine + 1;
                Number    sum =
                  prolongation_matrix_scratch(base_kernel) * tmp1(base_tmp1);

                for (int k = 1; k < Nk; ++k)
                  sum += prolongation_matrix_scratch(base_kernel +
                                                     k * stride_kernel) *
                         tmp1(base_tmp1 + k * stride_tmp1);

                const int index_tmp2 = i + (j * Nm + m) * Ni;
                tmp2(index_tmp2)     = sum;
              });
          }

          team_member.team_barrier();

          {
            const int Ni = degree_coarse + 1;
            const int Nj = degree_coarse + 1;
            const int Nm = degree_coarse + 1;
            const int Nk = degree_fine + 1;

            auto thread_policy =
              Kokkos::TeamThreadMDRange<Kokkos::Rank<3>, TeamHandle>(
                team_member, Ni, Nj, Nm);

            Kokkos::parallel_for(
              thread_policy, [&](const int i, const int j, const int m) {
                const int base_kernel   = m * (degree_fine + 1);
                const int stride_kernel = 1;

                const int base_tmp2   = (i * Nj + j) * Nk;
                const int stride_tmp2 = 1;

                Number sum =
                  prolongation_matrix_scratch(base_kernel) * tmp2(base_tmp2);

                for (int k = 1; k < Nk; ++k)
                  sum += prolongation_matrix_scratch(base_kernel +
                                                     k * stride_kernel) *
                         tmp2(base_tmp2 + k * stride_tmp2);

                const int index_coarse      = (i * Nj + j) * Nm + m;
                values_coarse(index_coarse) = sum;
              });
          }

          team_member.team_barrier();
        }

#else
      AssertThrow(false,
                  ExcMessage(
                    "This feature requires Kokkos version 4.0.0 or newer."));
#endif

      // distribute coarse dofs values
      Kokkos::parallel_for(Kokkos::TeamThreadRange(
                             team_member, transfer_data.n_dofs_per_cell_coarse),
                           [&](const int &i) {
                             const unsigned int dof_index =
                               dof_indices_coarse(i, coarse_cell_index);
                             if (dof_index != numbers::invalid_unsigned_int)
                               Kokkos::atomic_add(&dst[dof_index],
                                                  values_coarse(i));
                           });
      team_member.team_barrier();
    }

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

        // 1) compute weights globally
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
        Assert((mg_level_fine == numbers::invalid_unsigned_int &&
                mg_level_coarse == numbers::invalid_unsigned_int),
               ExcNotImplemented());

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

        // nothing to do since for identity prolongation matrices a short-cut
        // code path is used during prolongation/restriction

        // Hasn't been tested with local mesh refinement yet!
        Assert(transfer.schemes[0].n_coarse_cells == 0, ExcNotImplemented());

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
        // }

        // ------------------- fill dof indices-------------------

        {
          unsigned int scheme_counter = 0;
          unsigned int cell_counter   = 0;
          for (auto &scheme : transfer.schemes)
            {
              if (scheme.n_coarse_cells == 0)
                continue;

              scheme.dof_indices_coarse =
                Kokkos::View<unsigned int **,
                             MemorySpace::Default::kokkos_space>(
                  Kokkos::view_alloc("h_transfer_dof_indices_coarse_scheme_" +
                                       std::to_string(scheme_counter),
                                     Kokkos::WithoutInitializing),
                  scheme.n_dofs_per_cell_coarse,
                  scheme.n_coarse_cells);

              scheme.dof_indices_fine =
                Kokkos::View<unsigned int **,
                             MemorySpace::Default::kokkos_space>(
                  Kokkos::view_alloc("h_transfer_dof_indices_fine_scheme_" +
                                       std::to_string(scheme_counter),
                                     Kokkos::WithoutInitializing),
                  scheme.n_dofs_per_cell_fine,
                  scheme.n_coarse_cells);

              auto dof_indices_coarse_host =
                Kokkos::create_mirror_view(scheme.dof_indices_coarse);

              auto dof_indices_fine_host =
                Kokkos::create_mirror_view(scheme.dof_indices_fine);

              const unsigned int first_cell = cell_counter;
              for (unsigned int cell = 0; cell < scheme.n_coarse_cells; ++cell)
                {
                  const unsigned int coarse_cell_index = first_cell + cell;

                  // fill coarse indices
                  {
                    const unsigned int *dof_indices_coarse =
                      transfer_cpu.constraint_info_coarse.dof_indices.data() +
                      transfer_cpu.constraint_info_coarse
                        .row_starts[coarse_cell_index]
                        .first;
                    unsigned int index_indicators =
                      transfer_cpu.constraint_info_coarse
                        .row_starts[coarse_cell_index]
                        .second;
                    unsigned int next_index_indicators =
                      transfer_cpu.constraint_info_coarse
                        .row_starts[coarse_cell_index + 1]
                        .second;

                    unsigned int ind_local = 0;
                    for (; index_indicators != next_index_indicators;
                         ++index_indicators)
                      {
                        const std::pair<unsigned short, unsigned short>
                          indicator = transfer_cpu.constraint_info_coarse
                                        .constraint_indicator[index_indicators];

                        for (unsigned int j = 0; j < indicator.first; ++j)
                          {
                            dof_indices_coarse_host(ind_local + j, cell) =
                              dof_indices_coarse[j];
                          }
                        ind_local += indicator.first;
                        dof_indices_coarse += indicator.first;

                        dof_indices_coarse_host(ind_local, cell) =
                          numbers::invalid_unsigned_int;

                        ++ind_local;
                      }

                    for (; ind_local < scheme.n_dofs_per_cell_coarse;
                         ++dof_indices_coarse, ++ind_local)
                      dof_indices_coarse_host(ind_local, cell) =
                        *dof_indices_coarse;
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

                  ++cell_counter;
                }

              // copy to the device memory
              {
                Kokkos::deep_copy(scheme.dof_indices_coarse,
                                  dof_indices_coarse_host);
                Kokkos::fence();

                Kokkos::deep_copy(scheme.dof_indices_fine,
                                  dof_indices_fine_host);
                Kokkos::fence();
              }

              ++scheme_counter;
            }
        }

        // ------------------- weights-------------------
        if (transfer.fine_element_is_continuous)
          setup_weights(constraints_fine, transfer, is_feq);
      }
    };

  } // namespace internal

  template <int dim, typename VectorType>
  void
  MGTwoLevelTransfer<dim, VectorType>::prolongate_and_add_internal(
    VectorType       &dst,
    const VectorType &src) const
  {
    using TeamPolicy =
      Kokkos::TeamPolicy<MemorySpace::Default::kokkos_space::execution_space>;

    using Functor = internal::CellProlongationKernel<dim, VectorType>;

    MemorySpace::Default::kokkos_space::execution_space exec;

    unsigned int scheme_index = 0;
    for (auto &scheme : schemes)
      {
        if (scheme.n_coarse_cells == 0)
          continue;

        Functor prolongator;

        auto team_policy =
          TeamPolicy(exec, scheme.n_coarse_cells, Kokkos::AUTO);

        internal::ApplyCellKernel<dim, VectorType, Functor> apply_prolongation(
          prolongator, scheme, src, dst);

        Kokkos::parallel_for("prolongate_and_add_h_transfer_scheme_" +
                               std::to_string(scheme_index),
                             team_policy,
                             apply_prolongation);
        ++scheme_index;
      }
  }

  template <int dim, typename VectorType>
  void
  MGTwoLevelTransfer<dim, VectorType>::restrict_and_add_internal(
    VectorType       &dst,
    const VectorType &src) const
  {
    using TeamPolicy =
      Kokkos::TeamPolicy<MemorySpace::Default::kokkos_space::execution_space>;

    using Functor = internal::CellRestrictionKernel<dim, VectorType>;

    MemorySpace::Default::kokkos_space::execution_space exec;

    unsigned int scheme_index = 0;
    for (auto &scheme : schemes)
      {
        if (scheme.n_coarse_cells == 0)
          continue;

        Functor restrictor;

        auto team_policy =
          TeamPolicy(exec, scheme.n_coarse_cells, Kokkos::AUTO);

        internal::ApplyCellKernel<dim, VectorType, Functor> apply_restriction(
          restrictor, scheme, src, dst);

        Kokkos::parallel_for("prolongate_and_add_h_transfer_scheme_" +
                               std::to_string(scheme_index),
                             team_policy,
                             apply_restriction);
        ++scheme_index;
      }
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
