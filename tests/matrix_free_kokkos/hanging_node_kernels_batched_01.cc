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



// Compare the result of Portable::internal::batched::resolve_hanging_nodes()
// with FEEvaluationImplHangingNodes. Same as
// tests/matrix_free/hanging_node_kernels_01.cc, but all constraint masks of
// one dimension and degree are resolved on the device in a single kernel
// launch, with the masks distributed as cells over the batches of several
// teams, once with batches of four cells and once with a single cell per
// team. The printed values are the ones of the batches of four cells; the
// ones of a single cell per team must be identical. The interpolation
// weights are the ones Portable::MatrixFree uses, staged in team scratch
// memory.

#include <deal.II/base/memory_space.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/matrix_free/evaluation_kernels_hanging_nodes.h>
#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/portable_hanging_nodes_internal.h>

#include <deal.II/numerics/vector_tools.h>

#include <Kokkos_Core.hpp>

#include <vector>

#include "../tests.h"


using ConstraintKinds = dealii::internal::MatrixFreeFunctions::ConstraintKinds;

using TeamHandle = Kokkos::TeamPolicy<
  MemorySpace::Default::kokkos_space::execution_space>::member_type;

template <typename T>
using DeviceView = Kokkos::View<T, MemorySpace::Default::kokkos_space>;

using ScratchView = Kokkos::View<
  double *,
  MemorySpace::Default::kokkos_space::execution_space::scratch_memory_space,
  Kokkos::MemoryTraits<Kokkos::Unmanaged>>;



// One team resolves the hanging nodes of the cells [first, first + n) of its
// batch, staged cell after cell in shared memory. The weights are copied to
// shared memory together with the values, before the same barrier.
template <int dim, int fe_degree, bool transpose>
struct ResolveBatch
{
  static constexpr int n_dofs = Utilities::pow(fe_degree + 1, dim);

  DeviceView<double *>          weights;
  DeviceView<ConstraintKinds *> masks;
  DeviceView<double *>          values;
  int                           n_cells;
  int                           cells_per_team;

  DEAL_II_HOST_DEVICE void
  operator()(const TeamHandle &team_member) const
  {
    const int first_cell       = team_member.league_rank() * cells_per_team;
    const int n_cells_in_batch = (n_cells - first_cell < cells_per_team) ?
                                   (n_cells - first_cell) :
                                   cells_per_team;
    const int n_values         = n_cells_in_batch * n_dofs;
    const int offset           = first_cell * n_dofs;

    ScratchView shared_values(team_member.team_shmem(),
                              cells_per_team * n_dofs);
    const Portable::internal::batched::ConstraintWeightsView<double>
      shared_weights(team_member.team_shmem(), weights.size());

    Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, n_values),
                         [&](const int i) {
                           shared_values(i) = values(offset + i);
                         });
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, weights.size()),
                         [&](const int i) { shared_weights(i) = weights(i); });
    team_member.team_barrier();

    Portable::internal::batched::
      resolve_hanging_nodes<dim, fe_degree, transpose>(
        team_member,
        shared_weights,
        Kokkos::subview(
          masks, Kokkos::make_pair(first_cell, first_cell + n_cells_in_batch)),
        shared_values,
        n_cells_in_batch);

    Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, n_values),
                         [&](const int i) {
                           values(offset + i) = shared_values(i);
                         });
  }
};



// Resolve the hanging nodes of all cells on the device, every cell starting
// from the values 0, 1, 2, ... like in hanging_node_kernels_01.
template <int dim, int fe_degree, bool transpose>
std::vector<double>
run_batched(const AlignedVector<double>        &subface_interpolation_matrix,
            const std::vector<ConstraintKinds> &cell_masks,
            const int                           cells_per_team)
{
  constexpr int n_dofs  = Utilities::pow(fe_degree + 1, dim);
  const int     n_cells = cell_masks.size();

  DeviceView<double *> weights("weights", subface_interpolation_matrix.size());
  DeviceView<ConstraintKinds *> masks("masks", n_cells);
  DeviceView<double *>          values("values", n_cells * n_dofs);

  auto weights_host = Kokkos::create_mirror_view(weights);
  for (unsigned int i = 0; i < subface_interpolation_matrix.size(); ++i)
    weights_host(i) = subface_interpolation_matrix[i];
  Kokkos::deep_copy(weights, weights_host);

  auto masks_host = Kokkos::create_mirror_view(masks);
  for (int cell = 0; cell < n_cells; ++cell)
    masks_host(cell) = cell_masks[cell];
  Kokkos::deep_copy(masks, masks_host);

  auto values_host = Kokkos::create_mirror_view(values);
  for (int cell = 0; cell < n_cells; ++cell)
    for (int i = 0; i < n_dofs; ++i)
      values_host(cell * n_dofs + i) = i;
  Kokkos::deep_copy(values, values_host);

  const int n_teams = (n_cells + cells_per_team - 1) / cells_per_team;

  MemorySpace::Default::kokkos_space::execution_space exec;
  Kokkos::TeamPolicy<MemorySpace::Default::kokkos_space::execution_space>
    team_policy(exec, n_teams, Kokkos::AUTO);
  team_policy.set_scratch_size(
    0,
    Kokkos::PerTeam(ScratchView::shmem_size(cells_per_team * n_dofs) +
                    ScratchView::shmem_size(weights.size())));

  Kokkos::parallel_for(team_policy,
                       ResolveBatch<dim, fe_degree, transpose>{
                         weights, masks, values, n_cells, cells_per_team});

  Kokkos::deep_copy(values_host, values);

  return std::vector<double>(values_host.data(),
                             values_host.data() + values_host.size());
}



// The masks come in groups; in 2d, an empty line follows each group like in
// hanging_node_kernels_01.
template <int dim, int fe_degree>
void
test(const std::vector<std::vector<ConstraintKinds>> &groups)
{
  Triangulation<dim> tria;
  GridGenerator::subdivided_hyper_cube(tria, 2);
  tria.begin()->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  QGauss<dim>    quad(fe_degree + 1);
  FE_Q<dim>      fe(fe_degree);
  MappingQ1<dim> mapping;

  DoFHandler<dim> dof_handler;
  dof_handler.reinit(tria);
  dof_handler.distribute_dofs(fe);

  AffineConstraints<double> constraints;

  typename MatrixFree<dim, double, VectorizedArray<double>>::AdditionalData
    additional_data;
  additional_data.mapping_update_flags = update_values | update_gradients |
                                         update_JxW_values |
                                         dealii::update_quadrature_points;

  MatrixFree<dim, double, VectorizedArray<double>> matrix_free;
  matrix_free.reinit(mapping, dof_handler, constraints, quad, additional_data);

  FEEvaluation<dim, -1, 0, 1, double> eval(matrix_free);
  eval.reinit(0);

  // the weights Portable::MatrixFree copies to the device
  const AlignedVector<double> &subface_interpolation_matrix =
    eval.get_shape_info().data.front().subface_interpolation_matrices[0];

  std::vector<ConstraintKinds> cell_masks;
  for (const auto &group : groups)
    cell_masks.insert(cell_masks.end(), group.begin(), group.end());

  // batches of 4 cells, where the last team gets a shorter batch, and of a
  // single cell, which is how the cell-wise callers would use the kernel
  const auto run = [&](const int cells_per_team) {
    return std::array<std::vector<double>, 2>{
      {run_batched<dim, fe_degree, false>(subface_interpolation_matrix,
                                          cell_masks,
                                          cells_per_team),
       run_batched<dim, fe_degree, true>(subface_interpolation_matrix,
                                         cell_masks,
                                         cells_per_team)}};
  };
  const auto batched     = run(4);
  const auto single_cell = run(1);

  const unsigned int n_dofs = fe.n_dofs_per_cell();
  unsigned int       cell   = 0;
  for (const auto &group : groups)
    {
      for (const auto mask : group)
        {
          std::array<
            dealii::internal::MatrixFreeFunctions::compressed_constraint_kind,
            VectorizedArray<double>::size()>
            cmask;
          std::fill(cmask.begin(),
                    cmask.end(),
                    dealii::internal::MatrixFreeFunctions::
                      unconstrained_compressed_constraint_kind);
          cmask[0] = dealii::internal::MatrixFreeFunctions::compress(mask, dim);

          for (unsigned int b = 0; b < 2; ++b)
            {
              AlignedVector<VectorizedArray<double>> values(n_dofs);
              for (unsigned int i = 0; i < n_dofs; ++i)
                values[i][0] = i;

              for (const auto i : values)
                deallog << i[0] << ' ';
              deallog << std::endl;

              internal::FEEvaluationImplHangingNodes<dim,
                                                     VectorizedArray<double>>::
                template run<-1>(
                  1, eval.get_shape_info(), b == 1, cmask, values.data());

              for (const auto i : values)
                deallog << i[0] << ' ';
              deallog << std::endl;

              for (unsigned int i = 0; i < n_dofs; ++i)
                deallog << batched[b][cell * n_dofs + i] << ' ';
              deallog << std::endl;
              deallog << std::endl;

              for (unsigned int i = 0; i < n_dofs; ++i)
                {
                  AssertThrow(std::abs(values[i][0] -
                                       batched[b][cell * n_dofs + i]) < 1e-10,
                              ExcInternalError());
                  AssertThrow(single_cell[b][cell * n_dofs + i] ==
                                batched[b][cell * n_dofs + i],
                              ExcInternalError());
                }
            }
          ++cell;
        }
      if (dim == 2)
        deallog << std::endl;
    }
}



template <int fe_degree>
void
test_2d()
{
  test<2, fe_degree>(
    {{ConstraintKinds::unconstrained},
     {ConstraintKinds::face_x | ConstraintKinds::subcell_x |
        ConstraintKinds::subcell_y,                         // face 0/0
      ConstraintKinds::face_x | ConstraintKinds::subcell_x, // face 0/1
      ConstraintKinds::face_x | ConstraintKinds::subcell_y, // face 1/0
      ConstraintKinds::face_x},                             // face 1/1
     {ConstraintKinds::face_y | ConstraintKinds::subcell_y |
        ConstraintKinds::subcell_x,                         // face 2/0
      ConstraintKinds::face_y | ConstraintKinds::subcell_y, // face 2/1
      ConstraintKinds::face_y | ConstraintKinds::subcell_x, // face 3/0
      ConstraintKinds::face_y}});                           // face 3/1
}



template <int fe_degree>
void
test_3d()
{
  using CK = ConstraintKinds;
  test<3, fe_degree>(
    {{// edge 2
      CK::edge_x | CK::subcell_y | CK::subcell_z,
      CK::edge_x | CK::subcell_y | CK::subcell_z | CK::subcell_x,
      // edge 3
      CK::edge_x | CK::subcell_z,
      CK::edge_x | CK::subcell_z | CK::subcell_x,
      // edge 6
      CK::edge_x | CK::subcell_y,
      CK::edge_x | CK::subcell_y | CK::subcell_x,
      // edge 7
      CK::edge_x,
      CK::edge_x | CK::subcell_x,

      // edge 0
      CK::edge_y | CK::subcell_x | CK::subcell_z,
      CK::edge_y | CK::subcell_x | CK::subcell_z | CK::subcell_y,
      // edge 1
      CK::edge_y | CK::subcell_z,
      CK::edge_y | CK::subcell_z | CK::subcell_y,
      // edge 4
      CK::edge_y | CK::subcell_x,
      CK::edge_y | CK::subcell_x | CK::subcell_y,
      // edge 5
      CK::edge_y,
      CK::edge_y | CK::subcell_y,

      // edge 8
      CK::edge_z | CK::subcell_x | CK::subcell_y,
      CK::edge_z | CK::subcell_x | CK::subcell_y | CK::subcell_z,
      // edge 9
      CK::edge_z | CK::subcell_y,
      CK::edge_z | CK::subcell_y | CK::subcell_z,
      // edge 10
      CK::edge_z | CK::subcell_x,
      CK::edge_z | CK::subcell_x | CK::subcell_z,
      // edge 11
      CK::edge_z,
      CK::edge_z | CK::subcell_z,

      // face 0
      CK::face_x | CK::subcell_x,
      // face 1
      CK::face_x,
      // face 2
      CK::face_y | CK::subcell_y,
      // face 3
      CK::face_y,
      // face 4
      CK::face_z | CK::subcell_z,
      // face 5
      CK::face_z}});
}



int
main()
{
  initlog();

  Kokkos::initialize();
  {
    test_2d<1>();
    test_2d<2>();
    test_2d<3>();

    test_3d<1>();
    test_3d<2>();
    test_3d<3>();
  }
  Kokkos::finalize();

  return 0;
}
