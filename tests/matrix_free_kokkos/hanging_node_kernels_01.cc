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



// Compare the result of Portable::internal::resolve_hanging_nodes() with
// FEEvaluationImplHangingNodes. Same as
// tests/matrix_free/hanging_node_kernels_01.cc, but all constraint masks of
// one dimension and degree are resolved on the device in a single kernel
// launch with one cell per team, like Portable::FEEvaluation does. The
// interpolation weights are the ones Portable::MatrixFree uses.

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



// One team resolves the hanging nodes of one cell, staged in shared memory,
// like Portable::FEEvaluation::read_dof_values() and
// distribute_local_to_global() do.
template <int dim, int fe_degree, bool transpose>
struct ResolveCell
{
  static constexpr int n_dofs = Utilities::pow(fe_degree + 1, dim);

  DeviceView<double *>          weights;
  DeviceView<ConstraintKinds *> masks;
  DeviceView<double *>          values;

  DEAL_II_HOST_DEVICE void
  operator()(const TeamHandle &team_member) const
  {
    const int cell   = team_member.league_rank();
    const int offset = cell * n_dofs;

    ScratchView shared_values(team_member.team_shmem(), n_dofs);

    Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, n_dofs),
                         [&](const int i) {
                           shared_values(i) = values(offset + i);
                         });
    team_member.team_barrier();

    if (masks(cell) != ConstraintKinds::unconstrained)
      Portable::internal::resolve_hanging_nodes<dim, fe_degree, transpose>(
        team_member, weights, masks(cell), shared_values);

    Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, n_dofs),
                         [&](const int i) {
                           values(offset + i) = shared_values(i);
                         });
  }
};



// Resolve the hanging nodes of all cells on the device, every cell starting
// from the values 0, 1, 2, ... like in hanging_node_kernels_01.
template <int dim, int fe_degree, bool transpose>
std::vector<double>
run_on_device(const AlignedVector<double>        &subface_interpolation_matrix,
              const std::vector<ConstraintKinds> &cell_masks)
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

  MemorySpace::Default::kokkos_space::execution_space exec;
  Kokkos::TeamPolicy<MemorySpace::Default::kokkos_space::execution_space>
    team_policy(exec, n_cells, Kokkos::AUTO);
  team_policy.set_scratch_size(
    0, Kokkos::PerTeam(ScratchView::shmem_size(n_dofs)));

  Kokkos::parallel_for(team_policy,
                       ResolveCell<dim, fe_degree, transpose>{weights,
                                                              masks,
                                                              values});

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

  const std::array<std::vector<double>, 2> device = {
    {run_on_device<dim, fe_degree, false>(subface_interpolation_matrix,
                                          cell_masks),
     run_on_device<dim, fe_degree, true>(subface_interpolation_matrix,
                                         cell_masks)}};

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
                deallog << device[b][cell * n_dofs + i] << ' ';
              deallog << std::endl;
              deallog << std::endl;

              for (unsigned int i = 0; i < n_dofs; ++i)
                AssertThrow(std::abs(values[i][0] -
                                     device[b][cell * n_dofs + i]) < 1e-10,
                            ExcInternalError());
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
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  test_2d<1>();
  test_2d<2>();
  test_2d<3>();

  test_3d<1>();
  test_3d<2>();
  test_3d<3>();
}
