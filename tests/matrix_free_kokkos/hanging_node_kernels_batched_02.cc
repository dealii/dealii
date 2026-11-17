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



// Test Portable::internal::batched::resolve_hanging_nodes() against the
// cell-wise kernels with one thread per degree of freedom, copied into this
// file as namespace reference like in hanging_node_kernels_02. Every
// combination of the ConstraintKinds bits (subcell, face and, in 3d, edge
// flags) is one cell with random values, the 1d interpolation weights are
// random as well. The reference kernels resolve one cell per team, like in
// Portable::FEEvaluation; the batched kernel resolves batches of 1, 3 and 8
// cells per team, where batches of 3 and 8 cells mix constrained and
// unconstrained cells and the last batch is shorter. The checksums are the
// same as in hanging_node_kernels_02.

#include <deal.II/base/memory_space.h>
#include <deal.II/base/utilities.h>

#include <deal.II/matrix_free/hanging_nodes_internal.h>
#include <deal.II/matrix_free/portable_hanging_nodes_internal.h>

#include <Kokkos_Core.hpp>

#include <fstream>
#include <iostream>
#include <vector>

#include "../tests.h"


using ConstraintKinds = dealii::internal::MatrixFreeFunctions::ConstraintKinds;



// The kernels Portable::internal::resolve_hanging_nodes() used before: per
// direction, one thread per degree of freedom computes the constrained values
// into a temporary array, and a second pass copies them back.
namespace reference
{
  using Portable::internal::index2;
  using Portable::internal::index3;
  using Portable::internal::is_constrained_dof_2d;
  using Portable::internal::is_constrained_dof_3d;

  template <unsigned int fe_degree,
            unsigned int direction,
            bool         transpose,
            typename Number,
            typename ViewType>
  DEAL_II_HOST_DEVICE inline void
  interpolate_boundary_2d(
    const Kokkos::TeamPolicy<
      MemorySpace::Default::kokkos_space::execution_space>::member_type
      &team_member,
    Kokkos::View<Number *, MemorySpace::Default::kokkos_space>
      constraint_weights,
    const dealii::internal::MatrixFreeFunctions::ConstraintKinds
            &constraint_mask,
    ViewType values)
  {
    constexpr unsigned int n_q_points_1d = fe_degree + 1;
    constexpr unsigned int n_q_points    = Utilities::pow(n_q_points_1d, 2);

    // Flag is true if dof is constrained for the given direction and the
    // given face.
    const bool constrained_face =
      (constraint_mask &
       (((direction == 0) ?
           dealii::internal::MatrixFreeFunctions::ConstraintKinds::face_y :
           dealii::internal::MatrixFreeFunctions::ConstraintKinds::
             unconstrained) |
        ((direction == 1) ?
           dealii::internal::MatrixFreeFunctions::ConstraintKinds::face_x :
           dealii::internal::MatrixFreeFunctions::ConstraintKinds::
             unconstrained))) !=
      dealii::internal::MatrixFreeFunctions::ConstraintKinds::unconstrained;

    Number tmp[n_q_points] = {};
    Kokkos::parallel_for(
      Kokkos::TeamThreadRange(team_member, n_q_points),
      [&](const int &q_point) {
        const unsigned int x_idx = q_point % n_q_points_1d;
        const unsigned int y_idx = q_point / n_q_points_1d;

        const auto this_type =
          (direction == 0) ?
            dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_x :
            dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_y;

        const unsigned int interp_idx = (direction == 0) ? x_idx : y_idx;

        // Flag is true if for the given direction, the dof is constrained
        // with the right type and is on the correct side (left (= 0) or right
        // (= fe_degree))
        const bool constrained_dof =
          is_constrained_dof_2d<fe_degree, direction>(constraint_mask,
                                                      x_idx,
                                                      y_idx);

        if (constrained_face && constrained_dof)
          {
            Number     sum  = 0.;
            const bool type = (constraint_mask & this_type) !=
                              dealii::internal::MatrixFreeFunctions::
                                ConstraintKinds::unconstrained;

            if (type)
              {
                for (unsigned int i = 0; i <= fe_degree; ++i)
                  {
                    const unsigned int real_idx =
                      (direction == 0) ? index2<n_q_points_1d>(i, y_idx) :
                                         index2<n_q_points_1d>(x_idx, i);

                    const Number w =
                      transpose ?
                        constraint_weights[i * n_q_points_1d + interp_idx] :
                        constraint_weights[interp_idx * n_q_points_1d + i];
                    sum += w * values[real_idx];
                  }
              }
            else
              {
                for (unsigned int i = 0; i <= fe_degree; ++i)
                  {
                    const unsigned int real_idx =
                      (direction == 0) ? index2<n_q_points_1d>(i, y_idx) :
                                         index2<n_q_points_1d>(x_idx, i);

                    const Number w =
                      transpose ?
                        constraint_weights[(fe_degree - i) * n_q_points_1d +
                                           fe_degree - interp_idx] :
                        constraint_weights[(fe_degree - interp_idx) *
                                             n_q_points_1d +
                                           fe_degree - i];
                    sum += w * values[real_idx];
                  }
              }
            tmp[q_point] = sum;
          }
      });

    // The synchronization is done for all the threads in one team with
    // each team being assigned to one element.
    team_member.team_barrier();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, n_q_points),
                         [&](const int &q_point) {
                           const unsigned int x_idx = q_point % n_q_points_1d;
                           const unsigned int y_idx = q_point / n_q_points_1d;
                           const bool         constrained_dof =
                             is_constrained_dof_2d<fe_degree, direction>(
                               constraint_mask, x_idx, y_idx);
                           if (constrained_face && constrained_dof)
                             values[index2<fe_degree + 1>(x_idx, y_idx)] =
                               tmp[q_point];
                         });

    team_member.team_barrier();
  }



  template <unsigned int fe_degree,
            unsigned int direction,
            bool         transpose,
            typename Number,
            typename ViewType>
  DEAL_II_HOST_DEVICE inline void
  interpolate_boundary_3d(
    const Kokkos::TeamPolicy<
      MemorySpace::Default::kokkos_space::execution_space>::member_type
      &team_member,
    Kokkos::View<Number *, MemorySpace::Default::kokkos_space>
      constraint_weights,
    const dealii::internal::MatrixFreeFunctions::ConstraintKinds
             constraint_mask,
    ViewType values)
  {
    constexpr unsigned int n_q_points_1d = fe_degree + 1;
    constexpr unsigned int n_q_points    = Utilities::pow(n_q_points_1d, 3);

    const auto this_type =
      (direction == 0) ?
        dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_x :
      (direction == 1) ?
        dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_y :
        dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_z;
    const auto face1_type =
      (direction == 0) ?
        dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_y :
      (direction == 1) ?
        dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_z :
        dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_x;
    const auto face2_type =
      (direction == 0) ?
        dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_z :
      (direction == 1) ?
        dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_x :
        dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_y;

    // If computing in x-direction, need to match against face_y or
    // face_z
    const auto face1 =
      (direction == 0) ?
        dealii::internal::MatrixFreeFunctions::ConstraintKinds::face_y :
      (direction == 1) ?
        dealii::internal::MatrixFreeFunctions::ConstraintKinds::face_z :
        dealii::internal::MatrixFreeFunctions::ConstraintKinds::face_x;
    const auto face2 =
      (direction == 0) ?
        dealii::internal::MatrixFreeFunctions::ConstraintKinds::face_z :
      (direction == 1) ?
        dealii::internal::MatrixFreeFunctions::ConstraintKinds::face_x :
        dealii::internal::MatrixFreeFunctions::ConstraintKinds::face_y;
    const auto edge =
      (direction == 0) ?
        dealii::internal::MatrixFreeFunctions::ConstraintKinds::edge_x :
      (direction == 1) ?
        dealii::internal::MatrixFreeFunctions::ConstraintKinds::edge_y :
        dealii::internal::MatrixFreeFunctions::ConstraintKinds::edge_z;
    const auto constrained_face = constraint_mask & (face1 | face2 | edge);

    Number tmp[n_q_points] = {};
    Kokkos::parallel_for(
      Kokkos::TeamThreadRange(team_member, n_q_points),
      [&](const int &q_point) {
        const unsigned int x_idx = q_point % n_q_points_1d;
        const unsigned int y_idx = (q_point / n_q_points_1d) % n_q_points_1d;
        const unsigned int z_idx = q_point / (n_q_points_1d * n_q_points_1d);

        const unsigned int interp_idx = (direction == 0) ? x_idx :
                                        (direction == 1) ? y_idx :
                                                           z_idx;
        const bool         constrained_dof =
          is_constrained_dof_3d<fe_degree, direction>(constraint_mask,
                                                      x_idx,
                                                      y_idx,
                                                      z_idx,
                                                      face1_type,
                                                      face2_type,
                                                      face1,
                                                      face2,
                                                      edge);
        if ((constrained_face != dealii::internal::MatrixFreeFunctions::
                                   ConstraintKinds::unconstrained) &&
            constrained_dof)
          {
            Number     sum  = 0.;
            const bool type = (constraint_mask & this_type) !=
                              dealii::internal::MatrixFreeFunctions::
                                ConstraintKinds::unconstrained;
            if (type)
              {
                for (unsigned int i = 0; i <= fe_degree; ++i)
                  {
                    const unsigned int real_idx =
                      (direction == 0) ?
                        index3<fe_degree + 1>(i, y_idx, z_idx) :
                      (direction == 1) ?
                        index3<fe_degree + 1>(x_idx, i, z_idx) :
                        index3<fe_degree + 1>(x_idx, y_idx, i);

                    const Number w =
                      transpose ?
                        constraint_weights[i * n_q_points_1d + interp_idx] :
                        constraint_weights[interp_idx * n_q_points_1d + i];
                    sum += w * values[real_idx];
                  }
              }
            else
              {
                for (unsigned int i = 0; i <= fe_degree; ++i)
                  {
                    const unsigned int real_idx =
                      (direction == 0) ?
                        index3<n_q_points_1d>(i, y_idx, z_idx) :
                      (direction == 1) ?
                        index3<n_q_points_1d>(x_idx, i, z_idx) :
                        index3<n_q_points_1d>(x_idx, y_idx, i);

                    const Number w =
                      transpose ?
                        constraint_weights[(fe_degree - i) * n_q_points_1d +
                                           fe_degree - interp_idx] :
                        constraint_weights[(fe_degree - interp_idx) *
                                             n_q_points_1d +
                                           fe_degree - i];
                    sum += w * values[real_idx];
                  }
              }
            tmp[q_point] = sum;
          }
      });

    // The synchronization is done for all the threads in one team with
    // each team being assigned to one element.
    team_member.team_barrier();

    Kokkos::parallel_for(
      Kokkos::TeamThreadRange(team_member, n_q_points),
      [&](const int &q_point) {
        const unsigned int x_idx = q_point % n_q_points_1d;
        const unsigned int y_idx = (q_point / n_q_points_1d) % n_q_points_1d;
        const unsigned int z_idx = q_point / (n_q_points_1d * n_q_points_1d);
        const bool         constrained_dof =
          is_constrained_dof_3d<fe_degree, direction>(constraint_mask,
                                                      x_idx,
                                                      y_idx,
                                                      z_idx,
                                                      face1_type,
                                                      face2_type,
                                                      face1,
                                                      face2,
                                                      edge);
        if ((constrained_face != dealii::internal::MatrixFreeFunctions::
                                   ConstraintKinds::unconstrained) &&
            constrained_dof)
          values[index3<fe_degree + 1>(x_idx, y_idx, z_idx)] = tmp[q_point];
      });

    team_member.team_barrier();
  }



  /**
   * This function resolves the hanging nodes using tensor product.
   *
   * The implementation of this class is explained in Section 3 of
   * @cite ljungkvist2017matrix and in Section 3.4 of
   * @cite kronbichler2019multigrid.
   */
  template <int  dim,
            int  fe_degree,
            bool transpose,
            typename Number,
            typename ViewType>
  DEAL_II_HOST_DEVICE void
  resolve_hanging_nodes(
    const Kokkos::TeamPolicy<
      MemorySpace::Default::kokkos_space::execution_space>::member_type
      &team_member,
    Kokkos::View<Number *, MemorySpace::Default::kokkos_space>
      constraint_weights,
    const dealii::internal::MatrixFreeFunctions::ConstraintKinds
             constraint_mask,
    ViewType values)
  {
    if constexpr (dim == 2)
      {
        interpolate_boundary_2d<fe_degree, 0, transpose>(team_member,
                                                         constraint_weights,
                                                         constraint_mask,
                                                         values);

        interpolate_boundary_2d<fe_degree, 1, transpose>(team_member,
                                                         constraint_weights,
                                                         constraint_mask,
                                                         values);
      }
    else if constexpr (dim == 3)
      {
        // Interpolate y and z faces (x-direction)
        interpolate_boundary_3d<fe_degree, 0, transpose>(team_member,
                                                         constraint_weights,
                                                         constraint_mask,
                                                         values);
        // Interpolate x and z faces (y-direction)
        interpolate_boundary_3d<fe_degree, 1, transpose>(team_member,
                                                         constraint_weights,
                                                         constraint_mask,
                                                         values);
        // Interpolate x and y faces (z-direction)
        interpolate_boundary_3d<fe_degree, 2, transpose>(team_member,
                                                         constraint_weights,
                                                         constraint_mask,
                                                         values);
      }
  }
} // namespace reference



using TeamHandle = Kokkos::TeamPolicy<
  MemorySpace::Default::kokkos_space::execution_space>::member_type;

template <typename T>
using DeviceView = Kokkos::View<T, MemorySpace::Default::kokkos_space>;

using ScratchView = Kokkos::View<
  double *,
  MemorySpace::Default::kokkos_space::execution_space::scratch_memory_space,
  Kokkos::MemoryTraits<Kokkos::Unmanaged>>;



// One team resolves the hanging nodes of the cells [first_cell, first_cell +
// n_cells_in_batch), staged cell after cell in shared memory, either with the
// batched kernel, which gets the masks of the batch as a subview of masks,
// or, for cells_per_team == 1 only, with the reference kernels.
// Both are called also for unconstrained cells, like in
// MatrixFreeTools::compute_diagonal().
template <int dim, int fe_degree, bool transpose, bool batched>
struct ResolveBatch
{
  static constexpr int n_dofs = Utilities::pow(fe_degree + 1, dim);

  DeviceView<double *>          weights;
  DeviceView<ConstraintKinds *> masks;
  DeviceView<double *>          src;
  DeviceView<double *>          dst;
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

    ScratchView values(team_member.team_shmem(), cells_per_team * n_dofs);

    Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, n_values),
                         [&](const int i) { values(i) = src(offset + i); });
    team_member.team_barrier();

    if constexpr (batched)
      Portable::internal::batched::resolve_hanging_nodes<dim,
                                                         fe_degree,
                                                         transpose>(
        team_member,
        weights,
        Kokkos::subview(
          masks, Kokkos::make_pair(first_cell, first_cell + n_cells_in_batch)),
        values,
        n_cells_in_batch);
    else
      reference::resolve_hanging_nodes<dim, fe_degree, transpose>(
        team_member, weights, masks(first_cell), values);

    Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, n_values),
                         [&](const int i) { dst(offset + i) = values(i); });
  }
};



template <int dim, int fe_degree, bool transpose, bool batched>
Kokkos::View<double *, Kokkos::HostSpace>
run(const DeviceView<double *>          &weights,
    const DeviceView<ConstraintKinds *> &masks,
    const DeviceView<double *>          &src,
    const int                            cells_per_team)
{
  constexpr int n_dofs  = Utilities::pow(fe_degree + 1, dim);
  const int     n_cells = masks.size();
  const int     n_teams = (n_cells + cells_per_team - 1) / cells_per_team;

  DeviceView<double *> dst("dst", n_cells * n_dofs);

  MemorySpace::Default::kokkos_space::execution_space exec;
  Kokkos::TeamPolicy<MemorySpace::Default::kokkos_space::execution_space>
    team_policy(exec, n_teams, Kokkos::AUTO);
  team_policy.set_scratch_size(
    0, Kokkos::PerTeam(ScratchView::shmem_size(cells_per_team * n_dofs)));

  Kokkos::parallel_for(team_policy,
                       ResolveBatch<dim, fe_degree, transpose, batched>{
                         weights, masks, src, dst, n_cells, cells_per_team});

  return Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, dst);
}



template <int dim, int fe_degree, bool transpose>
void
test()
{
  constexpr int n_dofs_1d = fe_degree + 1;
  constexpr int n_dofs    = Utilities::pow(n_dofs_1d, dim);

  // all combinations of the subcell, face and (3d) edge bits
  constexpr unsigned int n_bits    = (dim == 2) ? 4 : 9;
  constexpr unsigned int n_cells   = 1u << n_bits;
  const auto             make_mask = [](const unsigned int bits) {
    unsigned int mask = 0;
    if (dim == 2)
      {
        // bits: subcell_x, subcell_y, face_x, face_y
        mask |= (bits & 1) ?
                  static_cast<unsigned int>(ConstraintKinds::subcell_x) :
                  0;
        mask |= (bits & 2) ?
                  static_cast<unsigned int>(ConstraintKinds::subcell_y) :
                  0;
        mask |=
          (bits & 4) ? static_cast<unsigned int>(ConstraintKinds::face_x) : 0;
        mask |=
          (bits & 8) ? static_cast<unsigned int>(ConstraintKinds::face_y) : 0;
      }
    else
      // the nine bits of ConstraintKinds in order
      mask = bits;
    return static_cast<ConstraintKinds>(mask);
  };

  deallog << "dim=" << dim << " degree=" << fe_degree
          << " transpose=" << transpose << ": n_cells=" << n_cells << std::endl;

  DeviceView<double *>          weights("weights", n_dofs_1d * n_dofs_1d);
  DeviceView<ConstraintKinds *> masks("masks", n_cells);
  DeviceView<double *>          src("src", n_cells * n_dofs);
  {
    auto weights_host = Kokkos::create_mirror_view(weights);
    for (unsigned int i = 0; i < weights_host.size(); ++i)
      weights_host(i) =
        -1. + 2. * static_cast<double>(Testing::rand()) / RAND_MAX;
    Kokkos::deep_copy(weights, weights_host);

    auto masks_host = Kokkos::create_mirror_view(masks);
    for (unsigned int cell = 0; cell < n_cells; ++cell)
      masks_host(cell) = make_mask(cell);
    Kokkos::deep_copy(masks, masks_host);

    auto src_host = Kokkos::create_mirror_view(src);
    for (unsigned int i = 0; i < src_host.size(); ++i)
      src_host(i) = -1. + 2. * static_cast<double>(Testing::rand()) / RAND_MAX;
    Kokkos::deep_copy(src, src_host);
  }

  const auto reference_values =
    run<dim, fe_degree, transpose, false>(weights, masks, src, 1);

  double sum = 0., sum_of_squares = 0.;
  for (unsigned int i = 0; i < n_cells * n_dofs; ++i)
    {
      sum += reference_values(i);
      sum_of_squares += reference_values(i) * reference_values(i);
    }
  deallog << "  checksum reference: " << sum << " " << sum_of_squares
          << std::endl;

  for (const int cells_per_team : {1, 3, 8})
    {
      const auto batched = run<dim, fe_degree, transpose, true>(weights,
                                                                masks,
                                                                src,
                                                                cells_per_team);

      double max_difference = 0.;
      for (unsigned int i = 0; i < n_cells * n_dofs; ++i)
        max_difference =
          std::max(max_difference, std::abs(batched(i) - reference_values(i)));
      deallog << "  batches of " << cells_per_team
              << " cells, max |batched - reference|: " << max_difference
              << std::endl;

      // the smallest case in full
      if (dim == 2 && fe_degree == 1 && cells_per_team == 3)
        for (unsigned int cell = 0; cell < n_cells; ++cell)
          {
            deallog << "  mask " << static_cast<unsigned int>(make_mask(cell))
                    << ":";
            for (unsigned int i = 0; i < n_dofs; ++i)
              deallog << " " << batched(cell * n_dofs + i);
            deallog << std::endl;
          }
    }
}



int
main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog << std::setprecision(12);

  Kokkos::initialize();
  {
    test<2, 1, false>();
    test<2, 1, true>();
    test<2, 2, false>();
    test<2, 2, true>();
    test<2, 3, false>();
    test<2, 3, true>();

    test<3, 1, false>();
    test<3, 1, true>();
    test<3, 2, false>();
    test<3, 2, true>();
    test<3, 3, false>();
    test<3, 3, true>();
  }
  Kokkos::finalize();

  return 0;
}
