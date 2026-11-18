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


// Test Portable::internal::resolve_hanging_nodes(), which resolves the
// hanging-node constraints with one thread per constrained line of degrees
// of freedom, against the kernels it replaced (one thread per degree of
// freedom, copied into this file as namespace reference). Every combination
// of the ConstraintKinds bits (subcell, face and, in 3d, edge flags) is one
// cell with random values, resolved by one team like in
// Portable::FEEvaluation; the 1d interpolation weights are random as well.
// The expected output was recorded with the reference kernels.

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



// One team resolves the hanging nodes of one cell twice, with
// reference::resolve_hanging_nodes() and with
// Portable::internal::resolve_hanging_nodes(). Both are called also for
// unconstrained cells, like in MatrixFreeTools::compute_diagonal().
template <int dim, int fe_degree, bool transpose>
struct ResolveCell
{
  static constexpr int n_dofs = Utilities::pow(fe_degree + 1, dim);

  DeviceView<double *>          weights;
  DeviceView<ConstraintKinds *> masks;
  DeviceView<double *>          src;
  DeviceView<double *>          dst_reference;
  DeviceView<double *>          dst;

  DEAL_II_HOST_DEVICE void
  operator()(const TeamHandle &team_member) const
  {
    const int cell   = team_member.league_rank();
    const int offset = cell * n_dofs;

    ScratchView values(team_member.team_shmem(), n_dofs);

    const auto copy = [&](const auto from,
                          auto       to,
                          const int  offset_from,
                          const int  offset_to) {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, n_dofs),
                           [&](const int i) {
                             to(offset_to + i) = from(offset_from + i);
                           });
      team_member.team_barrier();
    };

    copy(src, values, offset, 0);
    reference::resolve_hanging_nodes<dim, fe_degree, transpose>(team_member,
                                                                weights,
                                                                masks(cell),
                                                                values);
    copy(values, dst_reference, 0, offset);

    copy(src, values, offset, 0);
    Portable::internal::resolve_hanging_nodes<dim, fe_degree, transpose>(
      team_member, weights, masks(cell), values);
    copy(values, dst, 0, offset);
  }
};



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
  DeviceView<double *> dst_reference("dst_reference", n_cells * n_dofs);
  DeviceView<double *> dst("dst", n_cells * n_dofs);

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

  MemorySpace::Default::kokkos_space::execution_space exec;
  Kokkos::TeamPolicy<MemorySpace::Default::kokkos_space::execution_space>
    team_policy(exec, n_cells, Kokkos::AUTO);
  team_policy.set_scratch_size(
    0, Kokkos::PerTeam(ScratchView::shmem_size(n_dofs)));

  Kokkos::parallel_for(team_policy,
                       ResolveCell<dim, fe_degree, transpose>{
                         weights, masks, src, dst_reference, dst});
  Kokkos::fence();

  auto reference_host =
    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, dst_reference);
  auto host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, dst);

  double sum = 0., sum_of_squares = 0., max_difference = 0.;
  for (unsigned int i = 0; i < n_cells * n_dofs; ++i)
    {
      sum += reference_host(i);
      sum_of_squares += reference_host(i) * reference_host(i);
      max_difference =
        std::max(max_difference, std::abs(host(i) - reference_host(i)));
    }

  deallog << "  checksum reference: " << sum << " " << sum_of_squares
          << std::endl;
  deallog << "  max |new - reference|: " << max_difference << std::endl;

  // the smallest case in full
  if (dim == 2 && fe_degree == 1)
    for (unsigned int cell = 0; cell < n_cells; ++cell)
      {
        deallog << "  mask " << static_cast<unsigned int>(make_mask(cell))
                << ":";
        for (unsigned int i = 0; i < n_dofs; ++i)
          deallog << " " << host(cell * n_dofs + i);
        deallog << std::endl;
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
