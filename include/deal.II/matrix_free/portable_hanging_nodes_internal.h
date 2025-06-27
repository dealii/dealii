// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_portable_hanging_nodes_internal_h
#define dealii_portable_hanging_nodes_internal_h

#include <deal.II/base/config.h>

#include <deal.II/matrix_free/hanging_nodes_internal.h>

#include <Kokkos_Core.hpp>


DEAL_II_NAMESPACE_OPEN
namespace Portable
{
  namespace internal
  {
    //------------------------------------------------------------------------//
    // Functions for resolving the hanging node constraints on the GPU        //
    //------------------------------------------------------------------------//
    template <unsigned int size>
    DEAL_II_HOST_DEVICE inline unsigned int
    index2(unsigned int i, unsigned int j)
    {
      return i + size * j;
    }



    template <unsigned int size>
    DEAL_II_HOST_DEVICE inline unsigned int
    index3(unsigned int i, unsigned int j, unsigned int k)
    {
      return i + size * j + size * size * k;
    }



    template <unsigned int fe_degree, unsigned int direction>
    DEAL_II_HOST_DEVICE inline bool
    is_constrained_dof_2d(
      const dealii::internal::MatrixFreeFunctions::ConstraintKinds
                        &constraint_mask,
      const unsigned int x_idx,
      const unsigned int y_idx)
    {
      return ((direction == 0) &&
              (((constraint_mask & dealii::internal::MatrixFreeFunctions::
                                     ConstraintKinds::subcell_y) !=
                dealii::internal::MatrixFreeFunctions::ConstraintKinds::
                  unconstrained) ?
                 (y_idx == 0) :
                 (y_idx == fe_degree))) ||
             ((direction == 1) &&
              (((constraint_mask & dealii::internal::MatrixFreeFunctions::
                                     ConstraintKinds::subcell_x) !=
                dealii::internal::MatrixFreeFunctions::ConstraintKinds::
                  unconstrained) ?
                 (x_idx == 0) :
                 (x_idx == fe_degree)));
    }

    template <unsigned int fe_degree, unsigned int direction>
    DEAL_II_HOST_DEVICE inline bool
    is_constrained_dof_3d(
      const dealii::internal::MatrixFreeFunctions::ConstraintKinds
                        &constraint_mask,
      const unsigned int x_idx,
      const unsigned int y_idx,
      const unsigned int z_idx,
      const dealii::internal::MatrixFreeFunctions::ConstraintKinds face1_type,
      const dealii::internal::MatrixFreeFunctions::ConstraintKinds face2_type,
      const dealii::internal::MatrixFreeFunctions::ConstraintKinds face1,
      const dealii::internal::MatrixFreeFunctions::ConstraintKinds face2,
      const dealii::internal::MatrixFreeFunctions::ConstraintKinds edge)
    {
      const unsigned int face1_idx = (direction == 0) ? y_idx :
                                     (direction == 1) ? z_idx :
                                                        x_idx;
      const unsigned int face2_idx = (direction == 0) ? z_idx :
                                     (direction == 1) ? x_idx :
                                                        y_idx;

      const bool on_face1 = ((constraint_mask & face1_type) !=
                             dealii::internal::MatrixFreeFunctions::
                               ConstraintKinds::unconstrained) ?
                              (face1_idx == 0) :
                              (face1_idx == fe_degree);
      const bool on_face2 = ((constraint_mask & face2_type) !=
                             dealii::internal::MatrixFreeFunctions::
                               ConstraintKinds::unconstrained) ?
                              (face2_idx == 0) :
                              (face2_idx == fe_degree);
      return (
        (((constraint_mask & face1) != dealii::internal::MatrixFreeFunctions::
                                         ConstraintKinds::unconstrained) &&
         on_face1) ||
        (((constraint_mask & face2) != dealii::internal::MatrixFreeFunctions::
                                         ConstraintKinds::unconstrained) &&
         on_face2) ||
        (((constraint_mask & edge) != dealii::internal::MatrixFreeFunctions::
                                        ConstraintKinds::unconstrained) &&
         on_face1 && on_face2));
    }



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

      Number tmp[n_q_points];
      Kokkos::parallel_for(
        Kokkos::TeamThreadRange(team_member, n_q_points),
        [&](const int &q_point) {
          const unsigned int x_idx = q_point % n_q_points_1d;
          const unsigned int y_idx = q_point / n_q_points_1d;

          const auto this_type =
            (direction == 0) ?
              dealii::internal::MatrixFreeFunctions::ConstraintKinds::
                subcell_x :
              dealii::internal::MatrixFreeFunctions::ConstraintKinds::subcell_y;

          const unsigned int interp_idx = (direction == 0) ? x_idx : y_idx;
          tmp[q_point]                  = 0;

          // Flag is true if for the given direction, the dof is constrained
          // with the right type and is on the correct side (left (= 0) or right
          // (= fe_degree))
          const bool constrained_dof =
            is_constrained_dof_2d<fe_degree, direction>(constraint_mask,
                                                        x_idx,
                                                        y_idx);

          if (constrained_face && constrained_dof)
            {
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
                      tmp[q_point] += w * values[real_idx];
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
                      tmp[q_point] += w * values[real_idx];
                    }
                }
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

      Number tmp[n_q_points];
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
          tmp[q_point] = 0;
          if ((constrained_face != dealii::internal::MatrixFreeFunctions::
                                     ConstraintKinds::unconstrained) &&
              constrained_dof)
            {
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
                      tmp[q_point] += w * values[real_idx];
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
                      tmp[q_point] += w * values[real_idx];
                    }
                }
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
  } // namespace internal
} // namespace Portable

DEAL_II_NAMESPACE_CLOSE
#endif
