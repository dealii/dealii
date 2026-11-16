// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2021 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

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



    /**
     * Apply the hanging-node interpolation along one @p direction to the
     * values of a cell, stored in @p values with `(fe_degree + 1)^dim`
     * entries in lexicographic order.
     *
     * For a given direction, the constrained degrees of freedom of a cell
     * are complete lines along that direction: whether a degree of freedom
     * is constrained depends only on its indices orthogonal to the
     * direction (see is_constrained_dof_2d() and is_constrained_dof_3d()),
     * and the interpolation replaces the values on a line by a linear
     * combination of the values on the same line. Each thread therefore
     * handles one line: it loads the line into registers, applies the 1d
     * weights and writes the line back. Since the lines are disjoint, no
     * temporary storage and only a single team barrier are needed.
     */
    template <int  dim,
              int  fe_degree,
              int  direction,
              bool transpose,
              typename Number,
              typename ViewType>
    DEAL_II_HOST_DEVICE inline void
    interpolate_boundary(
      const Kokkos::TeamPolicy<
        MemorySpace::Default::kokkos_space::execution_space>::member_type
        &team_member,
      Kokkos::View<Number *, MemorySpace::Default::kokkos_space>
        constraint_weights,
      const dealii::internal::MatrixFreeFunctions::ConstraintKinds
               constraint_mask,
      ViewType values)
    {
      static_assert(dim == 2 || dim == 3, "Only implemented for 2d and 3d");
      static_assert(direction >= 0 && direction < dim,
                    "direction must be in [0, dim)");

      using ConstraintKinds =
        dealii::internal::MatrixFreeFunctions::ConstraintKinds;

      constexpr unsigned int n_dofs_1d = fe_degree + 1;

      // lines along the direction per cell, and the stride along a line
      constexpr unsigned int n_lines = Utilities::pow(n_dofs_1d, dim - 1);
      constexpr unsigned int stride  = Utilities::pow(n_dofs_1d, direction);

      // The subcell bit of this direction selects between the weights and
      // their point-reflected variant, the face (and edge) bits of the other
      // directions whether the lines on those faces are constrained at all.
      // face1_type and face2_type are only needed in 3d.
      const auto this_type = (direction == 0) ? ConstraintKinds::subcell_x :
                             (direction == 1) ? ConstraintKinds::subcell_y :
                                                ConstraintKinds::subcell_z;
      [[maybe_unused]] const auto face1_type =
        (direction == 0) ? ConstraintKinds::subcell_y :
        (direction == 1) ? ConstraintKinds::subcell_z :
                           ConstraintKinds::subcell_x;
      [[maybe_unused]] const auto face2_type =
        (direction == 0) ? ConstraintKinds::subcell_z :
        (direction == 1) ? ConstraintKinds::subcell_x :
                           ConstraintKinds::subcell_y;

      const auto face1 = (direction == 0) ? ConstraintKinds::face_y :
                         (direction == 1) ? ConstraintKinds::face_z :
                                            ConstraintKinds::face_x;
      const auto face2 = (direction == 0) ? ConstraintKinds::face_z :
                         (direction == 1) ? ConstraintKinds::face_x :
                                            ConstraintKinds::face_y;
      const auto edge  = (direction == 0) ? ConstraintKinds::edge_x :
                         (direction == 1) ? ConstraintKinds::edge_y :
                                            ConstraintKinds::edge_z;
      // in 2d the lines along x lie on the y faces and vice versa
      const auto constrained_faces =
        (dim == 2) ? ((direction == 0) ? ConstraintKinds::face_y :
                                         ConstraintKinds::face_x) :
                     (face1 | face2 | edge);

      // the same for all threads of the team, so they all skip the barrier
      if ((constraint_mask & constrained_faces) ==
          ConstraintKinds::unconstrained)
        return;

      Kokkos::parallel_for(
        Kokkos::TeamThreadRange(team_member, n_lines),
        [&](const unsigned int line) {
          // indices of the line: zero along the direction, the others
          // follow from the line index in increasing coordinate order
          unsigned int x_idx = 0, y_idx = 0, z_idx = 0;
          if constexpr (dim == 2)
            {
              if (direction == 0)
                y_idx = line;
              else
                x_idx = line;
            }
          else
            {
              const unsigned int i1 = line % n_dofs_1d;
              const unsigned int i2 = line / n_dofs_1d;
              if (direction == 0)
                {
                  y_idx = i1;
                  z_idx = i2;
                }
              else if (direction == 1)
                {
                  x_idx = i1;
                  z_idx = i2;
                }
              else
                {
                  x_idx = i1;
                  y_idx = i2;
                }
            }

          bool constrained_line;
          if constexpr (dim == 2)
            constrained_line =
              is_constrained_dof_2d<fe_degree, direction>(constraint_mask,
                                                          x_idx,
                                                          y_idx);
          else
            constrained_line =
              is_constrained_dof_3d<fe_degree, direction>(constraint_mask,
                                                          x_idx,
                                                          y_idx,
                                                          z_idx,
                                                          face1_type,
                                                          face2_type,
                                                          face1,
                                                          face2,
                                                          edge);
          if (!constrained_line)
            return;

          const unsigned int offset = (dim == 2) ?
                                        index2<n_dofs_1d>(x_idx, y_idx) :
                                        index3<n_dofs_1d>(x_idx, y_idx, z_idx);

          // the line is read completely before it is written
          Number line_values[n_dofs_1d];
          for (unsigned int i = 0; i < n_dofs_1d; ++i)
            line_values[i] = values(offset + i * stride);

          const bool type =
            (constraint_mask & this_type) != ConstraintKinds::unconstrained;

          for (unsigned int j = 0; j < n_dofs_1d; ++j)
            {
              Number sum = 0.;
              if (type)
                for (unsigned int i = 0; i < n_dofs_1d; ++i)
                  {
                    const Number w = transpose ?
                                       constraint_weights[i * n_dofs_1d + j] :
                                       constraint_weights[j * n_dofs_1d + i];
                    sum += w * line_values[i];
                  }
              else
                for (unsigned int i = 0; i < n_dofs_1d; ++i)
                  {
                    const Number w =
                      transpose ?
                        constraint_weights[(fe_degree - i) * n_dofs_1d +
                                           fe_degree - j] :
                        constraint_weights[(fe_degree - j) * n_dofs_1d +
                                           fe_degree - i];
                    sum += w * line_values[i];
                  }
              values(offset + j * stride) = sum;
            }
        });

      team_member.team_barrier();
    }



    /**
     * This function resolves the hanging nodes using tensor product, i.e.,
     * it applies interpolate_boundary() in all directions one after the
     * other.
     *
     * The implementation of this function is explained in Section 3 of
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
          interpolate_boundary<2, fe_degree, 0, transpose>(team_member,
                                                           constraint_weights,
                                                           constraint_mask,
                                                           values);
          interpolate_boundary<2, fe_degree, 1, transpose>(team_member,
                                                           constraint_weights,
                                                           constraint_mask,
                                                           values);
        }
      else if constexpr (dim == 3)
        {
          // Interpolate y and z faces (x-direction)
          interpolate_boundary<3, fe_degree, 0, transpose>(team_member,
                                                           constraint_weights,
                                                           constraint_mask,
                                                           values);
          // Interpolate x and z faces (y-direction)
          interpolate_boundary<3, fe_degree, 1, transpose>(team_member,
                                                           constraint_weights,
                                                           constraint_mask,
                                                           values);
          // Interpolate x and y faces (z-direction)
          interpolate_boundary<3, fe_degree, 2, transpose>(team_member,
                                                           constraint_weights,
                                                           constraint_mask,
                                                           values);
        }
    }
  } // namespace internal
} // namespace Portable

DEAL_II_NAMESPACE_CLOSE
#endif
