// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


#ifndef dealii_matrix_free_evaluation_kernels_hanging_nodes_h
#define dealii_matrix_free_evaluation_kernels_hanging_nodes_h

#include <deal.II/base/config.h>

#include <deal.II/base/vectorization.h>

#include <deal.II/matrix_free/hanging_nodes_internal.h>


DEAL_II_NAMESPACE_OPEN


// forward declaration
template <int, typename, bool, typename>
class FEEvaluationBaseData;



namespace internal
{
  template <int dim, typename Number, bool is_face>
  struct FEEvaluationImplHangingNodes
  {
    template <int fe_degree, int n_q_points_1d>
    static bool
    run(const unsigned int                  n_desired_components,
        const FEEvaluationBaseData<dim,
                                   typename Number::value_type,
                                   is_face,
                                   Number> &fe_eval,
        const bool                          transpose,
        const std::array<MatrixFreeFunctions::ConstraintKinds, Number::size()>
          &     c_mask,
        Number *values)
    {
      if (transpose)
        run_internal<fe_degree, true>(n_desired_components,
                                      fe_eval,
                                      c_mask,
                                      values);
      else
        run_internal<fe_degree, false>(n_desired_components,
                                       fe_eval,
                                       c_mask,
                                       values);

      return false;
    }

  private:
    template <int fe_degree, unsigned int side, bool transpose>
    static void
    interpolate_2D(const unsigned int given_degree,
                   const unsigned int v,
                   const Number *     weight,
                   Number *           values)
    {
      typename Number::value_type temp[40];

      const unsigned int points =
        (fe_degree != -1 ? fe_degree : given_degree) + 1;

      AssertIndexRange(points, 40);

      const unsigned int d = side / 2; // direction
      const unsigned int s = side % 2; // left or right surface

      const unsigned int offset = dealii::Utilities::pow(points, d + 1);
      const unsigned int stride =
        (s == 0 ? 0 : (points - 1)) * dealii::Utilities::pow(points, d);

      const unsigned int r1 = dealii::Utilities::pow(points, dim - d - 1);
      const unsigned int r2 = dealii::Utilities::pow(points, d);

      // copy result back
      for (unsigned int i = 0, k = 0; i < r1; ++i)
        for (unsigned int j = 0; j < r2; ++j, ++k)
          temp[k] = values[i * offset + stride + j][v];

      // perform interpolation point by point (note: r1 * r2 == points^(dim-1))
      for (unsigned int i = 0, k = 0; i < r1; ++i)
        for (unsigned int j = 0; j < r2; ++j, ++k)
          {
            typename Number::value_type sum = 0.0;
            for (unsigned int h = 0; h < points; ++h)
              sum += weight[(transpose ? 1 : points) * k +
                            (transpose ? points : 1) * h][v] *
                     temp[h];
            values[i * offset + stride + j][v] = sum;
          }
    }

    template <int          fe_degree,
              unsigned int direction,
              unsigned int d,
              bool         transpose>
    static void
    interpolate_3D_face(const unsigned int dof_offset,
                        const unsigned int given_degree,
                        const unsigned int v,
                        const Number *     weight,
                        Number *           values)
    {
      typename Number::value_type temp[40];

      const unsigned int points =
        (fe_degree != -1 ? fe_degree : given_degree) + 1;

      AssertIndexRange(points, 40);

      const unsigned int stride = Utilities::pow(points, direction);

      // direction   side0   side1   side2
      // 0             -      p^2      p
      // 1            p^2      -       1
      // 2             p       -       1
      const unsigned int stride2 =
        ((direction == 0 && d == 1) || (direction == 1 && d == 0)) ?
          (points * points) :
          (((direction == 0 && d == 2) || (direction == 2 && d == 0)) ? points :
                                                                        1);

      for (unsigned int g = 1; g < points - 1; ++g)
        {
          // copy result back
          for (unsigned int k = 0; k < points; ++k)
            temp[k] = values[dof_offset + k * stride + stride2 * g][v];

          // perform interpolation point by point
          for (unsigned int k = 0; k < points; ++k)
            {
              typename Number::value_type sum = 0.0;
              for (unsigned int h = 0; h < points; ++h)
                sum += weight[(transpose ? 1 : points) * k +
                              (transpose ? points : 1) * h][v] *
                       temp[h];
              values[dof_offset + k * stride + stride2 * g][v] = sum;
            }
        }
    }

    template <int fe_degree, unsigned int direction, bool transpose>
    static void
    interpolate_3D_edge(const unsigned int p,
                        const unsigned int given_degree,
                        const unsigned int v,
                        const Number *     weight,
                        Number *           values)
    {
      typename Number::value_type temp[40];

      const unsigned int points =
        (fe_degree != -1 ? fe_degree : given_degree) + 1;

      AssertIndexRange(points, 40);

      const unsigned int stride = Utilities::pow(points, direction);

      // copy result back
      for (unsigned int k = 0; k < points; ++k)
        temp[k] = values[p + k * stride][v];

      // perform interpolation point by point
      for (unsigned int k = 0; k < points; ++k)
        {
          typename Number::value_type sum = 0.0;
          for (unsigned int h = 0; h < points; ++h)
            sum += weight[(transpose ? 1 : points) * k +
                          (transpose ? points : 1) * h][v] *
                   temp[h];
          values[p + k * stride][v] = sum;
        }
    }

    template <int fe_degree, bool transpose>
    static void
    run_internal(const unsigned int                  n_desired_components,
                 const FEEvaluationBaseData<dim,
                                            typename Number::value_type,
                                            is_face,
                                            Number> &fe_eval,
                 const std::array<MatrixFreeFunctions::ConstraintKinds,
                                  Number::size()> &  constraint_mask,
                 Number *                            values)
    {
      using Kinds = MatrixFreeFunctions::ConstraintKinds;
      const unsigned int given_degree =
        fe_degree != -1 ? fe_degree :
                          fe_eval.get_shape_info().data.front().fe_degree;

      const auto &interpolation_matrices =
        fe_eval.get_shape_info().data.front().subface_interpolation_matrices;

      const auto is_set = [](const auto a, const auto b) -> bool {
        return (a & b) == b;
      };

      const auto not_set = [](const auto a, const auto b) -> bool {
        return (a & b) == Kinds::unconstrained;
      };

      const unsigned int points = given_degree + 1;

      for (unsigned int c = 0; c < n_desired_components; ++c)
        {
          for (unsigned int v = 0; v < Number::size(); ++v)
            {
              const auto mask = constraint_mask[v];

              if (mask == Kinds::unconstrained)
                continue;

              if (dim == 2) // 2D: only faces
                {
                  // direction 0:
                  if ((mask & Kinds::face_y) != Kinds::unconstrained)
                    {
                      const bool is_subface_0 =
                        (mask & Kinds::subcell_x) == Kinds::unconstrained;

                      const Number *weights =
                        interpolation_matrices[is_subface_0].data();

                      if (is_set(mask, Kinds::subcell_y))
                        interpolate_2D<fe_degree, 2, transpose>(
                          given_degree,
                          v,
                          weights,
                          values); // face 2
                      else
                        interpolate_2D<fe_degree, 3, transpose>(
                          given_degree,
                          v,
                          weights,
                          values); // face 3
                    }

                  // direction 1:
                  if ((mask & Kinds::face_x) != Kinds::unconstrained)
                    {
                      const bool is_subface_0 =
                        (mask & Kinds::subcell_y) == Kinds::unconstrained;

                      const Number *weights =
                        interpolation_matrices[is_subface_0].data();

                      if (is_set(mask, Kinds::subcell_x))
                        interpolate_2D<fe_degree, 0, transpose>(
                          given_degree,
                          v,
                          weights,
                          values); // face 0
                      else
                        interpolate_2D<fe_degree, 1, transpose>(
                          given_degree,
                          v,
                          weights,
                          values); // face 1
                    }
                }
              else if (dim == 3) // 3D faces and edges
                {
                  const unsigned int p0 = 0;
                  const unsigned int p1 = points - 1;
                  const unsigned int p2 = points * points - points;
                  const unsigned int p3 = points * points - 1;
                  const unsigned int p4 =
                    points * points * points - points * points;
                  const unsigned int p5 =
                    points * points * points - points * points + points - 1;
                  const unsigned int p6 = points * points * points - points;

                  std::array<std::array<char, 4>, 3> process_edge = {};
                  std::array<std::array<char, 4>, 3> process_face = {};

                  if (is_set(mask, Kinds::face_x))
                    {
                      const unsigned int side = not_set(mask, Kinds::subcell_x);
                      process_face[1][side]   = 1;
                      process_edge[1][side]   = 1;
                      process_edge[1][2 + side] = 1;
                      process_face[2][side]     = 1;
                      process_edge[2][side]     = 1;
                      process_edge[2][2 + side] = 1;
                    }
                  if (is_set(mask, Kinds::face_y))
                    {
                      const unsigned int side = not_set(mask, Kinds::subcell_y);
                      process_face[0][side]   = 1;
                      process_edge[0][side]   = 1;
                      process_edge[0][2 + side]     = 1;
                      process_face[2][2 + side]     = 1;
                      process_edge[2][2 * side]     = 1;
                      process_edge[2][2 * side + 1] = 1;
                    }
                  if (is_set(mask, Kinds::face_z))
                    {
                      const unsigned int side = not_set(mask, Kinds::subcell_z);
                      process_face[0][2 + side]     = 1;
                      process_edge[0][2 * side]     = 1;
                      process_edge[0][2 * side + 1] = 1;
                      process_face[1][2 + side]     = 1;
                      process_edge[1][2 * side]     = 1;
                      process_edge[1][2 * side + 1] = 1;
                    }
                  if (is_set(mask, Kinds::edge_x))
                    process_edge[0][not_set(mask, Kinds::subcell_z) * 2 +
                                    not_set(mask, Kinds::subcell_y)] = 1;
                  if (is_set(mask, Kinds::edge_y))
                    process_edge[1][not_set(mask, Kinds::subcell_z) * 2 +
                                    not_set(mask, Kinds::subcell_x)] = 1;
                  if (is_set(mask, Kinds::edge_z))
                    process_edge[2][not_set(mask, Kinds::subcell_y) * 2 +
                                    not_set(mask, Kinds::subcell_x)] = 1;

                  // direction 0:
                  {
                    const bool is_subface_0 =
                      (mask & Kinds::subcell_x) == Kinds::unconstrained;

                    const Number *weights =
                      interpolation_matrices[is_subface_0].data();

                    unsigned int face_offsets[4] = {p0, p2, p0, p4};
                    // face 2, 3
                    for (unsigned int face = 0; face < 2; ++face)
                      if (process_face[0][face])
                        interpolate_3D_face<fe_degree, 0, 1, transpose>(
                          face_offsets[face], given_degree, v, weights, values);
                    // face 4, 5
                    for (unsigned int face = 2; face < 4; ++face)
                      if (process_face[0][face])
                        interpolate_3D_face<fe_degree, 0, 2, transpose>(
                          face_offsets[face], given_degree, v, weights, values);

                    // edges
                    unsigned int edge_offsets[4] = {p0, p2, p4, p6};
                    for (unsigned int edge = 0; edge < 4; ++edge)
                      if (process_edge[0][edge])
                        interpolate_3D_edge<fe_degree, 0, transpose>(
                          edge_offsets[edge], given_degree, v, weights, values);
                  }

                  // direction 1:
                  {
                    const bool is_subface_0 =
                      (mask & Kinds::subcell_y) == Kinds::unconstrained;

                    const Number *weights =
                      interpolation_matrices[is_subface_0].data();

                    unsigned int face_offsets[4] = {p0, p1, p0, p4};
                    // face 0, 1
                    for (unsigned int face = 0; face < 2; ++face)
                      if (process_face[1][face])
                        interpolate_3D_face<fe_degree, 1, 0, transpose>(
                          face_offsets[face], given_degree, v, weights, values);
                    // face 4, 5
                    for (unsigned int face = 2; face < 4; ++face)
                      if (process_face[1][face])
                        interpolate_3D_face<fe_degree, 1, 2, transpose>(
                          face_offsets[face], given_degree, v, weights, values);

                    // edges
                    unsigned int edge_offsets[4] = {p0, p1, p4, p5};
                    for (unsigned int edge = 0; edge < 4; ++edge)
                      if (process_edge[1][edge])
                        interpolate_3D_edge<fe_degree, 1, transpose>(
                          edge_offsets[edge], given_degree, v, weights, values);
                  }

                  // direction 2:
                  {
                    const bool is_subface_0 =
                      (mask & Kinds::subcell_z) == Kinds::unconstrained;

                    const Number *weights =
                      interpolation_matrices[is_subface_0].data();

                    unsigned int face_offsets[4] = {p0, p1, p0, p2};
                    // face 0, 1
                    for (unsigned int face = 0; face < 2; ++face)
                      if (process_face[2][face])
                        interpolate_3D_face<fe_degree, 2, 0, transpose>(
                          face_offsets[face], given_degree, v, weights, values);
                    // face 2, 3
                    for (unsigned int face = 2; face < 4; ++face)
                      if (process_face[2][face])
                        interpolate_3D_face<fe_degree, 2, 1, transpose>(
                          face_offsets[face], given_degree, v, weights, values);

                    // edges
                    unsigned int edge_offsets[4] = {p0, p1, p2, p3};
                    for (unsigned int edge = 0; edge < 4; ++edge)
                      if (process_edge[2][edge])
                        interpolate_3D_edge<fe_degree, 2, transpose>(
                          edge_offsets[edge], given_degree, v, weights, values);
                  }
                }
              else
                {
                  Assert(false, ExcNotImplemented());
                }
            }

          values += fe_eval.get_shape_info().dofs_per_component_on_cell;
        }
    }
  };


} // end of namespace internal


DEAL_II_NAMESPACE_CLOSE

#endif
