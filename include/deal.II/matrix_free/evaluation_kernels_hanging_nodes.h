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
      const unsigned int given_degree =
        fe_degree != -1 ? fe_degree :
                          fe_eval.get_shape_info().data.front().fe_degree;

      const auto &interpolation_matrices =
        fe_eval.get_shape_info().data.front().subface_interpolation_matrices;

      const auto is_set = [](const auto a, const auto b) -> bool {
        return (a & b) == b;
      };

      const auto not_set = [](const auto a, const auto b) -> bool {
        return (a & b) == MatrixFreeFunctions::ConstraintKinds::unconstrained;
      };

      const unsigned int points = given_degree + 1;

      for (unsigned int c = 0; c < n_desired_components; ++c)
        {
          for (unsigned int v = 0; v < Number::size(); ++v)
            {
              const auto mask = constraint_mask[v];

              if (mask == MatrixFreeFunctions::ConstraintKinds::unconstrained)
                continue;

              if (dim == 2) // 2D: only faces
                {
                  // direction 0:
                  if ((mask & MatrixFreeFunctions::ConstraintKinds::face_y) !=
                      MatrixFreeFunctions::ConstraintKinds::unconstrained)
                    {
                      const bool is_subface_0 =
                        (mask & MatrixFreeFunctions::ConstraintKinds::type_x) ==
                        MatrixFreeFunctions::ConstraintKinds::unconstrained;

                      const Number *weights =
                        interpolation_matrices[is_subface_0].data();

                      if (is_set(mask,
                                 MatrixFreeFunctions::ConstraintKinds::type_y))
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
                  if ((mask & MatrixFreeFunctions::ConstraintKinds::face_x) !=
                      MatrixFreeFunctions::ConstraintKinds::unconstrained)
                    {
                      const bool is_subface_0 =
                        (mask & MatrixFreeFunctions::ConstraintKinds::type_y) ==
                        MatrixFreeFunctions::ConstraintKinds::unconstrained;

                      const Number *weights =
                        interpolation_matrices[is_subface_0].data();

                      if (is_set(mask,
                                 MatrixFreeFunctions::ConstraintKinds::type_x))
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

                  const bool is_face_0 =
                    is_set(mask,
                           MatrixFreeFunctions::ConstraintKinds::face_x) &&
                    is_set(mask, MatrixFreeFunctions::ConstraintKinds::type_x);
                  const bool is_face_1 =
                    is_set(mask,
                           MatrixFreeFunctions::ConstraintKinds::face_x) &&
                    not_set(mask, MatrixFreeFunctions::ConstraintKinds::type_x);
                  const bool is_face_2 =
                    is_set(mask,
                           MatrixFreeFunctions::ConstraintKinds::face_y) &&
                    is_set(mask, MatrixFreeFunctions::ConstraintKinds::type_y);
                  const bool is_face_3 =
                    is_set(mask,
                           MatrixFreeFunctions::ConstraintKinds::face_y) &&
                    not_set(mask, MatrixFreeFunctions::ConstraintKinds::type_y);
                  const bool is_face_4 =
                    is_set(mask,
                           MatrixFreeFunctions::ConstraintKinds::face_z) &&
                    is_set(mask, MatrixFreeFunctions::ConstraintKinds::type_z);
                  const bool is_face_5 =
                    is_set(mask,
                           MatrixFreeFunctions::ConstraintKinds::face_z) &&
                    not_set(mask, MatrixFreeFunctions::ConstraintKinds::type_z);

                  const bool is_edge_2 =
                    is_set(mask,
                           MatrixFreeFunctions::ConstraintKinds::edge_yz) &&
                    is_set(mask,
                           MatrixFreeFunctions::ConstraintKinds::type_y) &&
                    is_set(mask, MatrixFreeFunctions::ConstraintKinds::type_z);
                  const bool is_edge_3 =
                    is_set(mask,
                           MatrixFreeFunctions::ConstraintKinds::edge_yz) &&
                    not_set(mask,
                            MatrixFreeFunctions::ConstraintKinds::type_y) &&
                    is_set(mask, MatrixFreeFunctions::ConstraintKinds::type_z);
                  const bool is_edge_6 =
                    is_set(mask,
                           MatrixFreeFunctions::ConstraintKinds::edge_yz) &&
                    is_set(mask,
                           MatrixFreeFunctions::ConstraintKinds::type_y) &&
                    not_set(mask, MatrixFreeFunctions::ConstraintKinds::type_z);
                  const bool is_edge_7 =
                    is_set(mask,
                           MatrixFreeFunctions::ConstraintKinds::edge_yz) &&
                    not_set(mask,
                            MatrixFreeFunctions::ConstraintKinds::type_y) &&
                    not_set(mask, MatrixFreeFunctions::ConstraintKinds::type_z);

                  const bool is_edge_0 =
                    is_set(mask,
                           MatrixFreeFunctions::ConstraintKinds::edge_zx) &&
                    is_set(mask,
                           MatrixFreeFunctions::ConstraintKinds::type_x) &&
                    is_set(mask, MatrixFreeFunctions::ConstraintKinds::type_z);
                  const bool is_edge_1 =
                    is_set(mask,
                           MatrixFreeFunctions::ConstraintKinds::edge_zx) &&
                    not_set(mask,
                            MatrixFreeFunctions::ConstraintKinds::type_x) &&
                    is_set(mask, MatrixFreeFunctions::ConstraintKinds::type_z);
                  const bool is_edge_4 =
                    is_set(mask,
                           MatrixFreeFunctions::ConstraintKinds::edge_zx) &&
                    is_set(mask,
                           MatrixFreeFunctions::ConstraintKinds::type_x) &&
                    not_set(mask, MatrixFreeFunctions::ConstraintKinds::type_z);
                  const bool is_edge_5 =
                    is_set(mask,
                           MatrixFreeFunctions::ConstraintKinds::edge_zx) &&
                    not_set(mask,
                            MatrixFreeFunctions::ConstraintKinds::type_x) &&
                    not_set(mask, MatrixFreeFunctions::ConstraintKinds::type_z);

                  const bool is_edge_8 =
                    is_set(mask,
                           MatrixFreeFunctions::ConstraintKinds::edge_xy) &&
                    is_set(mask,
                           MatrixFreeFunctions::ConstraintKinds::type_x) &&
                    is_set(mask, MatrixFreeFunctions::ConstraintKinds::type_y);
                  const bool is_edge_9 =
                    is_set(mask,
                           MatrixFreeFunctions::ConstraintKinds::edge_xy) &&
                    not_set(mask,
                            MatrixFreeFunctions::ConstraintKinds::type_x) &&
                    is_set(mask, MatrixFreeFunctions::ConstraintKinds::type_y);
                  const bool is_edge_10 =
                    is_set(mask,
                           MatrixFreeFunctions::ConstraintKinds::edge_xy) &&
                    is_set(mask,
                           MatrixFreeFunctions::ConstraintKinds::type_x) &&
                    not_set(mask, MatrixFreeFunctions::ConstraintKinds::type_y);
                  const bool is_edge_11 =
                    is_set(mask,
                           MatrixFreeFunctions::ConstraintKinds::edge_xy) &&
                    not_set(mask,
                            MatrixFreeFunctions::ConstraintKinds::type_x) &&
                    not_set(mask, MatrixFreeFunctions::ConstraintKinds::type_y);

                  // direction 0:
                  {
                    const bool is_subface_0 =
                      (mask & MatrixFreeFunctions::ConstraintKinds::type_x) ==
                      MatrixFreeFunctions::ConstraintKinds::unconstrained;

                    const Number *weights =
                      interpolation_matrices[is_subface_0].data();

                    // ... faces
                    if (is_face_2)
                      interpolate_3D_face<fe_degree, 0, 1, transpose>(
                        p0,
                        given_degree,
                        v,
                        weights,
                        values); // face 2
                    else if (is_face_3)
                      interpolate_3D_face<fe_degree, 0, 1, transpose>(
                        p2,
                        given_degree,
                        v,
                        weights,
                        values); // face 3
                    if (is_face_4)
                      interpolate_3D_face<fe_degree, 0, 2, transpose>(
                        p0,
                        given_degree,
                        v,
                        weights,
                        values); // face 4
                    else if (is_face_5)
                      interpolate_3D_face<fe_degree, 0, 2, transpose>(
                        p4,
                        given_degree,
                        v,
                        weights,
                        values); // face 5

                    // ... edges
                    if (is_face_2 || is_face_4 || is_edge_2)
                      interpolate_3D_edge<fe_degree, 0, transpose>(
                        p0,
                        given_degree,
                        v,
                        weights,
                        values); // edge 2
                    if (is_face_3 || is_face_4 || is_edge_3)
                      interpolate_3D_edge<fe_degree, 0, transpose>(
                        p2,
                        given_degree,
                        v,
                        weights,
                        values); // edge 3
                    if (is_face_2 || is_face_5 || is_edge_6)
                      interpolate_3D_edge<fe_degree, 0, transpose>(
                        p4,
                        given_degree,
                        v,
                        weights,
                        values); // edge 6
                    if (is_face_3 || is_face_5 || is_edge_7)
                      interpolate_3D_edge<fe_degree, 0, transpose>(
                        p6,
                        given_degree,
                        v,
                        weights,
                        values); // edge 7
                  }

                  // direction 1:
                  {
                    const bool is_subface_0 =
                      (mask & MatrixFreeFunctions::ConstraintKinds::type_y) ==
                      MatrixFreeFunctions::ConstraintKinds::unconstrained;

                    const Number *weights =
                      interpolation_matrices[is_subface_0].data();

                    // ... faces
                    if (is_face_0)
                      interpolate_3D_face<fe_degree, 1, 0, transpose>(
                        p0,
                        given_degree,
                        v,
                        weights,
                        values); // face 0
                    else if (is_face_1)
                      interpolate_3D_face<fe_degree, 1, 0, transpose>(
                        p1,
                        given_degree,
                        v,
                        weights,
                        values); // face 1
                    if (is_face_4)
                      interpolate_3D_face<fe_degree, 1, 2, transpose>(
                        p0,
                        given_degree,
                        v,
                        weights,
                        values); // face 4
                    else if (is_face_5)
                      interpolate_3D_face<fe_degree, 1, 2, transpose>(
                        p4,
                        given_degree,
                        v,
                        weights,
                        values); // face 5

                    // ... edges
                    if (is_face_0 || is_face_4 || is_edge_0)
                      interpolate_3D_edge<fe_degree, 1, transpose>(
                        p0,
                        given_degree,
                        v,
                        weights,
                        values); // edge 0
                    if (is_face_1 || is_face_4 || is_edge_1)
                      interpolate_3D_edge<fe_degree, 1, transpose>(
                        p1,
                        given_degree,
                        v,
                        weights,
                        values); // edge 1
                    if (is_face_0 || is_face_5 || is_edge_4)
                      interpolate_3D_edge<fe_degree, 1, transpose>(
                        p4,
                        given_degree,
                        v,
                        weights,
                        values); // edge 4
                    if (is_face_1 || is_face_5 || is_edge_5)
                      interpolate_3D_edge<fe_degree, 1, transpose>(
                        p5,
                        given_degree,
                        v,
                        weights,
                        values); // edge 5
                  }

                  // direction 2:
                  {
                    const bool is_subface_0 =
                      (mask & MatrixFreeFunctions::ConstraintKinds::type_z) ==
                      MatrixFreeFunctions::ConstraintKinds::unconstrained;

                    const Number *weights =
                      interpolation_matrices[is_subface_0].data();

                    // ... faces
                    if (is_face_0)
                      interpolate_3D_face<fe_degree, 2, 0, transpose>(
                        p0,
                        given_degree,
                        v,
                        weights,
                        values); // face 0
                    else if (is_face_1)
                      interpolate_3D_face<fe_degree, 2, 0, transpose>(
                        p1,
                        given_degree,
                        v,
                        weights,
                        values); // face 1
                    if (is_face_2)
                      interpolate_3D_face<fe_degree, 2, 1, transpose>(
                        p0,
                        given_degree,
                        v,
                        weights,
                        values); // face 2
                    else if (is_face_3)
                      interpolate_3D_face<fe_degree, 2, 1, transpose>(
                        p2,
                        given_degree,
                        v,
                        weights,
                        values); // face 3

                    // ... edges
                    if (is_face_0 || is_face_2 || is_edge_8)
                      interpolate_3D_edge<fe_degree, 2, transpose>(
                        p0,
                        given_degree,
                        v,
                        weights,
                        values); // edge 8
                    if (is_face_1 || is_face_2 || is_edge_9)
                      interpolate_3D_edge<fe_degree, 2, transpose>(
                        p1,
                        given_degree,
                        v,
                        weights,
                        values); // edge 9
                    if (is_face_0 || is_face_3 || is_edge_10)
                      interpolate_3D_edge<fe_degree, 2, transpose>(
                        p2,
                        given_degree,
                        v,
                        weights,
                        values); // edge 10
                    if (is_face_1 || is_face_3 || is_edge_11)
                      interpolate_3D_edge<fe_degree, 2, transpose>(
                        p3,
                        given_degree,
                        v,
                        weights,
                        values); // edge 11
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
