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

#include <deal.II/base/ndarray.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/matrix_free/hanging_nodes_internal.h>


DEAL_II_NAMESPACE_OPEN


// forward declaration
template <int, typename, bool, typename>
class FEEvaluationBaseData;



namespace internal
{
  enum class FEEvaluationImplHangingNodesRunnerTypes
  {
    vectorized
  };



  template <FEEvaluationImplHangingNodesRunnerTypes,
            int dim,
            int fe_degree,
            typename Number,
            bool is_face>
  class FEEvaluationImplHangingNodesRunner;



  template <int dim, int fe_degree, typename Number, bool is_face>
  class FEEvaluationImplHangingNodesRunner<
    FEEvaluationImplHangingNodesRunnerTypes::vectorized,
    dim,
    fe_degree,
    Number,
    is_face>
  {
  private:
    template <int structdim, unsigned int direction, bool transpose>
    static void
    interpolate(const unsigned int             offset,
                const unsigned int             outer_stride,
                const unsigned int             given_degree,
                const Number                   mask_weight,
                const Number                   mask_write,
                const Number *DEAL_II_RESTRICT weights,
                Number *DEAL_II_RESTRICT       values)
    {
      static_assert(structdim == 1 || structdim == 2,
                    "Only 1D and 2D interpolation implemented");
      Number temp[fe_degree != -1 ? fe_degree + 1 : 40];

      const unsigned int points =
        (fe_degree != -1 ? fe_degree : given_degree) + 1;

      AssertIndexRange(points, 40);

      const unsigned int stride = Utilities::pow(points, direction);

      const unsigned int end_of_outer_loop = structdim == 1 ? 2 : points - 1;
      for (unsigned int g = 1; g < end_of_outer_loop; ++g)
        {
          const unsigned int my_offset =
            offset + (structdim > 1 ? g * outer_stride : 0);

          // extract values to interpolate
          for (unsigned int k = 0; k < points; ++k)
            temp[k] = values[my_offset + k * stride];

          // perform interpolation point by point and write back
          for (unsigned int k = 0; k < points / 2; ++k)
            {
              const unsigned int kmirror = points - 1 - k;
              Number sum0 = Number(), sum1 = Number(), sum2 = Number(),
                     sum3 = Number();
              for (unsigned int h = 0; h < points; ++h)
                {
                  const unsigned int hmirror = points - 1 - h;
                  // load from both sides of the interpolation matrix to
                  // reflect symmetry between the two subfaces along that
                  // direction
                  const Number w0 = weights[(transpose ? 1 : points) * kmirror +
                                            (transpose ? points : 1) * hmirror];
                  const Number w1 = weights[(transpose ? 1 : points) * k +
                                            (transpose ? points : 1) * h];
                  sum0 += temp[h] * w0;
                  sum1 += temp[h] * w1;
                  sum2 += temp[hmirror] * w1;
                  sum3 += temp[hmirror] * w0;
                }
              values[my_offset + k * stride] =
                temp[k] +
                mask_write * (sum0 + mask_weight * (sum1 - sum0) - temp[k]);
              values[my_offset + kmirror * stride] =
                temp[kmirror] +
                mask_write *
                  (sum2 + mask_weight * (sum3 - sum2) - temp[kmirror]);
            }

          // cleanup case
          if (points % 2)
            {
              const unsigned int k = points / 2;
              Number sum0 = temp[k] * weights[(transpose ? 1 : points) * k +
                                              (transpose ? points : 1) * k],
                     sum1 = sum0;
              for (unsigned int h = 0; h < points / 2; ++h)
                {
                  const unsigned int hmirror = points - 1 - h;
                  const Number       w0 = weights[(transpose ? 1 : points) * k +
                                            (transpose ? points : 1) * hmirror];
                  const Number       w1 = weights[(transpose ? 1 : points) * k +
                                            (transpose ? points : 1) * h];
                  sum0 += temp[h] * w0;
                  sum0 += temp[hmirror] * w1;
                  sum1 += temp[h] * w1;
                  sum1 += temp[hmirror] * w0;
                }
              values[my_offset + k * stride] =
                temp[k] +
                mask_write * (sum0 + mask_weight * (sum1 - sum0) - temp[k]);
            }
        }
    }

  public:
    template <bool transpose>
    static void
    run_internal(const unsigned int                  n_components,
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

      const Number *DEAL_II_RESTRICT weights =
        fe_eval.get_shape_info()
          .data.front()
          .subface_interpolation_matrices[0]
          .data();

      const unsigned int points = given_degree + 1;
      const unsigned int n_dofs =
        fe_eval.get_shape_info().dofs_per_component_on_cell;

      if (dim == 2)
        {
          dealii::ndarray<Number, 2>    mask_weights = {};
          dealii::ndarray<Number, 2, 2> mask_write   = {};
          dealii::ndarray<bool, 2, 2>   do_face      = {};

          for (unsigned int v = 0; v < Number::size(); ++v)
            {
              const auto kind = static_cast<std::uint16_t>(constraint_mask[v]);
              const bool subcell_x = (kind >> 0) & 1;
              const bool subcell_y = (kind >> 1) & 1;
              const bool face_x    = (kind >> 3) & 1;
              const bool face_y    = (kind >> 4) & 1;

              if (face_y)
                {
                  const unsigned int side = !subcell_y;
                  mask_write[0][side][v]  = 1;
                  do_face[0][side]        = true;
                  mask_weights[0][v]      = subcell_x;
                }

              if (face_x)
                {
                  const unsigned int side = !subcell_x;
                  mask_write[1][side][v]  = 1;
                  do_face[1][side]        = true;
                  mask_weights[1][v]      = subcell_y;
                }
            }

          // x direction
          {
            const std::array<unsigned int, 2> offsets = {
              {0, (points - 1) * points}};
            for (unsigned int c = 0; c < n_components; ++c)
              for (unsigned int face = 0; face < 2; ++face)
                if (do_face[0][face])
                  interpolate<1, 0, transpose>(offsets[face],
                                               0,
                                               given_degree,
                                               mask_weights[0],
                                               mask_write[0][face],
                                               weights,
                                               values + c * n_dofs);
          }

          // y direction
          {
            const std::array<unsigned int, 2> offsets = {{0, points - 1}};
            for (unsigned int c = 0; c < n_components; ++c)
              for (unsigned int face = 0; face < 2; ++face)
                if (do_face[1][face])
                  interpolate<1, 1, transpose>(offsets[face],
                                               0,
                                               given_degree,
                                               mask_weights[1],
                                               mask_write[1][face],
                                               weights,
                                               values + c * n_dofs);
          }
        }
      else if (dim == 3)
        {
          const unsigned int p0 = 0;
          const unsigned int p1 = points - 1;
          const unsigned int p2 = points * points - points;
          const unsigned int p3 = points * points - 1;
          const unsigned int p4 = points * points * points - points * points;
          const unsigned int p5 =
            points * points * points - points * points + points - 1;
          const unsigned int p6 = points * points * points - points;

          dealii::ndarray<bool, 3, 4>   process_edge = {};
          dealii::ndarray<bool, 3, 4>   process_face = {};
          dealii::ndarray<Number, 3, 4> mask_edge    = {};
          dealii::ndarray<Number, 3, 4> mask_face    = {};
          dealii::ndarray<Number, 3>    mask_weights = {};

          for (unsigned int v = 0; v < Number::size(); ++v)
            {
              const auto kind = static_cast<std::uint16_t>(constraint_mask[v]);
              const bool subcell_x = (kind >> 0) & 1;
              const bool subcell_y = (kind >> 1) & 1;
              const bool subcell_z = (kind >> 2) & 1;
              const bool face_x    = (kind >> 3) & 1;
              const bool face_y    = (kind >> 4) & 1;
              const bool face_z    = (kind >> 5) & 1;
              const bool edge_x    = (kind >> 6) & 1;
              const bool edge_y    = (kind >> 7) & 1;
              const bool edge_z    = (kind >> 8) & 1;

              if (subcell_x)
                mask_weights[0][v] = 1;
              if (subcell_y)
                mask_weights[1][v] = 1;
              if (subcell_z)
                mask_weights[2][v] = 1;

              if (face_x)
                {
                  const unsigned int side = !subcell_x;

                  mask_face[1][side][v] = process_face[1][side] = true;
                  mask_edge[1][side][v] = process_edge[1][side] = true;
                  mask_edge[1][2 + side][v] = process_edge[1][2 + side] = true;
                  mask_face[2][side][v] = process_face[2][side] = true;
                  mask_edge[2][side][v] = process_edge[2][side] = true;
                  mask_edge[2][2 + side][v] = process_edge[2][2 + side] = true;
                }
              if (face_y)
                {
                  const unsigned int side = !subcell_y;

                  mask_face[0][side][v] = process_face[0][side] = true;
                  mask_edge[0][side][v] = process_edge[0][side] = true;
                  mask_edge[0][2 + side][v] = process_edge[0][2 + side] = true;
                  mask_face[2][2 + side][v] = process_face[2][2 + side] = true;
                  mask_edge[2][2 * side][v] = process_edge[2][2 * side] = true;
                  mask_edge[2][2 * side + 1][v] =
                    process_edge[2][2 * side + 1] = true;
                }
              if (face_z)
                {
                  const unsigned int side = !subcell_z;

                  mask_face[0][2 + side][v] = process_face[0][2 + side] = true;
                  mask_edge[0][2 * side][v] = process_edge[0][2 * side] = true;
                  mask_edge[0][2 * side + 1][v] =
                    process_edge[0][2 * side + 1] = true;
                  mask_face[1][2 + side][v] = process_face[1][2 + side] = true;
                  mask_edge[1][2 * side][v] = process_edge[1][2 * side] = true;
                  mask_edge[1][2 * side + 1][v] =
                    process_edge[1][2 * side + 1] = true;
                }
              if (edge_x)
                {
                  const unsigned int index = (!subcell_z) * 2 + (!subcell_y);
                  mask_edge[0][index][v] = process_edge[0][index] = true;
                }
              if (edge_y)
                {
                  const unsigned int index = (!subcell_z) * 2 + (!subcell_x);
                  mask_edge[1][index][v] = process_edge[1][index] = true;
                }
              if (edge_z)
                {
                  const unsigned int index = (!subcell_y) * 2 + (!subcell_x);
                  mask_edge[2][index][v] = process_edge[2][index] = true;
                }
            }

          // direction 0:
          if (given_degree > 1)
            {
              const std::array<unsigned int, 4> face_offsets = {
                {p0, p2, p0, p4}};
              const std::array<unsigned int, 2> outer_strides = {
                {points * points, points}};
              for (unsigned int c = 0; c < n_components; ++c)
                for (unsigned int face = 0; face < 4; ++face)
                  if (process_face[0][face])
                    interpolate<2, 0, transpose>(face_offsets[face],
                                                 outer_strides[face / 2],
                                                 given_degree,
                                                 mask_weights[0],
                                                 mask_face[0][face],
                                                 weights,
                                                 values + c * n_dofs);
            }
          {
            const std::array<unsigned int, 4> edge_offsets = {{p0, p2, p4, p6}};
            for (unsigned int c = 0; c < n_components; ++c)
              for (unsigned int edge = 0; edge < 4; ++edge)
                if (process_edge[0][edge])
                  interpolate<1, 0, transpose>(edge_offsets[edge],
                                               0,
                                               given_degree,
                                               mask_weights[0],
                                               mask_edge[0][edge],
                                               weights,
                                               values + c * n_dofs);
          }

          // direction 1:
          if (given_degree > 1)
            {
              const std::array<unsigned int, 4> face_offsets = {
                {p0, p1, p0, p4}};
              const std::array<unsigned int, 2> outer_strides = {
                {points * points, 1}};
              for (unsigned int c = 0; c < n_components; ++c)
                for (unsigned int face = 0; face < 4; ++face)
                  if (process_face[1][face])
                    interpolate<2, 1, transpose>(face_offsets[face],
                                                 outer_strides[face / 2],
                                                 given_degree,
                                                 mask_weights[1],
                                                 mask_face[1][face],
                                                 weights,
                                                 values + c * n_dofs);
            }

          {
            const std::array<unsigned int, 4> edge_offsets = {{p0, p1, p4, p5}};
            for (unsigned int c = 0; c < n_components; ++c)
              for (unsigned int edge = 0; edge < 4; ++edge)
                if (process_edge[1][edge])
                  interpolate<1, 1, transpose>(edge_offsets[edge],
                                               0,
                                               given_degree,
                                               mask_weights[1],
                                               mask_edge[1][edge],
                                               weights,
                                               values + c * n_dofs);
          }

          // direction 2:
          if (given_degree > 1)
            {
              const std::array<unsigned int, 4> face_offsets = {
                {p0, p1, p0, p2}};
              const std::array<unsigned int, 2> outer_strides = {{points, 1}};
              for (unsigned int c = 0; c < n_components; ++c)
                for (unsigned int face = 0; face < 4; ++face)
                  if (process_face[2][face])
                    interpolate<2, 2, transpose>(face_offsets[face],
                                                 outer_strides[face / 2],
                                                 given_degree,
                                                 mask_weights[2],
                                                 mask_face[2][face],
                                                 weights,
                                                 values + c * n_dofs);
            }

          {
            const std::array<unsigned int, 4> edge_offsets = {{p0, p1, p2, p3}};
            for (unsigned int c = 0; c < n_components; ++c)
              for (unsigned int edge = 0; edge < 4; ++edge)
                if (process_edge[2][edge])
                  interpolate<1, 2, transpose>(edge_offsets[edge],
                                               0,
                                               given_degree,
                                               mask_weights[2],
                                               mask_edge[2][edge],
                                               weights,
                                               values + c * n_dofs);
          }
        }
      else
        {
          Assert(false, ExcNotImplemented());
        }
    }
  };



  template <int dim, typename Number, bool is_face>
  struct FEEvaluationImplHangingNodes
  {
  public:
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
      using RunnerType = FEEvaluationImplHangingNodesRunner<
        FEEvaluationImplHangingNodesRunnerTypes::vectorized,
        dim,
        fe_degree,
        Number,
        is_face>;

      if (transpose)
        RunnerType::template run_internal<true>(n_desired_components,
                                                fe_eval,
                                                c_mask,
                                                values);
      else
        RunnerType::template run_internal<false>(n_desired_components,
                                                 fe_eval,
                                                 c_mask,
                                                 values);

      return false;
    }
  };


} // end of namespace internal


DEAL_II_NAMESPACE_CLOSE

#endif
