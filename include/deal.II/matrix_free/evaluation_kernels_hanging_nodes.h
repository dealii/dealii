// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#ifndef dealii_matrix_free_evaluation_kernels_hanging_nodes_h
#define dealii_matrix_free_evaluation_kernels_hanging_nodes_h

#include <deal.II/base/config.h>

#include <deal.II/base/ndarray.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/matrix_free/hanging_nodes_internal.h>
#include <deal.II/matrix_free/shape_info.h>


DEAL_II_NAMESPACE_OPEN

#ifdef DEBUG
#  define DEAL_II_ALWAYS_INLINE_RELEASE
#else
#  define DEAL_II_ALWAYS_INLINE_RELEASE DEAL_II_ALWAYS_INLINE
#endif



namespace internal
{
  enum class FEEvaluationImplHangingNodesRunnerTypes
  {
    scalar,
    vectorized
  };



  /**
   * Helper enum to specify the type of vectorization for
   * FEEvaluationImplHangingNodesRunnerTypes::scalar.
   */
  enum class VectorizationTypes
  {
    /**
     * Process cell by cell.
     */
    index,
    /**
     * Process cells with the same refinement configuration together.
     */
    group,
    /**
     * Like index but without access to individual lanes and instead use masks
     * with a single entry with the value one.
     */
    mask,
    /**
     * Assume that all lanes have the same refinement configuration.
     */
    sorted
  };



  template <FEEvaluationImplHangingNodesRunnerTypes,
            int dim,
            int fe_degree,
            typename Number>
  class FEEvaluationImplHangingNodesRunner;



  template <int dim, int fe_degree, typename Number>
  class FEEvaluationImplHangingNodesRunner<
    FEEvaluationImplHangingNodesRunnerTypes::vectorized,
    dim,
    fe_degree,
    Number>
  {
  private:
    template <int          structdim,
              unsigned int direction,
              bool         transpose,
              typename Number2>
    static void
    interpolate(const unsigned int              offset,
                const unsigned int              outer_stride,
                const unsigned int              given_degree,
                const Number                    mask_weight,
                const Number                    mask_write,
                const Number2 *DEAL_II_RESTRICT weights,
                Number *DEAL_II_RESTRICT        values)
    {
      static constexpr unsigned int max_n_points_1D = 40;

      static_assert(structdim == 1 || structdim == 2,
                    "Only 1D and 2d interpolation implemented");
      Number temp[fe_degree != -1 ? fe_degree + 1 : max_n_points_1D];

      const unsigned int points =
        (fe_degree != -1 ? fe_degree : given_degree) + 1;

      AssertIndexRange(points, max_n_points_1D);

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
    template <bool transpose, typename Number2>
    static void
    run_internal(
      const unsigned int                             n_components,
      const MatrixFreeFunctions::ShapeInfo<Number2> &shape_info,
      const std::array<MatrixFreeFunctions::compressed_constraint_kind,
                       Number::size()>              &constraint_mask,
      Number                                        *values)
    {
      const unsigned int given_degree =
        fe_degree != -1 ? fe_degree : shape_info.data.front().fe_degree;

      const Number2 *DEAL_II_RESTRICT weights =
        shape_info.data.front().subface_interpolation_matrices[0].data();

      const unsigned int points = given_degree + 1;
      const unsigned int n_dofs = shape_info.dofs_per_component_on_cell;

      if (dim == 2)
        {
          dealii::ndarray<Number, 2>    mask_weights = {};
          dealii::ndarray<Number, 2, 2> mask_write   = {};
          dealii::ndarray<bool, 2, 2>   do_face      = {};

          for (unsigned int v = 0; v < Number::size(); ++v)
            {
              const auto kind      = constraint_mask[v];
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
              const auto kind = constraint_mask[v];

              const bool subcell_x = (kind >> 0) & 1;
              const bool subcell_y = (kind >> 1) & 1;
              const bool subcell_z = (kind >> 2) & 1;
              const bool face_x    = ((kind >> 3) & 1) ? (kind >> 5) & 1 : 0;
              const bool face_y    = ((kind >> 3) & 1) ? (kind >> 6) & 1 : 0;
              const bool face_z    = ((kind >> 3) & 1) ? (kind >> 7) & 1 : 0;
              const bool edge_x    = ((kind >> 4) & 1) ? (kind >> 5) & 1 : 0;
              const bool edge_y    = ((kind >> 4) & 1) ? (kind >> 6) & 1 : 0;
              const bool edge_z    = ((kind >> 4) & 1) ? (kind >> 7) & 1 : 0;

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
          DEAL_II_NOT_IMPLEMENTED();
        }
    }
  };

  template <typename T1, VectorizationTypes VT>
  struct Trait;

  template <typename T1>
  struct Trait<T1, VectorizationTypes::index>
  {
    using value_type         = typename T1::value_type;
    using index_type         = unsigned int;
    using interpolation_type = value_type;

    template <typename T>
    static inline const std::array<AlignedVector<interpolation_type>, 2> &
    get_interpolation_matrix(const T &shape_info)
    {
      return shape_info.data.front().subface_interpolation_matrices_scalar;
    }

    static inline DEAL_II_ALWAYS_INLINE_RELEASE unsigned int
    create(const std::array<MatrixFreeFunctions::compressed_constraint_kind,
                            T1::size()> mask,
           const std::array<MatrixFreeFunctions::compressed_constraint_kind,
                            T1::size()> mask_new,
           const unsigned int           v)
    {
      (void)mask;
      (void)mask_new;
      return v;
    }

    static inline DEAL_II_ALWAYS_INLINE_RELEASE bool
    do_break(unsigned int                                           v,
             const MatrixFreeFunctions::compressed_constraint_kind &kind)
    {
      (void)v;
      (void)kind;
      return false;
    }

    static inline DEAL_II_ALWAYS_INLINE_RELEASE bool
    do_continue(unsigned int                                           v,
                const MatrixFreeFunctions::compressed_constraint_kind &kind)
    {
      (void)v;
      return kind ==
             MatrixFreeFunctions::unconstrained_compressed_constraint_kind;
    }

    static inline DEAL_II_ALWAYS_INLINE_RELEASE
      std::array<MatrixFreeFunctions::compressed_constraint_kind, T1::size()>
      create_mask(
        const std::array<MatrixFreeFunctions::compressed_constraint_kind,
                         T1::size()> mask)
    {
      return mask;
    }

    static inline DEAL_II_ALWAYS_INLINE_RELEASE typename T1::value_type
    get_value(const typename T1::value_type &value, const index_type &i)
    {
      (void)i;
      return value;
    }

    static inline DEAL_II_ALWAYS_INLINE_RELEASE typename T1::value_type
    get_value(const T1 &value, const index_type &i)
    {
      return value[i];
    }

    static inline DEAL_II_ALWAYS_INLINE_RELEASE void
    set_value(T1                            &result,
              const typename T1::value_type &value,
              const index_type              &i)
    {
      result[i] = value;
    }
  };

  template <typename T1>
  struct Trait<T1, VectorizationTypes::mask>
  {
    using value_type         = T1;
    using index_type         = std::pair<T1, T1>;
    using interpolation_type = T1;

    template <typename T>
    static inline const std::array<AlignedVector<T1>, 2> &
    get_interpolation_matrix(const T &shape_info)
    {
      return shape_info.data.front().subface_interpolation_matrices;
    }

    static inline DEAL_II_ALWAYS_INLINE_RELEASE bool
    do_break(unsigned int                                           v,
             const MatrixFreeFunctions::compressed_constraint_kind &kind)
    {
      (void)v;
      (void)kind;
      return false;
    }

    static inline DEAL_II_ALWAYS_INLINE_RELEASE bool
    do_continue(unsigned int                                           v,
                const MatrixFreeFunctions::compressed_constraint_kind &kind)
    {
      (void)v;
      return kind ==
             MatrixFreeFunctions::unconstrained_compressed_constraint_kind;
    }

    static inline DEAL_II_ALWAYS_INLINE_RELEASE index_type
    create(const std::array<MatrixFreeFunctions::compressed_constraint_kind,
                            T1::size()> mask,
           const std::array<MatrixFreeFunctions::compressed_constraint_kind,
                            T1::size()> mask_new,
           const unsigned int           v)
    {
      (void)mask;
      (void)mask_new;
      T1 result = 0.0;
      result[v] = 1.0;
      return {result, T1(1.0) - result};
    }

    static inline DEAL_II_ALWAYS_INLINE_RELEASE
      std::array<MatrixFreeFunctions::compressed_constraint_kind, T1::size()>
      create_mask(
        const std::array<MatrixFreeFunctions::compressed_constraint_kind,
                         T1::size()> mask)
    {
      return mask;
    }

    static inline DEAL_II_ALWAYS_INLINE_RELEASE T1
    get_value(const T1 &value, const index_type &)
    {
      return value;
    }

    static inline DEAL_II_ALWAYS_INLINE_RELEASE void
    set_value(T1 &result, const T1 &value, const index_type &i)
    {
      result = result * i.second + value * i.first;
    }
  };

  template <typename T1>
  struct Trait<T1, VectorizationTypes::group>
  {
    using value_type         = T1;
    using index_type         = std::pair<T1, T1>;
    using interpolation_type = T1;

    template <typename T>
    static inline const std::array<AlignedVector<T1>, 2> &
    get_interpolation_matrix(const T &shape_info)
    {
      return shape_info.data.front().subface_interpolation_matrices;
    }

    static inline DEAL_II_ALWAYS_INLINE_RELEASE bool
    do_break(unsigned int                                           v,
             const MatrixFreeFunctions::compressed_constraint_kind &kind)
    {
      (void)v;
      return kind ==
             MatrixFreeFunctions::unconstrained_compressed_constraint_kind;
    }

    static inline DEAL_II_ALWAYS_INLINE_RELEASE bool
    do_continue(unsigned int                                           v,
                const MatrixFreeFunctions::compressed_constraint_kind &kind)
    {
      (void)v;
      return kind ==
             MatrixFreeFunctions::unconstrained_compressed_constraint_kind;
    }

    static inline DEAL_II_ALWAYS_INLINE_RELEASE index_type
    create(const std::array<MatrixFreeFunctions::compressed_constraint_kind,
                            T1::size()> mask,
           const std::array<MatrixFreeFunctions::compressed_constraint_kind,
                            T1::size()> mask_new,
           const unsigned int           v)
    {
      T1 result;

      for (unsigned int i = 0; i < T1::size(); ++i)
        result[i] = mask_new[v] == mask[i];

      return {result, T1(1.0) - result};
    }

    static inline DEAL_II_ALWAYS_INLINE_RELEASE
      std::array<MatrixFreeFunctions::compressed_constraint_kind, T1::size()>
      create_mask(
        const std::array<MatrixFreeFunctions::compressed_constraint_kind,
                         T1::size()> mask)
    {
      auto new_mask = mask;

      std::sort(new_mask.begin(), new_mask.end());
      std::fill(std::unique(new_mask.begin(), new_mask.end()),
                new_mask.end(),
                MatrixFreeFunctions::unconstrained_compressed_constraint_kind);

      return new_mask;
    }

    static inline DEAL_II_ALWAYS_INLINE_RELEASE T1
    get_value(const T1 &value, const index_type &)
    {
      return value;
    }

    static inline DEAL_II_ALWAYS_INLINE_RELEASE void
    set_value(T1 &result, const T1 &value, const index_type &i)
    {
      result = result * i.second + value * i.first;
    }
  };

  template <typename T1>
  struct Trait<T1, VectorizationTypes::sorted>
  {
    using value_type         = T1;
    using index_type         = T1;
    using interpolation_type = T1;

    template <typename T>
    static inline const std::array<AlignedVector<T1>, 2> &
    get_interpolation_matrix(const T &shape_info)
    {
      return shape_info.data.front().subface_interpolation_matrices;
    }

    static inline DEAL_II_ALWAYS_INLINE_RELEASE bool
    do_break(unsigned int                                           v,
             const MatrixFreeFunctions::compressed_constraint_kind &kind)
    {
      (void)kind;
      return v > 0;
    }

    static inline DEAL_II_ALWAYS_INLINE_RELEASE bool
    do_continue(unsigned int                                           v,
                const MatrixFreeFunctions::compressed_constraint_kind &kind)
    {
      (void)kind;

      DEAL_II_ASSERT_UNREACHABLE();

      return v > 0; // should not be called
    }

    static inline DEAL_II_ALWAYS_INLINE_RELEASE T1
    create(const std::array<MatrixFreeFunctions::compressed_constraint_kind,
                            T1::size()> mask,
           const std::array<MatrixFreeFunctions::compressed_constraint_kind,
                            T1::size()> mask_new,
           const unsigned int           v)
    {
      (void)mask;
      (void)mask_new;
      (void)v;
      return 1.0; // return something since not used
    }

    static inline DEAL_II_ALWAYS_INLINE_RELEASE
      std::array<MatrixFreeFunctions::compressed_constraint_kind, T1::size()>
      create_mask(
        const std::array<MatrixFreeFunctions::compressed_constraint_kind,
                         T1::size()> mask)
    {
      return mask;
    }

    static inline DEAL_II_ALWAYS_INLINE_RELEASE T1
    get_value(const T1 &value, const index_type &)
    {
      return value;
    }

    static inline DEAL_II_ALWAYS_INLINE_RELEASE void
    set_value(T1 &result, const T1 &value, const index_type &i)
    {
      (void)i;
      result = value;
    }
  };



  template <typename T,
            typename Number,
            VectorizationTypes VectorizationType,
            int                fe_degree,
            bool               transpose>
  class HelperBase
  {
  public:
    inline DEAL_II_ALWAYS_INLINE_RELEASE
    HelperBase(
      const T                                                     &t,
      const unsigned int                                          &given_degree,
      const bool                                                  &type_x,
      const bool                                                  &type_y,
      const bool                                                  &type_z,
      const typename Trait<Number, VectorizationType>::index_type &v,
      const std::array<
        AlignedVector<
          typename Trait<Number, VectorizationType>::interpolation_type>,
        2>   &interpolation_matrices,
      Number *values)
      : t(t)
      , given_degree(given_degree)
      , type_x(type_x)
      , type_y(type_y)
      , type_z(type_z)
      , v(v)
      , interpolation_matrices(interpolation_matrices)
      , values(values)
    {}

    template <unsigned int direction, unsigned int d, bool skip_borders>
    static inline DEAL_II_ALWAYS_INLINE_RELEASE void
    interpolate_3D_face(
      const unsigned int                                          dof_offset,
      const unsigned int                                          given_degree,
      const typename Trait<Number, VectorizationType>::index_type v,
      const typename Trait<Number, VectorizationType>::interpolation_type
        *DEAL_II_RESTRICT      weight,
      Number *DEAL_II_RESTRICT values)
    {
      static constexpr unsigned int max_n_points_1D = 40;

      typename Trait<Number, VectorizationType>::value_type
        temp[fe_degree != -1 ? (fe_degree + 1) : max_n_points_1D];

      const unsigned int points =
        (fe_degree != -1 ? fe_degree : given_degree) + 1;

      AssertIndexRange(given_degree, max_n_points_1D);

      const unsigned int stride = fe_degree != -1 ?
                                    Utilities::pow(fe_degree + 1, direction) :
                                    Utilities::pow(given_degree + 1, direction);

      // direction   side0   side1   side2
      // 0             -      p^2      p
      // 1            p^2      -       1
      // 2             p       -       1
      const unsigned int stride2 =
        ((direction == 0 && d == 1) || (direction == 1 && d == 0)) ?
          (points * points) :
          (((direction == 0 && d == 2) || (direction == 2 && d == 0)) ? points :
                                                                        1);

      for (unsigned int g = (skip_borders ? 1 : 0);
           g < points - (skip_borders ? 1 : 0);
           ++g)
        {
          // copy result back
          for (unsigned int k = 0; k < points; ++k)
            temp[k] = Trait<Number, VectorizationType>::get_value(
              values[dof_offset + k * stride + stride2 * g], v);

          // perform interpolation point by point
          for (unsigned int k = 0; k < points; ++k)
            {
              auto sum = Trait<Number, VectorizationType>::get_value(
                           weight[(transpose ? 1 : points) * k], v) *
                         temp[0];
              for (unsigned int h = 1; h < points; ++h)
                sum += Trait<Number, VectorizationType>::get_value(
                         weight[(transpose ? 1 : points) * k +
                                (transpose ? points : 1) * h],
                         v) *
                       temp[h];
              Trait<Number, VectorizationType>::set_value(
                values[dof_offset + k * stride + stride2 * g], sum, v);
            }
        }
    }

    template <unsigned int direction>
    static inline DEAL_II_ALWAYS_INLINE_RELEASE void
    interpolate_3D_edge(
      const unsigned int                                          p,
      const unsigned int                                          given_degree,
      const typename Trait<Number, VectorizationType>::index_type v,
      const typename Trait<Number, VectorizationType>::interpolation_type
        *DEAL_II_RESTRICT      weight,
      Number *DEAL_II_RESTRICT values)
    {
      static constexpr unsigned int max_n_points_1D = 40;

      typename Trait<Number, VectorizationType>::value_type
        temp[fe_degree != -1 ? (fe_degree + 1) : max_n_points_1D];

      const unsigned int points =
        (fe_degree != -1 ? fe_degree : given_degree) + 1;

      AssertIndexRange(given_degree, max_n_points_1D);

      const unsigned int stride = fe_degree != -1 ?
                                    Utilities::pow(fe_degree + 1, direction) :
                                    Utilities::pow(given_degree + 1, direction);

      // copy result back
      for (unsigned int k = 0; k < points; ++k)
        temp[k] =
          Trait<Number, VectorizationType>::get_value(values[p + k * stride],
                                                      v);

      // perform interpolation point by point
      for (unsigned int k = 0; k < points; ++k)
        {
          auto sum = Trait<Number, VectorizationType>::get_value(
                       weight[(transpose ? 1 : points) * k], v) *
                     temp[0];
          for (unsigned int h = 1; h < points; ++h)
            sum += Trait<Number, VectorizationType>::get_value(
                     weight[(transpose ? 1 : points) * k +
                            (transpose ? points : 1) * h],
                     v) *
                   temp[h];
          Trait<Number, VectorizationType>::set_value(values[p + k * stride],
                                                      sum,
                                                      v);
        }
    }

    template <bool do_x, bool do_y, bool do_z>
    inline DEAL_II_ALWAYS_INLINE_RELEASE void
    process_edge() const
    {
      if (do_x)
        interpolate_3D_edge<0>(t.line(0, type_y, type_z),
                               given_degree,
                               v,
                               interpolation_matrices[!type_x].data(),
                               values);

      if (do_y)
        interpolate_3D_edge<1>(t.line(1, type_x, type_z),
                               given_degree,
                               v,
                               interpolation_matrices[!type_y].data(),
                               values);

      if (do_z)
        interpolate_3D_edge<2>(t.line(2, type_x, type_y),
                               given_degree,
                               v,
                               interpolation_matrices[!type_z].data(),
                               values);
    }

    template <bool do_x, bool do_y, bool do_z>
    inline DEAL_II_ALWAYS_INLINE_RELEASE void
    process_faces_fast() const
    {
      static_assert((do_x && !do_y && !do_z) || (!do_x && do_y && !do_z) ||
                      (!do_x && !do_y && do_z),
                    "Only one face can be chosen.");

      static const unsigned int direction = do_x ? 0 : (do_y ? 1 : 2);
      const bool                type = do_x ? type_x : (do_y ? type_y : type_z);

      if (!do_x)
        interpolate_3D_face<0, direction, false>(
          t.face(direction, type),
          given_degree,
          v,
          interpolation_matrices[!type_x].data(),
          values);

      if (!do_y)
        interpolate_3D_face<1, direction, false>(
          t.face(direction, type),
          given_degree,
          v,
          interpolation_matrices[!type_y].data(),
          values);

      if (!do_z)
        interpolate_3D_face<2, direction, false>(
          t.face(direction, type),
          given_degree,
          v,
          interpolation_matrices[!type_z].data(),
          values);
    }

    template <bool do_x, bool do_y, bool do_z>
    inline DEAL_II_ALWAYS_INLINE_RELEASE void
    process_faces() const
    {
      static_assert(((do_x && !do_y && !do_z) || (!do_x && do_y && !do_z) ||
                     (!do_x && !do_y && do_z)) == false,
                    "Only one face can be chosen.");

      // direction 0
      {
        const auto inpterolation_matrix =
          interpolation_matrices[!type_x].data();

        // faces
        if (do_y && given_degree > 1)
          interpolate_3D_face<0, 1, true>(
            t.face(1, type_y), given_degree, v, inpterolation_matrix, values);

        if (do_z && given_degree > 1)
          interpolate_3D_face<0, 2, true>(
            t.face(2, type_z), given_degree, v, inpterolation_matrix, values);

        // direction 0 -> edges
        interpolate_3D_edge<0>((do_x && do_y && !do_z) ?
                                 (t.lines_plane(0, type_x, type_y, 0)) :
                                 ((do_x && !do_y && do_z) ?
                                    (t.lines_plane(1, type_x, type_z, 0)) :
                                    (t.lines(0, type_y, type_z, 0))),
                               given_degree,
                               v,
                               inpterolation_matrix,
                               values);


        interpolate_3D_edge<0>((do_x && do_y && !do_z) ?
                                 (t.lines_plane(0, type_x, type_y, 1)) :
                                 ((do_x && !do_y && do_z) ?
                                    (t.lines_plane(1, type_x, type_z, 1)) :
                                    (t.lines(0, type_y, type_z, 1))),
                               given_degree,
                               v,
                               inpterolation_matrix,
                               values);

        if (do_y && do_z)
          interpolate_3D_edge<0>(t.lines(0, type_y, type_z, 2),
                                 given_degree,
                                 v,
                                 inpterolation_matrix,
                                 values);
      }

      // direction 1
      {
        const auto inpterolation_matrix =
          interpolation_matrices[!type_y].data();

        // faces
        if (do_x && given_degree > 1)
          interpolate_3D_face<1, 0, true>(
            t.face(0, type_x), given_degree, v, inpterolation_matrix, values);

        if (do_z && given_degree > 1)
          interpolate_3D_face<1, 2, true>(
            t.face(2, type_z), given_degree, v, inpterolation_matrix, values);

        // lines
        interpolate_3D_edge<1>((do_x && do_y && !do_z) ?
                                 (t.lines_plane(0, type_x, type_y, 2)) :
                                 ((!do_x && do_y && do_z) ?
                                    (t.lines_plane(2, type_y, type_z, 0)) :
                                    (t.lines(1, type_x, type_z, 0))),
                               given_degree,
                               v,
                               inpterolation_matrix,
                               values);

        interpolate_3D_edge<1>((do_x && do_y && !do_z) ?
                                 (t.lines_plane(0, type_x, type_y, 3)) :
                                 ((!do_x && do_y && do_z) ?
                                    (t.lines_plane(2, type_y, type_z, 1)) :
                                    (t.lines(1, type_x, type_z, 1))),
                               given_degree,
                               v,
                               inpterolation_matrix,
                               values);

        if (do_x && do_z)
          interpolate_3D_edge<1>(t.lines(1, type_x, type_z, 2),
                                 given_degree,
                                 v,
                                 inpterolation_matrix,
                                 values);
      }

      // direction 2 -> faces
      {
        const auto inpterolation_matrix =
          interpolation_matrices[!type_z].data();

        if (do_x && given_degree > 1)
          interpolate_3D_face<2, 0, true>(
            t.face(0, type_x), given_degree, v, inpterolation_matrix, values);

        if (do_y && given_degree > 1)
          interpolate_3D_face<2, 1, true>(
            t.face(1, type_y), given_degree, v, inpterolation_matrix, values);

        // direction 2 -> edges
        interpolate_3D_edge<2>((do_x && !do_y && do_z) ?
                                 (t.lines_plane(1, type_x, type_z, 2)) :
                                 ((!do_x && do_y && do_z) ?
                                    (t.lines_plane(2, type_y, type_z, 2)) :
                                    (t.lines(2, type_x, type_y, 0))),
                               given_degree,
                               v,
                               inpterolation_matrix,
                               values);

        interpolate_3D_edge<2>((do_x && !do_y && do_z) ?
                                 (t.lines_plane(1, type_x, type_z, 3)) :
                                 ((!do_x && do_y && do_z) ?
                                    (t.lines_plane(2, type_y, type_z, 3)) :
                                    (t.lines(2, type_x, type_y, 1))),
                               given_degree,
                               v,
                               inpterolation_matrix,
                               values);

        if (do_x && do_y)
          interpolate_3D_edge<2>(t.lines(2, type_x, type_y, 2),
                                 given_degree,
                                 v,
                                 inpterolation_matrix,
                                 values);
      }
    }

  private:
    const T                                                     &t;
    const unsigned int                                          &given_degree;
    const bool                                                  &type_x;
    const bool                                                  &type_y;
    const bool                                                  &type_z;
    const typename Trait<Number, VectorizationType>::index_type &v;
    const std::array<
      AlignedVector<
        typename Trait<Number, VectorizationType>::interpolation_type>,
      2>   &interpolation_matrices;
    Number *values;
  };

  /**
   * Helper enum to specify which Helper implementation should be used.
   */
  enum class HelperType
  {
    /**
     * Compute the start indices of faces and edges based on the template
     * argument fe_degree.
     */
    constant,
    /**
     * Compute the start indices of faces and edges based on the fe_degree
     * passed to the constructor (to be used if the template argument is -1).
     */
    dynamic
  };

  template <HelperType helper_type,
            typename Number,
            VectorizationTypes VectorizationType,
            int                fe_degree,
            bool               transpose>
  class Helper;

  template <typename Number,
            VectorizationTypes VectorizationType,
            int                fe_degree,
            bool               transpose>
  class Helper<HelperType::dynamic,
               Number,
               VectorizationType,
               fe_degree,
               transpose> : public HelperBase<Helper<HelperType::dynamic,
                                                     Number,
                                                     VectorizationType,
                                                     fe_degree,
                                                     transpose>,
                                              Number,
                                              VectorizationType,
                                              fe_degree,
                                              transpose>
  {
  public:
    inline DEAL_II_ALWAYS_INLINE_RELEASE
    Helper(const unsigned int &given_degree,
           const bool         &type_x,
           const bool         &type_y,
           const bool         &type_z,
           const typename Trait<Number, VectorizationType>::index_type &v,
           const std::array<
             AlignedVector<
               typename Trait<Number, VectorizationType>::interpolation_type>,
             2>   &interpolation_matrices,
           Number *values)
      : HelperBase<Helper<HelperType::dynamic,
                          Number,
                          VectorizationType,
                          fe_degree,
                          transpose>,
                   Number,
                   VectorizationType,
                   fe_degree,
                   transpose>(*this,
                              given_degree,
                              type_x,
                              type_y,
                              type_z,
                              v,
                              interpolation_matrices,
                              values)
      , points(given_degree + 1)
    {
      static_assert(fe_degree == -1, "Only working for fe_degree = -1.");
    }

    const unsigned int points;

    inline DEAL_II_ALWAYS_INLINE_RELEASE unsigned int
    line(unsigned int i, unsigned int j, unsigned int k) const
    {
      return line_array[i][j][k];
    }

    inline DEAL_II_ALWAYS_INLINE_RELEASE unsigned int
    face(unsigned int i, unsigned int j) const
    {
      return face_array[i][j];
    }

    inline DEAL_II_ALWAYS_INLINE_RELEASE unsigned int
    lines_plane(unsigned int i,
                unsigned int j,
                unsigned int k,
                unsigned int l) const
    {
      return lines_plane_array[i][j][k][l];
    }

    inline DEAL_II_ALWAYS_INLINE_RELEASE unsigned int
    lines(unsigned int i, unsigned int j, unsigned int k, unsigned int l) const
    {
      return lines_array[i][j][k][l];
    }

  private:
    const dealii::ndarray<unsigned int, 3, 2, 2> line_array = {
      {{{{{points * points * points - points, points *points - points}},
         {{points * points * points - points * points, 0}}}},
       {{{{points * points * points - points * points + points - 1,
           points - 1}},
         {{points * points * points - points * points, 0}}}},
       {{{{points * points - 1, points - 1}},
         {{points * points - points, 0}}}}}};

    const dealii::ndarray<unsigned int, 3, 2> face_array = {
      {{{points - 1, 0}},
       {{points * points - points, 0}},
       {{points * points * points - points * points, 0}}}};

    const dealii::ndarray<unsigned int, 3, 2, 2, 4> lines_plane_array = {
      {{{{{{{points * points - points,
             points *points *points - points,
             points - 1,
             points *points *points - points *points + points - 1}},
           {{0,
             points *points *points - points *points,
             points - 1,
             points *points *points - points *points + points - 1}}}},
         {{{{points * points - points,
             points *points *points - points,
             0,
             points *points *points - points *points}},
           {{0,
             points *points *points - points *points,
             0,
             points *points *points - points *points}}}}}},
       {{{{{{points * points * points - points * points,
             points *points *points - points,
             points - 1,
             points *points - 1}},
           {{0, points *points - points, points - 1, points *points - 1}}}},
         {{{{points * points * points - points * points,
             points *points *points - points,
             0,
             points *points - points}},
           {{0, points *points - points, 0, points *points - points}}}}}},
       {{{{{{points * points * points - points * points,
             points *points *points - points *points + points - 1,
             points                          *points - points,
             points                          *points - 1}},
           {{0, points - 1, points *points - points, points *points - 1}}}},
         {{{{points * points * points - points * points,
             points *points *points - points *points + points - 1,
             0,
             points - 1}},
           {{0, points - 1, 0, points - 1}}}}}}}};

    const dealii::ndarray<unsigned int, 3, 2, 2, 3> lines_array = {
      {{{{{{{points * points - points,
             points *points *points - points *points,
             points *points                  *points - points}},
           {{0, points *points - points, points *points *points - points}}}},
         {{{{0,
             points *points *points - points *points,
             points *points                  *points - points}},
           {{0,
             points         *points - points,
             points *points *points - points *points}}}}}},
       {{{{{{points - 1,
             points *points *points - points *points,
             points *points *points - points *points + points - 1}},
           {{0,
             points - 1,
             points *points *points - points *points + points - 1}}}},
         {{{{0,
             points *points *points - points *points,
             points *points *points - points *points + points - 1}},
           {{0, points - 1, points *points *points - points *points}}}}}},
       {{{{{{points - 1, points *points - points, points *points - 1}},
           {{0, points - 1, points *points - 1}}}},
         {{{{0, points *points - points, points *points - 1}},
           {{0, points - 1, points *points - points}}}}}}}};
  };

  template <typename Number,
            VectorizationTypes VectorizationType,
            int                fe_degree,
            bool               transpose>
  class Helper<HelperType::constant,
               Number,
               VectorizationType,
               fe_degree,
               transpose> : public HelperBase<Helper<HelperType::constant,
                                                     Number,
                                                     VectorizationType,
                                                     fe_degree,
                                                     transpose>,
                                              Number,
                                              VectorizationType,
                                              fe_degree,
                                              transpose>
  {
  public:
    inline DEAL_II_ALWAYS_INLINE_RELEASE
    Helper(const unsigned int &given_degree,
           const bool         &type_x,
           const bool         &type_y,
           const bool         &type_z,
           const typename Trait<Number, VectorizationType>::index_type &v,
           const std::array<
             AlignedVector<
               typename Trait<Number, VectorizationType>::interpolation_type>,
             2>   &interpolation_matrices,
           Number *values)
      : HelperBase<Helper<HelperType::constant,
                          Number,
                          VectorizationType,
                          fe_degree,
                          transpose>,
                   Number,
                   VectorizationType,
                   fe_degree,
                   transpose>(*this,
                              given_degree,
                              type_x,
                              type_y,
                              type_z,
                              v,
                              interpolation_matrices,
                              values)
    {
      static_assert(fe_degree != -1, "Only working for fe_degree != -1.");
    }


    inline DEAL_II_ALWAYS_INLINE_RELEASE unsigned int
    line(unsigned int i, unsigned int j, unsigned int k) const
    {
      static constexpr unsigned int points = fe_degree + 1;

      static constexpr dealii::ndarray<unsigned int, 3, 2, 2> line_array = {
        {{{{{points * points * points - points, points * points - points}},
           {{points * points * points - points * points, 0}}}},
         {{{{points * points * points - points * points + points - 1,
             points - 1}},
           {{points * points * points - points * points, 0}}}},
         {{{{points * points - 1, points - 1}},
           {{points * points - points, 0}}}}}};

      return line_array[i][j][k];
    }

    inline DEAL_II_ALWAYS_INLINE_RELEASE unsigned int
    face(unsigned int i, unsigned int j) const
    {
      static constexpr unsigned int points = fe_degree + 1;

      static constexpr dealii::ndarray<unsigned int, 3, 2> face_array = {
        {{{points - 1, 0}},
         {{points * points - points, 0}},
         {{points * points * points - points * points, 0}}}};

      return face_array[i][j];
    }

    inline DEAL_II_ALWAYS_INLINE_RELEASE unsigned int
    lines_plane(unsigned int i,
                unsigned int j,
                unsigned int k,
                unsigned int l) const
    {
      static constexpr unsigned int points = fe_degree + 1;

      static constexpr dealii::ndarray<unsigned int, 3, 2, 2, 4>
        lines_plane_array = {
          {{{{{{{points * points - points,
                 points * points * points - points,
                 points - 1,
                 points * points * points - points * points + points - 1}},
               {{0,
                 points * points * points - points * points,
                 points - 1,
                 points * points * points - points * points + points - 1}}}},
             {{{{points * points - points,
                 points * points * points - points,
                 0,
                 points * points * points - points * points}},
               {{0,
                 points * points * points - points * points,
                 0,
                 points * points * points - points * points}}}}}},
           {{{{{{points * points * points - points * points,
                 points * points * points - points,
                 points - 1,
                 points * points - 1}},
               {{0,
                 points * points - points,
                 points - 1,
                 points * points - 1}}}},
             {{{{points * points * points - points * points,
                 points * points * points - points,
                 0,
                 points * points - points}},
               {{0, points * points - points, 0, points * points - points}}}}}},
           {{{{{{points * points * points - points * points,
                 points * points * points - points * points + points - 1,
                 points * points - points,
                 points * points - 1}},
               {{0,
                 points - 1,
                 points * points - points,
                 points * points - 1}}}},
             {{{{points * points * points - points * points,
                 points * points * points - points * points + points - 1,
                 0,
                 points - 1}},
               {{0, points - 1, 0, points - 1}}}}}}}};

      return lines_plane_array[i][j][k][l];
    }

    inline DEAL_II_ALWAYS_INLINE_RELEASE unsigned int
    lines(unsigned int i, unsigned int j, unsigned int k, unsigned int l) const
    {
      static constexpr unsigned int points = fe_degree + 1;

      static constexpr dealii::ndarray<unsigned int, 3, 2, 2, 3> lines_array = {
        {{{{{{{points * points - points,
               points * points * points - points * points,
               points * points * points - points}},
             {{0,
               points * points - points,
               points * points * points - points}}}},
           {{{{0,
               points * points * points - points * points,
               points * points * points - points}},
             {{0,
               points * points - points,
               points * points * points - points * points}}}}}},
         {{{{{{points - 1,
               points * points * points - points * points,
               points * points * points - points * points + points - 1}},
             {{0,
               points - 1,
               points * points * points - points * points + points - 1}}}},
           {{{{0,
               points * points * points - points * points,
               points * points * points - points * points + points - 1}},
             {{0, points - 1, points * points * points - points * points}}}}}},
         {{{{{{points - 1, points * points - points, points * points - 1}},
             {{0, points - 1, points * points - 1}}}},
           {{{{0, points * points - points, points * points - 1}},
             {{0, points - 1, points * points - points}}}}}}}};

      return lines_array[i][j][k][l];
    }
  };


  template <int dim, int fe_degree, typename Number>
  class FEEvaluationImplHangingNodesRunner<
    FEEvaluationImplHangingNodesRunnerTypes::scalar,
    dim,
    fe_degree,
    Number>
  {
  public:
    static const VectorizationTypes VectorizationType =
      VectorizationTypes::index;

  private:
    template <unsigned int side, bool transpose>
    static inline DEAL_II_ALWAYS_INLINE_RELEASE void
    interpolate_2D(
      const unsigned int                                          given_degree,
      const typename Trait<Number, VectorizationType>::index_type v,
      const typename Trait<Number, VectorizationType>::interpolation_type
        *DEAL_II_RESTRICT      weight,
      Number *DEAL_II_RESTRICT values)
    {
      static constexpr unsigned int max_n_points_1D = 40;

      typename Trait<Number, VectorizationType>::value_type
        temp[fe_degree != -1 ? (fe_degree + 1) : max_n_points_1D];

      const unsigned int points =
        (fe_degree != -1 ? fe_degree : given_degree) + 1;

      AssertIndexRange(given_degree, max_n_points_1D);

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
          temp[k] = Trait<Number, VectorizationType>::get_value(
            values[i * offset + stride + j], v);

      // perform interpolation point by point (note: r1 * r2 ==
      // points^(dim-1))
      for (unsigned int i = 0, k = 0; i < r1; ++i)
        for (unsigned int j = 0; j < r2; ++j, ++k)
          {
            typename Trait<Number, VectorizationType>::value_type sum = 0.0;
            for (unsigned int h = 0; h < points; ++h)
              sum += Trait<Number, VectorizationType>::get_value(
                       weight[(transpose ? 1 : points) * k +
                              (transpose ? points : 1) * h],
                       v) *
                     temp[h];
            Trait<Number, VectorizationType>::set_value(
              values[i * offset + stride + j], sum, v);
          }
    }

  public:
    template <bool transpose, typename Number2>
    static void
    run_internal(
      const unsigned int                             n_desired_components,
      const MatrixFreeFunctions::ShapeInfo<Number2> &shape_info,
      const std::array<MatrixFreeFunctions::compressed_constraint_kind,
                       Number::size()>              &constraint_mask,
      Number                                        *values)
    {
      const unsigned int given_degree =
        fe_degree != -1 ? fe_degree : shape_info.data.front().fe_degree;

      const auto &interpolation_matrices =
        Trait<Number, VectorizationType>::get_interpolation_matrix(shape_info);

      const auto constraint_mask_sorted =
        Trait<Number, VectorizationType>::create_mask(constraint_mask);

      for (unsigned int c = 0; c < n_desired_components; ++c)
        {
          for (unsigned int v = 0; v < Number::size(); ++v)
            {
              const auto mask = constraint_mask_sorted[v];

              if (Trait<Number, VectorizationType>::do_break(v, mask))
                break;

              if (Trait<Number, VectorizationType>::do_continue(v, mask))
                continue;

              const auto vv =
                Trait<Number, VectorizationType>::create(constraint_mask,
                                                         constraint_mask_sorted,
                                                         v);

              if (dim == 2) // 2d: only faces
                {
                  const bool subcell_x = (mask >> 0) & 1;
                  const bool subcell_y = (mask >> 1) & 1;
                  const bool face_x    = (mask >> 3) & 1;
                  const bool face_y    = (mask >> 4) & 1;

                  // direction 0:
                  if (face_y)
                    {
                      const auto *weights =
                        interpolation_matrices[!subcell_x].data();

                      if (subcell_y)
                        interpolate_2D<2, transpose>(given_degree,
                                                     vv,
                                                     weights,
                                                     values); // face 2
                      else
                        interpolate_2D<3, transpose>(given_degree,
                                                     vv,
                                                     weights,
                                                     values); // face 3
                    }

                  // direction 1:
                  if (face_x)
                    {
                      const auto *weights =
                        interpolation_matrices[!subcell_y].data();

                      if (subcell_x)
                        interpolate_2D<0, transpose>(given_degree,
                                                     vv,
                                                     weights,
                                                     values); // face 0
                      else
                        interpolate_2D<1, transpose>(given_degree,
                                                     vv,
                                                     weights,
                                                     values); // face 1
                    }
                }
              else if (dim == 3) // 3d faces and edges
                {
                  const bool type_x = (mask >> 0) & 1;
                  const bool type_y = (mask >> 1) & 1;
                  const bool type_z = (mask >> 2) & 1;

                  const auto flag_0 = (mask >> 3) & 3;
                  const auto flag_1 = (mask >> 5) & 7;
                  const auto faces  = (flag_0 & 0b01) ? flag_1 : 0;
                  const auto edges  = (flag_0 & 0b10) ? flag_1 : 0;

                  Helper<fe_degree == -1 ? HelperType::dynamic :
                                           HelperType::constant,
                         Number,
                         VectorizationType,
                         fe_degree,
                         transpose>
                    helper(given_degree,
                           type_x,
                           type_y,
                           type_z,
                           vv,
                           interpolation_matrices,
                           values);

                  if (faces > 0)
                    switch (faces)
                      {
                        case 0:
                          break;
                        case 1:
                          helper
                            .template process_faces_fast<true, false, false>();
                          break;
                        case 2:
                          helper
                            .template process_faces_fast<false, true, false>();
                          break;
                        case 3:
                          helper.template process_faces<true, true, false>();
                          break;
                        case 4:
                          helper
                            .template process_faces_fast<false, false, true>();
                          break;
                        case 5:
                          helper.template process_faces<true, false, true>();
                          break;
                        case 6:
                          helper.template process_faces<false, true, true>();
                          break;
                        case 7:
                          helper.template process_faces<true, true, true>();
                          break;
                      }

                  if (edges > 0)
                    switch (edges)
                      {
                        case 0:
                          break;
                        case 1:
                          helper.template process_edge<true, false, false>();
                          break;
                        case 2:
                          helper.template process_edge<false, true, false>();
                          break;
                        case 3:
                          helper.template process_edge<true, true, false>();
                          break;
                        case 4:
                          helper.template process_edge<false, false, true>();
                          break;
                        case 5:
                          helper.template process_edge<true, false, true>();
                          break;
                        case 6:
                          helper.template process_edge<false, true, true>();
                          break;
                        case 7:
                          helper.template process_edge<true, true, true>();
                          break;
                      }
                }
              else
                {
                  DEAL_II_NOT_IMPLEMENTED();
                }
            }

          values += shape_info.dofs_per_component_on_cell;
        }
    }
  };



  template <int dim, typename Number>
  struct FEEvaluationImplHangingNodes
  {
  public:
    template <int fe_degree, typename Number2>
    static bool
    run(const unsigned int                             n_desired_components,
        const MatrixFreeFunctions::ShapeInfo<Number2> &shape_info,
        const bool                                     transpose,
        const std::array<MatrixFreeFunctions::compressed_constraint_kind,
                         Number::size()>              &c_mask,
        Number                                        *values)
    {
      using RunnerType =
        FEEvaluationImplHangingNodesRunner<used_runner_type<fe_degree>(),
                                           dim,
                                           fe_degree,
                                           Number>;

      if (transpose)
        RunnerType::template run_internal<true>(n_desired_components,
                                                shape_info,
                                                c_mask,
                                                values);
      else
        RunnerType::template run_internal<false>(n_desired_components,
                                                 shape_info,
                                                 c_mask,
                                                 values);

      return false;
    }

    template <int fe_degree>
    static constexpr FEEvaluationImplHangingNodesRunnerTypes
    used_runner_type()
    {
      return ((Number::size() > 2) && (fe_degree == -1 || fe_degree > 2)) ?
               FEEvaluationImplHangingNodesRunnerTypes::vectorized :
               FEEvaluationImplHangingNodesRunnerTypes::scalar;
    }
  };


} // end of namespace internal

#undef DEAL_II_ALWAYS_INLINE_RELEASE


DEAL_II_NAMESPACE_CLOSE

#endif
