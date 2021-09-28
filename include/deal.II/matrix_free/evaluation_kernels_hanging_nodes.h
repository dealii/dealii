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
    // using VectorizationType = unsigned int;
    using VectorizationType = Number;

    static const unsigned int max_n_points_1D = 40;

    template <typename T1, typename T2>
    struct Trait;

    template <typename T1>
    struct Trait<T1, unsigned int>
    {
      using value_type = typename T1::value_type;

      static inline DEAL_II_ALWAYS_INLINE unsigned int
      create(const unsigned int v)
      {
        return v;
      }
    };

    template <typename T1>
    struct Trait<T1, T1>
    {
      using value_type = T1;

      static inline DEAL_II_ALWAYS_INLINE Number
      create(const unsigned int v)
      {
        Number result = 0.0;
        result[v]     = 1.0;
        return result;
      }
    };

    static inline DEAL_II_ALWAYS_INLINE typename Number::value_type
    get_value(const Number &value, const unsigned int &i)
    {
      return value[i];
    }

    static inline DEAL_II_ALWAYS_INLINE Number
    get_value(const Number &value, const Number &)
    {
      return value;
    }

    static inline DEAL_II_ALWAYS_INLINE void
    set_value(Number &                           result,
              const typename Number::value_type &value,
              const unsigned int &               i)
    {
      result[i] = value;
    }

    static inline DEAL_II_ALWAYS_INLINE void
    set_value(Number &result, const Number &value, const Number &i)
    {
      result = result * (Number(1.0) - i) + value * i;
    }

    template <int fe_degree, unsigned int side, bool transpose>
    static inline DEAL_II_ALWAYS_INLINE void
    interpolate_2D(const unsigned int      given_degree,
                   const VectorizationType v,
                   const Number *          weight,
                   Number *                values)
    {
      typename Trait<Number, VectorizationType>::value_type
        temp[max_n_points_1D];

      const unsigned int points =
        (fe_degree != -1 ? fe_degree : given_degree) + 1;

      AssertIndexRange(points, max_n_points_1D);

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
          temp[k] = get_value(values[i * offset + stride + j], v);

      // perform interpolation point by point (note: r1 * r2 == points^(dim-1))
      for (unsigned int i = 0, k = 0; i < r1; ++i)
        for (unsigned int j = 0; j < r2; ++j, ++k)
          {
            typename Trait<Number, VectorizationType>::value_type sum = 0.0;
            for (unsigned int h = 0; h < points; ++h)
              sum += get_value(weight[(transpose ? 1 : points) * k +
                                      (transpose ? points : 1) * h],
                               v) *
                     temp[h];
            set_value(values[i * offset + stride + j], sum, v);
          }
    }

    template <int          fe_degree,
              unsigned int direction,
              unsigned int d,
              bool         transpose,
              bool         skip_borders>
    static inline DEAL_II_ALWAYS_INLINE void
    interpolate_3D_face(const unsigned int      dof_offset,
                        const unsigned int      given_degree,
                        const VectorizationType v,
                        const Number *          weight,
                        Number *                values)
    {
      typename Trait<Number, VectorizationType>::value_type
        temp[max_n_points_1D];

      const unsigned int points =
        (fe_degree != -1 ? fe_degree : given_degree) + 1;

      AssertIndexRange(points, max_n_points_1D);

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

      for (unsigned int g = (skip_borders ? 1 : 0);
           g < points - (skip_borders ? 1 : 0);
           ++g)
        {
          // copy result back
          for (unsigned int k = 0; k < points; ++k)
            temp[k] =
              get_value(values[dof_offset + k * stride + stride2 * g], v);

          // perform interpolation point by point
          for (unsigned int k = 0; k < points; ++k)
            {
              typename Trait<Number, VectorizationType>::value_type sum = 0.0;
              for (unsigned int h = 0; h < points; ++h)
                sum += get_value(weight[(transpose ? 1 : points) * k +
                                        (transpose ? points : 1) * h],
                                 v) *
                       temp[h];
              set_value(values[dof_offset + k * stride + stride2 * g], sum, v);
            }
        }
    }

    template <int fe_degree, unsigned int direction, bool transpose>
    static inline DEAL_II_ALWAYS_INLINE void
    interpolate_3D_edge(const unsigned int      p,
                        const unsigned int      given_degree,
                        const VectorizationType v,
                        const Number *          weight,
                        Number *                values)
    {
      typename Trait<Number, VectorizationType>::value_type
        temp[max_n_points_1D];

      const unsigned int points =
        (fe_degree != -1 ? fe_degree : given_degree) + 1;

      AssertIndexRange(points, max_n_points_1D);

      const unsigned int stride = Utilities::pow(points, direction);

      // copy result back
      for (unsigned int k = 0; k < points; ++k)
        temp[k] = get_value(values[p + k * stride], v);

      // perform interpolation point by point
      for (unsigned int k = 0; k < points; ++k)
        {
          typename Trait<Number, VectorizationType>::value_type sum = 0.0;
          for (unsigned int h = 0; h < points; ++h)
            sum += get_value(weight[(transpose ? 1 : points) * k +
                                    (transpose ? points : 1) * h],
                             v) *
                   temp[h];
          set_value(values[p + k * stride], sum, v);
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

      const unsigned int points = given_degree + 1;

      const unsigned int p0 = 0;
      const unsigned int p1 = points - 1;
      const unsigned int p2 = points * points - points;
      const unsigned int p3 = points * points - 1;
      const unsigned int p4 = points * points * points - points * points;
      const unsigned int p5 =
        points * points * points - points * points + points - 1;
      const unsigned int p6 = points * points * points - points;

      static const std::array<unsigned int, 12> line_to_point = {
        {p0, p1, p0, p2, p4, p5, p4, p6, p0, p1, p2, p3}};

      static const std::array<unsigned int, 6> face_to_point{
        {p0, p1, p0, p2, p0, p4}};

      for (unsigned int c = 0; c < n_desired_components; ++c)
        {
          for (unsigned int v = 0; v < Number::size(); ++v)
            {
              const auto mask = constraint_mask[v];

              if (mask == MatrixFreeFunctions::ConstraintKinds::unconstrained)
                continue;

              const auto vv = Trait<Number, VectorizationType>::create(v);

              if (dim == 2) // 2D: only faces
                {
                  const auto is_set = [](const auto a, const auto b) -> bool {
                    return (a & b) == b;
                  };

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
                          vv,
                          weights,
                          values); // face 2
                      else
                        interpolate_2D<fe_degree, 3, transpose>(
                          given_degree,
                          vv,
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
                          vv,
                          weights,
                          values); // face 0
                      else
                        interpolate_2D<fe_degree, 1, transpose>(
                          given_degree,
                          vv,
                          weights,
                          values); // face 1
                    }
                }
              else if (dim == 3) // 3D faces and edges
                {
                  const auto m = static_cast<std::uint16_t>(mask);

                  const bool type_x = (m >> 0) & 1;
                  const bool type_y = (m >> 1) & 1;
                  const bool type_z = (m >> 2) & 1;
                  const auto faces  = (m >> 3) & 7;
                  const auto edges  = (m >> 6);

                  Helper<fe_degree, transpose> helper(given_degree,
                                                      type_x,
                                                      type_y,
                                                      type_z,
                                                      vv,
                                                      interpolation_matrices,
                                                      line_to_point,
                                                      face_to_point,
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
                          helper.template process_edge<false, false, true>();
                          break;
                        case 2:
                          helper.template process_edge<true, false, false>();
                          break;
                        case 3:
                          helper.template process_edge<true, false, true>();
                          break;
                        case 4:
                          helper.template process_edge<false, true, false>();
                          break;
                        case 5:
                          helper.template process_edge<false, true, true>();
                          break;
                        case 6:
                          helper.template process_edge<true, true, false>();
                          break;
                        case 7:
                          helper.template process_edge<true, true, true>();
                          break;
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

    template <int fe_degree, bool transpose>
    class Helper
    {
    public:
      inline DEAL_II_ALWAYS_INLINE
      Helper(const unsigned int &                        given_degree,
             const bool &                                type_x,
             const bool &                                type_y,
             const bool &                                type_z,
             const VectorizationType &                   v,
             const std::array<AlignedVector<Number>, 2> &interpolation_matrices,
             const std::array<unsigned int, 12> &        line_to_point,
             const std::array<unsigned int, 6> &         face_to_point,
             Number *                                    values)
        : given_degree(given_degree)
        , type_x(type_x)
        , type_y(type_y)
        , type_z(type_z)
        , v(v)
        , interpolation_matrices(interpolation_matrices)
        , line_to_point(line_to_point)
        , face_to_point(face_to_point)
        , values(values)
      {}

      const unsigned int &                        given_degree;
      const bool &                                type_x;
      const bool &                                type_y;
      const bool &                                type_z;
      const VectorizationType &                   v;
      const std::array<AlignedVector<Number>, 2> &interpolation_matrices;
      const std::array<unsigned int, 12> &        line_to_point;
      const std::array<unsigned int, 6> &         face_to_point;
      Number *                                    values;

      template <bool do_x, bool do_y, bool do_z>
      inline DEAL_II_ALWAYS_INLINE void
      process_edge() const
      {
        if (do_x)
          interpolate_3D_edge<fe_degree, 0, transpose>(
            line_to_point[line[0][type_y][type_z]],
            given_degree,
            v,
            interpolation_matrices[!type_x].data(),
            values);

        if (do_y)
          interpolate_3D_edge<fe_degree, 1, transpose>(
            line_to_point[line[1][type_x][type_z]],
            given_degree,
            v,
            interpolation_matrices[!type_y].data(),
            values);

        if (do_z)
          interpolate_3D_edge<fe_degree, 2, transpose>(
            line_to_point[line[2][type_x][type_y]],
            given_degree,
            v,
            interpolation_matrices[!type_z].data(),
            values);
      }

      template <bool do_x, bool do_y, bool do_z>
      inline DEAL_II_ALWAYS_INLINE void
      process_faces_fast() const
      {
        static_assert((do_x && !do_y && !do_z) || (!do_x && do_y && !do_z) ||
                        (!do_x && !do_y && do_z),
                      "Only one face can be chosen.");

        static const unsigned int direction = do_x ? 0 : (do_y ? 1 : 2);
        const bool type = do_x ? type_x : (do_y ? type_y : type_z);

        if (!do_x)
          interpolate_3D_face<fe_degree, 0, direction, transpose, false>(
            face_to_point[face[direction][type]],
            given_degree,
            v,
            interpolation_matrices[!type_x].data(),
            values);

        if (!do_y)
          interpolate_3D_face<fe_degree, 1, direction, transpose, false>(
            face_to_point[face[direction][type]],
            given_degree,
            v,
            interpolation_matrices[!type_y].data(),
            values);

        if (!do_z)
          interpolate_3D_face<fe_degree, 2, direction, transpose, false>(
            face_to_point[face[direction][type]],
            given_degree,
            v,
            interpolation_matrices[!type_z].data(),
            values);
      }

      template <bool do_x, bool do_y, bool do_z>
      inline DEAL_II_ALWAYS_INLINE void
      process_faces() const
      {
        static_assert(((do_x && !do_y && !do_z) || (!do_x && do_y && !do_z) ||
                       (!do_x && !do_y && do_z)) == false,
                      "Only one face can be chosen.");

        // clang-format off
          
        // direction 0 -> faces
        if(do_y && given_degree > 1)
            interpolate_3D_face<fe_degree, 0, 1, transpose, true>(
              face_to_point[face[1][type_y]], given_degree, v, interpolation_matrices[!type_x].data(), values);

        if(do_z && given_degree > 1)
            interpolate_3D_face<fe_degree, 0, 2, transpose, true>(
              face_to_point[face[2][type_z]], given_degree, v, interpolation_matrices[!type_x].data(), values);

        // direction 0 -> edges
        interpolate_3D_edge<fe_degree, 0, transpose>(
          (do_x && do_y && !do_z) ? line_to_point[lines_plane[0][type_x][type_y][0]] : 
         ((do_x && !do_y && do_z) ? line_to_point[lines_plane[1][type_x][type_z][0]] : 
                                    line_to_point[lines[0][type_y][type_z][0]]) , given_degree, v, interpolation_matrices[!type_x].data(), values);


        interpolate_3D_edge<fe_degree, 0, transpose>(
          (do_x && do_y && !do_z) ? line_to_point[lines_plane[0][type_x][type_y][1]] :
         ((do_x && !do_y && do_z) ? line_to_point[lines_plane[1][type_x][type_z][1]] : 
                                    line_to_point[lines[0][type_y][type_z][1]]), given_degree, v, interpolation_matrices[!type_x].data(), values);

        if(do_y && do_z)
            interpolate_3D_edge<fe_degree, 0, transpose>(
              line_to_point[lines[0][type_y][type_z][2]], given_degree, v, interpolation_matrices[!type_x].data(), values);

        // direction 1 -> faces
        if(do_x && given_degree > 1)
            interpolate_3D_face<fe_degree, 1, 0, transpose, true>(
              face_to_point[face[0][type_x]], given_degree, v, interpolation_matrices[!type_y].data(), values);

        if(do_z && given_degree > 1)
            interpolate_3D_face<fe_degree, 1, 2, transpose, true>(
              face_to_point[face[2][type_z]], given_degree, v, interpolation_matrices[!type_y].data(), values);

        // direction 1 -> lines
        interpolate_3D_edge<fe_degree, 1, transpose>(
          (do_x && do_y && !do_z) ? line_to_point[lines_plane[0][type_x][type_y][2]] : 
         ((!do_x && do_y && do_z) ? line_to_point[lines_plane[2][type_y][type_z][0]] :  
                                    line_to_point[lines[1][type_x][type_z][0]]),
                given_degree, v, interpolation_matrices[!type_y].data(), values);

        interpolate_3D_edge<fe_degree, 1, transpose>(
          (do_x && do_y && !do_z) ? line_to_point[lines_plane[0][type_x][type_y][3]] : 
         ((!do_x && do_y && do_z) ? line_to_point[lines_plane[2][type_y][type_z][1]] :
                                    line_to_point[lines[1][type_x][type_z][1]]),
                given_degree, v, interpolation_matrices[!type_y].data(), values);

        if(do_x && do_z)
            interpolate_3D_edge<fe_degree, 1, transpose>(
              line_to_point[lines[1][type_x][type_z][2]], given_degree, v, interpolation_matrices[!type_y].data(), values);

        // direction 2 -> faces
        if(do_x && given_degree > 1)
            interpolate_3D_face<fe_degree, 2, 0, transpose, true>(
              face_to_point[face[0][type_x]], given_degree, v, interpolation_matrices[!type_z].data(), values);

        if(do_y && given_degree > 1)
            interpolate_3D_face<fe_degree, 2, 1, transpose, true>(
              face_to_point[face[1][type_y]], given_degree, v, interpolation_matrices[!type_z].data(), values);

        // direction 2 -> edges
        interpolate_3D_edge<fe_degree, 2, transpose>(
          (do_x && !do_y && do_z) ? line_to_point[lines_plane[1][type_x][type_z][2]] : 
         ((!do_x && do_y && do_z) ? line_to_point[lines_plane[2][type_y][type_z][2]] :
                                    line_to_point[lines[2][type_x][type_y][0]]),
            given_degree, v, interpolation_matrices[!type_z].data(), values);

        interpolate_3D_edge<fe_degree, 2, transpose>(
          (do_x && !do_y && do_z) ? line_to_point[lines_plane[1][type_x][type_z][3]] :
         ((!do_x && do_y && do_z) ? line_to_point[lines_plane[2][type_y][type_z][3]] : 
                                    line_to_point[lines[2][type_x][type_y][1]]), 
            given_degree, v, interpolation_matrices[!type_z].data(), values);

        if(do_x && do_y)
            interpolate_3D_edge<fe_degree, 2, transpose>(
              line_to_point[lines[2][type_x][type_y][2]], given_degree, v, interpolation_matrices[!type_z].data(), values);

        // clang-format on
      }

    private:
      static constexpr std::array<std::array<std::array<unsigned int, 2>, 2>, 3>
        line = {{{{{{7, 3}}, {{6, 2}}}},
                 {{{{5, 1}}, {{4, 0}}}},
                 {{{{11, 9}}, {{10, 8}}}}}};



      static constexpr std::array<std::array<unsigned int, 2>, 3> face = {
        {{{1, 0}}, {{3, 2}}, {{5, 4}}}};

      static constexpr std::
        array<std::array<std::array<std::array<unsigned int, 4>, 2>, 2>, 3>
          lines_plane{{{{{{{{3, 7, 1, 5}}, {{2, 6, 1, 5}}}},
                         {{{{3, 7, 0, 4}}, {{2, 6, 0, 4}}}}}},
                       {{{{{{6, 7, 9, 11}}, {{2, 3, 9, 11}}}},
                         {{{{6, 7, 8, 10}}, {{2, 3, 8, 10}}}}}},
                       {{{{{{4, 5, 10, 11}}, {{0, 1, 10, 11}}}},
                         {{{{4, 5, 8, 9}}, {{0, 1, 8, 9}}}}}}}};

      static constexpr std::
        array<std::array<std::array<std::array<unsigned int, 3>, 2>, 2>, 3>
          lines{{{{{{{{3, 6, 7}}, {{2, 3, 7}}}}, {{{{2, 6, 7}}, {{2, 3, 6}}}}}},
                 {{{{{{1, 4, 5}}, {{0, 1, 5}}}}, {{{{0, 4, 5}}, {{0, 1, 4}}}}}},
                 {{{{{{9, 10, 11}}, {{8, 9, 11}}}},
                   {{{{8, 10, 11}}, {{8, 9, 10}}}}}}}};
    };
  };


} // end of namespace internal


DEAL_II_NAMESPACE_CLOSE

#endif
