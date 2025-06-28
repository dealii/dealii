// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/derivative_form.h>
#include <deal.II/base/geometry_info.h>
#include <deal.II/base/qprojector.h>

#include <deal.II/grid/reference_cell.h>
#include <deal.II/grid/tria_orientation.h>

#include <boost/container/small_vector.hpp>


DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace QProjector
  {
    namespace
    {
      // Reflect points across the y = x line.
      std::vector<Point<2>>
      reflect(const std::vector<Point<2>> &points)
      {
        std::vector<Point<2>> q_points;
        q_points.reserve(points.size());
        for (const Point<2> &p : points)
          q_points.emplace_back(p[1], p[0]);

        return q_points;
      }


      // Rotate points in the plane around the positive z axis (i.e.,
      // counter-clockwise).
      std::vector<Point<2>>
      rotate(const std::vector<Point<2>> &points, const unsigned int n_times)
      {
        std::vector<Point<2>> q_points;
        q_points.reserve(points.size());
        switch (n_times % 4)
          {
            case 0:
              // 0 degree. the point remains as it is.
              for (const Point<2> &p : points)
                q_points.push_back(p);
              break;
            case 1:
              // 90 degree counterclockwise
              for (const Point<2> &p : points)
                q_points.emplace_back(1.0 - p[1], p[0]);
              break;
            case 2:
              // 180 degree counterclockwise
              for (const Point<2> &p : points)
                q_points.emplace_back(1.0 - p[0], 1.0 - p[1]);
              break;
            case 3:
              // 270 degree counterclockwise
              for (const Point<2> &p : points)
                q_points.emplace_back(p[1], 1.0 - p[0]);
              break;
          }

        return q_points;
      }

      /**
       * Internal function to translate a 2-dimensional quadrature formula to
       * a 3-dimensional quadrature formula on hex elements, addressing both
       * standard faces (FEFaceValues) and subfaces (FESubfaceValues),
       * depending on the given refinement case and subface.
       */
      void
      project_to_hex_face_and_append(
        const std::vector<Point<2>> &points,
        const unsigned int           face_no,
        std::vector<Point<3>>       &q_points,
        const RefinementCase<2> &ref_case   = RefinementCase<2>::no_refinement,
        const unsigned int       subface_no = 0)
      {
        // one coordinate is at a const value. for faces 0, 2 and 4 this value
        // is 0.0, for faces 1, 3 and 5 it is 1.0
        const double const_value = face_no % 2;

        // local 2d coordinates are xi and eta, global 3d coordinates are x, y
        // and z. those have to be mapped. the following indices tell, which
        // global coordinate (0->x, 1->y, 2->z) corresponds to which local one
        const unsigned int xi_index    = (1 + face_no / 2) % 3,
                           eta_index   = (2 + face_no / 2) % 3,
                           const_index = face_no / 2;

        // for a standard face (no refinement), we use the default values of
        // the xi and eta scales and translations, otherwise the xi and eta
        // values will be scaled (by factor 0.5 or factor 1.0) depending on
        // the refinement case and translated (by 0.0 or 0.5) depending on the
        // refinement case and subface_no
        double xi_scale = 1.0, eta_scale = 1.0, xi_translation = 0.0,
               eta_translation = 0.0;

        // set the scale and translation parameter for individual subfaces
        switch (ref_case)
          {
            case RefinementCase<2>::no_refinement:
              break;
            case RefinementCase<2>::cut_x:
              xi_scale       = 0.5;
              xi_translation = subface_no % 2 * 0.5;
              break;
            case RefinementCase<2>::cut_y:
              eta_scale       = 0.5;
              eta_translation = subface_no % 2 * 0.5;
              break;
            case RefinementCase<2>::cut_xy:
              xi_scale        = 0.5;
              eta_scale       = 0.5;
              xi_translation  = int(subface_no % 2) * 0.5;
              eta_translation = int(subface_no / 2) * 0.5;
              break;
            default:
              DEAL_II_ASSERT_UNREACHABLE();
              break;
          }

        // finally, compute the scaled, translated, projected quadrature
        // points
        for (const Point<2> &p : points)
          {
            Point<3> cell_point;
            cell_point[xi_index]    = xi_scale * p[0] + xi_translation;
            cell_point[eta_index]   = eta_scale * p[1] + eta_translation;
            cell_point[const_index] = const_value;
            q_points.push_back(cell_point);
          }
      }

      std::vector<Point<2>>
      mutate_points_with_offset(
        const std::vector<Point<2>>       &points,
        const types::geometric_orientation combined_orientation)
      {
        // These rotations are backwards (relative to the standard notion of,
        // e.g., what rotation index 7 means) since they are rotations about the
        // positive z axis in 2d: i.e., they are done from the perspective of
        // 'inside' a cell instead of the perspective of an abutting cell.
        //
        // For example: consider points on face 4 of a hexahedron with
        // orientation 3. In 2d, rotating such points clockwise is the same as
        // rotating them counter-clockwise from the perspective of the abutting
        // face. Hence, such points must be rotated 90 degrees
        // counter-clockwise.
        switch (combined_orientation)
          {
            case 1:
              return reflect(points);
            case 3:
              return rotate(reflect(points), 3);
            case 5:
              return rotate(reflect(points), 2);
            case 7:
              return rotate(reflect(points), 1);
            case 0:
              return points;
            case 2:
              return rotate(points, 1);
            case 4:
              return rotate(points, 2);
            case 6:
              return rotate(points, 3);
            default:
              DEAL_II_ASSERT_UNREACHABLE();
          }
        return {};
      }



      /**
       * Append the points and weights of a quadrature rule projected onto a
       * subobject (either face or subface) to the end of two vectors. These
       * vectors should ultimately be indexed with a
       * QProjector::DataSetDescriptor.
       *
       * The goal of this function (as used by QProjector and then FEFaceValues)
       * is to compute identical sets of quadrature points on the common face of
       * two abutting cells. Our orientation convention is that, given such a
       * pair of abutting cells:
       *
       * 1. The shared face, from the perspective of the first cell, is
       *    in the default orientation.
       * 2. The shared face, from the perspective of the second cell, has
       *    its orientation computed relative to the first cell: i.e.,
       *    'orientation' is the vertex permutation applied to the first
       *    cell's face to get the second cell's face.
       *
       * The first case is trivial since points do not need to be
       * oriented. However, in the second case, we need to use the
       * *reverse* of the stored orientation (i.e., the permutation
       * applied to the second cell's face which yields the first cell's
       * face) so that we get identical quadrature points.
       *
       * For more information see connectivity.h.
       */
      template <int dim>
      void
      append_subobject_rule(
        const ReferenceCell               &face_reference_cell,
        const Quadrature<dim - 1>         &quadrature,
        const std::vector<Point<dim>>     &vertices,
        const double                       measure,
        const types::geometric_orientation combined_orientation,
        std::vector<Point<dim>>           &points,
        std::vector<double>               &weights)
      {
        const auto support_points =
          face_reference_cell.permute_by_combined_orientation(
            make_const_array_view(vertices),
            face_reference_cell.get_inverse_combined_orientation(
              combined_orientation));

        for (unsigned int j = 0; j < quadrature.size(); ++j)
          {
            Point<dim> mapped_point;

            // map reference quadrature point
            for (const unsigned int vertex_no :
                 face_reference_cell.vertex_indices())
              mapped_point +=
                support_points[vertex_no] *
                face_reference_cell.d_linear_shape_function(quadrature.point(j),
                                                            vertex_no);

            points.push_back(mapped_point);

            // rescale quadrature weights so that the sum of the weights on
            // each face equals the measure of that face.
            weights.push_back(quadrature.weight(j) * measure /
                              face_reference_cell.volume());
          }
      }
    } // namespace
  }   // namespace QProjector
} // namespace internal



template <>
void
QProjector<1>::project_to_face(const ReferenceCell   &reference_cell,
                               const Quadrature<0>   &quadrature,
                               const unsigned int     face_no,
                               std::vector<Point<1>> &q_points)
{
  AssertDimension(quadrature.size(), q_points.size());
  const auto face_quadrature =
    QProjector<1>::project_to_face(reference_cell,
                                   quadrature,
                                   face_no,
                                   numbers::default_geometric_orientation);
  q_points = face_quadrature.get_points();
}



template <>
void
QProjector<2>::project_to_face(const ReferenceCell   &reference_cell,
                               const Quadrature<1>   &quadrature,
                               const unsigned int     face_no,
                               std::vector<Point<2>> &q_points)
{
  AssertDimension(quadrature.size(), q_points.size());
  const auto face_quadrature =
    QProjector<2>::project_to_face(reference_cell,
                                   quadrature,
                                   face_no,
                                   numbers::default_geometric_orientation);
  q_points = face_quadrature.get_points();
}



template <>
void
QProjector<3>::project_to_face(const ReferenceCell   &reference_cell,
                               const Quadrature<2>   &quadrature,
                               const unsigned int     face_no,
                               std::vector<Point<3>> &q_points)
{
  Assert(reference_cell == ReferenceCells::Hexahedron, ExcNotImplemented());
  (void)reference_cell;

  AssertIndexRange(face_no, GeometryInfo<3>::faces_per_cell);
  Assert(q_points.size() == quadrature.size(),
         ExcDimensionMismatch(q_points.size(), quadrature.size()));
  q_points.clear();
  internal::QProjector::project_to_hex_face_and_append(quadrature.get_points(),
                                                       face_no,
                                                       q_points);
}



template <int dim>
Quadrature<dim>
QProjector<dim>::project_to_oriented_face(const ReferenceCell &reference_cell,
                                          const Quadrature<dim - 1> &quadrature,
                                          const unsigned int         face_no,
                                          const bool face_orientation,
                                          const bool face_flip,
                                          const bool face_rotation)
{
  return QProjector<dim>::project_to_face(
    reference_cell,
    quadrature,
    face_no,
    internal::combined_face_orientation(face_orientation,
                                        face_rotation,
                                        face_flip));
}



template <int dim>
Quadrature<dim>
QProjector<dim>::project_to_face(
  const ReferenceCell               &reference_cell,
  const Quadrature<dim - 1>         &quadrature,
  const unsigned int                 face_no,
  const types::geometric_orientation combined_orientation)
{
  AssertIndexRange(face_no, reference_cell.n_faces());
  AssertIndexRange(combined_orientation,
                   reference_cell.n_face_orientations(face_no));
  AssertDimension(reference_cell.get_dimension(), dim);

  std::vector<Point<dim>> points;
  std::vector<double>     weights;

  const ReferenceCell face_reference_cell =
    reference_cell.face_reference_cell(face_no);
  std::vector<Point<dim>> face_vertices(face_reference_cell.n_vertices());
  for (const unsigned int vertex_no : face_reference_cell.vertex_indices())
    face_vertices[vertex_no] =
      reference_cell.face_vertex_location<dim>(face_no, vertex_no);
  internal::QProjector::append_subobject_rule(face_reference_cell,
                                              quadrature,
                                              face_vertices,
                                              reference_cell.face_measure(
                                                face_no),
                                              combined_orientation,
                                              points,
                                              weights);

  return Quadrature<dim>(std::move(points), std::move(weights));
}



template <>
void
QProjector<1>::project_to_subface(const ReferenceCell &reference_cell,
                                  const Quadrature<0> &,
                                  const unsigned int face_no,
                                  const unsigned int,
                                  std::vector<Point<1>> &q_points,
                                  const RefinementCase<0> &)
{
  Assert(reference_cell == ReferenceCells::Line, ExcNotImplemented());
  (void)reference_cell;

  const unsigned int dim = 1;
  AssertIndexRange(face_no, GeometryInfo<dim>::faces_per_cell);
  AssertDimension(q_points.size(), 1);

  q_points[0] = Point<dim>(static_cast<double>(face_no));
}



template <>
void
QProjector<2>::project_to_subface(const ReferenceCell     &reference_cell,
                                  const Quadrature<1>     &quadrature,
                                  const unsigned int       face_no,
                                  const unsigned int       subface_no,
                                  std::vector<Point<2>>   &q_points,
                                  const RefinementCase<1> &ref_case)
{
  AssertDimension(quadrature.size(), q_points.size());
  const auto face_quadrature =
    project_to_subface(reference_cell,
                       quadrature,
                       face_no,
                       subface_no,
                       numbers::default_geometric_orientation,
                       ref_case);
  q_points = face_quadrature.get_points();
}



template <>
void
QProjector<3>::project_to_subface(const ReferenceCell     &reference_cell,
                                  const Quadrature<2>     &quadrature,
                                  const unsigned int       face_no,
                                  const unsigned int       subface_no,
                                  std::vector<Point<3>>   &q_points,
                                  const RefinementCase<2> &ref_case)
{
  Assert(reference_cell == ReferenceCells::Hexahedron, ExcNotImplemented());
  (void)reference_cell;
  AssertDimension(quadrature.size(), q_points.size());
  const auto face_quadrature =
    project_to_subface(reference_cell,
                       quadrature,
                       face_no,
                       subface_no,
                       numbers::default_geometric_orientation,
                       ref_case);
  q_points = face_quadrature.get_points();
}



template <int dim>
Quadrature<dim>
QProjector<dim>::project_to_oriented_subface(
  const ReferenceCell       &reference_cell,
  const Quadrature<dim - 1> &quadrature,
  const unsigned int         face_no,
  const unsigned int         subface_no,
  const bool,
  const bool,
  const bool,
  const internal::SubfaceCase<dim>)
{
  return QProjector<dim>::project_to_subface(
    reference_cell,
    quadrature,
    face_no,
    subface_no,
    RefinementCase<dim - 1>::isotropic_refinement);
}



template <int dim>
Quadrature<dim>
QProjector<dim>::project_to_subface(
  const ReferenceCell               &reference_cell,
  const SubQuadrature               &quadrature,
  const unsigned int                 face_no,
  const unsigned int                 subface_no,
  const types::geometric_orientation combined_orientation,
  const RefinementCase<dim - 1>     &ref_case)
{
  AssertIndexRange(face_no, reference_cell.n_faces());
  AssertIndexRange(combined_orientation,
                   reference_cell.n_face_orientations(face_no));
  AssertDimension(reference_cell.get_dimension(), dim);
  AssertIndexRange(subface_no,
                   reference_cell.face_reference_cell(face_no)
                     .template n_children<dim - 1>(ref_case));

  std::vector<Point<dim>> q_points;
  std::vector<double>     q_weights = quadrature.get_weights();
  q_points.reserve(quadrature.size());

  if constexpr (dim == 1)
    {
      AssertDimension(quadrature.size(), 1);
      q_points.emplace_back(static_cast<double>(face_no));
    }
  else if constexpr (dim == 2)
    {
      if (reference_cell == ReferenceCells::Triangle)
        // use linear polynomial to map the reference quadrature points
        // correctly on faces, i.e., BarycentricPolynomials<1>(1)
        for (unsigned int p = 0; p < quadrature.size(); ++p)
          {
            if (face_no == 0)
              {
                if (subface_no == 0)
                  q_points.emplace_back(quadrature.point(p)[0] / 2, 0);
                else
                  q_points.emplace_back(0.5 + quadrature.point(p)[0] / 2, 0);
              }
            else if (face_no == 1)
              {
                if (subface_no == 0)
                  q_points.emplace_back(1 - quadrature.point(p)[0] / 2,
                                        quadrature.point(p)[0] / 2);
                else
                  q_points.emplace_back(0.5 - quadrature.point(p)[0] / 2,
                                        0.5 + quadrature.point(p)[0] / 2);
              }
            else if (face_no == 2)
              {
                if (subface_no == 0)
                  q_points.emplace_back(0, 1 - quadrature.point(p)[0] / 2);
                else
                  q_points.emplace_back(0, 0.5 - quadrature.point(p)[0] / 2);
              }
            else
              DEAL_II_ASSERT_UNREACHABLE();
          }
      else if (reference_cell == ReferenceCells::Quadrilateral)
        for (unsigned int p = 0; p < quadrature.size(); ++p)
          {
            if (face_no == 0)
              {
                if (subface_no == 0)
                  q_points.emplace_back(0, quadrature.point(p)[0] / 2);
                else
                  q_points.emplace_back(0, quadrature.point(p)[0] / 2 + 0.5);
              }
            else if (face_no == 1)
              {
                if (subface_no == 0)
                  q_points.emplace_back(1, quadrature.point(p)[0] / 2);
                else
                  q_points.emplace_back(1, quadrature.point(p)[0] / 2 + 0.5);
              }
            else if (face_no == 2)
              {
                if (subface_no == 0)
                  q_points.emplace_back(quadrature.point(p)[0] / 2, 0);
                else
                  q_points.emplace_back(quadrature.point(p)[0] / 2 + 0.5, 0);
              }
            else if (face_no == 3)
              {
                if (subface_no == 0)
                  q_points.emplace_back(quadrature.point(p)[0] / 2, 1);
                else
                  q_points.emplace_back(quadrature.point(p)[0] / 2 + 0.5, 1);
              }
            else
              DEAL_II_ASSERT_UNREACHABLE();
          }
      else
        DEAL_II_ASSERT_UNREACHABLE();

      if (combined_orientation == numbers::reverse_line_orientation)
        {
          std::reverse(q_points.begin(), q_points.end());
          std::reverse(q_weights.begin(), q_weights.end());
        }
      for (auto &w : q_weights)
        w *= reference_cell.face_measure(face_no);
    }
  else if constexpr (dim == 3)
    {
      Assert(reference_cell == ReferenceCells::Hexahedron, ExcNotImplemented());
      internal::QProjector::project_to_hex_face_and_append(
        quadrature.get_points(), face_no, q_points, ref_case, subface_no);
    }
  else
    {
      DEAL_II_ASSERT_UNREACHABLE();
    }

  return Quadrature<dim>(std::move(q_points), std::move(q_weights));
}



template <int dim>
Quadrature<dim>
QProjector<dim>::project_to_all_faces(
  const ReferenceCell            &reference_cell,
  const hp::QCollection<dim - 1> &quadrature)
{
  std::vector<Point<dim>> points;
  std::vector<double>     weights;

  for (const unsigned int face_no : reference_cell.face_indices())
    {
      const ReferenceCell face_reference_cell =
        reference_cell.face_reference_cell(face_no);
      std::vector<Point<dim>> face_vertices(face_reference_cell.n_vertices());
      for (const unsigned int vertex_no : face_reference_cell.vertex_indices())
        face_vertices[vertex_no] =
          reference_cell.face_vertex_location<dim>(face_no, vertex_no);

      for (types::geometric_orientation combined_orientation = 0;
           combined_orientation < reference_cell.n_face_orientations(face_no);
           ++combined_orientation)
        internal::QProjector::append_subobject_rule(
          face_reference_cell,
          quadrature[quadrature.size() == 1 ? 0 : face_no],
          face_vertices,
          reference_cell.face_measure(face_no),
          combined_orientation,
          points,
          weights);
    }

  return Quadrature<dim>(std::move(points), std::move(weights));
}



template <>
Quadrature<1>
QProjector<1>::project_to_all_subfaces(const ReferenceCell &reference_cell,
                                       const Quadrature<0> &quadrature)
{
  Assert(reference_cell == ReferenceCells::Line, ExcNotImplemented());
  (void)reference_cell;

  const unsigned int dim = 1;

  const unsigned int n_points = 1, n_faces = GeometryInfo<dim>::faces_per_cell,
                     subfaces_per_face =
                       GeometryInfo<dim>::max_children_per_face;

  // first fix quadrature points
  std::vector<Point<dim>> q_points;
  q_points.reserve(n_points * n_faces * subfaces_per_face);
  std::vector<Point<dim>> help(n_points);

  // project to each face and copy
  // results
  for (unsigned int face = 0; face < n_faces; ++face)
    for (unsigned int subface = 0; subface < subfaces_per_face; ++subface)
      {
        project_to_subface(reference_cell, quadrature, face, subface, help);
        std::copy(help.begin(), help.end(), std::back_inserter(q_points));
      }

  // next copy over weights
  std::vector<double> weights;
  weights.reserve(n_points * n_faces * subfaces_per_face);
  for (unsigned int face = 0; face < n_faces; ++face)
    for (unsigned int subface = 0; subface < subfaces_per_face; ++subface)
      std::copy(quadrature.get_weights().begin(),
                quadrature.get_weights().end(),
                std::back_inserter(weights));

  Assert(q_points.size() == n_points * n_faces * subfaces_per_face,
         ExcInternalError());
  Assert(weights.size() == n_points * n_faces * subfaces_per_face,
         ExcInternalError());

  return Quadrature<dim>(std::move(q_points), std::move(weights));
}



template <>
Quadrature<2>
QProjector<2>::project_to_all_subfaces(const ReferenceCell &reference_cell,
                                       const SubQuadrature &quadrature)
{
  Assert(reference_cell == ReferenceCells::Quadrilateral ||
           reference_cell == ReferenceCells::Triangle,
         ExcNotImplemented());

  const unsigned int dim = 2;

  std::vector<Point<dim>> q_points;
  std::vector<double>     weights;

  // project to each face and copy
  // results
  for (unsigned int face = 0; face < reference_cell.n_faces(); ++face)
    for (types::geometric_orientation orientation = 0;
         orientation < reference_cell.n_face_orientations(face);
         ++orientation)
      for (unsigned int subface = 0;
           subface <
           reference_cell.face_reference_cell(face).n_isotropic_children();
           ++subface)
        {
          const unsigned int local_subface =
            orientation == numbers::reverse_line_orientation ? 1 - subface :
                                                               subface;
          const auto sub_quadrature =
            project_to_subface(reference_cell,
                               quadrature,
                               face,
                               local_subface,
                               orientation,
                               RefinementCase<dim - 1>::isotropic_refinement);
          q_points.insert(q_points.end(),
                          sub_quadrature.get_points().begin(),
                          sub_quadrature.get_points().end());
          weights.insert(weights.end(),
                         sub_quadrature.get_weights().begin(),
                         sub_quadrature.get_weights().end());
        }

  return Quadrature<dim>(std::move(q_points), std::move(weights));
}



template <>
Quadrature<3>
QProjector<3>::project_to_all_subfaces(const ReferenceCell &reference_cell,
                                       const SubQuadrature &quadrature)
{
  if (reference_cell == ReferenceCells::Triangle ||
      reference_cell == ReferenceCells::Tetrahedron)
    return Quadrature<3>(); // nothing to do

  Assert(reference_cell == ReferenceCells::Hexahedron, ExcNotImplemented());

  const unsigned int dim = 3;

  const unsigned int n_points = quadrature.size(),
                     n_faces  = GeometryInfo<dim>::faces_per_cell,
                     total_subfaces_per_face = 2 + 2 + 4;

  // first fix quadrature points
  std::vector<Point<dim>> q_points;
  q_points.reserve(n_points * n_faces * total_subfaces_per_face * 8);

  std::vector<double> weights;
  weights.reserve(n_points * n_faces * total_subfaces_per_face * 8);

  // do the following for all possible mutations of a face (mutation==0
  // corresponds to a face with standard orientation, no flip and no rotation)
  for (unsigned char offset = 0; offset < 8; ++offset)
    {
      const auto mutation =
        internal::QProjector::mutate_points_with_offset(quadrature.get_points(),
                                                        offset);

      // project to each face and copy results
      for (unsigned int face = 0; face < n_faces; ++face)
        for (unsigned int ref_case = RefinementCase<dim - 1>::cut_xy;
             ref_case >= RefinementCase<dim - 1>::cut_x;
             --ref_case)
          for (unsigned int subface = 0;
               subface < GeometryInfo<dim - 1>::n_children(
                           RefinementCase<dim - 1>(ref_case));
               ++subface)
            {
              internal::QProjector::project_to_hex_face_and_append(
                mutation,
                face,
                q_points,
                RefinementCase<dim - 1>(ref_case),
                subface);

              // next copy over weights
              std::copy(quadrature.get_weights().begin(),
                        quadrature.get_weights().end(),
                        std::back_inserter(weights));
            }
    }

  Assert(q_points.size() == n_points * n_faces * total_subfaces_per_face * 8,
         ExcInternalError());
  Assert(weights.size() == n_points * n_faces * total_subfaces_per_face * 8,
         ExcInternalError());

  return Quadrature<dim>(std::move(q_points), std::move(weights));
}



template <int dim>
Quadrature<dim>
QProjector<dim>::project_to_child(const ReferenceCell   &reference_cell,
                                  const Quadrature<dim> &quadrature,
                                  const unsigned int     child_no)
{
  Assert(reference_cell == ReferenceCells::get_hypercube<dim>(),
         ExcNotImplemented());
  (void)reference_cell;

  AssertIndexRange(child_no, GeometryInfo<dim>::max_children_per_cell);

  const unsigned int n_q_points = quadrature.size();

  std::vector<Point<dim>> q_points(n_q_points);
  for (unsigned int i = 0; i < n_q_points; ++i)
    q_points[i] =
      GeometryInfo<dim>::child_to_cell_coordinates(quadrature.point(i),
                                                   child_no);

  // for the weights, things are
  // equally simple: copy them and
  // scale them
  std::vector<double> weights = quadrature.get_weights();
  for (unsigned int i = 0; i < n_q_points; ++i)
    weights[i] *= (1. / GeometryInfo<dim>::max_children_per_cell);

  return Quadrature<dim>(q_points, weights);
}



template <int dim>
Quadrature<dim>
QProjector<dim>::project_to_all_children(const ReferenceCell   &reference_cell,
                                         const Quadrature<dim> &quadrature)
{
  Assert(reference_cell == ReferenceCells::get_hypercube<dim>(),
         ExcNotImplemented());
  (void)reference_cell;

  const unsigned int n_points   = quadrature.size(),
                     n_children = GeometryInfo<dim>::max_children_per_cell;

  std::vector<Point<dim>> q_points(n_points * n_children);
  std::vector<double>     weights(n_points * n_children);

  // project to each child and copy
  // results
  for (unsigned int child = 0; child < n_children; ++child)
    {
      Quadrature<dim> help =
        project_to_child(reference_cell, quadrature, child);
      for (unsigned int i = 0; i < n_points; ++i)
        {
          q_points[child * n_points + i] = help.point(i);
          weights[child * n_points + i]  = help.weight(i);
        }
    }
  return Quadrature<dim>(q_points, weights);
}



template <int dim>
Quadrature<dim>
QProjector<dim>::project_to_line(const ReferenceCell &reference_cell,
                                 const Quadrature<1> &quadrature,
                                 const Point<dim>    &p1,
                                 const Point<dim>    &p2)
{
  Assert(reference_cell == ReferenceCells::get_hypercube<dim>(),
         ExcNotImplemented());
  (void)reference_cell;

  const unsigned int      n = quadrature.size();
  std::vector<Point<dim>> points(n);
  std::vector<double>     weights(n);
  const double            length = p1.distance(p2);

  for (unsigned int k = 0; k < n; ++k)
    {
      const double alpha = quadrature.point(k)[0];
      points[k]          = alpha * p2;
      points[k] += (1. - alpha) * p1;
      weights[k] = length * quadrature.weight(k);
    }
  return Quadrature<dim>(points, weights);
}



template <int dim>
typename QProjector<dim>::DataSetDescriptor
QProjector<dim>::DataSetDescriptor::face(const ReferenceCell &reference_cell,
                                         const unsigned int   face_no,
                                         const bool           face_orientation,
                                         const bool           face_flip,
                                         const bool           face_rotation,
                                         const unsigned int n_quadrature_points)
{
  return face(reference_cell,
              face_no,
              internal::combined_face_orientation(face_orientation,
                                                  face_rotation,
                                                  face_flip),
              n_quadrature_points);
}



template <int dim>
typename QProjector<dim>::DataSetDescriptor
QProjector<dim>::DataSetDescriptor::face(
  const ReferenceCell               &reference_cell,
  const unsigned int                 face_no,
  const types::geometric_orientation combined_orientation,
  const unsigned int                 n_quadrature_points)
{
  AssertIndexRange(face_no, reference_cell.n_faces());
  AssertIndexRange(combined_orientation,
                   reference_cell.n_face_orientations(face_no));
  AssertDimension(reference_cell.get_dimension(), dim);


  return {(reference_cell.n_face_orientations(face_no) * face_no +
           combined_orientation) *
          n_quadrature_points};
}



template <int dim>
typename QProjector<dim>::DataSetDescriptor
QProjector<dim>::DataSetDescriptor::face(
  const ReferenceCell            &reference_cell,
  const unsigned int              face_no,
  const bool                      face_orientation,
  const bool                      face_flip,
  const bool                      face_rotation,
  const hp::QCollection<dim - 1> &quadrature)
{
  return face(reference_cell,
              face_no,
              internal::combined_face_orientation(face_orientation,
                                                  face_rotation,
                                                  face_flip),
              quadrature);
}



template <int dim>
typename QProjector<dim>::DataSetDescriptor
QProjector<dim>::DataSetDescriptor::face(
  const ReferenceCell               &reference_cell,
  const unsigned int                 face_no,
  const types::geometric_orientation combined_orientation,
  const hp::QCollection<dim - 1>    &quadrature)
{
  AssertIndexRange(face_no, reference_cell.n_faces());
  AssertIndexRange(combined_orientation,
                   reference_cell.n_face_orientations(face_no));
  AssertDimension(reference_cell.get_dimension(), dim);

  unsigned int offset = 0;
  if (quadrature.size() == 1)
    offset =
      reference_cell.n_face_orientations(0) * quadrature[0].size() * face_no;
  else
    for (unsigned int i = 0; i < face_no; ++i)
      offset += reference_cell.n_face_orientations(i) * quadrature[i].size();

  return {offset + combined_orientation *
                     quadrature[quadrature.size() == 1 ? 0 : face_no].size()};
}



template <int dim>
typename QProjector<dim>::DataSetDescriptor
QProjector<dim>::DataSetDescriptor::subface(
  const ReferenceCell             &reference_cell,
  const unsigned int               face_no,
  const unsigned int               subface_no,
  const bool                       face_orientation,
  const bool                       face_flip,
  const bool                       face_rotation,
  const unsigned int               n_quadrature_points,
  const internal::SubfaceCase<dim> ref_case)
{
  return QProjector<dim>::DataSetDescriptor::subface(
    reference_cell,
    face_no,
    subface_no,
    internal::combined_face_orientation(face_orientation,
                                        face_rotation,
                                        face_flip),
    n_quadrature_points,
    ref_case);
}



template <>
QProjector<1>::DataSetDescriptor
QProjector<1>::DataSetDescriptor::subface(
  const ReferenceCell &reference_cell,
  const unsigned int   face_no,
  const unsigned int   subface_no,
  const types::geometric_orientation /*combined_orientation*/,
  const unsigned int n_quadrature_points,
  const internal::SubfaceCase<1>)
{
  Assert(reference_cell == ReferenceCells::Line, ExcNotImplemented());
  (void)reference_cell;

  Assert(face_no < GeometryInfo<1>::faces_per_cell, ExcInternalError());
  Assert(subface_no < GeometryInfo<1>::max_children_per_face,
         ExcInternalError());

  return ((face_no * GeometryInfo<1>::max_children_per_face + subface_no) *
          n_quadrature_points);
}



template <>
QProjector<2>::DataSetDescriptor
QProjector<2>::DataSetDescriptor::subface(
  const ReferenceCell               &reference_cell,
  const unsigned int                 face_no,
  const unsigned int                 subface_no,
  const types::geometric_orientation combined_orientation,
  const unsigned int                 n_quadrature_points,
  const internal::SubfaceCase<2>)
{
  Assert(reference_cell == ReferenceCells::Quadrilateral ||
           reference_cell == ReferenceCells::Triangle,
         ExcNotImplemented());

  const unsigned int n_faces = reference_cell.n_faces();
  const unsigned int n_subfaces =
    reference_cell.face_reference_cell(face_no).n_isotropic_children();
  const unsigned int n_orientations =
    reference_cell.n_face_orientations(face_no);

  AssertIndexRange(face_no, n_faces);
  AssertIndexRange(subface_no, n_subfaces);
  AssertIndexRange(combined_orientation, n_orientations);

  return ((face_no * n_orientations * n_subfaces +
           combined_orientation * n_subfaces + subface_no) *
          n_quadrature_points);
}



template <>
QProjector<3>::DataSetDescriptor
QProjector<3>::DataSetDescriptor::subface(
  const ReferenceCell               &reference_cell,
  const unsigned int                 face_no,
  const unsigned int                 subface_no,
  const types::geometric_orientation combined_orientation,
  const unsigned int                 n_quadrature_points,
  const internal::SubfaceCase<3>     ref_case)
{
  const unsigned int dim = 3;

  Assert(reference_cell == ReferenceCells::Hexahedron, ExcNotImplemented());
  (void)reference_cell;

  Assert(face_no < GeometryInfo<dim>::faces_per_cell, ExcInternalError());
  Assert(subface_no < GeometryInfo<dim>::max_children_per_face,
         ExcInternalError());

  // in 3d, we have to account for faces that
  // have non-standard orientation. thus, we
  // have to store _eight_ data sets per face
  // or subface already for the isotropic
  // case. Additionally, we have three
  // different refinement cases, resulting in
  // <tt>4 + 2 + 2 = 8</tt> different subfaces
  // for each face.
  const unsigned int total_subfaces_per_face = 8;

  // set up a table with the offsets for a
  // given refinement case respecting the
  // corresponding number of subfaces. the
  // index corresponds to (RefineCase::Type - 1)

  // note, that normally we should use the
  // obvious offsets 0,2,6. However, prior to
  // the implementation of anisotropic
  // refinement, in many places of the library
  // the convention was used, that the first
  // dataset with offset 0 corresponds to a
  // standard (isotropic) face
  // refinement. therefore we use the offsets
  // 6,4,0 here to stick to that (implicit)
  // convention
  static const unsigned int ref_case_offset[3] = {
    6, // cut_x
    4, // cut_y
    0  // cut_xy
  };

  const auto [final_subface_no, refinement_case] =
    reference_cell.equivalent_refinement_case(combined_orientation,
                                              ref_case,
                                              subface_no);

  return ((face_no * total_subfaces_per_face +
           ref_case_offset[refinement_case - 1] + final_subface_no) +
          GeometryInfo<dim>::faces_per_cell * total_subfaces_per_face *
            combined_orientation) *
         n_quadrature_points;
}



template <int dim>
Quadrature<dim>
QProjector<dim>::project_to_face(const ReferenceCell &reference_cell,
                                 const SubQuadrature &quadrature,
                                 const unsigned int   face_no)
{
  Assert(reference_cell == ReferenceCells::get_hypercube<dim>(),
         ExcNotImplemented());
  (void)reference_cell;

  std::vector<Point<dim>> points(quadrature.size());
  project_to_face(reference_cell, quadrature, face_no, points);
  return Quadrature<dim>(points, quadrature.get_weights());
}



template <int dim>
Quadrature<dim>
QProjector<dim>::project_to_subface(const ReferenceCell &reference_cell,
                                    const SubQuadrature &quadrature,
                                    const unsigned int   face_no,
                                    const unsigned int   subface_no,
                                    const RefinementCase<dim - 1> &ref_case)
{
  Assert(reference_cell == ReferenceCells::get_hypercube<dim>(),
         ExcNotImplemented());
  (void)reference_cell;

  std::vector<Point<dim>> points(quadrature.size());
  project_to_subface(
    reference_cell, quadrature, face_no, subface_no, points, ref_case);
  return Quadrature<dim>(points, quadrature.get_weights());
}


// explicit instantiations; note: we need them all for all dimensions
template class QProjector<1>;
template class QProjector<2>;
template class QProjector<3>;

DEAL_II_NAMESPACE_CLOSE
