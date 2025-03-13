// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2024 by the deal.II authors
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
#include <deal.II/base/polynomials_barycentric.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/base/tensor_product_polynomials.h>

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
            case 0:
              return reflect(points);
            case 2:
              return rotate(reflect(points), 3);
            case 4:
              return rotate(reflect(points), 2);
            case 6:
              return rotate(reflect(points), 1);
            case 1:
              return points;
            case 3:
              return rotate(points, 1);
            case 5:
              return rotate(points, 2);
            case 7:
              return rotate(points, 3);
            default:
              DEAL_II_ASSERT_UNREACHABLE();
          }
        return {};
      }

      Quadrature<2>
      mutate_quadrature(const Quadrature<2>               &quadrature,
                        const types::geometric_orientation combined_orientation)
      {
        return Quadrature<2>(mutate_points_with_offset(quadrature.get_points(),
                                                       combined_orientation),
                             quadrature.get_weights());
      }

      std::pair<unsigned int, RefinementCase<2>>
      select_subface_no_and_refinement_case(
        const unsigned int             subface_no,
        const bool                     face_orientation,
        const bool                     face_flip,
        const bool                     face_rotation,
        const internal::SubfaceCase<3> ref_case)
      {
        constexpr int dim = 3;
        // for each subface of a given FaceRefineCase
        // there is a corresponding equivalent
        // subface number of one of the "standard"
        // RefineCases (cut_x, cut_y, cut_xy). Map
        // the given values to those equivalent
        // ones.

        // first, define an invalid number
        static const unsigned int e = numbers::invalid_unsigned_int;

        static const RefinementCase<dim - 1>
          equivalent_refine_case[internal::SubfaceCase<dim>::case_isotropic + 1]
                                [GeometryInfo<3>::max_children_per_face] = {
                                  // case_none. there should be only
                                  // invalid values here. However, as
                                  // this function is also called (in
                                  // tests) for cells which have no
                                  // refined faces, use isotropic
                                  // refinement instead
                                  {RefinementCase<dim - 1>::cut_xy,
                                   RefinementCase<dim - 1>::cut_xy,
                                   RefinementCase<dim - 1>::cut_xy,
                                   RefinementCase<dim - 1>::cut_xy},
                                  // case_x
                                  {RefinementCase<dim - 1>::cut_x,
                                   RefinementCase<dim - 1>::cut_x,
                                   RefinementCase<dim - 1>::no_refinement,
                                   RefinementCase<dim - 1>::no_refinement},
                                  // case_x1y
                                  {RefinementCase<dim - 1>::cut_xy,
                                   RefinementCase<dim - 1>::cut_xy,
                                   RefinementCase<dim - 1>::cut_x,
                                   RefinementCase<dim - 1>::no_refinement},
                                  // case_x2y
                                  {RefinementCase<dim - 1>::cut_x,
                                   RefinementCase<dim - 1>::cut_xy,
                                   RefinementCase<dim - 1>::cut_xy,
                                   RefinementCase<dim - 1>::no_refinement},
                                  // case_x1y2y
                                  {RefinementCase<dim - 1>::cut_xy,
                                   RefinementCase<dim - 1>::cut_xy,
                                   RefinementCase<dim - 1>::cut_xy,
                                   RefinementCase<dim - 1>::cut_xy},
                                  // case_y
                                  {RefinementCase<dim - 1>::cut_y,
                                   RefinementCase<dim - 1>::cut_y,
                                   RefinementCase<dim - 1>::no_refinement,
                                   RefinementCase<dim - 1>::no_refinement},
                                  // case_y1x
                                  {RefinementCase<dim - 1>::cut_xy,
                                   RefinementCase<dim - 1>::cut_xy,
                                   RefinementCase<dim - 1>::cut_y,
                                   RefinementCase<dim - 1>::no_refinement},
                                  // case_y2x
                                  {RefinementCase<dim - 1>::cut_y,
                                   RefinementCase<dim - 1>::cut_xy,
                                   RefinementCase<dim - 1>::cut_xy,
                                   RefinementCase<dim - 1>::no_refinement},
                                  // case_y1x2x
                                  {RefinementCase<dim - 1>::cut_xy,
                                   RefinementCase<dim - 1>::cut_xy,
                                   RefinementCase<dim - 1>::cut_xy,
                                   RefinementCase<dim - 1>::cut_xy},
                                  // case_xy (case_isotropic)
                                  {RefinementCase<dim - 1>::cut_xy,
                                   RefinementCase<dim - 1>::cut_xy,
                                   RefinementCase<dim - 1>::cut_xy,
                                   RefinementCase<dim - 1>::cut_xy}};

        static const unsigned int
          equivalent_subface_number[internal::SubfaceCase<dim>::case_isotropic +
                                    1][GeometryInfo<3>::max_children_per_face] =
            {// case_none, see above
             {0, 1, 2, 3},
             // case_x
             {0, 1, e, e},
             // case_x1y
             {0, 2, 1, e},
             // case_x2y
             {0, 1, 3, e},
             // case_x1y2y
             {0, 2, 1, 3},
             // case_y
             {0, 1, e, e},
             // case_y1x
             {0, 1, 1, e},
             // case_y2x
             {0, 2, 3, e},
             // case_y1x2x
             {0, 1, 2, 3},
             // case_xy (case_isotropic)
             {0, 1, 2, 3}};

        // If face-orientation or face_rotation are
        // non-standard, cut_x and cut_y have to be
        // exchanged.
        static const RefinementCase<dim - 1> ref_case_permutation[4] = {
          RefinementCase<dim - 1>::no_refinement,
          RefinementCase<dim - 1>::cut_y,
          RefinementCase<dim - 1>::cut_x,
          RefinementCase<dim - 1>::cut_xy};

        // set a corresponding (equivalent)
        // RefineCase and subface number
        const RefinementCase<dim - 1> equ_ref_case =
          equivalent_refine_case[ref_case][subface_no];
        const unsigned int equ_subface_no =
          equivalent_subface_number[ref_case][subface_no];
        // make sure, that we got a valid subface and RefineCase
        Assert(equ_ref_case != RefinementCase<dim>::no_refinement,
               ExcInternalError());
        Assert(equ_subface_no != e, ExcInternalError());
        // now, finally respect non-standard faces
        const RefinementCase<dim - 1> final_ref_case =
          (face_orientation == face_rotation ?
             ref_case_permutation[equ_ref_case] :
             equ_ref_case);

        const unsigned int final_subface_no =
          GeometryInfo<dim>::child_cell_on_face(RefinementCase<dim>(
                                                  final_ref_case),
                                                4,
                                                equ_subface_no,
                                                face_orientation,
                                                face_flip,
                                                face_rotation,
                                                equ_ref_case);

        return std::make_pair(final_subface_no, final_ref_case);
      }


      Quadrature<2>
      project_to_all_faces_or_subfaces(const ReferenceCell      &reference_cell,
                                       const hp::QCollection<1> &quadrature,
                                       const bool project_to_subface)
      {
        Assert(reference_cell == ReferenceCells::Triangle ||
                 reference_cell == ReferenceCells::Quadrilateral,
               ExcInternalError());

        constexpr unsigned int dim = 2;

        const unsigned int n_faces = reference_cell.n_faces();
        const unsigned int num_subfaces_or_face =
          project_to_subface ?
            reference_cell.face_reference_cell(0).n_isotropic_children() :
            1;

        unsigned int n_points = 0;
        if (quadrature.size() == 1)
          n_points = quadrature[0].size() * n_faces;
        else
          {
            AssertDimension(quadrature.size(), n_faces);
            for (const auto &q : quadrature)
              n_points += q.size();
          }

        // new (projected) quadrature points and weights
        std::vector<Point<2>> q_points;
        q_points.reserve(2 * n_points * num_subfaces_or_face);
        std::vector<double> weights;
        weights.reserve(2 * n_points * num_subfaces_or_face);

        // loop over all faces ...
        for (unsigned int face = 0; face < n_faces; ++face)
          // ... and over all possible orientations ...
          for (unsigned int orientation = 0; orientation < 2; ++orientation)
            // ... and over all subfaces
            for (unsigned int subface = 0; subface < num_subfaces_or_face;
                 ++subface)
              {
                const unsigned int quadrature_index =
                  quadrature.size() == 1 ? 0 : face;
                const unsigned int n_q_points =
                  quadrature[quadrature_index].size();

                // handle the points
                std::vector<Point<dim>> points(n_q_points);

                // orientation numbers::default_geometric_orientation is the
                // standard orientation so invert lines for orientation
                // numbers::reverse_line_orientation
                if (project_to_subface)
                  {
                    const unsigned int local_subface =
                      orientation == numbers::reverse_line_orientation ?
                        num_subfaces_or_face - 1 - subface :
                        subface;

                    dealii::QProjector<dim>::project_to_subface(
                      reference_cell,
                      quadrature[quadrature_index],
                      face,
                      local_subface,
                      points);
                  }
                else
                  {
                    dealii::QProjector<dim>::project_to_face(
                      reference_cell,
                      quadrature[quadrature_index],
                      face,
                      points);
                  }

                if (orientation == numbers::reverse_line_orientation)
                  std::reverse(points.begin(), points.end());

                std::copy(points.begin(),
                          points.end(),
                          std::back_inserter(q_points));


                // handle the weights
                std::vector<double> scaled_weights =
                  quadrature[quadrature_index].get_weights();

                // triangular face 1 is given by (1,0) (0,1) so scale the
                // weights by the length
                if (reference_cell == ReferenceCells::Triangle && face == 1)
                  for (auto &w : scaled_weights)
                    w *= std::sqrt(2);

                if (orientation == numbers::reverse_line_orientation)
                  std::reverse(scaled_weights.begin(), scaled_weights.end());

                std::copy(scaled_weights.begin(),
                          scaled_weights.end(),
                          std::back_inserter(weights));
              }

        Assert(q_points.size() == 2 * n_points * num_subfaces_or_face,
               ExcInternalError());
        Assert(weights.size() == 2 * n_points * num_subfaces_or_face,
               ExcInternalError());

        // construct new quadrature rule
        return Quadrature<2>(std::move(q_points), std::move(weights));
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
QProjector<dim>::project_to_face(const ReferenceCell       &reference_cell,
                                 const Quadrature<dim - 1> &quadrature,
                                 const unsigned int         face_no,
                                 const types::geometric_orientation orientation)
{
  AssertIndexRange(face_no, reference_cell.n_faces());
  // TODO once the default orientation is zero we can remove this special case
  AssertIndexRange(dim == 1 ? 1 - orientation : orientation,
                   reference_cell.n_face_orientations(face_no));
  AssertDimension(reference_cell.get_dimension(), dim);

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
              q_points.emplace_back(quadrature.point(p)[0], 0);
            else if (face_no == 1)
              q_points.emplace_back(1.0 - quadrature.point(p)[0],
                                    quadrature.point(p)[0]);
            else if (face_no == 2)
              q_points.emplace_back(0, 1.0 - quadrature.point(p)[0]);
            else
              DEAL_II_ASSERT_UNREACHABLE();
          }
      else if (reference_cell == ReferenceCells::Quadrilateral)
        for (unsigned int p = 0; p < quadrature.size(); ++p)
          {
            if (face_no == 0)
              q_points.emplace_back(0, quadrature.point(p)[0]);
            else if (face_no == 1)
              q_points.emplace_back(1, quadrature.point(p)[0]);
            else if (face_no == 2)
              q_points.emplace_back(quadrature.point(p)[0], 0);
            else if (face_no == 3)
              q_points.emplace_back(quadrature.point(p)[0], 1);
            else
              DEAL_II_ASSERT_UNREACHABLE();
          }
      else
        DEAL_II_ASSERT_UNREACHABLE();

      if (orientation == numbers::reverse_line_orientation)
        {
          std::reverse(q_points.begin(), q_points.end());
          std::reverse(q_weights.begin(), q_weights.end());
        }
    }
  else if constexpr (dim == 3)
    {
      Assert(reference_cell == ReferenceCells::Hexahedron, ExcNotImplemented());
      const Quadrature<2> mutation =
        internal::QProjector::mutate_quadrature(quadrature, orientation);
      return QProjector<3>::project_to_face(reference_cell, mutation, face_no);
    }
  else
    {
      DEAL_II_ASSERT_UNREACHABLE();
    }

  return Quadrature<dim>(std::move(q_points), std::move(q_weights));
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
  const bool                 face_orientation,
  const bool,
  const bool,
  const internal::SubfaceCase<dim>)
{
  if (!face_orientation)
    DEAL_II_NOT_IMPLEMENTED();
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
  const types::geometric_orientation orientation,
  const RefinementCase<dim - 1>     &ref_case)
{
  AssertIndexRange(face_no, reference_cell.n_faces());
  // TODO once the default orientation is zero we can remove this special case
  AssertIndexRange(dim == 1 ? 1 - orientation : orientation,
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

      if (orientation == numbers::reverse_line_orientation)
        {
          std::reverse(q_points.begin(), q_points.end());
          std::reverse(q_weights.begin(), q_weights.end());
        }
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



template <>
Quadrature<1>
QProjector<1>::project_to_all_faces(const ReferenceCell      &reference_cell,
                                    const hp::QCollection<0> &quadrature)
{
  AssertDimension(quadrature.size(), 1);
  Assert(reference_cell == ReferenceCells::Line, ExcNotImplemented());
  (void)reference_cell;

  const unsigned int dim = 1;

  const unsigned int n_points = 1, n_faces = GeometryInfo<dim>::faces_per_cell;

  // first fix quadrature points
  std::vector<Point<dim>> q_points;
  q_points.reserve(n_points * n_faces);
  std::vector<Point<dim>> help(n_points);


  // project to each face and append
  // results
  for (unsigned int face = 0; face < n_faces; ++face)
    {
      project_to_face(reference_cell,
                      quadrature[quadrature.size() == 1 ? 0 : face],
                      face,
                      help);
      std::copy(help.begin(), help.end(), std::back_inserter(q_points));
    }

  // next copy over weights
  std::vector<double> weights;
  weights.reserve(n_points * n_faces);
  for (unsigned int face = 0; face < n_faces; ++face)
    std::copy(
      quadrature[quadrature.size() == 1 ? 0 : face].get_weights().begin(),
      quadrature[quadrature.size() == 1 ? 0 : face].get_weights().end(),
      std::back_inserter(weights));

  Assert(q_points.size() == n_points * n_faces, ExcInternalError());
  Assert(weights.size() == n_points * n_faces, ExcInternalError());

  return Quadrature<dim>(std::move(q_points), std::move(weights));
}



template <>
Quadrature<2>
QProjector<2>::project_to_all_faces(const ReferenceCell      &reference_cell,
                                    const hp::QCollection<1> &quadrature)
{
  return internal::QProjector::project_to_all_faces_or_subfaces(reference_cell,
                                                                quadrature,
                                                                false);
}



template <>
Quadrature<3>
QProjector<3>::project_to_all_faces(const ReferenceCell      &reference_cell,
                                    const hp::QCollection<2> &quadrature)
{
  const auto process = [&](const std::vector<std::vector<Point<3>>> &faces) {
    // new (projected) quadrature points and weights
    std::vector<Point<3>> points;
    std::vector<double>   weights;

    const auto poly_tri = BarycentricPolynomials<2>::get_fe_p_basis(1);
    const TensorProductPolynomials<2> poly_quad(
      Polynomials::generate_complete_Lagrange_basis(
        {Point<1>(0.0), Point<1>(1.0)}));

    // loop over all faces (triangles) ...
    for (unsigned int face_no = 0; face_no < faces.size(); ++face_no)
      {
        const ReferenceCell face_reference_cell =
          reference_cell.face_reference_cell(face_no);
        // We will use linear polynomials to map the reference quadrature
        // points correctly to on faces. There are as many linear shape
        // functions as there are vertices in the face.
        const unsigned int n_linear_shape_functions = faces[face_no].size();
        std::vector<Tensor<1, 2>> shape_derivatives;

        const auto &poly =
          (n_linear_shape_functions == 3 ?
             static_cast<const ScalarPolynomialsBase<2> &>(poly_tri) :
             static_cast<const ScalarPolynomialsBase<2> &>(poly_quad));

        // ... and over all possible orientations
        for (types::geometric_orientation orientation = 0;
             orientation < reference_cell.n_face_orientations(face_no);
             ++orientation)
          {
            const auto &face = faces[face_no];

            // The goal of this function is to compute identical sets of
            // quadrature points on the common face of two abutting cells. Our
            // orientation convention is that, given such a pair of abutting
            // cells:
            //
            // 1. The shared face, from the perspective of the first cell, is
            //    in the default orientation.
            // 2. The shared face, from the perspective of the second cell, has
            //    its orientation computed relative to the first cell: i.e.,
            //    'orientation' is the vertex permutation applied to the first
            //    cell's face to get the second cell's face.
            //
            // The first case is trivial since points do not need to be
            // oriented. However, in the second case, we need to use the
            // *reverse* of the stored orientation (i.e., the permutation
            // applied to the second cell's face which yields the first cell's
            // face) so that we get identical quadrature points.
            //
            // For more information see connectivity.h.
            const boost::container::small_vector<Point<3>, 8> support_points =
              face_reference_cell.permute_by_combined_orientation<Point<3>>(
                face,
                face_reference_cell.get_inverse_combined_orientation(
                  orientation));

            // the quadrature rule to be projected ...
            const auto &sub_quadrature_points =
              quadrature[quadrature.size() == 1 ? 0 : face_no].get_points();
            const auto &sub_quadrature_weights =
              quadrature[quadrature.size() == 1 ? 0 : face_no].get_weights();

            // loop over all quadrature points
            for (unsigned int j = 0; j < sub_quadrature_points.size(); ++j)
              {
                Point<3> mapped_point;

                // map reference quadrature point
                for (unsigned int i = 0; i < n_linear_shape_functions; ++i)
                  mapped_point +=
                    support_points[i] *
                    poly.compute_value(i, sub_quadrature_points[j]);

                points.push_back(mapped_point);

                // scale quadrature weight
                const double scaling = [&]() {
                  const unsigned int dim_     = 2;
                  const unsigned int spacedim = 3;

                  Tensor<1, dim_, Tensor<1, spacedim>> DX_t;

                  shape_derivatives.resize(n_linear_shape_functions);

                  for (unsigned int i = 0; i < n_linear_shape_functions; ++i)
                    shape_derivatives[i] =
                      poly.compute_1st_derivative(i, sub_quadrature_points[j]);

                  for (unsigned int k = 0; k < n_linear_shape_functions; ++k)
                    for (unsigned int i = 0; i < spacedim; ++i)
                      for (unsigned int j = 0; j < dim_; ++j)
                        DX_t[j][i] +=
                          shape_derivatives[k][j] * support_points[k][i];

                  Tensor<2, dim_> G;
                  for (unsigned int i = 0; i < dim_; ++i)
                    for (unsigned int j = 0; j < dim_; ++j)
                      G[i][j] = DX_t[i] * DX_t[j];

                  return std::sqrt(determinant(G));
                }();

                weights.push_back(sub_quadrature_weights[j] * scaling);
              }
          }
      }

    // construct new quadrature rule
    return Quadrature<3>(std::move(points), std::move(weights));
  };

  if ((reference_cell == ReferenceCells::Tetrahedron) ||
      (reference_cell == ReferenceCells::Wedge) ||
      (reference_cell == ReferenceCells::Pyramid))
    {
      std::vector<std::vector<Point<3>>> face_vertex_locations(
        reference_cell.n_faces());
      for (const unsigned int f : reference_cell.face_indices())
        {
          face_vertex_locations[f].resize(
            reference_cell.face_reference_cell(f).n_vertices());
          for (const unsigned int v :
               reference_cell.face_reference_cell(f).vertex_indices())
            face_vertex_locations[f][v] =
              reference_cell.face_vertex_location<3>(f, v);
        }

      return process(face_vertex_locations);
    }
  else
    {
      Assert(reference_cell == ReferenceCells::Hexahedron, ExcNotImplemented());

      const unsigned int dim = 3;

      unsigned int n_points_total = 0;

      if (quadrature.size() == 1)
        n_points_total =
          quadrature[0].size() * GeometryInfo<dim>::faces_per_cell;
      else
        {
          AssertDimension(quadrature.size(), GeometryInfo<dim>::faces_per_cell);
          for (const auto &q : quadrature)
            n_points_total += q.size();
        }

      n_points_total *= 8;

      // first fix quadrature points
      std::vector<Point<dim>> q_points;
      q_points.reserve(n_points_total);

      std::vector<double> weights;
      weights.reserve(n_points_total);

      for (unsigned char offset = 0; offset < 8; ++offset)
        {
          const auto mutation = internal::QProjector::mutate_points_with_offset(
            quadrature[0].get_points(), offset);
          // project to each face and append results
          for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell;
               ++face)
            {
              const unsigned int q_index = quadrature.size() == 1 ? 0 : face;

              internal::QProjector::project_to_hex_face_and_append(
                q_index > 0 ? internal::QProjector::mutate_points_with_offset(
                                quadrature[face].get_points(), offset) :
                              mutation,
                face,
                q_points);

              std::copy(quadrature[q_index].get_weights().begin(),
                        quadrature[q_index].get_weights().end(),
                        std::back_inserter(weights));
            }
        }

      Assert(q_points.size() == n_points_total, ExcInternalError());
      Assert(weights.size() == n_points_total, ExcInternalError());

      return Quadrature<dim>(std::move(q_points), std::move(weights));
    }
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
  return internal::QProjector::project_to_all_faces_or_subfaces(
    reference_cell, hp::QCollection<1>(quadrature), true);
}



template <>
Quadrature<3>
QProjector<3>::project_to_all_subfaces(const ReferenceCell &reference_cell,
                                       const SubQuadrature &quadrature)
{
  if (reference_cell == ReferenceCells::Tetrahedron)
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
  // TODO: once we move to representing the default orientation as 0 (instead of
  // 1) we can get rid of the dim = 1 check
  Assert(dim == 1 ||
           (combined_orientation < reference_cell.n_face_orientations(face_no)),
         ExcInternalError());

  Assert(reference_cell == ReferenceCells::get_hypercube<dim>() ||
           reference_cell == ReferenceCells::get_simplex<dim>(),
         ExcNotImplemented());

  Assert(face_no < reference_cell.n_faces(), ExcInternalError());

  switch (dim)
    {
      if (dim == 2)
        return {(2 * face_no + (combined_orientation ==
                                    numbers::default_geometric_orientation ?
                                  1 :
                                  0)) *
                n_quadrature_points};
      case 3:
        if (reference_cell == ReferenceCells::Tetrahedron)
          return {(reference_cell.n_face_orientations(face_no) * face_no +
                   combined_orientation) *
                  n_quadrature_points};
        return (face_no +
                GeometryInfo<dim>::faces_per_cell * combined_orientation) *
               n_quadrature_points;
      default:
        DEAL_II_ASSERT_UNREACHABLE();
    }
  return numbers::invalid_unsigned_int;
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
  if (reference_cell == ReferenceCells::Triangle ||
      reference_cell == ReferenceCells::Tetrahedron ||
      reference_cell == ReferenceCells::Wedge ||
      reference_cell == ReferenceCells::Pyramid ||
      reference_cell == ReferenceCells::Quadrilateral)
    {
      unsigned int offset = 0;
      Assert(combined_orientation < reference_cell.n_face_orientations(face_no),
             ExcInternalError());

      static const unsigned int X = numbers::invalid_unsigned_int;
      static const std::array<unsigned int, 5> scale_tri   = {{2, 2, 2, X, X}};
      static const std::array<unsigned int, 5> scale_quad  = {{2, 2, 2, 2, X}};
      static const std::array<unsigned int, 5> scale_tet   = {{6, 6, 6, 6, X}};
      static const std::array<unsigned int, 5> scale_wedge = {{6, 6, 8, 8, 8}};
      static const std::array<unsigned int, 5> scale_pyramid = {
        {8, 6, 6, 6, 6}};

      const auto &scale =
        (reference_cell == ReferenceCells::Quadrilateral) ?
          scale_quad :
          ((reference_cell == ReferenceCells::Triangle) ?
             scale_tri :
             ((reference_cell == ReferenceCells::Tetrahedron) ?
                scale_tet :
                ((reference_cell == ReferenceCells::Wedge) ? scale_wedge :
                                                             scale_pyramid)));

      if (quadrature.size() == 1)
        offset = scale[0] * quadrature[0].size() * face_no;
      else
        for (unsigned int i = 0; i < face_no; ++i)
          offset += scale[i] * quadrature[i].size();

      if (dim == 2)
        return {
          offset +
          (combined_orientation == numbers::default_geometric_orientation) *
            quadrature[quadrature.size() == 1 ? 0 : face_no].size()};
      else if (dim == 3)
        {
          return {offset +
                  combined_orientation *
                    quadrature[quadrature.size() == 1 ? 0 : face_no].size()};
        }
    }

  Assert(reference_cell == ReferenceCells::get_hypercube<dim>(),
         ExcNotImplemented());

  Assert(face_no < GeometryInfo<dim>::faces_per_cell, ExcInternalError());

  switch (dim)
    {
      case 1:
        {
          if (quadrature.size() == 1)
            return quadrature[0].size() * face_no;
          else
            {
              unsigned int result = 0;
              for (unsigned int i = 0; i < face_no; ++i)
                result += quadrature[i].size();
              return result;
            }
        }
      case 3:
        {
          if (quadrature.size() == 1)
            return (face_no +
                    combined_orientation * GeometryInfo<dim>::faces_per_cell) *
                   quadrature[0].size();
          else
            {
              unsigned int n_points_i = 0;
              for (unsigned int i = 0; i < face_no; ++i)
                n_points_i += quadrature[i].size();

              unsigned int n_points = 0;
              for (unsigned int i = 0; i < GeometryInfo<dim>::faces_per_cell;
                   ++i)
                n_points += quadrature[i].size();

              return n_points_i + combined_orientation * n_points;
            }
        }

      default:
        DEAL_II_ASSERT_UNREACHABLE();
    }
  return numbers::invalid_unsigned_int;
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
  Assert(face_no < reference_cell.n_faces(), ExcInternalError());
  Assert(subface_no <
           reference_cell.face_reference_cell(0).n_isotropic_children(),
         ExcInternalError());

  Assert(reference_cell == ReferenceCells::Triangle ||
           reference_cell == ReferenceCells::Quadrilateral,
         ExcInternalError());

  const unsigned int orientation_offset =
    combined_orientation == ReferenceCell::default_combined_face_orientation() ?
      1 :
      0;
  const unsigned int n_subfaces =
    reference_cell.face_reference_cell(0).n_isotropic_children();
  return (2 * face_no * n_subfaces + n_subfaces * orientation_offset +
          subface_no) *
         n_quadrature_points;
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

  const auto [face_orientation, face_rotation, face_flip] =
    internal::split_face_orientation(combined_orientation);

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

  const std::pair<unsigned int, RefinementCase<2>>
    final_subface_no_and_ref_case =
      internal::QProjector::select_subface_no_and_refinement_case(
        subface_no, face_orientation, face_flip, face_rotation, ref_case);

  return ((face_no * total_subfaces_per_face +
           ref_case_offset[final_subface_no_and_ref_case.second - 1] +
           final_subface_no_and_ref_case.first) +
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
