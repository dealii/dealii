// ---------------------------------------------------------------------
//
// Copyright (C) 2020 - 2022 by the deal.II authors
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

#include <deal.II/base/derivative_form.h>
#include <deal.II/base/geometry_info.h>
#include <deal.II/base/polynomials_barycentric.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/base/tensor_product_polynomials.h>

DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace QProjector
  {
    namespace
    {
      Quadrature<2>
      reflect(const Quadrature<2> &q)
      {
        // Take the points and reflect them by the diagonal
        std::vector<Point<2>> q_points(q.get_points());
        for (Point<2> &p : q_points)
          std::swap(p[0], p[1]);

        return Quadrature<2>(q_points, q.get_weights());
      }


      Quadrature<2>
      rotate(const Quadrature<2> &q, const unsigned int n_times)
      {
        std::vector<Point<2>> q_points(q.size());
        for (unsigned int i = 0; i < q.size(); ++i)
          {
            switch (n_times % 4)
              {
                case 0:
                  // 0 degree. the point remains as it is.
                  q_points[i] = q.point(i);
                  break;

                case 1:
                  // 90 degree counterclockwise
                  q_points[i][0] = 1.0 - q.point(i)[1];
                  q_points[i][1] = q.point(i)[0];
                  break;
                case 2:
                  // 180 degree counterclockwise
                  q_points[i][0] = 1.0 - q.point(i)[0];
                  q_points[i][1] = 1.0 - q.point(i)[1];
                  break;
                case 3:
                  // 270 degree counterclockwise
                  q_points[i][0] = q.point(i)[1];
                  q_points[i][1] = 1.0 - q.point(i)[0];
                  break;
              }
          }

        return Quadrature<2>(q_points, q.get_weights());
      }
    } // namespace
  }   // namespace QProjector
} // namespace internal



template <>
void
QProjector<1>::project_to_face(const Quadrature<0> &  quadrature,
                               const unsigned int     face_no,
                               std::vector<Point<1>> &q_points)
{
  project_to_face(ReferenceCells::Line, quadrature, face_no, q_points);
}



template <>
void
QProjector<1>::project_to_face(const ReferenceCell reference_cell,
                               const Quadrature<0> &,
                               const unsigned int     face_no,
                               std::vector<Point<1>> &q_points)
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
QProjector<2>::project_to_face(const Quadrature<1> &  quadrature,
                               const unsigned int     face_no,
                               std::vector<Point<2>> &q_points)
{
  project_to_face(ReferenceCells::Quadrilateral, quadrature, face_no, q_points);
}



template <>
void
QProjector<2>::project_to_face(const ReferenceCell    reference_cell,
                               const Quadrature<1> &  quadrature,
                               const unsigned int     face_no,
                               std::vector<Point<2>> &q_points)
{
  const unsigned int dim = 2;
  AssertIndexRange(face_no, GeometryInfo<dim>::faces_per_cell);
  Assert(q_points.size() == quadrature.size(),
         ExcDimensionMismatch(q_points.size(), quadrature.size()));

  if (reference_cell == ReferenceCells::Triangle)
    {
      // use linear polynomial to map the reference quadrature points correctly
      // on faces, i.e., BarycentricPolynomials<1>(1)
      for (unsigned int p = 0; p < quadrature.size(); ++p)
        switch (face_no)
          {
            case 0:
              q_points[p] = Point<dim>(quadrature.point(p)(0), 0);
              break;
            case 1:
              q_points[p] =
                Point<dim>(1 - quadrature.point(p)(0), quadrature.point(p)(0));
              break;
            case 2:
              q_points[p] = Point<dim>(0, 1 - quadrature.point(p)(0));
              break;
            default:
              Assert(false, ExcInternalError());
          }
    }
  else if (reference_cell == ReferenceCells::Quadrilateral)
    {
      for (unsigned int p = 0; p < quadrature.size(); ++p)
        switch (face_no)
          {
            case 0:
              q_points[p] = Point<dim>(0, quadrature.point(p)(0));
              break;
            case 1:
              q_points[p] = Point<dim>(1, quadrature.point(p)(0));
              break;
            case 2:
              q_points[p] = Point<dim>(quadrature.point(p)(0), 0);
              break;
            case 3:
              q_points[p] = Point<dim>(quadrature.point(p)(0), 1);
              break;
            default:
              Assert(false, ExcInternalError());
          }
    }
  else
    {
      Assert(false, ExcInternalError());
    }
}



template <>
void
QProjector<3>::project_to_face(const Quadrature<2> &  quadrature,
                               const unsigned int     face_no,
                               std::vector<Point<3>> &q_points)
{
  project_to_face(ReferenceCells::Hexahedron, quadrature, face_no, q_points);
}



template <>
void
QProjector<3>::project_to_face(const ReferenceCell    reference_cell,
                               const Quadrature<2> &  quadrature,
                               const unsigned int     face_no,
                               std::vector<Point<3>> &q_points)
{
  Assert(reference_cell == ReferenceCells::Hexahedron, ExcNotImplemented());
  (void)reference_cell;

  const unsigned int dim = 3;
  AssertIndexRange(face_no, GeometryInfo<dim>::faces_per_cell);
  Assert(q_points.size() == quadrature.size(),
         ExcDimensionMismatch(q_points.size(), quadrature.size()));

  for (unsigned int p = 0; p < quadrature.size(); ++p)
    switch (face_no)
      {
        case 0:
          q_points[p] =
            Point<dim>(0, quadrature.point(p)(0), quadrature.point(p)(1));
          break;
        case 1:
          q_points[p] =
            Point<dim>(1, quadrature.point(p)(0), quadrature.point(p)(1));
          break;
        case 2:
          q_points[p] =
            Point<dim>(quadrature.point(p)(1), 0, quadrature.point(p)(0));
          break;
        case 3:
          q_points[p] =
            Point<dim>(quadrature.point(p)(1), 1, quadrature.point(p)(0));
          break;
        case 4:
          q_points[p] =
            Point<dim>(quadrature.point(p)(0), quadrature.point(p)(1), 0);
          break;
        case 5:
          q_points[p] =
            Point<dim>(quadrature.point(p)(0), quadrature.point(p)(1), 1);
          break;

        default:
          Assert(false, ExcInternalError());
      }
}



template <>
void
QProjector<1>::project_to_subface(const Quadrature<0> &    quadrature,
                                  const unsigned int       face_no,
                                  const unsigned int       subface_no,
                                  std::vector<Point<1>> &  q_points,
                                  const RefinementCase<0> &ref_case)
{
  project_to_subface(
    ReferenceCells::Line, quadrature, face_no, subface_no, q_points, ref_case);
}



template <>
void
QProjector<1>::project_to_subface(const ReferenceCell reference_cell,
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
QProjector<2>::project_to_subface(const Quadrature<1> &    quadrature,
                                  const unsigned int       face_no,
                                  const unsigned int       subface_no,
                                  std::vector<Point<2>> &  q_points,
                                  const RefinementCase<1> &ref_case)
{
  project_to_subface(ReferenceCells::Quadrilateral,
                     quadrature,
                     face_no,
                     subface_no,
                     q_points,
                     ref_case);
}



template <>
void
QProjector<2>::project_to_subface(const ReferenceCell    reference_cell,
                                  const Quadrature<1> &  quadrature,
                                  const unsigned int     face_no,
                                  const unsigned int     subface_no,
                                  std::vector<Point<2>> &q_points,
                                  const RefinementCase<1> &)
{
  const unsigned int dim = 2;
  AssertIndexRange(face_no, GeometryInfo<dim>::faces_per_cell);
  AssertIndexRange(subface_no, GeometryInfo<dim>::max_children_per_face);

  Assert(q_points.size() == quadrature.size(),
         ExcDimensionMismatch(q_points.size(), quadrature.size()));

  if (reference_cell == ReferenceCells::Triangle)
    {
      // use linear polynomial to map the reference quadrature points correctly
      // on faces, i.e., BarycentricPolynomials<1>(1)
      for (unsigned int p = 0; p < quadrature.size(); ++p)
        switch (face_no)
          {
            case 0:
              switch (subface_no)
                {
                  case 0:
                    q_points[p] = Point<dim>(quadrature.point(p)(0) / 2, 0);
                    break;
                  case 1:
                    q_points[p] =
                      Point<dim>(0.5 + quadrature.point(p)(0) / 2, 0);
                    break;
                  default:
                    Assert(false, ExcInternalError());
                }
              break;
            case 1:
              switch (subface_no)
                {
                  case 0:
                    q_points[p] = Point<dim>(1 - quadrature.point(p)(0) / 2,
                                             quadrature.point(p)(0) / 2);
                    break;
                  case 1:
                    q_points[p] = Point<dim>(0.5 - quadrature.point(p)(0) / 2,
                                             0.5 + quadrature.point(p)(0) / 2);
                    break;
                  default:
                    Assert(false, ExcInternalError());
                }
              break;
            case 2:
              switch (subface_no)
                {
                  case 0:
                    q_points[p] = Point<dim>(0, 1 - quadrature.point(p)(0) / 2);
                    break;
                  case 1:
                    q_points[p] =
                      Point<dim>(0, 0.5 - quadrature.point(p)(0) / 2);
                    break;
                  default:
                    Assert(false, ExcInternalError());
                }
              break;
            default:
              Assert(false, ExcInternalError());
          }
    }
  else if (reference_cell == ReferenceCells::Quadrilateral)
    {
      for (unsigned int p = 0; p < quadrature.size(); ++p)
        switch (face_no)
          {
            case 0:
              switch (subface_no)
                {
                  case 0:
                    q_points[p] = Point<dim>(0, quadrature.point(p)(0) / 2);
                    break;
                  case 1:
                    q_points[p] =
                      Point<dim>(0, quadrature.point(p)(0) / 2 + 0.5);
                    break;
                  default:
                    Assert(false, ExcInternalError());
                }
              break;
            case 1:
              switch (subface_no)
                {
                  case 0:
                    q_points[p] = Point<dim>(1, quadrature.point(p)(0) / 2);
                    break;
                  case 1:
                    q_points[p] =
                      Point<dim>(1, quadrature.point(p)(0) / 2 + 0.5);
                    break;
                  default:
                    Assert(false, ExcInternalError());
                }
              break;
            case 2:
              switch (subface_no)
                {
                  case 0:
                    q_points[p] = Point<dim>(quadrature.point(p)(0) / 2, 0);
                    break;
                  case 1:
                    q_points[p] =
                      Point<dim>(quadrature.point(p)(0) / 2 + 0.5, 0);
                    break;
                  default:
                    Assert(false, ExcInternalError());
                }
              break;
            case 3:
              switch (subface_no)
                {
                  case 0:
                    q_points[p] = Point<dim>(quadrature.point(p)(0) / 2, 1);
                    break;
                  case 1:
                    q_points[p] =
                      Point<dim>(quadrature.point(p)(0) / 2 + 0.5, 1);
                    break;
                  default:
                    Assert(false, ExcInternalError());
                }
              break;

            default:
              Assert(false, ExcInternalError());
          }
    }
  else
    {
      Assert(false, ExcInternalError());
    }
}



template <>
void
QProjector<3>::project_to_subface(const Quadrature<2> &    quadrature,
                                  const unsigned int       face_no,
                                  const unsigned int       subface_no,
                                  std::vector<Point<3>> &  q_points,
                                  const RefinementCase<2> &ref_case)
{
  project_to_subface(ReferenceCells::Hexahedron,
                     quadrature,
                     face_no,
                     subface_no,
                     q_points,
                     ref_case);
}



template <>
void
QProjector<3>::project_to_subface(const ReferenceCell      reference_cell,
                                  const Quadrature<2> &    quadrature,
                                  const unsigned int       face_no,
                                  const unsigned int       subface_no,
                                  std::vector<Point<3>> &  q_points,
                                  const RefinementCase<2> &ref_case)
{
  Assert(reference_cell == ReferenceCells::Hexahedron, ExcNotImplemented());
  (void)reference_cell;

  const unsigned int dim = 3;
  AssertIndexRange(face_no, GeometryInfo<dim>::faces_per_cell);
  AssertIndexRange(subface_no, GeometryInfo<dim>::max_children_per_face);
  Assert(q_points.size() == quadrature.size(),
         ExcDimensionMismatch(q_points.size(), quadrature.size()));

  // one coordinate is at a const value. for
  // faces 0, 2 and 4 this value is 0.0, for
  // faces 1, 3 and 5 it is 1.0
  double const_value = face_no % 2;
  // local 2d coordinates are xi and eta,
  // global 3d coordinates are x, y and
  // z. those have to be mapped. the following
  // indices tell, which global coordinate
  // (0->x, 1->y, 2->z) corresponds to which
  // local one
  unsigned int xi_index    = numbers::invalid_unsigned_int,
               eta_index   = numbers::invalid_unsigned_int,
               const_index = face_no / 2;
  // the xi and eta values have to be scaled
  // (by factor 0.5 or factor 1.0) depending on
  // the refinement case and translated (by 0.0
  // or 0.5) depending on the refinement case
  // and subface_no.
  double xi_scale = 1.0, eta_scale = 1.0, xi_translation = 0.0,
         eta_translation = 0.0;
  // set the index mapping between local and
  // global coordinates
  switch (face_no / 2)
    {
      case 0:
        xi_index  = 1;
        eta_index = 2;
        break;
      case 1:
        xi_index  = 2;
        eta_index = 0;
        break;
      case 2:
        xi_index  = 0;
        eta_index = 1;
        break;
    }
  // set the scale and translation parameter
  // for individual subfaces
  switch (ref_case)
    {
      case RefinementCase<dim - 1>::cut_x:
        xi_scale       = 0.5;
        xi_translation = subface_no % 2 * 0.5;
        break;
      case RefinementCase<dim - 1>::cut_y:
        eta_scale       = 0.5;
        eta_translation = subface_no % 2 * 0.5;
        break;
      case RefinementCase<dim - 1>::cut_xy:
        xi_scale        = 0.5;
        eta_scale       = 0.5;
        xi_translation  = int(subface_no % 2) * 0.5;
        eta_translation = int(subface_no / 2) * 0.5;
        break;
      default:
        Assert(false, ExcInternalError());
        break;
    }
  // finally, compute the scaled, translated,
  // projected quadrature points
  for (unsigned int p = 0; p < quadrature.size(); ++p)
    {
      q_points[p][xi_index] =
        xi_scale * quadrature.point(p)(0) + xi_translation;
      q_points[p][eta_index] =
        eta_scale * quadrature.point(p)(1) + eta_translation;
      q_points[p][const_index] = const_value;
    }
}


template <>
Quadrature<1>
QProjector<1>::project_to_all_faces(const ReferenceCell       reference_cell,
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
      project_to_face(quadrature[quadrature.size() == 1 ? 0 : face],
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

  return Quadrature<dim>(q_points, weights);
}



template <>
Quadrature<2>
QProjector<2>::project_to_all_faces(const ReferenceCell       reference_cell,
                                    const hp::QCollection<1> &quadrature)
{
  if (reference_cell == ReferenceCells::Triangle)
    {
      const auto support_points_line =
        [](const auto &face, const auto &orientation) -> std::vector<Point<2>> {
        std::array<Point<2>, 2> vertices;
        std::copy_n(face.first.begin(), face.first.size(), vertices.begin());
        const auto temp =
          ReferenceCells::Line.permute_according_orientation(vertices,
                                                             orientation);
        return std::vector<Point<2>>(temp.begin(),
                                     temp.begin() + face.first.size());
      };

      // reference faces (defined by its support points and arc length)
      const std::array<std::pair<std::array<Point<2>, 2>, double>, 3> faces = {
        {{{{Point<2>(0.0, 0.0), Point<2>(1.0, 0.0)}}, 1.0},
         {{{Point<2>(1.0, 0.0), Point<2>(0.0, 1.0)}}, std::sqrt(2.0)},
         {{{Point<2>(0.0, 1.0), Point<2>(0.0, 0.0)}}, 1.0}}};

      // linear polynomial to map the reference quadrature points correctly
      // on faces
      const auto poly = BarycentricPolynomials<1>::get_fe_p_basis(1);

      // new (projected) quadrature points and weights
      std::vector<Point<2>> points;
      std::vector<double>   weights;

      // loop over all faces (lines) ...
      for (unsigned int face_no = 0; face_no < faces.size(); ++face_no)
        // ... and over all possible orientations
        for (unsigned int orientation = 0; orientation < 2; ++orientation)
          {
            const auto &face = faces[face_no];

            // determine support point of the current line with the correct
            // orientation
            std::vector<Point<2>> support_points =
              support_points_line(face, orientation);

            // the quadrature rule to be projected ...
            const auto &sub_quadrature_points =
              quadrature[quadrature.size() == 1 ? 0 : face_no].get_points();
            const auto &sub_quadrature_weights =
              quadrature[quadrature.size() == 1 ? 0 : face_no].get_weights();

            // loop over all quadrature points
            for (unsigned int j = 0; j < sub_quadrature_points.size(); ++j)
              {
                Point<2> mapped_point;

                // map reference quadrature point
                for (unsigned int i = 0; i < 2; ++i)
                  mapped_point +=
                    support_points[i] *
                    poly.compute_value(i, sub_quadrature_points[j]);

                points.emplace_back(mapped_point);

                // scale weight by arc length
                weights.emplace_back(sub_quadrature_weights[j] * face.second);
              }
          }

      // construct new quadrature rule
      return {points, weights};
    }

  Assert(reference_cell == ReferenceCells::Quadrilateral, ExcNotImplemented());

  const unsigned int dim = 2;

  const unsigned int n_faces = GeometryInfo<dim>::faces_per_cell;

  unsigned int n_points_total = 0;

  if (quadrature.size() == 1)
    n_points_total = quadrature[0].size() * GeometryInfo<dim>::faces_per_cell;
  else
    {
      AssertDimension(quadrature.size(), GeometryInfo<dim>::faces_per_cell);
      for (unsigned int i = 0; i < quadrature.size(); ++i)
        n_points_total += quadrature[i].size();
    }

  // first fix quadrature points
  std::vector<Point<dim>> q_points;
  q_points.reserve(n_points_total);
  std::vector<Point<dim>> help;
  help.reserve(quadrature.max_n_quadrature_points());

  // project to each face and append
  // results
  for (unsigned int face = 0; face < n_faces; ++face)
    {
      help.resize(quadrature[quadrature.size() == 1 ? 0 : face].size());
      project_to_face(quadrature[quadrature.size() == 1 ? 0 : face],
                      face,
                      help);
      std::copy(help.begin(), help.end(), std::back_inserter(q_points));
    }

  // next copy over weights
  std::vector<double> weights;
  weights.reserve(n_points_total);
  for (unsigned int face = 0; face < n_faces; ++face)
    std::copy(
      quadrature[quadrature.size() == 1 ? 0 : face].get_weights().begin(),
      quadrature[quadrature.size() == 1 ? 0 : face].get_weights().end(),
      std::back_inserter(weights));

  Assert(q_points.size() == n_points_total, ExcInternalError());
  Assert(weights.size() == n_points_total, ExcInternalError());

  return Quadrature<dim>(q_points, weights);
}



template <>
Quadrature<3>
QProjector<3>::project_to_all_faces(const ReferenceCell       reference_cell,
                                    const hp::QCollection<2> &quadrature)
{
  const auto support_points_tri =
    [](const auto &face, const auto &orientation) -> std::vector<Point<3>> {
    std::array<Point<3>, 3> vertices;
    std::copy_n(face.first.begin(), face.first.size(), vertices.begin());
    const auto temp =
      ReferenceCells::Triangle.permute_according_orientation(vertices,
                                                             orientation);
    return std::vector<Point<3>>(temp.begin(),
                                 temp.begin() + face.first.size());
  };

  const auto support_points_quad =
    [](const auto &face, const auto &orientation) -> std::vector<Point<3>> {
    std::array<Point<3>, 4> vertices;
    std::copy_n(face.first.begin(), face.first.size(), vertices.begin());
    const auto temp =
      ReferenceCells::Quadrilateral.permute_according_orientation(vertices,
                                                                  orientation);
    return std::vector<Point<3>>(temp.begin(),
                                 temp.begin() + face.first.size());
  };

  const auto process = [&](const auto &faces) {
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
        // linear polynomial to map the reference quadrature points correctly
        // on faces
        const unsigned int n_shape_functions = faces[face_no].first.size();

        const auto &poly =
          n_shape_functions == 3 ?
            static_cast<const ScalarPolynomialsBase<2> &>(poly_tri) :
            static_cast<const ScalarPolynomialsBase<2> &>(poly_quad);

        // ... and over all possible orientations
        for (unsigned int orientation = 0;
             orientation < (n_shape_functions * 2);
             ++orientation)
          {
            const auto &face = faces[face_no];

            const auto support_points =
              n_shape_functions == 3 ? support_points_tri(face, orientation) :
                                       support_points_quad(face, orientation);

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
                for (unsigned int i = 0; i < n_shape_functions; ++i)
                  mapped_point +=
                    support_points[i] *
                    poly.compute_value(i, sub_quadrature_points[j]);

                points.push_back(mapped_point);

                // scale quadrature weight
                const double scaling = [&]() {
                  const auto &       supp_pts = support_points;
                  const unsigned int dim_     = 2;
                  const unsigned int spacedim = 3;

                  double result[spacedim][dim_];

                  std::vector<Tensor<1, dim_>> shape_derivatives(
                    n_shape_functions);

                  for (unsigned int i = 0; i < n_shape_functions; ++i)
                    shape_derivatives[i] =
                      poly.compute_1st_derivative(i, sub_quadrature_points[j]);

                  for (unsigned int i = 0; i < spacedim; ++i)
                    for (unsigned int j = 0; j < dim_; ++j)
                      result[i][j] = shape_derivatives[0][j] * supp_pts[0][i];
                  for (unsigned int k = 1; k < n_shape_functions; ++k)
                    for (unsigned int i = 0; i < spacedim; ++i)
                      for (unsigned int j = 0; j < dim_; ++j)
                        result[i][j] +=
                          shape_derivatives[k][j] * supp_pts[k][i];

                  DerivativeForm<1, dim_, spacedim> contravariant;

                  for (unsigned int i = 0; i < spacedim; ++i)
                    for (unsigned int j = 0; j < dim_; ++j)
                      contravariant[i][j] = result[i][j];


                  Tensor<1, spacedim> DX_t[dim_];
                  for (unsigned int i = 0; i < spacedim; ++i)
                    for (unsigned int j = 0; j < dim_; ++j)
                      DX_t[j][i] = contravariant[i][j];

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
    return Quadrature<3>(points, weights);
  };

  if (reference_cell == ReferenceCells::Tetrahedron)
    {
      // reference faces (defined by its support points and its area)
      // note: the area is later not used as a scaling factor but recomputed
      const std::vector<std::pair<std::vector<Point<3>>, double>> faces = {
        {{{{Point<3>(0.0, 0.0, 0.0),
            Point<3>(1.0, 0.0, 0.0),
            Point<3>(0.0, 1.0, 0.0)}},
          0.5},
         {{{Point<3>(1.0, 0.0, 0.0),
            Point<3>(0.0, 0.0, 0.0),
            Point<3>(0.0, 0.0, 1.0)}},
          0.5},
         {{{Point<3>(0.0, 0.0, 0.0),
            Point<3>(0.0, 1.0, 0.0),
            Point<3>(0.0, 0.0, 1.0)}},
          0.5},
         {{{Point<3>(0.0, 1.0, 0.0),
            Point<3>(1.0, 0.0, 0.0),
            Point<3>(0.0, 0.0, 1.0)}},
          0.5 * sqrt(3.0) /*equilateral triangle*/}}};

      return process(faces);
    }
  else if (reference_cell == ReferenceCells::Wedge)
    {
      const std::vector<std::pair<std::vector<Point<3>>, double>> faces = {
        {{{{Point<3>(1.0, 0.0, 0.0),
            Point<3>(0.0, 0.0, 0.0),
            Point<3>(0.0, 1.0, 0.0)}},
          0.5},
         {{{Point<3>(0.0, 0.0, 1.0),
            Point<3>(1.0, 0.0, 1.0),
            Point<3>(0.0, 1.0, 1.0)}},
          0.5},
         {{{Point<3>(0.0, 0.0, 0.0),
            Point<3>(1.0, 0.0, 0.0),
            Point<3>(0.0, 0.0, 1.0),
            Point<3>(1.0, 0.0, 1.0)}},
          1.0},
         {{{Point<3>(1.0, 0.0, 0.0),
            Point<3>(0.0, 1.0, 0.0),
            Point<3>(1.0, 0.0, 1.0),
            Point<3>(0.0, 1.0, 1.0)}},
          std::sqrt(2.0)},
         {{{Point<3>(0.0, 1.0, 0.0),
            Point<3>(0.0, 0.0, 0.0),
            Point<3>(0.0, 1.0, 1.0),
            Point<3>(0.0, 0.0, 1.0)}},
          1.0}}};

      return process(faces);
    }
  else if (reference_cell == ReferenceCells::Pyramid)
    {
      const std::vector<std::pair<std::vector<Point<3>>, double>> faces = {
        {{{{Point<3>(-1.0, -1.0, 0.0),
            Point<3>(+1.0, -1.0, 0.0),
            Point<3>(-1.0, +1.0, 0.0),
            Point<3>(+1.0, +1.0, 0.0)}},
          4.0},
         {{{Point<3>(-1.0, -1.0, 0.0),
            Point<3>(-1.0, +1.0, 0.0),
            Point<3>(+0.0, +0.0, 1.0)}},
          std::sqrt(2.0)},
         {{{Point<3>(+1.0, +1.0, 0.0),
            Point<3>(+1.0, -1.0, 0.0),
            Point<3>(+0.0, +0.0, 1.0)}},
          std::sqrt(2.0)},
         {{{Point<3>(+1.0, -1.0, 0.0),
            Point<3>(-1.0, -1.0, 0.0),
            Point<3>(+0.0, +0.0, 1.0)}},
          std::sqrt(2.0)},
         {{{Point<3>(-1.0, +1.0, 0.0),
            Point<3>(+1.0, +1.0, 0.0),
            Point<3>(+0.0, +0.0, 1.0)}},
          std::sqrt(2.0)}}};

      return process(faces);
    }


  Assert(reference_cell == ReferenceCells::Hexahedron, ExcNotImplemented());

  const unsigned int dim = 3;

  unsigned int n_points_total = 0;

  if (quadrature.size() == 1)
    n_points_total = quadrature[0].size() * GeometryInfo<dim>::faces_per_cell;
  else
    {
      AssertDimension(quadrature.size(), GeometryInfo<dim>::faces_per_cell);
      for (unsigned int i = 0; i < quadrature.size(); ++i)
        n_points_total += quadrature[i].size();
    }

  n_points_total *= 8;

  // first fix quadrature points
  std::vector<Point<dim>> q_points;
  q_points.reserve(n_points_total);
  std::vector<Point<dim>> help;
  help.reserve(quadrature.max_n_quadrature_points());

  std::vector<double> weights;
  weights.reserve(n_points_total);

  // do the following for all possible
  // mutations of a face (mutation==0
  // corresponds to a face with standard
  // orientation, no flip and no rotation)
  for (unsigned int i = 0; i < 8; ++i)
    {
      // project to each face and append
      // results
      for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell;
           ++face)
        {
          SubQuadrature mutation;

          const auto quadrature_f =
            quadrature[quadrature.size() == 1 ? 0 : face];
          switch (i)
            {
              case 0:
                mutation = quadrature_f;
                break;
              case 1:
                mutation = internal::QProjector::rotate(quadrature_f, 1);
                break;
              case 2:
                mutation = internal::QProjector::rotate(quadrature_f, 2);
                break;
              case 3:
                mutation = internal::QProjector::rotate(quadrature_f, 3);
                break;
              case 4:
                mutation = internal::QProjector::reflect(quadrature_f);
                break;
              case 5:
                mutation = internal::QProjector::rotate(
                  internal::QProjector::reflect(quadrature_f), 3);
                break;
              case 6:
                mutation = internal::QProjector::rotate(
                  internal::QProjector::reflect(quadrature_f), 2);
                break;
              case 7:
                mutation = internal::QProjector::rotate(
                  internal::QProjector::reflect(quadrature_f), 1);
                break;
              default:
                Assert(false, ExcInternalError())
            }

          help.resize(quadrature[quadrature.size() == 1 ? 0 : face].size());
          project_to_face(mutation, face, help);
          std::copy(help.begin(), help.end(), std::back_inserter(q_points));

          std::copy(mutation.get_weights().begin(),
                    mutation.get_weights().end(),
                    std::back_inserter(weights));
        }
    }


  Assert(q_points.size() == n_points_total, ExcInternalError());
  Assert(weights.size() == n_points_total, ExcInternalError());

  return Quadrature<dim>(q_points, weights);
}



template <>
Quadrature<1>
QProjector<1>::project_to_all_subfaces(const Quadrature<0> &quadrature)
{
  return project_to_all_subfaces(ReferenceCells::Line, quadrature);
}



template <>
Quadrature<1>
QProjector<1>::project_to_all_subfaces(const ReferenceCell  reference_cell,
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
        project_to_subface(quadrature, face, subface, help);
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

  return Quadrature<dim>(q_points, weights);
}



template <>
Quadrature<2>
QProjector<2>::project_to_all_subfaces(const ReferenceCell  reference_cell,
                                       const SubQuadrature &quadrature)
{
  if (reference_cell == ReferenceCells::Triangle ||
      reference_cell == ReferenceCells::Tetrahedron)
    return Quadrature<2>(); // nothing to do

  Assert(reference_cell == ReferenceCells::Quadrilateral, ExcNotImplemented());

  const unsigned int dim = 2;

  const unsigned int n_points = quadrature.size(),
                     n_faces  = GeometryInfo<dim>::faces_per_cell,
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
        project_to_subface(quadrature, face, subface, help);
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

  return Quadrature<dim>(q_points, weights);
}



template <>
Quadrature<2>
QProjector<2>::project_to_all_subfaces(const SubQuadrature &quadrature)
{
  return project_to_all_subfaces(ReferenceCells::Quadrilateral, quadrature);
}



template <>
Quadrature<3>
QProjector<3>::project_to_all_subfaces(const ReferenceCell  reference_cell,
                                       const SubQuadrature &quadrature)
{
  if (reference_cell == ReferenceCells::Triangle ||
      reference_cell == ReferenceCells::Tetrahedron)
    return Quadrature<3>(); // nothing to do

  Assert(reference_cell == ReferenceCells::Hexahedron, ExcNotImplemented());

  const unsigned int dim         = 3;
  SubQuadrature      q_reflected = internal::QProjector::reflect(quadrature);
  SubQuadrature      q[8]        = {quadrature,
                        internal::QProjector::rotate(quadrature, 1),
                        internal::QProjector::rotate(quadrature, 2),
                        internal::QProjector::rotate(quadrature, 3),
                        q_reflected,
                        internal::QProjector::rotate(q_reflected, 3),
                        internal::QProjector::rotate(q_reflected, 2),
                        internal::QProjector::rotate(q_reflected, 1)};

  const unsigned int n_points = quadrature.size(),
                     n_faces  = GeometryInfo<dim>::faces_per_cell,
                     total_subfaces_per_face = 2 + 2 + 4;

  // first fix quadrature points
  std::vector<Point<dim>> q_points;
  q_points.reserve(n_points * n_faces * total_subfaces_per_face * 8);
  std::vector<Point<dim>> help(n_points);

  std::vector<double> weights;
  weights.reserve(n_points * n_faces * total_subfaces_per_face * 8);

  // do the following for all possible
  // mutations of a face (mutation==0
  // corresponds to a face with standard
  // orientation, no flip and no rotation)
  for (const auto &mutation : q)
    {
      // project to each face and copy
      // results
      for (unsigned int face = 0; face < n_faces; ++face)
        for (unsigned int ref_case = RefinementCase<dim - 1>::cut_xy;
             ref_case >= RefinementCase<dim - 1>::cut_x;
             --ref_case)
          for (unsigned int subface = 0;
               subface < GeometryInfo<dim - 1>::n_children(
                           RefinementCase<dim - 1>(ref_case));
               ++subface)
            {
              project_to_subface(mutation,
                                 face,
                                 subface,
                                 help,
                                 RefinementCase<dim - 1>(ref_case));
              std::copy(help.begin(), help.end(), std::back_inserter(q_points));
            }

      // next copy over weights
      for (unsigned int face = 0; face < n_faces; ++face)
        for (unsigned int ref_case = RefinementCase<dim - 1>::cut_xy;
             ref_case >= RefinementCase<dim - 1>::cut_x;
             --ref_case)
          for (unsigned int subface = 0;
               subface < GeometryInfo<dim - 1>::n_children(
                           RefinementCase<dim - 1>(ref_case));
               ++subface)
            std::copy(mutation.get_weights().begin(),
                      mutation.get_weights().end(),
                      std::back_inserter(weights));
    }

  Assert(q_points.size() == n_points * n_faces * total_subfaces_per_face * 8,
         ExcInternalError());
  Assert(weights.size() == n_points * n_faces * total_subfaces_per_face * 8,
         ExcInternalError());

  return Quadrature<dim>(q_points, weights);
}



template <>
Quadrature<3>
QProjector<3>::project_to_all_subfaces(const SubQuadrature &quadrature)
{
  return project_to_all_subfaces(ReferenceCells::Hexahedron, quadrature);
}



// This function is not used in the library
template <int dim>
Quadrature<dim>
QProjector<dim>::project_to_child(const Quadrature<dim> &quadrature,
                                  const unsigned int     child_no)
{
  return project_to_child(ReferenceCells::get_hypercube<dim>(),
                          quadrature,
                          child_no);
}



template <int dim>
Quadrature<dim>
QProjector<dim>::project_to_child(const ReferenceCell    reference_cell,
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
QProjector<dim>::project_to_all_children(const Quadrature<dim> &quadrature)
{
  return project_to_all_children(ReferenceCells::get_hypercube<dim>(),
                                 quadrature);
}



template <int dim>
Quadrature<dim>
QProjector<dim>::project_to_all_children(const ReferenceCell    reference_cell,
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
      Quadrature<dim> help = project_to_child(quadrature, child);
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
QProjector<dim>::project_to_line(const Quadrature<1> &quadrature,
                                 const Point<dim> &   p1,
                                 const Point<dim> &   p2)
{
  return project_to_line(ReferenceCells::get_hypercube<dim>(),
                         quadrature,
                         p1,
                         p2);
}



template <int dim>
Quadrature<dim>
QProjector<dim>::project_to_line(const ReferenceCell  reference_cell,
                                 const Quadrature<1> &quadrature,
                                 const Point<dim> &   p1,
                                 const Point<dim> &   p2)
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
      const double alpha = quadrature.point(k)(0);
      points[k]          = alpha * p2;
      points[k] += (1. - alpha) * p1;
      weights[k] = length * quadrature.weight(k);
    }
  return Quadrature<dim>(points, weights);
}



template <int dim>
typename QProjector<dim>::DataSetDescriptor
QProjector<dim>::DataSetDescriptor::face(const unsigned int face_no,
                                         const bool         face_orientation,
                                         const bool         face_flip,
                                         const bool         face_rotation,
                                         const unsigned int n_quadrature_points)
{
  return face(ReferenceCells::get_hypercube<dim>(),
              face_no,
              face_orientation,
              face_flip,
              face_rotation,
              n_quadrature_points);
}



template <int dim>
typename QProjector<dim>::DataSetDescriptor
QProjector<dim>::DataSetDescriptor::face(const ReferenceCell reference_cell,
                                         const unsigned int  face_no,
                                         const bool          face_orientation,
                                         const bool          face_flip,
                                         const bool          face_rotation,
                                         const unsigned int n_quadrature_points)
{
  if (reference_cell == ReferenceCells::Triangle ||
      reference_cell == ReferenceCells::Tetrahedron)
    {
      if (dim == 2)
        return {(2 * face_no + (face_orientation ? 1 : 0)) *
                n_quadrature_points};
      else if (dim == 3)
        {
          const unsigned int orientation = (face_flip ? 4 : 0) +
                                           (face_rotation ? 2 : 0) +
                                           (face_orientation ? 1 : 0);
          return {(6 * face_no + orientation) * n_quadrature_points};
        }
    }

  Assert(reference_cell == ReferenceCells::get_hypercube<dim>(),
         ExcNotImplemented());

  Assert(face_no < GeometryInfo<dim>::faces_per_cell, ExcInternalError());

  switch (dim)
    {
      case 1:
      case 2:
        return face_no * n_quadrature_points;


      case 3:
        {
          // in 3d, we have to account for faces that
          // have non-standard face orientation, flip
          // and rotation. thus, we have to store
          // _eight_ data sets per face or subface

          // set up a table with the according offsets
          // for non-standard orientation, first index:
          // face_orientation (standard true=1), second
          // index: face_flip (standard false=0), third
          // index: face_rotation (standard false=0)
          //
          // note, that normally we should use the
          // obvious offsets 0,1,2,3,4,5,6,7. However,
          // prior to the changes enabling flipped and
          // rotated faces, in many places of the
          // library the convention was used, that the
          // first dataset with offset 0 corresponds to
          // a face in standard orientation. therefore
          // we use the offsets 4,5,6,7,0,1,2,3 here to
          // stick to that (implicit) convention
          static const unsigned int offset[2][2][2] = {
            {{4 * GeometryInfo<dim>::faces_per_cell,
              5 * GeometryInfo<dim>::
                    faces_per_cell}, // face_orientation=false; face_flip=false;
                                     // face_rotation=false and true
             {6 * GeometryInfo<dim>::faces_per_cell,
              7 * GeometryInfo<dim>::
                    faces_per_cell}}, // face_orientation=false; face_flip=true;
                                      // face_rotation=false and true
            {{0 * GeometryInfo<dim>::faces_per_cell,
              1 * GeometryInfo<dim>::
                    faces_per_cell}, // face_orientation=true;  face_flip=false;
                                     // face_rotation=false and true
             {2 * GeometryInfo<dim>::faces_per_cell,
              3 * GeometryInfo<dim>::
                    faces_per_cell}}}; // face_orientation=true; face_flip=true;
                                       // face_rotation=false and true

          return (
            (face_no + offset[face_orientation][face_flip][face_rotation]) *
            n_quadrature_points);
        }

      default:
        Assert(false, ExcInternalError());
    }
  return numbers::invalid_unsigned_int;
}



template <int dim>
typename QProjector<dim>::DataSetDescriptor
QProjector<dim>::DataSetDescriptor::face(
  const ReferenceCell             reference_cell,
  const unsigned int              face_no,
  const bool                      face_orientation,
  const bool                      face_flip,
  const bool                      face_rotation,
  const hp::QCollection<dim - 1> &quadrature)
{
  if (reference_cell == ReferenceCells::Triangle ||
      reference_cell == ReferenceCells::Tetrahedron ||
      reference_cell == ReferenceCells::Wedge ||
      reference_cell == ReferenceCells::Pyramid)
    {
      unsigned int offset = 0;

      static const unsigned int X = numbers::invalid_unsigned_int;
      static const std::array<unsigned int, 5> scale_tri   = {{2, 2, 2, X, X}};
      static const std::array<unsigned int, 5> scale_tet   = {{6, 6, 6, 6, X}};
      static const std::array<unsigned int, 5> scale_wedge = {{6, 6, 8, 8, 8}};
      static const std::array<unsigned int, 5> scale_pyramid = {
        {8, 6, 6, 6, 6}};

      const auto &scale =
        (reference_cell == ReferenceCells::Triangle) ?
          scale_tri :
          ((reference_cell == ReferenceCells::Tetrahedron) ?
             scale_tet :
             ((reference_cell == ReferenceCells::Wedge) ? scale_wedge :
                                                          scale_pyramid));

      if (quadrature.size() == 1)
        offset = scale[0] * quadrature[0].size() * face_no;
      else
        for (unsigned int i = 0; i < face_no; ++i)
          offset += scale[i] * quadrature[i].size();

      if (dim == 2)
        return {offset +
                face_orientation *
                  quadrature[quadrature.size() == 1 ? 0 : face_no].size()};
      else if (dim == 3)
        {
          const unsigned int orientation = (face_flip ? 4 : 0) +
                                           (face_rotation ? 2 : 0) +
                                           (face_orientation ? 1 : 0);

          return {offset +
                  orientation *
                    quadrature[quadrature.size() == 1 ? 0 : face_no].size()};
        }
    }

  Assert(reference_cell == ReferenceCells::get_hypercube<dim>(),
         ExcNotImplemented());

  Assert(face_no < GeometryInfo<dim>::faces_per_cell, ExcInternalError());

  switch (dim)
    {
      case 1:
      case 2:
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
          // in 3d, we have to account for faces that
          // have non-standard face orientation, flip
          // and rotation. thus, we have to store
          // _eight_ data sets per face or subface

          // set up a table with the according offsets
          // for non-standard orientation, first index:
          // face_orientation (standard true=1), second
          // index: face_flip (standard false=0), third
          // index: face_rotation (standard false=0)
          //
          // note, that normally we should use the
          // obvious offsets 0,1,2,3,4,5,6,7. However,
          // prior to the changes enabling flipped and
          // rotated faces, in many places of the
          // library the convention was used, that the
          // first dataset with offset 0 corresponds to
          // a face in standard orientation. therefore
          // we use the offsets 4,5,6,7,0,1,2,3 here to
          // stick to that (implicit) convention
          static const unsigned int offset[2][2][2] = {
            {{4, 5},   // face_orientation=false; face_flip=false;
                       // face_rotation=false and true
             {6, 7}},  // face_orientation=false; face_flip=true;
                       // face_rotation=false and true
            {{0, 1},   // face_orientation=true;  face_flip=false;
                       // face_rotation=false and true
             {2, 3}}}; // face_orientation=true; face_flip=true;
                       // face_rotation=false and true


          if (quadrature.size() == 1)
            return (face_no +
                    offset[face_orientation][face_flip][face_rotation] *
                      GeometryInfo<dim>::faces_per_cell) *
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

              return (n_points_i +
                      offset[face_orientation][face_flip][face_rotation] *
                        n_points);
            }
        }

      default:
        Assert(false, ExcInternalError());
    }
  return numbers::invalid_unsigned_int;
}



template <>
QProjector<1>::DataSetDescriptor
QProjector<1>::DataSetDescriptor::subface(
  const ReferenceCell reference_cell,
  const unsigned int  face_no,
  const unsigned int  subface_no,
  const bool,
  const bool,
  const bool,
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
QProjector<1>::DataSetDescriptor
QProjector<1>::DataSetDescriptor::subface(
  const unsigned int             face_no,
  const unsigned int             subface_no,
  const bool                     face_orientation,
  const bool                     face_flip,
  const bool                     face_rotation,
  const unsigned int             n_quadrature_points,
  const internal::SubfaceCase<1> ref_case)
{
  return subface(ReferenceCells::Line,
                 face_no,
                 subface_no,
                 face_orientation,
                 face_flip,
                 face_rotation,
                 n_quadrature_points,
                 ref_case);
}



template <>
QProjector<2>::DataSetDescriptor
QProjector<2>::DataSetDescriptor::subface(
  const ReferenceCell reference_cell,
  const unsigned int  face_no,
  const unsigned int  subface_no,
  const bool,
  const bool,
  const bool,
  const unsigned int n_quadrature_points,
  const internal::SubfaceCase<2>)
{
  Assert(reference_cell == ReferenceCells::Quadrilateral, ExcNotImplemented());
  (void)reference_cell;

  Assert(face_no < GeometryInfo<2>::faces_per_cell, ExcInternalError());
  Assert(subface_no < GeometryInfo<2>::max_children_per_face,
         ExcInternalError());

  return ((face_no * GeometryInfo<2>::max_children_per_face + subface_no) *
          n_quadrature_points);
}



template <>
QProjector<2>::DataSetDescriptor
QProjector<2>::DataSetDescriptor::subface(
  const unsigned int             face_no,
  const unsigned int             subface_no,
  const bool                     face_orientation,
  const bool                     face_flip,
  const bool                     face_rotation,
  const unsigned int             n_quadrature_points,
  const internal::SubfaceCase<2> ref_case)
{
  return subface(ReferenceCells::Quadrilateral,
                 face_no,
                 subface_no,
                 face_orientation,
                 face_flip,
                 face_rotation,
                 n_quadrature_points,
                 ref_case);
}


template <>
QProjector<3>::DataSetDescriptor
QProjector<3>::DataSetDescriptor::subface(
  const ReferenceCell            reference_cell,
  const unsigned int             face_no,
  const unsigned int             subface_no,
  const bool                     face_orientation,
  const bool                     face_flip,
  const bool                     face_rotation,
  const unsigned int             n_quadrature_points,
  const internal::SubfaceCase<3> ref_case)
{
  const unsigned int dim = 3;

  Assert(reference_cell == ReferenceCells::Hexahedron, ExcNotImplemented());
  (void)reference_cell;

  Assert(face_no < GeometryInfo<dim>::faces_per_cell, ExcInternalError());
  Assert(subface_no < GeometryInfo<dim>::max_children_per_face,
         ExcInternalError());

  // As the quadrature points created by
  // QProjector are on subfaces in their
  // "standard location" we have to use a
  // permutation of the equivalent subface
  // number in order to respect face
  // orientation, flip and rotation. The
  // information we need here is exactly the
  // same as the
  // GeometryInfo<3>::child_cell_on_face info
  // for the bottom face (face 4) of a hex, as
  // on this the RefineCase of the cell matches
  // that of the face and the subfaces are
  // numbered in the same way as the child
  // cells.

  // in 3d, we have to account for faces that
  // have non-standard face orientation, flip
  // and rotation. thus, we have to store
  // _eight_ data sets per face or subface
  // already for the isotropic
  // case. Additionally, we have three
  // different refinement cases, resulting in
  // <tt>4 + 2 + 2 = 8</tt> different subfaces
  // for each face.
  const unsigned int total_subfaces_per_face = 8;

  // set up a table with the according offsets
  // for non-standard orientation, first index:
  // face_orientation (standard true=1), second
  // index: face_flip (standard false=0), third
  // index: face_rotation (standard false=0)
  //
  // note, that normally we should use the
  // obvious offsets 0,1,2,3,4,5,6,7. However,
  // prior to the changes enabling flipped and
  // rotated faces, in many places of the
  // library the convention was used, that the
  // first dataset with offset 0 corresponds to
  // a face in standard orientation. therefore
  // we use the offsets 4,5,6,7,0,1,2,3 here to
  // stick to that (implicit) convention
  static const unsigned int orientation_offset[2][2][2] = {
    {// face_orientation=false; face_flip=false; face_rotation=false and true
     {4 * GeometryInfo<dim>::faces_per_cell * total_subfaces_per_face,
      5 * GeometryInfo<dim>::faces_per_cell * total_subfaces_per_face},
     // face_orientation=false; face_flip=true;  face_rotation=false and true
     {6 * GeometryInfo<dim>::faces_per_cell * total_subfaces_per_face,
      7 * GeometryInfo<dim>::faces_per_cell * total_subfaces_per_face}},
    {// face_orientation=true;  face_flip=false; face_rotation=false and true
     {0 * GeometryInfo<dim>::faces_per_cell * total_subfaces_per_face,
      1 * GeometryInfo<dim>::faces_per_cell * total_subfaces_per_face},
     // face_orientation=true;  face_flip=true;  face_rotation=false and true
     {2 * GeometryInfo<dim>::faces_per_cell * total_subfaces_per_face,
      3 * GeometryInfo<dim>::faces_per_cell * total_subfaces_per_face}}};

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
    equivalent_subface_number[internal::SubfaceCase<dim>::case_isotropic + 1]
                             [GeometryInfo<3>::max_children_per_face] = {
                               // case_none, see above
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
    (face_orientation == face_rotation ? ref_case_permutation[equ_ref_case] :
                                         equ_ref_case);

  // what we have now is the number of
  // the subface in the natural
  // orientation of the *face*. what we
  // need to know is the number of the
  // subface concerning the standard face
  // orientation as seen from the *cell*.

  // this mapping is not trivial, but we
  // have done exactly this stuff in the
  // child_cell_on_face function. in
  // order to reduce the amount of code
  // as well as to make maintaining the
  // functionality easier we want to
  // reuse that information. So we note
  // that on the bottom face (face 4) of
  // a hex cell the local x and y
  // coordinates of the face and the cell
  // coincide, thus also the refinement
  // case of the face corresponds to the
  // refinement case of the cell
  // (ignoring cell refinement along the
  // z direction). Using this knowledge
  // we can (ab)use the
  // child_cell_on_face function to do
  // exactly the transformation we are in
  // need of now
  const unsigned int final_subface_no =
    GeometryInfo<dim>::child_cell_on_face(RefinementCase<dim>(final_ref_case),
                                          4,
                                          equ_subface_no,
                                          face_orientation,
                                          face_flip,
                                          face_rotation,
                                          equ_ref_case);

  return (((face_no * total_subfaces_per_face +
            ref_case_offset[final_ref_case - 1] + final_subface_no) +
           orientation_offset[face_orientation][face_flip][face_rotation]) *
          n_quadrature_points);
}


template <>
QProjector<3>::DataSetDescriptor
QProjector<3>::DataSetDescriptor::subface(
  const unsigned int             face_no,
  const unsigned int             subface_no,
  const bool                     face_orientation,
  const bool                     face_flip,
  const bool                     face_rotation,
  const unsigned int             n_quadrature_points,
  const internal::SubfaceCase<3> ref_case)
{
  return subface(ReferenceCells::Hexahedron,
                 face_no,
                 subface_no,
                 face_orientation,
                 face_flip,
                 face_rotation,
                 n_quadrature_points,
                 ref_case);
}



template <int dim>
Quadrature<dim>
QProjector<dim>::project_to_face(const SubQuadrature &quadrature,
                                 const unsigned int   face_no)
{
  return project_to_face(ReferenceCells::get_hypercube<dim>(),
                         quadrature,
                         face_no);
}



template <int dim>
Quadrature<dim>
QProjector<dim>::project_to_face(const ReferenceCell  reference_cell,
                                 const SubQuadrature &quadrature,
                                 const unsigned int   face_no)
{
  Assert(reference_cell == ReferenceCells::get_hypercube<dim>(),
         ExcNotImplemented());
  (void)reference_cell;

  std::vector<Point<dim>> points(quadrature.size());
  project_to_face(quadrature, face_no, points);
  return Quadrature<dim>(points, quadrature.get_weights());
}



template <int dim>
Quadrature<dim>
QProjector<dim>::project_to_subface(const SubQuadrature &          quadrature,
                                    const unsigned int             face_no,
                                    const unsigned int             subface_no,
                                    const RefinementCase<dim - 1> &ref_case)
{
  return project_to_subface(ReferenceCells::get_hypercube<dim>(),
                            quadrature,
                            face_no,
                            subface_no,
                            ref_case);
}



template <int dim>
Quadrature<dim>
QProjector<dim>::project_to_subface(const ReferenceCell  reference_cell,
                                    const SubQuadrature &quadrature,
                                    const unsigned int   face_no,
                                    const unsigned int   subface_no,
                                    const RefinementCase<dim - 1> &ref_case)
{
  Assert(reference_cell == ReferenceCells::get_hypercube<dim>(),
         ExcNotImplemented());
  (void)reference_cell;

  std::vector<Point<dim>> points(quadrature.size());
  project_to_subface(quadrature, face_no, subface_no, points, ref_case);
  return Quadrature<dim>(points, quadrature.get_weights());
}


// explicit instantiations; note: we need them all for all dimensions
template class QProjector<1>;
template class QProjector<2>;
template class QProjector<3>;

DEAL_II_NAMESPACE_CLOSE
