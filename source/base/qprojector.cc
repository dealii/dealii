// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2020 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#include <deal.II/base/derivative_form.h>
#include <deal.II/base/geometry_info.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/base/std_cxx26/inplace_vector.h>

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
        const ReferenceCell<dim - 1> &face_reference_cell,
        const Quadrature<dim - 1>    &quadrature,
        const std_cxx26::inplace_vector<
          Point<dim>,
          ReferenceCells::max_n_vertices<dim - 1>()> &vertices,
        const double                                  measure,
        const types::geometric_orientation            combined_orientation,
        std::vector<Point<dim>>                      &points,
        std::vector<double>                          &weights)
      {
        AssertDimension(points.size(), weights.size());
        points.reserve(points.size() + quadrature.size());
        weights.reserve(weights.size() + quadrature.size());

        const auto support_points =
          face_reference_cell.permute_by_combined_orientation(
            ArrayView<const Point<dim>>(vertices),
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



template <int dim>
Quadrature<dim>
QProjector<dim>::project_to_face(
  const ReferenceCell<dim>          &reference_cell,
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
  std_cxx26::inplace_vector<Point<dim>,
                            ReferenceCells::max_n_vertices<dim - 1>()>
    face_vertices(face_reference_cell.n_vertices());
  for (const unsigned int vertex_no : face_reference_cell.vertex_indices())
    face_vertices[vertex_no] =
      reference_cell.face_vertex_location(face_no, vertex_no);
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



template <int dim>
Quadrature<dim>
QProjector<dim>::project_to_subface(
  const ReferenceCell<dim>          &reference_cell,
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
  if (dim > 1)
    AssertIndexRange(subface_no,
                     reference_cell.face_reference_cell(face_no).n_children(
                       ref_case));
  if (dim == 1)
    AssertDimension(quadrature.size(), 1);

  std::vector<Point<dim>> points;
  std::vector<double>     weights;
  std_cxx26::inplace_vector<Point<dim>,
                            ReferenceCells::max_n_vertices<dim - 1>()>
    vertices;
  for (const unsigned int subface_vertex_no :
       reference_cell.face_reference_cell(face_no).vertex_indices())
    vertices.push_back(reference_cell.subface_vertex_location(
      face_no, subface_no, subface_vertex_no, ref_case));
  internal::QProjector::append_subobject_rule(
    reference_cell.face_reference_cell(face_no),
    quadrature,
    vertices,
    reference_cell.face_measure(face_no),
    combined_orientation,
    points,
    weights);

  return Quadrature<dim>(std::move(points), std::move(weights));
}



template <int dim>
Quadrature<dim>
QProjector<dim>::project_to_all_faces(
  const ReferenceCell<dim>       &reference_cell,
  const hp::QCollection<dim - 1> &quadrature)
{
  std::size_t n_points = 0;
  for (const unsigned int face_no : reference_cell.face_indices())
    n_points += quadrature[quadrature.size() == 1 ? 0 : face_no].size() *
                reference_cell.n_face_orientations(face_no);

  std::vector<Point<dim>> points;
  std::vector<double>     weights;
  points.reserve(n_points);
  weights.reserve(n_points);

  for (const unsigned int face_no : reference_cell.face_indices())
    {
      const ReferenceCell face_reference_cell =
        reference_cell.face_reference_cell(face_no);
      std_cxx26::inplace_vector<Point<dim>,
                                ReferenceCells::max_n_vertices<dim - 1>()>
        face_vertices(face_reference_cell.n_vertices());
      for (const unsigned int vertex_no : face_reference_cell.vertex_indices())
        face_vertices[vertex_no] =
          reference_cell.face_vertex_location(face_no, vertex_no);

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



template <int dim>
Quadrature<dim>
QProjector<dim>::project_to_all_subfaces(
  const ReferenceCell<dim>  &reference_cell,
  const Quadrature<dim - 1> &quadrature)
{
  AssertDimension(reference_cell.get_dimension(), dim);
  if (dim == 1)
    AssertDimension(quadrature.size(), 1);

  std::size_t n_points = 0;
  for (const unsigned int face_no : reference_cell.face_indices())
    {
      const auto n_orientations = reference_cell.n_face_orientations(face_no);
      const auto face_reference_cell =
        reference_cell.face_reference_cell(face_no);
      for (const auto &refinement_case : face_reference_cell.refinement_cases())
        {
          const auto n_children =
            dim > 1 ? face_reference_cell.n_children(refinement_case) : 1;
          n_points += n_children * n_orientations * quadrature.size();
        }
    }

  std::vector<Point<dim>> points;
  std::vector<double>     weights;
  points.reserve(n_points);
  weights.reserve(n_points);

  // project to each face and copy results
  for (unsigned int face_no = 0; face_no < reference_cell.n_faces(); ++face_no)
    {
      const auto face_reference_cell =
        reference_cell.face_reference_cell(face_no);
      const auto &refinement_cases = face_reference_cell.refinement_cases();
      for (const auto &refinement_case : refinement_cases)
        {
          const auto n_children =
            dim > 1 ? face_reference_cell.n_children(refinement_case) : 1;
          for (types::geometric_orientation combined_orientation = 0;
               combined_orientation <
               reference_cell.n_face_orientations(face_no);
               ++combined_orientation)
            for (unsigned int subface_no = 0; subface_no < n_children;
                 ++subface_no)
              {
                std_cxx26::inplace_vector<
                  Point<dim>,
                  ReferenceCells::max_n_vertices<dim - 1>()>
                  vertices;
                for (const unsigned int subface_vertex_no :
                     reference_cell.face_reference_cell(face_no)
                       .vertex_indices())
                  vertices.push_back(reference_cell.subface_vertex_location(
                    face_no, subface_no, subface_vertex_no, refinement_case));
                internal::QProjector::append_subobject_rule(
                  reference_cell.face_reference_cell(face_no),
                  quadrature,
                  vertices,
                  reference_cell.face_measure(face_no),
                  combined_orientation,
                  points,
                  weights);
              }
        }
    }

  return Quadrature<dim>(std::move(points), std::move(weights));
}



template <int dim>
Quadrature<dim>
QProjector<dim>::project_to_child(const ReferenceCell<dim> &reference_cell,
                                  const Quadrature<dim>    &quadrature,
                                  const unsigned int        child_no)
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
QProjector<dim>::project_to_all_children(
  const ReferenceCell<dim> &reference_cell,
  const Quadrature<dim>    &quadrature)
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
QProjector<dim>::project_to_line(const ReferenceCell<dim> &reference_cell,
                                 const Quadrature<1>      &quadrature,
                                 const Point<dim>         &p1,
                                 const Point<dim>         &p2)
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
QProjector<dim>::DataSetDescriptor::face(
  const ReferenceCell<dim>          &reference_cell,
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
  const ReferenceCell<dim>          &reference_cell,
  const unsigned int                 face_no,
  const types::geometric_orientation combined_orientation,
  const hp::QCollection<dim - 1>    &quadrature)
{
  AssertIndexRange(face_no, reference_cell.n_faces());
  AssertIndexRange(combined_orientation,
                   reference_cell.n_face_orientations(face_no));
  AssertDimension(reference_cell.get_dimension(), dim);

  unsigned int offset = 0;
  for (unsigned int i = 0; i < face_no; ++i)
    offset += reference_cell.n_face_orientations(i) *
              quadrature[quadrature.size() == 1 ? 0 : i].size();

  return {offset + combined_orientation *
                     quadrature[quadrature.size() == 1 ? 0 : face_no].size()};
}



template <int dim>
typename QProjector<dim>::DataSetDescriptor
QProjector<dim>::DataSetDescriptor::subface(
  const ReferenceCell<dim>          &reference_cell,
  const unsigned int                 face_no,
  const unsigned int                 subface_no,
  const types::geometric_orientation combined_orientation,
  const unsigned int                 n_quadrature_points,
  const internal::SubfaceCase<dim>   ref_case)
{
  AssertDimension(reference_cell.get_dimension(), dim);
  AssertIndexRange(face_no, reference_cell.n_faces());
  AssertIndexRange(combined_orientation,
                   reference_cell.n_face_orientations(face_no));

  const auto [final_subface_no, final_refinement_case] =
    reference_cell.equivalent_refinement_case(combined_orientation,
                                              ref_case,
                                              subface_no);
  if (dim > 1)
    AssertIndexRange(final_subface_no,
                     reference_cell.face_reference_cell(face_no).n_children(
                       final_refinement_case));

  // This function can't work with mixed elements since, in general, those may
  // have a different number of quadrature points per face
  Assert(reference_cell != ReferenceCells::Pyramid &&
           reference_cell != ReferenceCells::Wedge,
         ExcNotImplemented());

  // Calculate the total number of points per face:
  const auto  face_reference_cell = reference_cell.face_reference_cell(face_no);
  const auto &refinement_cases    = face_reference_cell.refinement_cases();
  unsigned int points_per_face    = 0;
  for (const auto &refinement_case : refinement_cases)
    {
      const auto n_children =
        dim > 1 ? face_reference_cell.n_children(refinement_case) : 1;
      points_per_face += reference_cell.n_face_orientations(face_no) *
                         n_children * n_quadrature_points;
    }

  // Next, calculate where we are in the current face's enumeration of
  // quadrature rules:
  unsigned int index = points_per_face * face_no;
  for (const auto &refinement_case : refinement_cases)
    {
      const auto n_children =
        dim > 1 ? face_reference_cell.n_children(refinement_case) : 1;

      if (refinement_case == final_refinement_case)
        return index + (combined_orientation * n_children + final_subface_no) *
                         n_quadrature_points;
      else
        index += reference_cell.n_face_orientations(face_no) * n_children *
                 n_quadrature_points;
    }
  DEAL_II_ASSERT_UNREACHABLE();

  return index;
}


// explicit instantiations; note: we need them all for all dimensions
template class QProjector<1>;
template class QProjector<2>;
template class QProjector<3>;

DEAL_II_NAMESPACE_CLOSE
