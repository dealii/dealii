// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#include <deal.II/base/tensor.h>
#include <deal.II/base/table.h>
#include <deal.II/grid/manifold.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/fe/fe_q.h>

#include <cmath>

#include <boost/container/small_vector.hpp>

DEAL_II_NAMESPACE_OPEN

using namespace Manifolds;

/* -------------------------- Manifold --------------------- */
template <int dim, int spacedim>
Point<spacedim>
Manifold<dim, spacedim>::
project_to_manifold (const ArrayView<const Point<spacedim>> &,
                     const Point<spacedim> &) const
{
  Assert (false, ExcPureFunctionCalled());
  return Point<spacedim>();
}



template <int dim, int spacedim>
Point<spacedim>
Manifold<dim, spacedim>::
get_intermediate_point (const Point<spacedim> &p1,
                        const Point<spacedim> &p2,
                        const double w) const
{
  const std::array<Point<spacedim>, 2> vertices {{p1, p2}};
  return project_to_manifold(make_array_view(vertices.begin(), vertices.end()),
                             w * p2 + (1-w)*p1);
}



template <int dim, int spacedim>
Point<spacedim>
Manifold<dim, spacedim>::
get_new_point (const ArrayView<const Point<spacedim>> &surrounding_points,
               const ArrayView<const double>          &weights) const
{
  const double tol = 1e-10;
  const unsigned int n_points = surrounding_points.size();

  Assert(n_points > 0,
         ExcMessage("There should be at least one point."));

  Assert(n_points == weights.size(),
         ExcMessage("There should be as many surrounding points as weights given."));

  Assert(std::abs(std::accumulate(weights.begin(), weights.end(), 0.0)-1.0) < tol,
         ExcMessage("The weights for the individual points should sum to 1!"));

  // First sort points in the order of their weights. This is done to
  // produce unique points even if get_intermediate_points is not
  // associative (as for the SphericalManifold).
  boost::container::small_vector<unsigned int, 100> permutation(n_points);
  std::iota(permutation.begin(), permutation.end(), 0u);
  std::sort(permutation.begin(), permutation.end(),
            [&weights](const std::size_t a, const std::size_t b)
  {
    return weights[a] < weights[b];
  });

  // Now loop over points in the order of their associated weight
  Point<spacedim> p = surrounding_points[permutation[0]];
  double w = weights[permutation[0]];

  for (unsigned int i=1; i<n_points; ++i)
    {
      double weight = 0.0;
      if ( (weights[permutation[i]] + w) < tol )
        weight = 0.0;
      else
        weight =  w/(weights[permutation[i]] + w);

      if (std::abs(weight) > 1e-14)
        p = get_intermediate_point(p, surrounding_points[permutation[i]],1.0 - weight );
      w += weights[permutation[i]];
    }

  return p;
}



template <int dim, int spacedim>
void
Manifold<dim, spacedim>::
get_new_points (const ArrayView<const Point<spacedim>> &surrounding_points,
                const Table<2,double>                  &weights,
                ArrayView<Point<spacedim>>              new_points) const
{
  AssertDimension(surrounding_points.size(), weights.size(1));

  for (unsigned int row=0; row<weights.size(0); ++row)
    {
      new_points[row] = get_new_point(make_array_view(surrounding_points.begin(),
                                                      surrounding_points.end()),
                                      make_array_view(weights,
                                                      row));
    }
}



template <>
Tensor<1,2>
Manifold<2, 2>::
normal_vector (const Triangulation<2, 2>::face_iterator &face,
               const Point<2> &p) const
{
  const int spacedim=2;

  // get the tangent vector from the point 'p' in the direction of the further
  // one of the two vertices that make up the face of this 2d cell
  const Tensor<1,spacedim> tangent
    = ((p-face->vertex(0)).norm_square() > (p-face->vertex(1)).norm_square() ?
       -get_tangent_vector(p, face->vertex(0)) :
       get_tangent_vector(p, face->vertex(1)));

  // then rotate it by 90 degrees
  const Tensor<1,spacedim> normal = cross_product_2d(tangent);
  return normal/normal.norm();
}



template <>
Tensor<1,3>
Manifold<3, 3>::
normal_vector (const Triangulation<3, 3>::face_iterator &face,
               const Point<3> &p) const
{
  const int spacedim=3;
  Tensor<1,spacedim> t1,t2;

  // Take the difference between p and all four vertices
  unsigned int min_index=0;
  double       min_distance = (p-face->vertex(0)).norm_square();

  for (unsigned int i=1; i<4; ++i)
    {
      const Tensor<1,spacedim> dp = p-face->vertex(i);
      double distance = dp.norm_square();
      if (distance < min_distance)
        {
          min_index = i;
          min_distance = distance;
        }
    }
  // Verify we have a valid vertex index
  AssertIndexRange(min_index, 4);

  // Now figure out which vertices are best to compute tangent vectors.
  // We split the cell in a central diamond of points closer to the
  // center than to any of the vertices, and the 4 triangles in the
  // corner. The central diamond is split into its upper and lower
  // half. For each of these 6 cases, the following encodes a list
  // of two vertices each to which we compute the tangent vectors,
  // and then take the cross product. See the documentation of this
  // function for exact details.
  if ((p-face->center()).norm_square() < min_distance)
    {
      // we are close to the face center: pick two consecutive vertices,
      // but not the closest one. We make sure the direction is always
      // the same.
      if (min_index < 2)
        {
          t1 = get_tangent_vector(p, face->vertex(3));
          t2 = get_tangent_vector(p, face->vertex(2));
        }
      else
        {
          t1 = get_tangent_vector(p, face->vertex(0));
          t2 = get_tangent_vector(p, face->vertex(1));
        }
    }
  else
    {
      // we are closer to one of the vertices than to the
      // center of the face
      switch (min_index)
        {
        case 0:
        {
          t1 = get_tangent_vector(p, face->vertex(1));
          t2 = get_tangent_vector(p, face->vertex(2));
          break;
        }
        case 1:
        {
          t1 = get_tangent_vector(p, face->vertex(3));
          t2 = get_tangent_vector(p, face->vertex(0));
          break;
        }
        case 2:
        {
          t1 = get_tangent_vector(p, face->vertex(0));
          t2 = get_tangent_vector(p, face->vertex(3));
          break;
        }
        case 3:
        {
          t1 = get_tangent_vector(p, face->vertex(2));
          t2 = get_tangent_vector(p, face->vertex(1));
          break;
        }
        default:
          Assert(false, ExcInternalError());
          break;
        }
    }

  const Tensor<1,spacedim> normal = cross_product_3d(t1,t2);
  return normal/normal.norm();
}



template <int dim, int spacedim>
Tensor<1,spacedim>
Manifold<dim, spacedim>::
normal_vector (const typename Triangulation<dim, spacedim>::face_iterator &/*face*/,
               const Point<spacedim> &/*p*/) const
{
  Assert(false, ExcPureFunctionCalled());
  return Tensor<1,spacedim>();
}



template <>
void
Manifold<2, 2>::
get_normals_at_vertices (const Triangulation<2, 2>::face_iterator &face,
                         FaceVertexNormals &n) const
{
  n[0] = cross_product_2d(get_tangent_vector(face->vertex(0), face->vertex(1)));
  n[1] = -cross_product_2d(get_tangent_vector(face->vertex(1), face->vertex(0)));

  for (unsigned int i=0; i<2; ++i)
    {
      Assert(n[i].norm() != 0, ExcInternalError("Something went wrong. The "
                                                "computed normals have "
                                                "zero length."));
      n[i] /= n[i].norm();
    }
}



template <>
void
Manifold<3, 3>::
get_normals_at_vertices (const Triangulation<3, 3>::face_iterator &face,
                         FaceVertexNormals &n) const
{
  n[0] = cross_product_3d
         (get_tangent_vector(face->vertex(0), face->vertex(1)),
          get_tangent_vector(face->vertex(0), face->vertex(2)));

  n[1] = cross_product_3d
         (get_tangent_vector(face->vertex(1), face->vertex(3)),
          get_tangent_vector(face->vertex(1), face->vertex(0)));

  n[2] = cross_product_3d
         (get_tangent_vector(face->vertex(2), face->vertex(0)),
          get_tangent_vector(face->vertex(2), face->vertex(3)));

  n[3] = cross_product_3d
         (get_tangent_vector(face->vertex(3), face->vertex(2)),
          get_tangent_vector(face->vertex(3), face->vertex(1)));

  for (unsigned int i=0; i<4; ++i)
    {
      Assert(n[i].norm() != 0, ExcInternalError("Something went wrong. The "
                                                "computed normals have "
                                                "zero length."));
      n[i] /=n[i].norm();
    }
}



template <int dim, int spacedim>
void
Manifold<dim, spacedim>::
get_normals_at_vertices (const typename Triangulation<dim, spacedim>::face_iterator &face,
                         FaceVertexNormals &n) const
{
  for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_face; ++v)
    {
      n[v] = normal_vector(face, face->vertex(v));
      n[v] /= n[v].norm();
    }
}



template <int dim, int spacedim>
Point<spacedim>
Manifold<dim, spacedim>::
get_new_point_on_line (const typename Triangulation<dim, spacedim>::line_iterator &line) const
{
  const auto points_weights = get_default_points_and_weights(line);
  return get_new_point (make_array_view(points_weights.first.begin(),
                                        points_weights.first.end()),
                        make_array_view(points_weights.second.begin(),
                                        points_weights.second.end()));
}



template <int dim, int spacedim>
Point<spacedim>
Manifold<dim, spacedim>::
get_new_point_on_quad (const typename Triangulation<dim, spacedim>::quad_iterator &quad) const
{
  const auto points_weights = get_default_points_and_weights(quad);
  return get_new_point (make_array_view(points_weights.first.begin(),
                                        points_weights.first.end()),
                        make_array_view(points_weights.second.begin(),
                                        points_weights.second.end()));
}



template <int dim, int spacedim>
Point<spacedim>
Manifold<dim,spacedim>::
get_new_point_on_face (const typename Triangulation<dim,spacedim>::face_iterator &face) const
{
  Assert (dim>1, ExcImpossibleInDim(dim));

  switch (dim)
    {
    case 2:
      return get_new_point_on_line (face);
    case 3:
      return get_new_point_on_quad (face);
    }

  return Point<spacedim>();
}



template <int dim, int spacedim>
Point<spacedim>
Manifold<dim,spacedim>::
get_new_point_on_cell (const typename Triangulation<dim,spacedim>::cell_iterator &cell) const
{
  switch (dim)
    {
    case 1:
      return get_new_point_on_line (cell);
    case 2:
      return get_new_point_on_quad (cell);
    case 3:
      return get_new_point_on_hex (cell);
    }

  return Point<spacedim>();
}



template <>
Point<1>
Manifold<1,1>::
get_new_point_on_face (const Triangulation<1,1>::face_iterator &) const
{
  Assert (false, ExcImpossibleInDim(1));
  return Point<1>();
}



template <>
Point<2>
Manifold<1,2>::
get_new_point_on_face (const Triangulation<1,2>::face_iterator &) const
{
  Assert (false, ExcImpossibleInDim(1));
  return Point<2>();
}



template <>
Point<3>
Manifold<1,3>::
get_new_point_on_face (const Triangulation<1,3>::face_iterator &) const
{
  Assert (false, ExcImpossibleInDim(1));
  return Point<3>();
}



template <>
Point<1>
Manifold<1,1>::
get_new_point_on_quad (const Triangulation<1,1>::quad_iterator &) const
{
  Assert (false, ExcImpossibleInDim(1));
  return Point<1>();
}



template <>
Point<2>
Manifold<1,2>::
get_new_point_on_quad (const Triangulation<1,2>::quad_iterator &) const
{
  Assert (false, ExcImpossibleInDim(1));
  return Point<2>();
}



template <>
Point<3>
Manifold<1,3>::
get_new_point_on_quad (const Triangulation<1,3>::quad_iterator &) const
{
  Assert (false, ExcImpossibleInDim(1));
  return Point<3>();
}



template <int dim, int spacedim>
Point<spacedim>
Manifold<dim, spacedim>::
get_new_point_on_hex (const typename Triangulation<dim, spacedim>::hex_iterator &/*hex*/) const
{
  Assert (false, ExcImpossibleInDim(dim));
  return Point<spacedim>();
}



template <>
Point<3>
Manifold<3,3>::
get_new_point_on_hex (const Triangulation<3, 3>::hex_iterator &hex) const
{
  const auto points_weights = get_default_points_and_weights(hex, true);
  return get_new_point (make_array_view(points_weights.first.begin(),
                                        points_weights.first.end()),
                        make_array_view(points_weights.second.begin(),
                                        points_weights.second.end()));
}



template <int dim, int spacedim>
Tensor<1,spacedim>
Manifold<dim,spacedim>::get_tangent_vector(const Point<spacedim> &x1,
                                           const Point<spacedim> &x2) const
{
  const double epsilon = 1e-8;

  const std::array<Point<spacedim>, 2> points {{x1, x2}};
  const std::array<double, 2> weights {{epsilon, 1.0 - epsilon}};
  const Point<spacedim> neighbor_point = get_new_point (make_array_view(points.begin(),
                                                        points.end()),
                                                        make_array_view(weights.begin(),
                                                            weights.end()));
  return (neighbor_point-x1)/epsilon;
}

/* -------------------------- FlatManifold --------------------- */


template <int dim, int spacedim>
FlatManifold<dim,spacedim>::FlatManifold (const Tensor<1,spacedim> &periodicity,
                                          const double tolerance)
  :
  periodicity(periodicity),
  tolerance(tolerance)
{}



template <int dim, int spacedim>
Point<spacedim>
FlatManifold<dim, spacedim>::
get_new_point (const ArrayView<const Point<spacedim>> &surrounding_points,
               const ArrayView<const double>          &weights) const
{
  Assert(std::abs(std::accumulate(weights.begin(), weights.end(), 0.0)-1.0) < 1e-10,
         ExcMessage("The weights for the individual points should sum to 1!"));

  Point<spacedim> p;

  // if there is no periodicity, use a shortcut
  if (periodicity == Tensor<1,spacedim>())
    {
      for (unsigned int i=0; i<surrounding_points.size(); ++i)
        p += surrounding_points[i] * weights[i];
    }
  else
    {
      Tensor<1,spacedim> minP = periodicity;

      for (unsigned int d=0; d<spacedim; ++d)
        if (periodicity[d] > 0)
          for (unsigned int i=0; i<surrounding_points.size(); ++i)
            {
              minP[d] = std::min(minP[d], surrounding_points[i][d]);
              Assert( (surrounding_points[i][d] < periodicity[d]+tolerance*periodicity[d]) ||
                      (surrounding_points[i][d] >= -tolerance*periodicity[d]),
                      ExcPeriodicBox(d, surrounding_points[i], periodicity[i]));
            }

      // compute the weighted average point, possibly taking into account periodicity
      for (unsigned int i=0; i<surrounding_points.size(); ++i)
        {
          Point<spacedim> dp;
          for (unsigned int d=0; d<spacedim; ++d)
            if (periodicity[d] > 0)
              dp[d] = ( (surrounding_points[i][d]-minP[d]) > periodicity[d]/2.0 ?
                        -periodicity[d] : 0.0 );

          p += (surrounding_points[i]+dp)*weights[i];
        }

      // if necessary, also adjust the weighted point by the periodicity
      for (unsigned int d=0; d<spacedim; ++d)
        if (periodicity[d] > 0)
          if (p[d] < 0)
            p[d] += periodicity[d];
    }

  return project_to_manifold(surrounding_points, p);
}



template <int dim, int spacedim>
void
FlatManifold<dim, spacedim>::
get_new_points (const ArrayView<const Point<spacedim>> &surrounding_points,
                const Table<2,double>                  &weights,
                ArrayView<Point<spacedim>>              new_points) const
{
  AssertDimension(surrounding_points.size(), weights.size(1));
  if (weights.size(0) == 0)
    return;

  const std::size_t n_points = surrounding_points.size();

  Tensor<1,spacedim> minP = periodicity;
  for (unsigned int d=0; d<spacedim; ++d)
    if (periodicity[d] > 0)
      for (unsigned int i=0; i<n_points; ++i)
        {
          minP[d] = std::min(minP[d], surrounding_points[i][d]);
          Assert( (surrounding_points[i][d] < periodicity[d]+tolerance*periodicity[d]) ||
                  (surrounding_points[i][d] >= -tolerance*periodicity[d]),
                  ExcPeriodicBox(d, surrounding_points[i], periodicity[i]));
        }

  // check whether periodicity shifts some of the points. Only do this if
  // necessary to avoid memory allocation
  const Point<spacedim> *surrounding_points_start = &surrounding_points[0];

  boost::container::small_vector<Point<spacedim>, 200> modified_points;
  bool adjust_periodicity = false;
  for (unsigned int d=0; d<spacedim; ++d)
    if (periodicity[d] > 0)
      for (unsigned int i=0; i<n_points; ++i)
        if ((surrounding_points[i][d]-minP[d]) > periodicity[d]/2.0)
          {
            adjust_periodicity = true;
            break;
          }
  if (adjust_periodicity == true)
    {
      modified_points.resize(surrounding_points.size());
      std::copy(surrounding_points.begin(), surrounding_points.end(),
                modified_points.begin());
      for (unsigned int d=0; d<spacedim; ++d)
        if (periodicity[d] > 0)
          for (unsigned int i=0; i<n_points; ++i)
            if ((surrounding_points[i][d]-minP[d]) > periodicity[d]/2.0)
              modified_points[i][d] -= periodicity[d];
      surrounding_points_start = &modified_points[0];
    }

  // Now perform the interpolation
  for (unsigned int row=0; row<weights.size(0); ++row)
    {
      Assert(std::abs(std::accumulate(&weights(row,0), &weights(row,0)+n_points, 0.0)-1.0) < 1e-10,
             ExcMessage("The weights for the individual points should sum to 1!"));
      Point<spacedim> new_point;
      for (unsigned int p=0; p<n_points; ++p)
        new_point += surrounding_points_start[p] * weights(row,p);

      // if necessary, also adjust the weighted point by the periodicity
      for (unsigned int d=0; d<spacedim; ++d)
        if (periodicity[d] > 0)
          if (new_point[d] < 0)
            new_point[d] += periodicity[d];

      // TODO should this use surrounding_points_start or surrounding_points?
      // The older version used surrounding_points
      new_points[row] = project_to_manifold(make_array_view(surrounding_points.begin(),
                                                            surrounding_points.end()),
                                            new_point);
    }
}


template <int dim, int spacedim>
Point<spacedim>
FlatManifold<dim, spacedim>::project_to_manifold
(const ArrayView<const Point<spacedim>> &/*vertices*/,
 const Point<spacedim>                  &candidate) const
{
  return candidate;
}



template <int dim, int spacedim>
const Tensor<1,spacedim> &
FlatManifold<dim, spacedim>::get_periodicity() const
{
  return periodicity;
}



template <int dim, int spacedim>
Tensor<1,spacedim>
FlatManifold<dim, spacedim>::get_tangent_vector (const Point<spacedim> &x1,
                                                 const Point<spacedim> &x2) const
{
  Tensor<1,spacedim> direction = x2-x1;

  // see if we have to take into account periodicity. if so, we need
  // to make sure that if a distance in one coordinate direction
  // is larger than half of the box length, then go the other way
  // around (i.e., via the periodic box)
  for (unsigned int d=0; d<spacedim; ++d)
    if (periodicity[d] > tolerance)
      {
        if (direction[d] < -periodicity[d]/2)
          direction[d] += periodicity[d];
        else if (direction[d] > periodicity[d]/2)
          direction[d] -= periodicity[d];
      }

  return direction;
}



/* -------------------------- ChartManifold --------------------- */
template <int dim, int spacedim, int chartdim>
ChartManifold<dim,spacedim,chartdim>::ChartManifold (const Tensor<1,chartdim> &periodicity)
  :
  sub_manifold(periodicity)
{}



template <int dim, int spacedim, int chartdim>
Point<spacedim>
ChartManifold<dim,spacedim,chartdim>::
get_new_point (const ArrayView<const Point<spacedim>> &surrounding_points,
               const ArrayView<const double>          &weights) const
{
  const std::size_t n_points = surrounding_points.size();

  boost::container::small_vector<Point<chartdim>, 200> chart_points(n_points);

  for (unsigned int i=0; i<n_points; ++i)
    chart_points[i] = pull_back(surrounding_points[i]);

  const Point<chartdim> p_chart = sub_manifold.get_new_point
                                  (make_array_view(chart_points.begin(), chart_points.end()),
                                   weights);

  return push_forward(p_chart);
}



template <int dim, int spacedim, int chartdim>
void
ChartManifold<dim,spacedim,chartdim>::
get_new_points (const ArrayView<const Point<spacedim>> &surrounding_points,
                const Table<2,double>                  &weights,
                ArrayView<Point<spacedim>>              new_points) const
{
  Assert(weights.size(0) > 0, ExcEmptyObject());
  AssertDimension(surrounding_points.size(), weights.size(1));

  const std::size_t n_points = surrounding_points.size();

  boost::container::small_vector<Point<chartdim>, 200> chart_points(n_points);
  for (std::size_t i=0; i<n_points; ++i)
    chart_points[i] = pull_back(surrounding_points[i]);

  boost::container::small_vector<Point<chartdim>, 200> new_points_on_chart(weights.size(0));
  sub_manifold.get_new_points(make_array_view(chart_points.begin(),
                                              chart_points.end()),
                              weights,
                              make_array_view(new_points_on_chart.begin(),
                                              new_points_on_chart.end()));

  for (std::size_t row=0; row<weights.size(0); ++row)
    new_points[row] = push_forward(new_points_on_chart[row]);
}



template <int dim, int spacedim, int chartdim>
DerivativeForm<1,chartdim,spacedim>
ChartManifold<dim,spacedim,chartdim>::
push_forward_gradient(const Point<chartdim> &) const
{
  // function must be implemented in a derived class to be usable,
  // as discussed in this function's documentation
  Assert (false, ExcPureFunctionCalled());
  return DerivativeForm<1,chartdim,spacedim>();
}



template <int dim, int spacedim, int chartdim>
Tensor<1,spacedim>
ChartManifold<dim,spacedim,chartdim>::
get_tangent_vector (const Point<spacedim> &x1,
                    const Point<spacedim> &x2) const
{
  const DerivativeForm<1,chartdim,spacedim> F_prime = push_forward_gradient(pull_back(x1));

  // ensure that the chart is not singular by asserting that its
  // derivative has a positive determinant. we need to make this
  // comparison relative to the size of the derivative. since the
  // determinant is the product of chartdim factors, take the
  // chartdim-th root of it in comparing against the size of the
  // derivative
  Assert (std::pow(std::abs(F_prime.determinant()), 1./chartdim) >= 1e-12 * F_prime.norm(),
          ExcMessage("The derivative of a chart function must not be singular."));

  const Tensor<1,chartdim>                  delta   = sub_manifold.get_tangent_vector(pull_back(x1),
                                                      pull_back(x2));

  Tensor<1,spacedim> result;
  for (unsigned int i=0; i<spacedim; ++i)
    result[i] += F_prime[i] * delta;

  return result;
}



template <int dim, int spacedim, int chartdim>
const Tensor<1, chartdim> &
ChartManifold<dim, spacedim, chartdim>::get_periodicity() const
{
  return sub_manifold.get_periodicity();
}

// explicit instantiations
#include "manifold.inst"

DEAL_II_NAMESPACE_CLOSE
