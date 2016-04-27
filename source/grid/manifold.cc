// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2016 by the deal.II authors
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
#include <deal.II/grid/manifold.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/fe/fe_q.h>
#include <cmath>

DEAL_II_NAMESPACE_OPEN

using namespace Manifolds;

/* -------------------------- Manifold --------------------- */


template <int dim, int spacedim>
Manifold<dim, spacedim>::~Manifold ()
{}



template <int dim, int spacedim>
Point<spacedim>
Manifold<dim, spacedim>::
project_to_manifold (const std::vector<Point<spacedim> > &,
                     const Point<spacedim> &) const
{
  Assert (false, ExcPureFunctionCalled());
  return Point<spacedim>();
}



template <int dim, int spacedim>
Point<spacedim>
Manifold<dim, spacedim>::
get_new_point (const Quadrature<spacedim> &quad) const
{
  Assert(std::abs(std::accumulate(quad.get_weights().begin(), quad.get_weights().end(), 0.0)-1.0) < 1e-10,
         ExcMessage("The weights for the individual points should sum to 1!"));

  Point<spacedim> p;
  for (unsigned int i=0; i<quad.size(); ++i)
    p += quad.weight(i) * quad.point(i);

  return project_to_manifold(quad.get_points(), p);
}



template <>
Tensor<1,2>
Manifold<2, 2>::
normal_vector (const Triangulation<2, 2>::face_iterator &face,
               const Point<2> &p) const
{
  const int spacedim=2;

  Tensor<1,spacedim> tangent = ((p-face->vertex(0)).norm_square() > (p-face->vertex(1)).norm_square() ?
                                -get_tangent_vector(p, face->vertex(0)) :
                                get_tangent_vector(p, face->vertex(1)));
  Tensor<1,spacedim> normal = cross_product_2d(tangent);
  return normal/normal.norm();
}

template<>
Tensor<1,3>
Manifold<3, 3>::
normal_vector (const Triangulation<3, 3>::face_iterator &face,
               const Point<3> &p) const
{
  const int spacedim=3;
  Tensor<1,spacedim> t1,t2;

  // Take the difference between p and all four vertices
  int min_index=0;
  Tensor<1,spacedim> dp = p-face->vertex(0);
  double min_distance = dp.norm_square();

  for (unsigned int i=1; i<4; ++i)
    {
      dp = p-face->vertex(i);
      double distance = dp.norm_square();
      if (distance < min_distance)
        {
          min_index = i;
          min_distance = distance;
        }
    }
  // Verify we have a valid vertex index
  AssertIndexRange(min_index, 4);

  // Now figure out which vertices are better to compute tangent vectors
  // we split the cell in 4 quadrants, and use the most orthogonal vertices
  // to the closest vertex if we ar far from the center, othewise we use
  // the two consecutive vertices, on the opposite side with respect to
  // the face center.
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
  Tensor<1,spacedim> normal = cross_product_3d(t1,t2);
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
  return get_new_point (get_default_quadrature(line));
}



template <int dim, int spacedim>
Point<spacedim>
Manifold<dim, spacedim>::
get_new_point_on_quad (const typename Triangulation<dim, spacedim>::quad_iterator &quad) const
{
  return get_new_point (get_default_quadrature(quad));
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
  return get_new_point(get_default_quadrature(hex, true));
}



template <int dim, int spacedim>
Tensor<1,spacedim>
Manifold<dim,spacedim>::get_tangent_vector(const Point<spacedim> &x1,
                                           const Point<spacedim> &x2) const
{
  const double epsilon = 1e-8;

  std::vector<Point<spacedim> > q;
  q.push_back(x1);
  q.push_back(x2);

  std::vector<double> w;
  w.push_back(epsilon);
  w.push_back(1.0-epsilon);

  const Tensor<1,spacedim> neighbor_point = get_new_point (Quadrature<spacedim>(q, w));
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
get_new_point (const Quadrature<spacedim> &quad) const
{
  Assert(std::abs(std::accumulate(quad.get_weights().begin(), quad.get_weights().end(), 0.0)-1.0) < 1e-10,
         ExcMessage("The weights for the individual points should sum to 1!"));

  const std::vector<Point<spacedim> > &surrounding_points = quad.get_points();
  const std::vector<double> &weights = quad.get_weights();

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
  Point<spacedim> p;
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

  return project_to_manifold(surrounding_points, p);
}



template <int dim, int spacedim>
Point<spacedim>
FlatManifold<dim, spacedim>::project_to_manifold (const std::vector<Point<spacedim> > &/*vertices*/,
                                                  const Point<spacedim> &candidate) const
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
ChartManifold<dim,spacedim,chartdim>::~ChartManifold ()
{}



template <int dim, int spacedim, int chartdim>
ChartManifold<dim,spacedim,chartdim>::ChartManifold (const Tensor<1,chartdim> &periodicity)
  :
  sub_manifold(periodicity)
{}



template <int dim, int spacedim, int chartdim>
Point<spacedim>
ChartManifold<dim,spacedim,chartdim>::
get_new_point (const Quadrature<spacedim> &quad) const
{
  const std::vector<Point<spacedim> > &surrounding_points = quad.get_points();
  const std::vector<double> &weights = quad.get_weights();
  std::vector<Point<chartdim> > chart_points(surrounding_points.size());

  for (unsigned int i=0; i<surrounding_points.size(); ++i)
    chart_points[i] = pull_back(surrounding_points[i]);

  const Quadrature<chartdim> chart_quad(chart_points, weights);
  const Point<chartdim> p_chart = sub_manifold.get_new_point(chart_quad);

  return push_forward(p_chart);
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

