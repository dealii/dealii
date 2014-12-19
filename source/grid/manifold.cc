// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2013 by the deal.II authors
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
  const std::vector<Point<spacedim> > &surrounding_points = quad.get_points();
  const std::vector<double> &weights = quad.get_weights();
  Point<spacedim> p;

#ifdef DEBUG
  double sum=0;
  for (unsigned int i=0; i<weights.size(); ++i)
    sum+= weights[i];
  Assert(std::abs(sum-1.0) < 1e-10, ExcMessage("Weights should sum to 1!"));
#endif

  for (unsigned int i=0; i<surrounding_points.size(); ++i)
    p += surrounding_points[i]*weights[i];

  return project_to_manifold(surrounding_points, p);
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
get_new_point_on_hex (const typename Triangulation<dim, spacedim>::hex_iterator &hex) const
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


/* -------------------------- FlatManifold --------------------- */


template <int dim, int spacedim>
FlatManifold<dim,spacedim>::FlatManifold (const Point<spacedim> periodicity,
                                          const double tolerance) :
  periodicity(periodicity),
  tolerance(tolerance)
{}

template <int dim, int spacedim>
Point<spacedim>
FlatManifold<dim, spacedim>::
get_new_point (const Quadrature<spacedim> &quad) const
{
  const std::vector<Point<spacedim> > &surrounding_points = quad.get_points();
  const std::vector<double> &weights = quad.get_weights();

#ifdef DEBUG
  double sum=0;
  for (unsigned int i=0; i<weights.size(); ++i)
    sum+= weights[i];
  // Here it is correct to use tolerance as an absolute one, since
  // this should be relative to unity.
  Assert(std::abs(sum-1.0) < tolerance, ExcMessage("Weights should sum to 1!"));
#endif


  Point<spacedim> p;
  Point<spacedim> dp;
  Point<spacedim> minP = periodicity;
  const bool check_period = (periodicity.norm() > tolerance);
  if (check_period)
    for (unsigned int i=0; i<surrounding_points.size(); ++i)
      for (unsigned int d=0; d<spacedim; ++d)
        {
          minP[d] = std::min(minP[d], surrounding_points[i][d]);
          if (periodicity[d] > 0)
            Assert( (surrounding_points[i][d] < periodicity[d]+tolerance*periodicity.norm()) ||
                    (surrounding_points[i][d] >= -tolerance*periodicity.norm()),
                    ExcPeriodicBox(d, surrounding_points[i], periodicity, tolerance*periodicity.norm()));
        }

  for (unsigned int i=0; i<surrounding_points.size(); ++i)
    {
      dp = Point<spacedim>();
      if (check_period)
        {
          for (unsigned int d=0; d<spacedim; ++d)
            if (periodicity[d] > 0)
              dp[d] = ( (surrounding_points[i][d]-minP[d]) > periodicity[d]/2.0 ?
                        -periodicity[d] : 0.0 );
        }
      p += (surrounding_points[i]+dp)*weights[i];
    }
  if (check_period)
    for (unsigned int d=0; d<spacedim; ++d)
      if (periodicity[d] > 0)
        p[d] = (p[d] < 0 ? p[d] + periodicity[d] : p[d]);

  return project_to_manifold(surrounding_points, p);
}

template <int dim, int spacedim>
Point<spacedim>
FlatManifold<dim, spacedim>::project_to_manifold (const std::vector<Point<spacedim> > &vertices,
                                                  const Point<spacedim> &candidate) const
{
  return candidate;
}


/* -------------------------- ChartManifold --------------------- */

template <int dim, int spacedim, int chartdim>
ChartManifold<dim,spacedim,chartdim>::~ChartManifold ()
{}

template <int dim, int spacedim, int chartdim>
ChartManifold<dim,spacedim,chartdim>::ChartManifold (const Point<chartdim> periodicity):
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






// explicit instantiations
#include "manifold.inst"

DEAL_II_NAMESPACE_CLOSE

