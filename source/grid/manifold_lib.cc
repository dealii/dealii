// ---------------------------------------------------------------------
// $Id: manifold_lib.cc 30130 2013-07-23 13:01:18Z heltai $
//
// Copyright (C) 2013 by the deal.II authors
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

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/base/tensor.h>
#include <cmath>

DEAL_II_NAMESPACE_OPEN

template <int dim, int spacedim>
SphericalManifold<dim,spacedim>::SphericalManifold(const Point<spacedim> center):
  ManifoldChart<dim,spacedim,spacedim>(SphericalManifold<dim,spacedim>::get_periodicity()),
  center(center)
{
  Assert(spacedim != 1, ExcImpossibleInDim(1));
}



template <int dim, int spacedim>
Point<spacedim>
SphericalManifold<dim,spacedim>::get_periodicity()  {
  Point<spacedim> periodicity;
  periodicity[spacedim-1] = 2*numbers::PI; // theta and phi period.
  return periodicity;
}


template <int dim, int spacedim>
Point<spacedim>
SphericalManifold<dim,spacedim>::get_new_point(const Quadrature<spacedim> &quad) const {
  if(spacedim == 2)
    return ManifoldChart<dim,spacedim,spacedim>::get_new_point(quad);
  else {
    double rho_average = 0;
    Point<spacedim> mid_point;
    for(unsigned int i=0; i<quad.size(); ++i) {
      rho_average += quad.weight(i)*(quad.point(i)-center).norm();
      mid_point += quad.weight(i)*quad.point(i);
    }
    // Project the mid_pont back to the right location 
    Point<spacedim> R = mid_point-center;
    // Scale it to have radius rho_average
    R *= rho_average/R.norm();
    // And return it.
    return center+R;
  }
}



template <int dim, int spacedim>
Point<spacedim>
SphericalManifold<dim,spacedim>::push_forward(const Point<spacedim> &spherical_point) const {
  Assert(spherical_point[0] >=0.0, 
	 ExcMessage("Negative radius for given point."));
  const double rho = spherical_point[0];
  const double theta = spherical_point[1];
  
  Point<spacedim> p;
  if(rho > 1e-10) 
    switch(spacedim) {
    case 2: 
      p[0] = rho*cos(theta);
      p[1] = rho*sin(theta);
      break;
    case 3:
      {
	const double &phi= spherical_point[2];
	p[0] = rho*sin(theta)*cos(phi);
	p[1] = rho*sin(theta)*sin(phi);
	p[2] = rho*cos(theta);
      }
      break;
    default:
      Assert(false, ExcInternalError());
  }
  return p+center;
}

template <int dim, int spacedim>
Point<spacedim>
SphericalManifold<dim,spacedim>::pull_back(const Point<spacedim> &space_point) const {
  const Point<spacedim> R = space_point-center;
  const double rho = R.norm();
  const double x = R[0];
  const double y = R[1];
  
  Point<spacedim> p;
  p[0] = rho;
  
  switch(spacedim) {
  case 2: 
    p[1] = atan2(y,x);
    if(p[1] < 0)
      p[1] += 2*numbers::PI;
    break;
  case 3:
    {
      const double z = R[2];
      p[2] = atan2(y,x); // phi
      if(p[2] < 0)
	p[2] += 2*numbers::PI; // phi is periodic
      p[1] = atan2(sqrt(x*x+y*y),z);  // theta
    }
    break;
  default:
    Assert(false, ExcInternalError());
  }
  return p;
}


// ============================================================
// CylindricalManifold
// ============================================================

template <int dim, int spacedim>
CylindricalManifold<dim,spacedim>::CylindricalManifold(const unsigned int axis, 
						       const double tolerance) :
  direction (Point<spacedim>::unit_vector(axis)),
  point_on_axis (Point<spacedim>()),
  tolerance(tolerance)
{
  Assert(spacedim > 1, ExcImpossibleInDim(1));
}


template <int dim, int spacedim>
CylindricalManifold<dim,spacedim>::CylindricalManifold(const Point<spacedim> &direction,
						       const Point<spacedim> &point_on_axis, 
						       const double tolerance) :
  direction (direction),
  point_on_axis (point_on_axis),
  tolerance(tolerance)
{
  Assert(spacedim > 2, ExcImpossibleInDim(spacedim));
}




template <int dim, int spacedim>
Point<spacedim>
CylindricalManifold<dim,spacedim>::
get_new_point (const Quadrature<spacedim> &quad) const
{
  const std::vector<Point<spacedim> > &surrounding_points = quad.get_points();
  const std::vector<double> &weights = quad.get_weights();

  // compute a proposed new point  
  Point<spacedim> middle = flat_manifold.get_new_point(quad);

  double radius = 0;
  Point<spacedim> on_plane;
  
  for(unsigned int i=0; i<surrounding_points.size(); ++i)
    {
      on_plane = surrounding_points[i]-point_on_axis;
      on_plane = on_plane - (on_plane*direction) * direction;
      radius += weights[i]*on_plane.norm();
    }
  
  // we then have to project this point out to the given radius from
  // the axis. to this end, we have to take into account the offset
  // point_on_axis and the direction of the axis
  const Point<spacedim> vector_from_axis = (middle-point_on_axis) -
					   ((middle-point_on_axis) * direction) * direction;

  // scale it to the desired length and put everything back together,
  // unless we have a point on the axis
  if (vector_from_axis.norm() <= tolerance * middle.norm())
    return middle;

  else
    return (vector_from_axis / vector_from_axis.norm() * radius +
	    ((middle-point_on_axis) * direction) * direction +
	    point_on_axis);
}
		
// explicit instantiations
#include "manifold_lib.inst"

DEAL_II_NAMESPACE_CLOSE
