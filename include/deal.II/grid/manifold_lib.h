// ---------------------------------------------------------------------
// $Id: tria_boundary_lib.h 30130 2013-07-23 13:01:18Z heltai $
//
// Copyright (C) 1999 - 2013 by the deal.II authors
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

#ifndef __deal2__manifold_lib_h
#define __deal2__manifold_lib_h


#include <deal.II/base/config.h>
#include <deal.II/grid/manifold.h>

DEAL_II_NAMESPACE_OPEN

/**
 * Manifold description for a spherical space coordinate system. 
 *
 * You can use this Manifold object to describe any sphere, circle,
 * hypersphere or hyperdisc in two or three dimensions, both as a
 * co-dimension one manifold descriptor or as co-dimension zero
 * manifold descriptor.
 *
 * The two template arguments match the meaning of the two template
 * arguments in Triangulation<dim, spacedim>, however this Manifold
 * can be used to describe both thin and thick objects, and the
 * behavior is identical when dim <= spacedim, i.e., the functionality
 * of SphericalManifold<2,3> is identical to SphericalManifold<3,3>.
 *
 * @ingroup manifold
 *
 * @author Luca Heltai, 2014
 */
template <int dim, int spacedim>
class SphericalManifold : public ManifoldChart<dim, spacedim, spacedim> 
{
public:
  /**
   * The Constructor takes the center of the spherical coordinates
   * system. This class uses the pull_back and push_forward mechanism
   * to transform from cartesian to spherical coordinate systems,
   * taking into account the periodicity of base Manifold.
   */
  SphericalManifold(const Point<spacedim> center = Point<spacedim>());

  /**
   * Pull back the given point from the Euclidean space. Will return
   * the polar coordinates associated with the point @p space_point.
   */
  virtual Point<spacedim>
  pull_back(const Point<spacedim> &space_point) const;

  /**
   * Given a point in the spherical coordinate system, this method
   * returns the Euclidean coordinates associated to the polar
   * coordinates @p chart_point.
   */
  virtual Point<spacedim>
  push_forward(const Point<spacedim> &chart_point) const;
  
  
  /**
   * The center of the spherical coordinate system.
   */
  const Point<spacedim> center;
private:
  
  /** Helper function which returns the periodicity associated with
      this coordinate system, according to dim, chartdim, and
      spacedim. */
  static Point<spacedim> get_periodicity();
};


DEAL_II_NAMESPACE_CLOSE

#endif
