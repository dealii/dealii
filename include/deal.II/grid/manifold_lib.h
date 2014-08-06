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
 * The two dimensional implementation of this class works by
 * transforming points to spherical coordinates, taking the average in
 * that coordinate system, and then transforming back the point to
 * Cartesian coordinates. For the three dimensional case, we use a
 * simpler approach: we take the average of the norm of the points,
 * and use this value to shift the average point along the radial
 * direction. In order for this manifold to work correctly, it cannot
 * be attached to cells containing the center of the coordinate
 * system. This point is a singular point of the coordinate
 * transformation, and there taking averages does not make any sense.
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
   * to transform from Cartesian to spherical coordinate systems,
   * taking into account the periodicity of base Manifold in two
   * dimensions, while in three dimensions it takes the middle point,
   * and project it along the radius using the average radius of the
   * surrounding points.
   */
  SphericalManifold(const Point<spacedim> center = Point<spacedim>());

  /**
   * Pull back the given point from the Euclidean space. Will return
   * the polar coordinates associated with the point @p
   * space_point. Only used when spacedim = 2.
   */
  virtual Point<spacedim>
  pull_back(const Point<spacedim> &space_point) const;

  /**
   * Given a point in the spherical coordinate system, this method
   * returns the Euclidean coordinates associated to the polar
   * coordinates @p chart_point. Only used when spacedim = 3.
   */
  virtual Point<spacedim>
  push_forward(const Point<spacedim> &chart_point) const;

  /**
   * Let the new point be the average sum of surrounding vertices.
   *
   * In the two dimensional implementation, we use the pull_back and
   * push_forward mechanism. For three dimensions, this does not work
   * well, so we overload the get_new_point function directly.
   */
  virtual Point<spacedim>
  get_new_point(const Quadrature<spacedim> &quad) const;
  
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


/**
 * Cylindrical Manifold description.  In three dimensions, points are
 * transformed using a cylindrical coordinate system along the
 * <tt>x-</tt>, <tt>y-</tt> or <tt>z</tt>-axis (when using the first
 * constructor of this class), or an arbitrarily oriented cylinder
 * described by the direction of its axis and a point located on the
 * axis.
 *
 * This class was developed to be used in conjunction with the @p
 * cylinder or @p cylinder_shell functions of GridGenerator. This
 * function will throw an exception whenever spacedim is not equal to
 * three. 
 *
 * @ingroup manifold
 *
 * @author Luca Heltai, 2014
 */
template <int dim, int spacedim = dim>
class CylindricalManifold : public Manifold<dim,spacedim>
{
public:
  /**
   * Constructor. Using default values for the constructor arguments
   * yields a cylinder along the x-axis (<tt>axis=0</tt>). Choose
   * <tt>axis=1</tt> or <tt>axis=2</tt> for a tube along the y- or
   * z-axis, respectively. The tolerance value is used to determine
   * if a point is on the axis.
   */
  CylindricalManifold (const unsigned int axis = 0, 
		       const double tolerance = 1e-10);

  /**
   * Constructor. If constructed with this constructor, the manifold
   * described is a cylinder with an axis that points in direction
   * #direction and goes through the given #point_on_axis. The
   * direction may be arbitrarily scaled, and the given point may be
   * any point on the axis. The tolerance value is used to determine
   * if a point is on the axis. 
   */
  CylindricalManifold (const Point<spacedim> &direction,
		       const Point<spacedim> &point_on_axis, 
		       const double tolerance = 1e-10);

 /**
   * Compute new points on the CylindricalManifold. See the documentation
   * of the base class for a detailed description of what this
   * function does.
   */
  virtual Point<spacedim>
  get_new_point(const Quadrature<spacedim> &quad) const;

protected:
  /**
   * The direction vector of the axis.
   */
  const Point<spacedim> direction;

  /**
   * An arbitrary point on the axis.
   */
  const Point<spacedim> point_on_axis;

private:
  /**
   * Helper FlatManifold to compute temptative midpoints.
   */
  FlatManifold<dim,spacedim> flat_manifold;
  
  /**
   * Relative tolerance to measure zero distances.
   */
  double tolerance;
};

DEAL_II_NAMESPACE_CLOSE

#endif
