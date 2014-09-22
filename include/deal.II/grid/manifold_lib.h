// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2014 by the deal.II authors
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
#include <deal.II/base/function.h>
#include <deal.II/base/function_parser.h>

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
template <int dim, int spacedim = dim>
class SphericalManifold : public ChartManifold<dim, spacedim, spacedim>
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


/**
 * Manifold description derived from ChartManifold, based on explicit
 * Function<spacedim> and Function<chartdim> objects describing the
 * push_forward() and pull_back() functions.
 *
 * You can use this Manifold object to describe any arbitray shape
 * domain, as long as you can express it in terms of an invertible
 * map, for which you provide both the forward expression, and the
 * inverse expression.
 *
 * In debug mode, a check is performed to verify that the
 * tranformations are actually one the inverse of the other.
 *
 * @ingroup manifold
 *
 * @author Luca Heltai, 2014
 */
template <int dim, int spacedim=dim, int chartdim=dim>
class FunctionManifold : public ChartManifold<dim, spacedim, chartdim>
{
public:
  /**
   * Explicit functions constructor. Takes a push_forward function of
   * spacedim components, and a pull_back function of chartdim
   * components. See the documentation of the base class ChartManifold
   * for the meaning of the optional #periodicity argument.
   *
   * The tolerance argument is used in debug mode to actually check
   * that the two functions are one the inverse of the other.
   */
  FunctionManifold(const Function<chartdim> &push_forward_function,
                   const Function<spacedim> &pull_back_function,
                   const Point<chartdim> periodicity=Point<chartdim>(),
                   const double tolerance=1e-10);

  /**
   * Expressions constructor. Takes the expressions of the
   * push_forward function of spacedim components, and of the
   * pull_back function of chartdim components. See the documentation
   * of the base class ChartManifold for the meaning of the optional
   * #periodicity argument.
   *
   * The strings should be the readable by the default constructor of
   * the FunctionParser classes. You can specify custom variable
   * expressions with the last two optional arguments. If you don't,
   * the default names are used, i.e., "x,y,z".
   *
   * The tolerance argument is used in debug mode to actually check
   * that the two functions are one the inverse of the other.
   */
  FunctionManifold(const std::string push_forward_expression,
                   const std::string pull_back_expression,
                   const Point<chartdim> periodicity=Point<chartdim>(),
                   const typename FunctionParser<spacedim>::ConstMap = typename FunctionParser<spacedim>::ConstMap(),
                   const std::string chart_vars=FunctionParser<chartdim>::default_variable_names(),
                   const std::string space_vars=FunctionParser<spacedim>::default_variable_names(),
                   const double tolerance=1e-10);

  /**
   * If needed, we delete the pointers we own.
   */
  ~FunctionManifold();

  /**
   * Given a point in the chartdim coordinate system, uses the
   * push_forward_function to compute the push_forward of points in
   * #chardim space dimensions to #spacedim space dimensions.
   */
  virtual Point<spacedim>
  push_forward(const Point<chartdim> &chart_point) const;

  /**
   * Given a point in the spacedim coordinate system, uses the
   * pull_back_function to compute the pull_back of points in
   * #spacedim space dimensions to #chartdim space dimensions.
   */
  virtual Point<chartdim>
  pull_back(const Point<spacedim> &space_point) const;

private:
  /**
   * Constants for the FunctionParser classes.
   */
  const typename FunctionParser<spacedim>::ConstMap const_map;

  /**
   * Pointer to the push_forward function.
   */
  SmartPointer<const Function<chartdim>,
               FunctionManifold<dim,spacedim,chartdim> > push_forward_function;

  /**
   * Pointer to the pull_back function.
   */
  SmartPointer<const Function<spacedim>,
               FunctionManifold<dim,spacedim,chartdim> > pull_back_function;

  /**
   * Relative tolerance. In debug mode, we check that the two
   * functions provided at construction time are actually one the
   * inverse of the other. This value is used as relative tolerance in
   * this check.
   */
  const double tolerance;

  /**
   * Check ownership of the smart pointers. Indicates whether this
   * class is the owner of the objects pointed to by the previous two
   * member variables.  This value is set in the constructor of the
   * class. If #true, then the destructor will delete the function
   * objects pointed to be the two pointers.
   */
  const bool owns_pointers;
};


DEAL_II_NAMESPACE_CLOSE

#endif
