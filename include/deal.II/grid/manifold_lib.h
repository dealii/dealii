// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_manifold_lib_h
#define dealii_manifold_lib_h


#include <deal.II/base/config.h>

#include <deal.II/base/function.h>
#include <deal.II/base/function_parser.h>

#include <deal.II/grid/manifold.h>

#include <boost/signals2/connection.hpp>


DEAL_II_NAMESPACE_OPEN

// forward declaration
namespace internal
{
  namespace MappingQImplementation
  {
    template <int, int>
    class InverseQuadraticApproximation;
  }
} // namespace internal


/**
 * Manifold description for a polar coordinate system.
 *
 * You can use this Manifold object to describe any sphere, circle,
 * hypersphere or hyperdisc in two or three dimensions, both as a
 * co-dimension one manifold descriptor or as co-dimension zero
 * manifold descriptor, provided that the north and south poles (in
 * three dimensions) and the center (in both two and three dimensions)
 * are excluded from the Manifold (as they are singular points of the
 * polar change of coordinates).
 *
 * The two template arguments match the meaning of the two template arguments
 * in Triangulation<dim, spacedim>; however, this Manifold can be used to
 * describe both thin and thick objects, and the behavior is identical when
 * dim <= spacedim, i.e., the functionality of PolarManifold<2,3> is
 * identical to PolarManifold<3,3>.
 *
 * This class works by transforming points from Cartesian coordinates into Polar
 * coordinates, doing computations (e.g., averages) in that coordinate system,
 * and then transforming the results back to Cartesian coordinates. For more
 * information on this approach see the documentation of the base class
 * (ChartManifold).
 *
 * In order for this manifold to work correctly, it cannot be attached to cells
 * containing the center of the coordinate system or, in 3d, the north and south
 * poles. Those points are singular points of the coordinate mapping;
 * consequently, pull_back() is not defined at those points and taking averages
 * across them does not make sense. We recommend setting the manifold_id of the
 * cell containing the center to numbers::flat_manifold_id and blending in one
 * of two ways:
 *
 * 1. Only set PolarManifold manifold ids on the curved boundary (i.e., set the
 *    manifold id of faces and not cells). deal.II will use a cell-by-cell
 *    version of transfinite interpolation to compute cell center locations for,
 *    e.g., mesh refinement.
 * 2. Set the PolarManifold manifold ids on cells and use
 *    TransfiniteInterpolationManifold.
 *
 * In 3d, we recommend using SphericalManifold to correctly handle the north and
 * south poles.
 *
 * @ingroup manifold
 */
template <int dim, int spacedim = dim>
class PolarManifold : public ChartManifold<dim, spacedim, spacedim>
{
public:
  /**
   * Constructor.
   *
   * @param[in] center The center of the coordinate system. Defaults to the
   * origin.
   */
  PolarManifold(const Point<spacedim> center = Point<spacedim>());

  /**
   * Make a clone of this Manifold object.
   */
  virtual std::unique_ptr<Manifold<dim, spacedim>>
  clone() const override;

  /**
   * Pull back the given point from the Euclidean space. Will return the polar
   * coordinates associated with the point @p space_point.
   */
  virtual Point<spacedim>
  pull_back(const Point<spacedim> &space_point) const override;

  /**
   * Given a point in the spherical coordinate system, this method returns the
   * Euclidean coordinates associated to the polar coordinates @p chart_point.
   */
  virtual Point<spacedim>
  push_forward(const Point<spacedim> &chart_point) const override;

  /**
   * Given a point in the spacedim dimensional Euclidean space, this
   * method returns the derivatives of the function $F$ that maps from
   * the polar coordinate system to the Euclidean coordinate
   * system. In other words, it is a matrix of size
   * $\text{spacedim}\times\text{spacedim}$.
   *
   * This function is used in the computations required by the
   * get_tangent_vector() function.
   *
   * Refer to the general documentation of this class for more information.
   */
  virtual DerivativeForm<1, spacedim, spacedim>
  push_forward_gradient(const Point<spacedim> &chart_point) const override;

  /**
   * @copydoc Manifold::normal_vector()
   */
  virtual Tensor<1, spacedim>
  normal_vector(
    const typename Triangulation<dim, spacedim>::face_iterator &face,
    const Point<spacedim> &p) const override;

  /**
   * Return the center of the spherical coordinate system.
   */
  const Point<spacedim> &
  get_center() const;

  /**
   * The center of the spherical coordinate system.
   *
   * @deprecated Use get_center() instead.
   */
  DEAL_II_DEPRECATED_WITH_COMMENT(
    "Access the center with get_center() instead.")
  const Point<spacedim> center;

private:
  /**
   * The center of the spherical coordinate system.
   *
   * @note This exists to avoid warnings when using center internally.
   */
  const Point<spacedim> p_center;

  /**
   * Helper function which returns the periodicity associated with this
   * coordinate system, according to dim, chartdim, and spacedim.
   */
  static Tensor<1, spacedim>
  get_periodicity();
};



/**
 * Manifold description for a spherical space coordinate system.
 *
 * You can use this Manifold object to describe any sphere, circle,
 * hypersphere or hyperdisc in two or three dimensions. This manifold
 * can be used as a co-dimension one manifold descriptor of a
 * spherical surface embedded in a higher dimensional space, or as a
 * co-dimension zero manifold descriptor for a body with positive
 * volume, provided that the center of the spherical space is excluded
 * from the domain. An example for the use of this function would be
 * in the description of a hyper-shell or hyper-ball geometry, for
 * example after creating a coarse mesh using GridGenerator::hyper_ball().
 * (However, it is worth mentioning that generating a good mesh for
 * a disk or ball is complicated and requires addition steps. See the
 * "Possibilities for extensions" section of step-6 for an extensive
 * discussion of how one would construct such meshes and what one
 * needs to do for it.)
 *
 * The two template arguments match the meaning of the two template arguments
 * in Triangulation<dim, spacedim>, however this Manifold can be used to
 * describe both thin and thick objects, and the behavior is identical when
 * dim <= spacedim, i.e., the functionality of SphericalManifold<2,3> is
 * identical to SphericalManifold<3,3>.
 *
 * While PolarManifold reflects the usual notion of polar coordinates,
 * it may not be suitable for domains that contain either the north or
 * south poles.  Consider for instance the pair of points
 * $x_1=(1,\pi/3,0)$ and $x_2=(1,\pi/3,\pi)$ in polar
 * coordinates (lying on the surface of a sphere with radius one, on
 * a parallel at height $\pi/3$). In this case connecting the points
 * with a straight line in polar coordinates would take the long road
 * around the globe, without passing through the north pole.
 *
 * These two points would be connected (using a PolarManifold) by the curve
 * @f{align*}{
 *   s: [0,1]  & \rightarrow &  \mathbb S^3 \\
 *           t & \mapsto     &  (1,\pi/3,0) + (0,0,t\pi)
 * @f}
 * This curve is not a geodesic on the sphere, and it is not how we
 * would connect those two points. A better curve, would be the one
 * passing through the North pole:
 * @f[
 *  s(t) = x_1 \cos(\alpha(t)) + \kappa \times x_1 \sin(\alpha(t)) +
 *  \kappa ( \kappa \cdot x_1) (1-\cos(\alpha(t))).
 * @f]
 * where $\kappa = \frac{x_1 \times x_2}{\Vert x_1 \times x_2 \Vert}$
 * and $\alpha(t) = t \cdot \arccos(x_1 \cdot x_2)$ for $t\in[0,1]$.
 * Indeed, this is a geodesic, and it is the natural choice when
 * connecting points on the surface of the sphere. In the examples above,
 * the PolarManifold class implements the first way of connecting two
 * points on the surface of a sphere, while SphericalManifold implements
 * the second way, i.e., this Manifold connects points using geodesics.
 * If more than two points are involved through a
 * SphericalManifold::get_new_points() call, a so-called spherical average is
 * used where the final point minimizes the weighted distance to all other
 * points via geodesics.
 *
 * In particular, this class implements a Manifold that joins any two
 * points in space by first projecting them onto the surface of a
 * sphere with unit radius, then connecting them with a geodesic, and
 * finally rescaling the final radius so that the resulting one is the
 * weighted average of the starting radii. This Manifold is identical
 * to PolarManifold in dimension two, while for dimension three it
 * returns points that are more uniformly distributed on the sphere,
 * and it is invariant with respect to rotations of the coordinate
 * system, therefore avoiding the problems that PolarManifold has at
 * the poles. Notice, in particular, that computing tangent vectors at
 * the poles with a PolarManifold is not well defined, while it is
 * perfectly fine with this class.
 *
 * For mathematical reasons, it is impossible to construct a unique
 * map of a sphere using only geodesic curves, and therefore, using
 * this class with MappingManifold is discouraged. If you use this
 * Manifold to describe the geometry of a sphere, you should use
 * MappingQ as the underlying mapping, and not MappingManifold.
 *
 * This Manifold can be used *only* on geometries where a ball with
 * finite radius is removed from the center. Indeed, the center is a
 * singular point for this manifold, and if you try to connect two
 * points across the center, they would travel on spherical
 * coordinates, avoiding the center.
 *
 * The ideal geometry for this Manifold is a hyper shell. If you plan to use
 * this Manifold on a hyper ball, you have to make sure you do not attach this
 * Manifold to the cell containing the center. It is advisable to combine this
 * class with TransfiniteInterpolationManifold to ensure a smooth transition
 * from a curved shape to the straight coordinate system in the center of the
 * ball. (See also the extensive discussion in step-65.)
 *
 * @ingroup manifold
 */
template <int dim, int spacedim = dim>
class SphericalManifold : public Manifold<dim, spacedim>
{
public:
  /**
   * Constructor.
   *
   * @param[in] center The center of the coordinate system. Defaults to the
   * origin.
   */
  SphericalManifold(const Point<spacedim> center = Point<spacedim>());

  /**
   * Make a clone of this Manifold object.
   */
  virtual std::unique_ptr<Manifold<dim, spacedim>>
  clone() const override;

  /**
   * Given any two points in space, first project them on the surface
   * of a sphere with unit radius, then connect them with a geodesic
   * and find the intermediate point, and finally rescale the final
   * radius so that the resulting one is the convex combination of the
   * starting radii.
   */
  virtual Point<spacedim>
  get_intermediate_point(const Point<spacedim> &p1,
                         const Point<spacedim> &p2,
                         const double           w) const override;

  /**
   * Compute the derivative of the get_intermediate_point() function
   * with parameter w equal to zero.
   */
  virtual Tensor<1, spacedim>
  get_tangent_vector(const Point<spacedim> &x1,
                     const Point<spacedim> &x2) const override;

  /**
   * @copydoc Manifold::normal_vector()
   */
  virtual Tensor<1, spacedim>
  normal_vector(
    const typename Triangulation<dim, spacedim>::face_iterator &face,
    const Point<spacedim> &p) const override;

  /**
   * Compute the normal vectors to the boundary at each vertex.
   */
  virtual void
  get_normals_at_vertices(
    const typename Triangulation<dim, spacedim>::face_iterator &face,
    typename Manifold<dim, spacedim>::FaceVertexNormals &face_vertex_normals)
    const override;

  /**
   * Compute a new set of points that interpolate between the given points @p
   * surrounding_points. @p weights is a table with as many columns as @p
   * surrounding_points.size(). The number of rows in @p weights must match
   * the length of @p new_points.
   *
   * This function is optimized to perform on a collection
   * of new points, by collecting operations that are not dependent on the
   * weights outside of the loop over all new points.
   *
   * The implementation does not allow for @p surrounding_points and
   * @p new_points to point to the same array, so make sure to pass different
   * objects into the function.
   */
  virtual void
  get_new_points(const ArrayView<const Point<spacedim>> &surrounding_points,
                 const Table<2, double>                 &weights,
                 ArrayView<Point<spacedim>> new_points) const override;

  /**
   * Return a point on the spherical manifold which is intermediate
   * with respect to the surrounding points.
   */
  virtual Point<spacedim>
  get_new_point(const ArrayView<const Point<spacedim>> &vertices,
                const ArrayView<const double>          &weights) const override;

  /**
   * Return the center of the spherical coordinate system.
   */
  const Point<spacedim> &
  get_center() const;

  /**
   * The center of the spherical coordinate system.
   *
   * @deprecated Use get_center() instead.
   */
  DEAL_II_DEPRECATED_WITH_COMMENT(
    "Access the center with get_center() instead.")
  const Point<spacedim> center;

private:
  /**
   * The center of the spherical coordinate system.
   *
   * @note This exists to avoid warnings when using center internally.
   */
  const Point<spacedim> p_center;

  /**
   * Return a point on the spherical manifold which is intermediate
   * with respect to the surrounding points. This function uses a linear
   * average of the directions to find an estimated point. It returns a pair
   * of radius and direction from the center point to the candidate point.
   */
  std::pair<double, Tensor<1, spacedim>>
  guess_new_point(const ArrayView<const Tensor<1, spacedim>> &directions,
                  const ArrayView<const double>              &distances,
                  const ArrayView<const double>              &weights) const;

  /**
   * This function provides an internal implementation of the get_new_points()
   * interface.
   *
   * It computes a new set of points that interpolate between the given points
   * @p
   * surrounding_points. @p weights is an array view with as many entries as @p
   * surrounding_points.size() times @p new_points.size().
   *
   * This function is optimized to perform on a collection
   * of new points, by collecting operations that are not dependent on the
   * weights outside of the loop over all new points.
   *
   * The implementation does not allow for @p surrounding_points and
   * @p new_points to point to the same array, so make sure to pass different
   * objects into the function.
   */
  void
  do_get_new_points(const ArrayView<const Point<spacedim>> &surrounding_points,
                    const ArrayView<const double>          &weights,
                    ArrayView<Point<spacedim>>              new_points) const;

  /**
   * A manifold description to be used for get_new_point in 2d.
   */
  const PolarManifold<spacedim> polar_manifold;
};

/**
 * Cylindrical Manifold description.  In three dimensions, points are
 * transformed using a cylindrical coordinate system along the <tt>x-</tt>,
 * <tt>y-</tt> or <tt>z</tt>-axis (when using the first constructor of this
 * class), or an arbitrarily oriented cylinder described by the direction of
 * its axis and a point located on the axis.
 *
 * This class was developed to be used in conjunction with the @p cylinder or
 * @p cylinder_shell functions of GridGenerator. This function will throw a
 * run time exception whenever spacedim is not equal to three.
 *
 * @ingroup manifold
 */
template <int dim, int spacedim = dim>
class CylindricalManifold : public ChartManifold<dim, spacedim, 3>
{
public:
  /**
   * Constructor. Using default values for the constructor arguments yields a
   * cylinder along the x-axis (<tt>axis=0</tt>). Choose <tt>axis=1</tt> or
   * <tt>axis=2</tt> for a tube along the y- or z-axis, respectively. The
   * tolerance value is used to determine if a point is on the axis.
   */
  CylindricalManifold(const unsigned int axis      = 0,
                      const double       tolerance = 1e-10);

  /**
   * Constructor. If constructed with this constructor, the manifold described
   * is a cylinder with an axis that points in direction #direction and goes
   * through the given #point_on_axis. The direction may be arbitrarily
   * scaled, and the given point may be any point on the axis. The tolerance
   * value is used to determine if a point is on the axis.
   */
  CylindricalManifold(const Tensor<1, spacedim> &direction,
                      const Point<spacedim>     &point_on_axis,
                      const double               tolerance = 1e-10);

  /**
   * Make a clone of this Manifold object.
   */
  virtual std::unique_ptr<Manifold<dim, spacedim>>
  clone() const override;

  /**
   * Compute the cylindrical coordinates $(r, \phi, \lambda)$ for the given
   * space point where $r$ denotes the distance from the axis,
   * $\phi$ the angle between the given point and the computed normal
   * direction, and $\lambda$ the axial position.
   */
  virtual Point<3>
  pull_back(const Point<spacedim> &space_point) const override;

  /**
   * Compute the Cartesian coordinates for a chart point given in cylindrical
   * coordinates $(r, \phi, \lambda)$, where $r$ denotes the distance from the
   * axis, $\phi$ the angle between the given point and the computed normal
   * direction, and $\lambda$ the axial position.
   */
  virtual Point<spacedim>
  push_forward(const Point<3> &chart_point) const override;

  /**
   * Compute the derivatives of the mapping from cylindrical coordinates
   * $(r, \phi, \lambda)$ to cartesian coordinates where $r$ denotes the
   * distance from the axis, $\phi$ the angle between the given point and the
   * computed normal direction, and $\lambda$ the axial position.
   */
  virtual DerivativeForm<1, 3, spacedim>
  push_forward_gradient(const Point<3> &chart_point) const override;

  /**
   * Compute new points on the CylindricalManifold. See the documentation of
   * the base class for a detailed description of what this function does.
   */
  virtual Point<spacedim>
  get_new_point(const ArrayView<const Point<spacedim>> &surrounding_points,
                const ArrayView<const double>          &weights) const override;

  /**
   * Get the Tensor parallel to the cylinder's axis.
   */
  const Tensor<1, spacedim> &
  get_direction() const;

  /**
   * Get a point on the Cylinder's axis.
   *
   * @note This argument, like get_direction() and get_tolerance(), just returns
   * the arguments set in the three-argument constructor.
   */
  const Point<spacedim> &
  get_point_on_axis() const;

  /**
   * Get the tolerance which determines if a point is on the Cylinder's axis.
   */
  double
  get_tolerance() const;

private:
  /**
   * A vector orthogonal to the normal direction.
   */
  const Tensor<1, spacedim> normal_direction;

  /**
   * The direction vector of the axis.
   */
  const Tensor<1, spacedim> direction;

  /**
   * An arbitrary point on the axis.
   */
  const Point<spacedim> point_on_axis;

  /**
   * Relative tolerance to measure zero distances.
   */
  const double tolerance;

  /**
   * The direction vector perpendicular to both direction and normal_direction.
   */
  const Tensor<1, spacedim> dxn;
};

/**
 * Elliptical manifold description derived from ChartManifold.
 * More information on the elliptical coordinate system can be found
 * at <a
 * href="https://en.wikipedia.org/wiki/Elliptic_coordinate_system">Wikipedia
 * </a>.
 *
 * This is based on the definition of elliptic coordinates $(u,v)$
 * @f[
 *  \left\lbrace\begin{aligned}
 *  x &=  x_0 + c \cosh(u) \cos(v) \\
 *  y &=  y_0 + c \sinh(u) \sin(v)
 *  \end{aligned}\right.
 * @f]
 * in which $(x_0,y_0)$ are coordinates of the center of the cartesian system.
 *
 * The current implementation uses coordinates $(c,v)$, instead of $(u,v)$, and
 * fixes $u$ according to a given eccentricity. Therefore, this choice
 * of coordinates generates an elliptical manifold characterized by a constant
 * eccentricity: $e=\frac{1}{\cosh(u)}$, with $e\in\left]0,1\right[$.
 *
 * The constructor of this class will throw an exception if both dim and
 * spacedim are different from two.
 *
 * This manifold can be used to produce hyper_shells with elliptical curvature.
 * As an example, the test <B>elliptical_manifold_01</B> produces the following
 * triangulation:
 * @image html elliptical_hyper_shell.png
 *
 * @ingroup manifold
 */
template <int dim, int spacedim = dim>
class EllipticalManifold : public ChartManifold<dim, spacedim, spacedim>
{
public:
  /**
   * Constructor that takes the center of the manifold system, the direction of
   * the major axis, and the manifold eccentricity.
   * The default major axis is the <tt>x</tt>-axis. The manifold is rotated in
   * order to align the major axis to the direction specified in input.
   * @param center Center of the manifold.
   * @param major_axis_direction Direction of the major axis of the
   * manifold.
   * @param eccentricity Eccentricity of the
   * manifold $e\in\left]0,1\right[$.
   */
  EllipticalManifold(const Point<spacedim>     &center,
                     const Tensor<1, spacedim> &major_axis_direction,
                     const double               eccentricity);

  virtual std::unique_ptr<Manifold<dim, spacedim>>
  clone() const override;

  /**
   * @copydoc ChartManifold::pull_back()
   */
  virtual Point<spacedim>
  pull_back(const Point<spacedim> &space_point) const override;

  /**
   * @copydoc ChartManifold::push_forward()
   */
  virtual Point<spacedim>
  push_forward(const Point<spacedim> &chart_point) const override;

  /**
   * @copydoc ChartManifold::push_forward_gradient()
   */
  virtual DerivativeForm<1, spacedim, spacedim>
  push_forward_gradient(const Point<spacedim> &chart_point) const override;

  /**
   * Get the Tensor parallel to the cylinder's major axis.
   */
  const Tensor<1, spacedim> &
  get_major_axis_direction() const;

  /**
   * Return the center of the elliptical coordinate system.
   */
  const Point<spacedim> &
  get_center() const;

  /**
   * Return the ellipse's eccentricity.
   */
  double
  get_eccentricity() const;

private:
  /**
   * The direction vector of the major axis.
   */
  const Tensor<1, spacedim> direction;

  /**
   * The center of the manifold.
   */
  const Point<spacedim> center;

  /**
   * The eccentricity.
   */
  const double eccentricity;

  /**
   * Parameters deriving from the eccentricity of the manifold.
   */
  const double cosh_u;
  const double sinh_u;

  /**
   * @copydoc ChartManifold::get_periodicity()
   *
   * For $\text{dim}=2$ and $\text{spacedim}=2$, the first coordinate is
   * non-periodic, while the second coordinate has a periodicity of $2\pi$.
   */
  static Tensor<1, spacedim>
  get_periodicity();
};


/**
 * Manifold description derived from ChartManifold, based on explicit
 * Function<spacedim> and Function<chartdim> objects describing the
 * push_forward() and pull_back() functions.
 *
 * You can use this Manifold object to describe any arbitrary shape domain, as
 * long as you can express it in terms of an invertible map, for which you
 * provide both the forward expression, and the inverse expression.
 *
 * In debug mode, a check is performed to verify that the transformations are
 * actually one the inverse of the other.
 *
 * @ingroup manifold
 */
template <int dim, int spacedim = dim, int chartdim = dim>
class FunctionManifold : public ChartManifold<dim, spacedim, chartdim>
{
public:
  /**
   * Explicit functions constructor. Takes a push_forward function of spacedim
   * components, and a pull_back function of @p chartdim components. See the
   * documentation of the base class ChartManifold for the meaning of the
   * optional @p periodicity argument.
   *
   * The tolerance argument is used in debug mode to actually check that the
   * two functions are one the inverse of the other.
   *
   * Note: the object constructed in this way stores pointers to the
   * push_forward and  pull_back functions. Therefore, one must guarantee that
   * the function objects are destroyed only after the constructed manifold.
   */
  FunctionManifold(
    const Function<chartdim>  &push_forward_function,
    const Function<spacedim>  &pull_back_function,
    const Tensor<1, chartdim> &periodicity = Tensor<1, chartdim>(),
    const double               tolerance   = 1e-10);

  /**
   * Same as previous, except this constructor takes ownership of the Function
   * objects passed as first and second argument, and is ultimately in charge
   * of deleting the pointers when the FunctionManifold object is destroyed.
   *
   * This constructor is useful because it allows creating function objects at
   * the place of calling the constructor without having to name and later
   * delete these objects. This allows the following idiom:
   * FunctionManifold<dim> manifold(std::make_unique<MyPushForward>(...),
   *                                std::make_unique<MyPullBack>(...));
   */
  FunctionManifold(
    std::unique_ptr<Function<chartdim>> push_forward,
    std::unique_ptr<Function<spacedim>> pull_back,
    const Tensor<1, chartdim>          &periodicity = Tensor<1, chartdim>(),
    const double                        tolerance   = 1e-10);

  /**
   * Expressions constructor. Takes the expressions of the push_forward
   * function of spacedim components, and of the pull_back function of @p
   * chartdim components. See the documentation of the base class
   * ChartManifold for the meaning of the optional @p periodicity argument.
   *
   * The strings should be the readable by the default constructor of the
   * FunctionParser classes. You can specify custom variable expressions with
   * the last two optional arguments. If you don't, the default names are
   * used, i.e., "x,y,z".
   *
   * The tolerance argument is used in debug mode to actually check that the
   * two functions are one the inverse of the other.
   */
  FunctionManifold(
    const std::string          push_forward_expression,
    const std::string          pull_back_expression,
    const Tensor<1, chartdim> &periodicity = Tensor<1, chartdim>(),
    const typename FunctionParser<spacedim>::ConstMap =
      typename FunctionParser<spacedim>::ConstMap(),
    const std::string chart_vars =
      FunctionParser<chartdim>::default_variable_names(),
    const std::string space_vars =
      FunctionParser<spacedim>::default_variable_names(),
    const double tolerance = 1e-10,
    const double h         = 1e-8);

  /**
   * If needed, we delete the pointers we own.
   */
  virtual ~FunctionManifold() override;

  /**
   * Make a clone of this Manifold object.
   */
  virtual std::unique_ptr<Manifold<dim, spacedim>>
  clone() const override;

  /**
   * Given a point in the @p chartdim coordinate system, uses the
   * push_forward_function to compute the push_forward of points in @p
   * chartdim space dimensions to @p spacedim space dimensions.
   */
  virtual Point<spacedim>
  push_forward(const Point<chartdim> &chart_point) const override;

  /**
   * Given a point in the chartdim dimensional Euclidean space, this
   * method returns the derivatives of the function $F$ that maps from
   * the sub_manifold coordinate system to the Euclidean coordinate
   * system. In other words, it is a matrix of size
   * $\text{spacedim}\times\text{chartdim}$.
   *
   * This function is used in the computations required by the
   * get_tangent_vector() function. The default implementation calls
   * the get_gradient() method of the
   * FunctionManifold::push_forward_function() member class. If you
   * construct this object using the constructor that takes two string
   * expression, then the default implementation of this method uses a
   * finite difference scheme to compute the gradients(see the
   * AutoDerivativeFunction() class for details), and you can specify
   * the size of the spatial step size at construction time with the
   * @p h parameter.
   *
   * Refer to the general documentation of this class for more information.
   */
  virtual DerivativeForm<1, chartdim, spacedim>
  push_forward_gradient(const Point<chartdim> &chart_point) const override;

  /**
   * Given a point in the spacedim coordinate system, uses the
   * pull_back_function to compute the pull_back of points in @p spacedim
   * space dimensions to @p chartdim space dimensions.
   */
  virtual Point<chartdim>
  pull_back(const Point<spacedim> &space_point) const override;

private:
  /**
   * Constants for the FunctionParser classes.
   */
  const typename FunctionParser<spacedim>::ConstMap const_map;

  /**
   * Pointer to the push_forward function.
   */
  ObserverPointer<const Function<chartdim>,
                  FunctionManifold<dim, spacedim, chartdim>>
    push_forward_function;

  /**
   * Pointer to the pull_back function.
   */
  ObserverPointer<const Function<spacedim>,
                  FunctionManifold<dim, spacedim, chartdim>>
    pull_back_function;

  /**
   * Relative tolerance. In debug mode, we check that the two functions
   * provided at construction time are actually one the inverse of the other.
   * This value is used as relative tolerance in this check.
   */
  const double tolerance;

  /**
   * Check ownership of the smart pointers. Indicates whether this class is
   * the owner of the objects pointed to by the previous two member variables.
   * This value is set in the constructor of the class. If @p true, then the
   * destructor will delete the function objects pointed to be the two
   * pointers.
   */
  bool owns_pointers;

  /**
   * The expression used to construct the push_forward function.
   */
  const std::string push_forward_expression;

  /**
   * The expression used to construct the pull_back function.
   */
  const std::string pull_back_expression;

  /**
   * Variable names in the chart domain.
   */
  const std::string chart_vars;

  /**
   * Variable names in the space domain.
   */
  const std::string space_vars;

  /**
   * The finite difference step to use internally.
   */
  const double finite_difference_step;
};



/**
 * Manifold description for the surface of a Torus in three dimensions. The
 * Torus is assumed to be in the x-z plane. The reference coordinate system
 * is given by the angle $phi$ around the y axis, the angle $theta$ around
 * the centerline of the torus, and the distance to the centerline $w$
 * (between 0 and 1).
 *
 * This class was developed to be used in conjunction with
 * GridGenerator::torus.
 *
 * @ingroup manifold
 */
template <int dim>
class TorusManifold : public ChartManifold<dim, 3, 3>
{
public:
  static const int chartdim = 3;
  static const int spacedim = 3;

  /**
   * Constructor. Specify the radius of the centerline @p centerline_radius and
   * the radius of the torus' inner circle (@p inner_radius). The variables have
   * the same meaning as the parameters in GridGenerator::torus().
   */
  TorusManifold(const double centerline_radius, const double inner_radius);

  /**
   * Make a clone of this Manifold object.
   */
  virtual std::unique_ptr<Manifold<dim, 3>>
  clone() const override;

  /**
   * Pull back operation.
   */
  virtual Point<3>
  pull_back(const Point<3> &p) const override;

  /**
   * Push forward operation.
   */
  virtual Point<3>
  push_forward(const Point<3> &chart_point) const override;

  /**
   * Gradient.
   */
  virtual DerivativeForm<1, 3, 3>
  push_forward_gradient(const Point<3> &chart_point) const override;

  /**
   * Get the radius of the centerline.
   */
  double
  get_centerline_radius() const;

  /**
   * Get the inner radius of the torus.
   */
  double
  get_inner_radius() const;

private:
  double centerline_radius;

  double inner_radius;
};



/**
 * A mapping class that extends curved boundary descriptions into the interior
 * of the computational domain. The outer curved boundary description is
 * assumed to be given by another manifold (e.g. a polar manifold on a circle).
 * The mechanism to extend the boundary information is a so-called transfinite
 * interpolation.
 * The use of this class is discussed extensively in step-65.
 *
 * The formula for extending such a description in 2d is, for example,
 * described on
 * <a href="https://en.wikipedia.org/wiki/Transfinite_interpolation">
 * Wikipedia</a>.  Given a point $(u,v)$ on the chart, the image of this point
 * in real space is given by
 * @f{align*}{
 * \mathbf S(u,v) &= (1-v)\mathbf c_0(u)+v \mathbf c_1(u) + (1-u)\mathbf c_2(v)
 * + u \mathbf c_3(v) \\
 * &\quad - \left[(1-u)(1-v) \mathbf x_0 + u(1-v) \mathbf x_1 + (1-u)v \mathbf
 * x_2 + uv \mathbf x_3 \right]
 * @f}
 * where $\bf x_0, \bf x_1, \bf x_2, \bf x_3$ denote the four bounding vertices
 * bounding the image space and $\bf c_0, \bf c_1, \bf c_2, \bf c_3$ are the
 * four curves describing the lines of the cell. If a curved manifold is
 * attached to any of these lines, the evaluation is done according to
 * Manifold::get_new_point() with the two end points of the line and
 * appropriate weight. In 3d, the generalization of this formula is
 * implemented, creating a weighted sum of the vertices (positive
 * contribution), the lines (negative), and the faces (positive contribution).
 *
 * This manifold is usually attached to a coarse mesh and then places new
 * points as a combination of the descriptions on the boundaries, weighted
 * appropriately according to the position of the point in the original chart
 * coordinates $(u,v)$. This manifold should be preferred over setting only a
 * curved manifold on the boundary of a mesh in most situations as it yields
 * more uniform mesh distributions as the mesh is refined because it switches
 * from a curved description to a straight description over all children of
 * the initial coarse cell this manifold was attached to. This way, the curved
 * nature of the manifold that is originally contained in one <i>coarse</i>
 * mesh layer will be applied to more than one <i>fine</i> mesh layer once the
 * mesh gets refined. Note that the mechanisms of
 * TransfiniteInterpolationManifold are also built into the MappingQ
 * class when only a surface of a cell is subject to a curved description,
 * ensuring that even the default case without this manifold gets optimal
 * convergence rates when applying curved boundary descriptions.
 *
 * If no curved boundaries surround a coarse cell, this class reduces to a flat
 * manifold description.
 *
 * To give an example of using this class, the following code attaches a
 * transfinite manifold to a circle:
 *
 * @code
 * PolarManifold<dim> polar_manifold;
 * TransfiniteInterpolationManifold<dim> inner_manifold;
 *
 * Triangulation<dim> triangulation;
 * GridGenerator::hyper_ball (triangulation);
 *
 * triangulation.set_all_manifold_ids(1);
 * triangulation.set_all_manifold_ids_on_boundary(0);
 * triangulation.set_manifold (0, polar_manifold);
 * inner_manifold.initialize(triangulation);
 * triangulation.set_manifold (1, inner_manifold);
 * triangulation.refine_global(4);
 * @endcode
 *
 * In this code, we first set all manifold ids to the id of the transfinite
 * interpolation, and then re-set the manifold ids on the boundary to identify
 * the curved boundary described by the polar manifold. With this code, one
 * gets a really nice mesh:
 *
 * <p ALIGN="center">
 * @image html circular_mesh_transfinite_interpolation.png
 * </p>
 *
 * which is obviously much nicer than the polar manifold applied to just the
 * boundary:
 *
 * <p ALIGN="center">
 * @image html circular_mesh_only_boundary_manifold.png
 * </p>
 *
 * This manifold is used in a few GridGenerator functions, including
 * GridGenerator::channel_with_cylinder.
 *
 * <h3>Implementation details</h3>
 *
 * In the implementation of this class, the manifolds surrounding a coarse
 * cell are queried repeatedly to compute points on their interior. For
 * optimal mesh quality, those manifolds should be compatible with a chart
 * notion. For example, computing a point that is 0.25 along the line between
 * two vertices using the weights 0.25 and 0.75 for the two vertices should
 * give the same result as first computing the mid point at 0.5 and then again
 * compute the midpoint between the first vertex and coarse mid point. This is
 * the case for most of the manifold classes provided by deal.II, such as
 * SphericalManifold or PolarManifold, but it might be violated by naive
 * implementations. In case the quality of the manifold is not good enough,
 * upon mesh refinement it may happen that the transformation to a chart
 * inside the get_new_point() or get_new_points() methods produces points that
 * are outside the unit cell. Then this class throws an exception of type
 * Mapping::ExcTransformationFailed. In that case, the mesh should be refined
 * before attaching this class, as done in the following example:
 *
 * @code
 * SphericalManifold<dim> spherical_manifold;
 * TransfiniteInterpolationManifold<dim> inner_manifold;
 * Triangulation<dim> triangulation;
 * GridGenerator::hyper_ball (triangulation);
 *
 * triangulation.set_all_manifold_ids(1);
 * triangulation.set_all_manifold_ids_on_boundary(0);
 * triangulation.set_manifold (0, spherical_manifold);
 * inner_manifold.initialize(triangulation);
 * triangulation.set_manifold (1, inner_manifold);
 * triangulation.refine_global(1);
 *
 * // initialize the transfinite manifold again
 * inner_manifold.initialize(triangulation);
 * triangulation.refine_global(4);
 * @endcode
 *
 * @note For performance and accuracy reasons, it is recommended to apply the
 * transfinite manifold to as coarse a mesh as possible. Regarding accuracy,
 * the curved description can only be applied to new points created from a
 * given neighborhood, and the grid quality is typically higher when extending
 * the curved description over as large a domain as possible. Regarding
 * performance, the identification of the correct coarse cell in the
 * get_new_point() method needs to pass all coarse cells, so expect a linear
 * complexity in the number of coarse cells for each single mapping operation,
 * i.e., at least quadratic in the number of coarse mesh cells for any global
 * operation on the whole mesh. Thus, the current implementation is only
 * economical when there are not more than a few hundreds of coarse cells. To
 * make performance better for larger numbers of cells, one could extend the
 * current implementation by a pre-identification of relevant cells with
 * axis-aligned bounding boxes.
 *
 * @ingroup manifold
 */
template <int dim, int spacedim = dim>
class TransfiniteInterpolationManifold : public Manifold<dim, spacedim>
{
public:
  /**
   * Constructor.
   */
  TransfiniteInterpolationManifold();

  /**
   * Destructor.
   */
  virtual ~TransfiniteInterpolationManifold() override;

  /**
   * Make a clone of this Manifold object.
   */
  virtual std::unique_ptr<Manifold<dim, spacedim>>
  clone() const override;

  /**
   * Initializes the manifold with a coarse mesh. The prerequisite for using
   * this class is that the input triangulation is uniformly refined and the
   * manifold is later attached to the same triangulation.
   *
   * Whenever the assignment of manifold ids changes on the level of the
   * triangulation which this class was initialized with, initialize() must be
   * called again to update the manifold ids connected to the coarse cells.
   *
   * @note The triangulation used to construct the manifold must not be
   * destroyed during the usage of this object.
   */
  void
  initialize(const Triangulation<dim, spacedim> &triangulation);

  /**
   * Return the point which shall become the new vertex surrounded by the
   * given points @p surrounding_points. @p weights contains appropriate
   * weights for the surrounding points according to which the manifold
   * determines the new point's position.
   *
   * The implementation in this class overrides the method in the base class
   * and computes the new point by a transfinite interpolation. The first step
   * in the implementation is to identify the coarse cell on which the
   * surrounding points are located. Then, the coordinates are transformed to
   * the unit coordinates on the coarse cell by a Newton iteration, where the
   * new point is then computed according to the weights. Finally, it is
   * pushed forward to the real space according to the transfinite
   * interpolation.
   */
  virtual Point<spacedim>
  get_new_point(const ArrayView<const Point<spacedim>> &surrounding_points,
                const ArrayView<const double>          &weights) const override;

  /**
   * Compute a new set of points that interpolate between the given points @p
   * surrounding_points. @p weights is a table with as many columns as @p
   * surrounding_points.size(). The number of columns in @p weights must match
   * the length of @p new_points.
   *
   * The implementation in this class overrides the method in the base class
   * and computes the new point by a transfinite interpolation. The first step
   * in the implementation is to identify the coarse cell on which the
   * surrounding points are located. Then, the coordinates are transformed to
   * the unit coordinates on the coarse cell by a Newton iteration, where the
   * new points are then computed according to the weights. Finally, the is
   * pushed forward to the real space according to the transfinite
   * interpolation.
   *
   * The implementation does not allow for @p surrounding_points and
   * @p new_points to point to the same vector, so make sure to pass different
   * objects into the function.
   */
  virtual void
  get_new_points(const ArrayView<const Point<spacedim>> &surrounding_points,
                 const Table<2, double>                 &weights,
                 ArrayView<Point<spacedim>> new_points) const override;

private:
  /**
   * Internal function to identify the most suitable cells (=charts) where the
   * given surrounding points are located. We use a cheap algorithm to
   * identify the cells and rank the cells by probability before we actually
   * do the search inside the relevant cells. The cells are sorted by the
   * distance of a Q1 approximation of the inverse mapping to the unit cell of
   * the surrounding points. We expect at most 20 cells (it should be up to 8
   * candidates on a 3d structured mesh and a bit more on unstructured ones,
   * typically we only get two or three), so get an array with 20 entries of a
   * the indices <tt>cell->index()</tt>.
   */
  std::array<unsigned int, 20>
  get_possible_cells_around_points(
    const ArrayView<const Point<spacedim>> &surrounding_points) const;

  /**
   * Finalizes the identification of the correct chart and populates @p
   * chart_points with the pullbacks of the surrounding points. This method
   * internally calls @p get_possible_cells_around_points().
   *
   * Return an iterator to the cell on which the chart is defined.
   */
  typename Triangulation<dim, spacedim>::cell_iterator
  compute_chart_points(
    const ArrayView<const Point<spacedim>> &surrounding_points,
    ArrayView<Point<dim>>                   chart_points) const;

  /**
   * Pull back operation into the unit coordinates on the given coarse cell.
   *
   * This method is currently based on a Newton-like iteration to find the
   * point in the origin. One may speed up the iteration by providing a good
   * initial guess as the third argument. If no better point is known, use
   * cell->real_to_unit_cell_affine_approximation(p)
   *
   * @note This internal function is currently not compatible with the
   * ChartManifold::pull_back() function because the given class represents an
   * atlas of charts, not a single chart. Thus, the pull_back() operation is
   * only valid with the additional information of the chart, given by a cell
   * on the coarse grid. An alternative implementation could shift the index
   * depending on the coarse cell for a 1-to-1 relation between the chart space
   * and the image space.
   */
  Point<dim>
  pull_back(const typename Triangulation<dim, spacedim>::cell_iterator &cell,
            const Point<spacedim>                                      &p,
            const Point<dim> &initial_guess) const;

  /**
   * Push forward operation.
   *
   * @note This internal function is currently not compatible with the
   * ChartManifold::push_forward() function because the given class represents
   * an atlas of charts, not a single chart. Thus, the push_forward()
   * operation is only valid with the additional information of the chart,
   * given by a cell on the coarse grid. An alternative implementation could
   * shift the index depending on the coarse cell for a 1-to-1 relation
   * between the chart space and the image space.
   */
  Point<spacedim>
  push_forward(const typename Triangulation<dim, spacedim>::cell_iterator &cell,
               const Point<dim> &chart_point) const;

  /**
   * Gradient of the push_forward method.
   *
   * @note This internal function is not compatible with the
   * ChartManifold::push_forward_gradient() function because the given class
   * represents an atlas of charts, not a single chart. Furthermore, this
   * private function also requires the user to provide the result of the
   * push_forward() call on the chart point for the single use case of this
   * function, namely inside a Newton iteration where the gradient is computed
   * by finite differences.
   */
  DerivativeForm<1, dim, spacedim>
  push_forward_gradient(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const Point<dim>                                           &chart_point,
    const Point<spacedim> &pushed_forward_chart_point) const;

  /**
   * The underlying triangulation.
   */
  const Triangulation<dim, spacedim> *triangulation;

  /**
   * The level of the mesh cells where the transfinite approximation is
   * applied, usually level 0.
   */
  int level_coarse;

  /**
   * In case there all surrounding manifolds are the transfinite manifold or
   * have default (invalid) manifold id, the manifold degenerates to a flat
   * manifold and we can choose cheaper algorithms for the push_forward method.
   */
  std::vector<bool> coarse_cell_is_flat;

  /**
   * A flat manifold used to compute new points in the chart space where we
   * use a FlatManifold description.
   */
  FlatManifold<dim> chart_manifold;

  /**
   * A vector of quadratic approximations to the inverse map from real points
   * to chart points for each of the coarse mesh cells.
   */
  std::vector<internal::MappingQImplementation::
                InverseQuadraticApproximation<dim, spacedim>>
    quadratic_approximation;

  /**
   * The connection to Triangulation::signals::clear that must be reset once
   * this class goes out of scope.
   */
  boost::signals2::connection clear_signal;
};

/*----------------------------- inline functions -----------------------------*/

template <int dim, int spacedim>
inline const Point<spacedim> &
PolarManifold<dim, spacedim>::get_center() const
{
  return p_center;
}



template <int dim, int spacedim>
inline const Point<spacedim> &
SphericalManifold<dim, spacedim>::get_center() const
{
  return p_center;
}



template <int dim, int spacedim>
inline const Tensor<1, spacedim> &
CylindricalManifold<dim, spacedim>::get_direction() const
{
  return direction;
}



template <int dim, int spacedim>
inline const Point<spacedim> &
CylindricalManifold<dim, spacedim>::get_point_on_axis() const
{
  return point_on_axis;
}



template <int dim, int spacedim>
inline double
CylindricalManifold<dim, spacedim>::get_tolerance() const
{
  return tolerance;
}



template <int dim, int spacedim>
inline const Tensor<1, spacedim> &
EllipticalManifold<dim, spacedim>::get_major_axis_direction() const
{
  return direction;
}



template <int dim, int spacedim>
inline const Point<spacedim> &
EllipticalManifold<dim, spacedim>::get_center() const
{
  return center;
}



template <int dim, int spacedim>
inline double
EllipticalManifold<dim, spacedim>::get_eccentricity() const
{
  return eccentricity;
}



template <int dim>
inline double
TorusManifold<dim>::get_centerline_radius() const
{
  return centerline_radius;
}



template <int dim>
inline double
TorusManifold<dim>::get_inner_radius() const
{
  return inner_radius;
}

DEAL_II_NAMESPACE_CLOSE

#endif
