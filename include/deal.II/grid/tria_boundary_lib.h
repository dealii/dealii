// ---------------------------------------------------------------------
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

#ifndef __deal2__tria_boundary_lib_h
#define __deal2__tria_boundary_lib_h


#include <deal.II/base/config.h>
#include <deal.II/grid/tria_boundary.h>

DEAL_II_NAMESPACE_OPEN

/**
 * Boundary object for the hull of a cylinder.  In three dimensions,
 * points are projected on a circular tube along the <tt>x-</tt>,
 * <tt>y-</tt> or <tt>z</tt>-axis (when using the first constructor of
 * this class), or an arbitrarily oriented cylinder described by the
 * direction of its axis and a point located on the axis. The radius
 * of the tube can be given independently. Similar to
 * HyperBallBoundary, new points are projected by dividing the
 * straight line between the old two points and adjusting the radius
 * from the axis.
 *
 * This class was developed to be used in conjunction with the
 * @p cylinder function of GridGenerator. It should be used for
 * the hull of the cylinder only (boundary indicator 0). Its use is
 * discussed in detail in the results section of step-49.
 *
 * This class is derived from StraightBoundary rather than from
 * Boundary, which would seem natural, since this way we can use the
 * StraightBoundary::in_between() function.
 *
 * @ingroup boundary
 *
 * @author Guido Kanschat, 2001, Wolfgang Bangerth, 2007
 */
template <int dim, int spacedim = dim>
class CylinderBoundary : public StraightBoundary<dim,spacedim>
{
public:
  /**
   * Constructor. Using default values for the constructor arguments yields a
   * circular tube along the x-axis (<tt>axis=0</tt>). Choose <tt>axis=1</tt>
   * or <tt>axis=2</tt> for a tube along the y- or z-axis, respectively.
   */
  CylinderBoundary (const double radius = 1.0,
                    const unsigned int axis = 0);

  /**
   * Constructor. If constructed
   * with this constructor, the
   * boundary described is a
   * cylinder with an axis that
   * points in direction #direction
   * and goes through the given
   * #point_on_axis. The direction
   * may be arbitrarily scaled, and
   * the given point may be any
   * point on the axis.
   */
  CylinderBoundary (const double           radius,
                    const Point<spacedim> &direction,
                    const Point<spacedim> &point_on_axis);

  /**
   * Refer to the general documentation of
   * this class and the documentation of the
   * base class.
   */
  virtual Point<spacedim>
  get_new_point_on_line (const typename Triangulation<dim,spacedim>::line_iterator &line) const;

  /**
   * Refer to the general documentation of
   * this class and the documentation of the
   * base class.
   */
  virtual Point<spacedim>
  get_new_point_on_quad (const typename Triangulation<dim,spacedim>::quad_iterator &quad) const;

  /**
   * Refer to the general
   * documentation of this class
   * and the documentation of the
   * base class.
   *
   * Calls
   * @p get_intermediate_points_between_points.
   */
  virtual void
  get_intermediate_points_on_line (const typename Triangulation<dim,spacedim>::line_iterator &line,
                                   std::vector<Point<spacedim> > &points) const;

  /**
   * Refer to the general
   * documentation of this class
   * and the documentation of the
   * base class.
   *
   * Only implemented for <tt>dim=3</tt>
   * and for <tt>points.size()==1</tt>.
   */
  virtual void
  get_intermediate_points_on_quad (const typename Triangulation<dim,spacedim>::quad_iterator &quad,
                                   std::vector<Point<spacedim> > &points) const;

  /**
   * Compute the normals to the
   * boundary at the vertices of
   * the given face.
   *
   * Refer to the general
   * documentation of this class
   * and the documentation of the
   * base class.
   */
  virtual void
  get_normals_at_vertices (const typename Triangulation<dim,spacedim>::face_iterator &face,
                           typename Boundary<dim,spacedim>::FaceVertexNormals &face_vertex_normals) const;

  /**
   * Return the radius of the cylinder.
   */
  double get_radius () const;

  /**
   * Exception. Thrown by the
   * @p get_radius if the
   * @p compute_radius_automatically,
   * see below, flag is set true.
   */
  DeclException0 (ExcRadiusNotSet);


protected:
  /**
   * Radius of the cylinder.
   */
  const double radius;

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
   * Called by
   * @p get_intermediate_points_on_line
   * and by
   * @p get_intermediate_points_on_quad.
   *
   * Refer to the general
   * documentation of
   * @p get_intermediate_points_on_line
   * in the documentation of the
   * base class.
   */
  void get_intermediate_points_between_points (const Point<spacedim> &p0, const Point<spacedim> &p1,
                                               std::vector<Point<spacedim> > &points) const;

  /**
   * Given a number for the axis,
   * return a vector that denotes
   * this direction.
   */
  static Point<spacedim> get_axis_vector (const unsigned int axis);
};



/**
 * Boundary object for the hull of a (truncated) cone with two
 * different radii at the two ends. If one radius is chosen to be 0
 * the object describes the boundary of a cone. In three dimensions,
 * points are projected on an arbitrarily oriented (truncated) cone
 * described by the two endpoints and the corresponding radii. Similar
 * to HyperBallBoundary, new points are projected by dividing the
 * straight line between the old two points and adjusting the radius
 * from the axis.
 *
 * This class is derived from StraightBoundary rather than from
 * Boundary, which would seem natural, since this way we can use the
 * StraightBoundary::in_between() function.
 *
 * As an example of use, consider the following code snippet:
 * @code
 *  Triangulation<dim> triangulation;
 *  GridGenerator::truncated_cone (triangulation);
 *  Point<dim> p1, p2;
 *  p1[0] = -1;
 *  p2[0] = 1;
 *  const ConeBoundary<dim> boundary (1, 0.5, p1, p2);
 *  triangulation.set_boundary (0, boundary);
 *  triangulation.refine_global (2);
 * @endcode
 * This will produce the following meshes after the two
 * refinements we perform, in 2d and 3d, respectively:
 *
 * @image html cone_2d.png
 * @image html cone_3d.png
 *
 * @author Markus B&uuml;rg, 2009
 */
template <int dim>
class ConeBoundary : public StraightBoundary<dim>
{
public:
  /**
   * Constructor. Here the boundary
   * object is constructed. The
   * points <tt>x_0</tt> and
   * <tt>x_1</tt> describe the
   * starting and ending points of
   * the axis of the (truncated)
   * cone. <tt>radius_0</tt>
   * denotes the radius
   * corresponding to <tt>x_0</tt>
   * and <tt>radius_1</tt> the one
   * corresponding to <tt>x_1</tt>.
   */
  ConeBoundary (const double radius_0,
                const double radius_1,
                const Point<dim> x_0,
                const Point<dim> x_1);

  /**
   * Return the radius of the
   * (truncated) cone at given point
   * <tt>x</tt> on the axis.
   */
  double get_radius (const Point<dim> x) const;

  /**
   * Refer to the general
   * documentation of this class
   * and the documentation of the
   * base class.
   */
  virtual
  Point<dim>
  get_new_point_on_line (const typename Triangulation<dim>::line_iterator &line) const;

  /**
   * Refer to the general
   * documentation of this class
   * and the documentation of the
   * base class.
   */
  virtual
  Point<dim>
  get_new_point_on_quad (const typename Triangulation<dim>::quad_iterator &quad) const;

  /**
   * Refer to the general
   * documentation of this class
   * and the documentation of the
   * base class.
   *
   * Calls @p
   * get_intermediate_points_between_points.
   */
  virtual
  void
  get_intermediate_points_on_line (const typename Triangulation<dim>::line_iterator &line,
                                   std::vector<Point<dim> > &points) const;

  /**
   * Refer to the general
   * documentation of this class
   * and the documentation of the
   * base class.
   *
   * Only implemented for
   * <tt>dim=3</tt> and for
   * <tt>points.size()==1</tt>.
   */
  virtual
  void
  get_intermediate_points_on_quad (const typename Triangulation<dim>::quad_iterator &quad,
                                   std::vector<Point<dim> > &points) const;

  /**
   * Compute the normals to the
   * boundary at the vertices of
   * the given face.
   *
   * Refer to the general
   * documentation of this class
   * and the documentation of the
   * base class.
   */
  virtual
  void
  get_normals_at_vertices (const typename Triangulation<dim>::face_iterator &face,
                           typename Boundary<dim>::FaceVertexNormals &face_vertex_normals) const;

protected:
  /**
   * First radius of the (truncated)
   * cone.
   */
  const double radius_0;

  /**
   * Second radius of the (truncated)
   * cone.
   */
  const double radius_1;

  /**
   * Starting point of the axis.
   */
  const Point<dim> x_0;

  /**
   * Ending point of the axis.
   */
  const Point<dim> x_1;

private:
  /**
   * Called by @p
   * get_intermediate_points_on_line
   * and by @p
   * get_intermediate_points_on_quad.
   *
   * Refer to the general
   * documentation of @p
   * get_intermediate_points_on_line
   * in the documentation of the
   * base class.
   */
  void
  get_intermediate_points_between_points (const Point<dim> &p0,
                                          const Point<dim> &p1,
                                          std::vector<Point<dim> > &points) const;
};



/**
 *   Specialization of Boundary<dim>, which places the new point on
 *   the boundary of a ball in arbitrary dimension. It works by projecting
 *   the point in the middle of the old points onto the ball. The middle is
 *   defined as the arithmetic mean of the points.
 *
 *   The center of the ball and its radius may be given upon construction of
 *   an object of this type. They default to the origin and a radius of 1.0.
 *
 *   This class is derived from StraightBoundary rather than from
 *   Boundary, which would seem natural, since this way we can use the
 *   StraightBoundary::in_between() function.
 *
 *   @ingroup boundary
 *
 *   @author Wolfgang Bangerth, 1998, Ralf Hartmann, 2001
 */
template <int dim, int spacedim=dim>
class HyperBallBoundary : public StraightBoundary<dim,spacedim>
{
public:
  /**
   * Constructor
   */
  HyperBallBoundary (const Point<spacedim> p      = Point<spacedim>(),
                     const double     radius = 1.0);

  /**
   * Refer to the general documentation of
   * this class and the documentation of the
   * base class.
   */
  virtual
  Point<spacedim>
  get_new_point_on_line (const typename Triangulation<dim,spacedim>::line_iterator &line) const;

  /**
   * Refer to the general documentation of
   * this class and the documentation of the
   * base class.
   */
  virtual
  Point<spacedim>
  get_new_point_on_quad (const typename Triangulation<dim,spacedim>::quad_iterator &quad) const;

  /**
   * Refer to the general
   * documentation of this class
   * and the documentation of the
   * base class.
   *
   * Calls
   * @p get_intermediate_points_between_points.
   */
  virtual
  void
  get_intermediate_points_on_line (const typename Triangulation<dim,spacedim>::line_iterator &line,
                                   std::vector<Point<spacedim> > &points) const;

  /**
   * Refer to the general
   * documentation of this class
   * and the documentation of the
   * base class.
   *
   * Only implemented for <tt>dim=3</tt>
   * and for <tt>points.size()==1</tt>.
   */
  virtual
  void
  get_intermediate_points_on_quad (const typename Triangulation<dim,spacedim>::quad_iterator &quad,
                                   std::vector<Point<spacedim> > &points) const;

  /**
   * Implementation of the function
   * declared in the base class.
   *
   * Refer to the general
   * documentation of this class
   * and the documentation of the
   * base class.
   */
  virtual
  Tensor<1,spacedim>
  normal_vector (const typename Triangulation<dim,spacedim>::face_iterator &face,
                 const Point<spacedim> &p) const;

  /**
   * Compute the normals to the
   * boundary at the vertices of
   * the given face.
   *
   * Refer to the general
   * documentation of this class
   * and the documentation of the
   * base class.
   */
  virtual
  void
  get_normals_at_vertices (const typename Triangulation<dim,spacedim>::face_iterator &face,
                           typename Boundary<dim,spacedim>::FaceVertexNormals &face_vertex_normals) const;

  /**
   * Return the center of the ball.
   */
  Point<spacedim>
  get_center () const;

  /**
   * Return the radius of the ball.
   */
  double
  get_radius () const;

  /**
   * Exception. Thrown by the
   * @p get_radius if the
   * @p compute_radius_automatically,
   * see below, flag is set true.
   */
  DeclException0 (ExcRadiusNotSet);


protected:

  /**
   * Center point of the hyperball.
   */
  const Point<spacedim> center;

  /**
   * Radius of the hyperball.
   */
  const double radius;

  /**
   * This flag is @p false for
   * this class and for all derived
   * classes that set the radius by
   * the constructor. For example
   * this flag is @p false for the
   * HalfHyperBallBoundary
   * class but it is @p true for
   * the HyperShellBoundary
   * class, for example.  The
   * latter class doesn't get its
   * radii by the constructor but
   * need to compute the radii
   * automatically each time one of
   * the virtual functions is
   * called.
   */
  bool compute_radius_automatically;

private:

  /**
   * Called by
   * @p get_intermediate_points_on_line
   * and by
   * @p get_intermediate_points_on_quad.
   *
   * Refer to the general
   * documentation of
   * @p get_intermediate_points_on_line
   * in the documentation of the
   * base class.
   */
  void get_intermediate_points_between_points (const Point<spacedim> &p0, const Point<spacedim> &p1,
                                               std::vector<Point<spacedim> > &points) const;
};



/**
 * Variant of HyperBallBoundary which denotes a half hyper ball
 * where the first coordinate is restricted to the range $x>=0$ (or
 * $x>=center(0)$). In two dimensions, this equals the right half
 * circle, in three space dimensions it is a half ball. This class
 * might be useful for computations with rotational symmetry, where
 * one dimension is the radius from the axis of rotation.
 *
 * @ingroup boundary
 *
 * @author Wolfgang Bangerth, 1999, 2001
 */
template <int dim>
class HalfHyperBallBoundary : public HyperBallBoundary<dim>
{
public:
  /**
   * Constructor
   */
  HalfHyperBallBoundary (const Point<dim> p      = Point<dim>(),
                         const double     radius = 1.0);

  /**
   * Check if on the line <tt>x==0</tt>,
   * otherwise pass to the base
   * class.
   */
  virtual Point<dim>
  get_new_point_on_line (const typename Triangulation<dim>::line_iterator &line) const;

  /**
   * Check if on the line <tt>x==0</tt>,
   * otherwise pass to the base
   * class.
   */
  virtual Point<dim>
  get_new_point_on_quad (const typename Triangulation<dim>::quad_iterator &quad) const;

  /**
   * Refer to the general
   * documentation of this class
   * and the documentation of the
   * base class.
   *
   * Calls
   * @p get_intermediate_points_between_points.
   */
  virtual void
  get_intermediate_points_on_line (const typename Triangulation<dim>::line_iterator &line,
                                   std::vector<Point<dim> > &points) const;

  /**
   * Refer to the general
   * documentation of this class
   * and the documentation of the
   * base class.
   *
   * Only implemented for <tt>dim=3</tt>
   * and for <tt>points.size()==1</tt>.
   */
  virtual void
  get_intermediate_points_on_quad (const typename Triangulation<dim>::quad_iterator &quad,
                                   std::vector<Point<dim> > &points) const;

  /**
   * Compute the normals to the
   * boundary at the vertices of
   * the given face.
   *
   * Refer to the general
   * documentation of this class
   * and the documentation of the
   * base class.
   */
  virtual void
  get_normals_at_vertices (const typename Triangulation<dim>::face_iterator &face,
                           typename Boundary<dim>::FaceVertexNormals &face_vertex_normals) const;
};



/**
 * Class describing the boundaries of a hyper shell. Only the center
 * of the two spheres needs to be given, the radii of inner and outer
 * sphere are computed automatically upon calling one of the virtual
 * functions.
 *
 * @ingroup boundary
 *
 * @author Wolfgang Bangerth, 1999
 */
template <int dim>
class HyperShellBoundary : public HyperBallBoundary<dim>
{
public:
  /**
   * Constructor. The center of the
   * spheres defaults to the
   * origin.
   *
   * Calls the constructor of its
   * base @p HyperBallBoundary
   * class with a dummy radius as
   * argument. This radius will be
   * ignored
   */
  HyperShellBoundary (const Point<dim> &center = Point<dim>());
};



/**
 * Variant of HyperShellBoundary which denotes a half hyper shell
 * where the first coordinate is restricted to the range $x>=0$ (or
 * $x>=center(0)$). In two dimensions, this equals the right half arc,
 * in three space dimensions it is a half shell. This class might be
 * useful for computations with rotational symmetry, where one
 * dimension is the radius from the axis of rotation.
 *
 * @ingroup boundary
 *
 * @author Wolfgang Bangerth, 2000, 2009
 */
template <int dim>
class HalfHyperShellBoundary : public HyperShellBoundary<dim>
{
public:
  /**
   * Constructor. The center of the
   * spheres defaults to the
   * origin.
   *
   * If the radii are not specified, the
   * class tries to infer them from the
   * location of points on the
   * boundary. This works in 2d, but not in
   * 3d. As a consequence, in 3d these
   * radii must be given.
   */
  HalfHyperShellBoundary (const Point<dim> &center = Point<dim>(),
                          const double inner_radius = -1,
                          const double outer_radius = -1);

  /**
   * Construct a new point on a line.
   */
  virtual Point<dim>
  get_new_point_on_line (const typename Triangulation<dim>::line_iterator &line) const;

  /**
   * Construct a new point on a quad.
   */
  virtual Point<dim>
  get_new_point_on_quad (const typename Triangulation<dim>::quad_iterator &quad) const;

  /**
   * Refer to the general
   * documentation of this class
   * and the documentation of the
   * base class.
   *
   * Calls
   * @p get_intermediate_points_between_points.
   */
  virtual void
  get_intermediate_points_on_line (const typename Triangulation<dim>::line_iterator &line,
                                   std::vector<Point<dim> > &points) const;

  /**
   * Refer to the general
   * documentation of this class
   * and the documentation of the
   * base class.
   *
   * Only implemented for <tt>dim=3</tt>
   * and for <tt>points.size()==1</tt>.
   */
  virtual void
  get_intermediate_points_on_quad (const typename Triangulation<dim>::quad_iterator &quad,
                                   std::vector<Point<dim> > &points) const;

  /**
   * Compute the normals to the
   * boundary at the vertices of
   * the given face.
   *
   * Refer to the general
   * documentation of this class
   * and the documentation of the
   * base class.
   */
  virtual void
  get_normals_at_vertices (const typename Triangulation<dim>::face_iterator &face,
                           typename Boundary<dim>::FaceVertexNormals &face_vertex_normals) const;

private:
  /**
   * Inner and outer radii of the shell.
   */
  const double inner_radius;
  const double outer_radius;
};


/**
 * Class describing the boundary of the torus. The axis of the torus is the
 * $y$-axis while the plane of the torus is the $x$-$z$ plane. A torus of this
 * kind can be generated by GridGenerator::torus.
 *
 * This class is only implemented for
 * <tt>dim=2</tt>,<tt>spacedim=3</tt>, that is, just the surface.
 */
template <int dim, int spacedim>
class TorusBoundary : public Boundary<dim,spacedim>
{
public:
  /**
   * Constructor.<tt>R</tt> has to be
   * greater than <tt>r</tt>.
   */
  TorusBoundary (const double R, const double r);

//Boundary Refinenment  Functions
  /**
   * Construct a new point on a line.
   */
  virtual Point<spacedim>
  get_new_point_on_line (const typename Triangulation<dim,spacedim>::line_iterator &line) const;

  /**
   * Construct a new point on a quad.
   */
  virtual Point<spacedim>
  get_new_point_on_quad (const typename Triangulation<dim,spacedim>::quad_iterator &quad) const;

  /**
   * Construct a new points on a line.
   */
  virtual void   get_intermediate_points_on_line (
    const typename Triangulation< dim, spacedim >::line_iterator   &line,
    std::vector< Point< spacedim > >         &points) const ;

  /**
   * Construct a new points on a quad.
   */
  virtual void  get_intermediate_points_on_quad (
    const typename Triangulation< dim, spacedim >::quad_iterator &quad,
    std::vector< Point< spacedim > >         &points ) const ;

  /**
   * Get the normal from cartesian
   * coordinates. This normal does not have
   * unit length.
   */
  virtual void get_normals_at_vertices (
    const typename Triangulation< dim, spacedim >::face_iterator &face,
    typename Boundary<dim,spacedim>::FaceVertexNormals &face_vertex_normals) const ;

private:
  //Handy functions
  /**
   * Function that corrects the value and
   * sign of angle, that is, given
   * <tt>angle=tan(abs(y/x))</tt>; if <tt>
   * (y > 0) && (x < 0) </tt> then
   * <tt>correct_angle = Pi - angle</tt>,
   * etc.
   */

  double           get_correct_angle(const double angle,const double x,const double y) const;

  /**
   * Get the cartesian coordinates of the
   * Torus, i.e., from <tt>(theta,phi)</tt>
   * to <tt>(x,y,z)</tt>.
   */
  Point<spacedim>  get_real_coord(const Point<dim> &surfP) const ;

  /**
   * Get the surface coordinates of the
   * Torus, i.e., from <tt>(x,y,z)</tt> to
   * <tt>(theta,phi)</tt>.
   */
  Point<dim>       get_surf_coord(const Point<spacedim> &p) const ;

  /**
   * Get the normal from surface
   * coordinates. This normal does not have
   * unit length.
   */
  Point<spacedim>  get_surf_norm_from_sp(const Point<dim> &surfP)      const ;

  /**
   * Get the normal from cartesian
   * coordinates. This normal does not have
   * unit length.
   */
  Point<spacedim>  get_surf_norm(const Point<spacedim> &p) const ;

  /**
   * Inner and outer radii of the shell.
   */
  const double R;
  const double r;
};



/* -------------- declaration of explicit specializations ------------- */

#ifndef DOXYGEN

template <>
Point<1>
HyperBallBoundary<1>::
get_new_point_on_quad (const Triangulation<1>::quad_iterator &) const;
template <>
void
HyperBallBoundary<1>::get_intermediate_points_on_line (
  const Triangulation<1>::line_iterator &,
  std::vector<Point<1> > &) const;
template <>
void
HyperBallBoundary<3>::get_intermediate_points_on_quad (
  const Triangulation<3>::quad_iterator &quad,
  std::vector<Point<3> > &points) const;
template <>
void
HyperBallBoundary<1>::
get_normals_at_vertices (const Triangulation<1,1>::face_iterator &,
                         Boundary<1,1>::FaceVertexNormals &) const;
template <>
Point<1>
HalfHyperBallBoundary<1>::
get_new_point_on_quad (const Triangulation<1>::quad_iterator &) const;
template <>
void
HalfHyperBallBoundary<1>::
get_intermediate_points_on_quad (const Triangulation<1>::quad_iterator &,
                                 std::vector<Point<1> > &) const;
template <>
void
HalfHyperBallBoundary<1>::
get_normals_at_vertices (const Triangulation<1,1>::face_iterator &,
                         Boundary<1,1>::FaceVertexNormals &) const;
template <>
Point<1>
HalfHyperShellBoundary<1>::
get_new_point_on_quad (const Triangulation<1>::quad_iterator &) const;
template <>
void
HalfHyperShellBoundary<1>::
get_intermediate_points_on_quad (const Triangulation<1>::quad_iterator &,
                                 std::vector<Point<1> > &) const;
template <>
void
HalfHyperShellBoundary<1>::
get_normals_at_vertices (const Triangulation<1,1>::face_iterator &,
                         Boundary<1,1>::FaceVertexNormals &) const;


#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif
