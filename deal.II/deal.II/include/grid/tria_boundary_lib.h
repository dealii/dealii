//----------------------------  tria_boundary_lib.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  tria_boundary_lib.h  ---------------------------
#ifndef __deal2__tria_boundary_lib_h
#define __deal2__tria_boundary_lib_h


#include <base/config.h>
#include <grid/tria_boundary.h>


/**
 * Boundary object for the hull of a cylinder.  In three dimensions,
 * points are projected on a circular tube along the @p{x}-axis. The
 * radius of the tube can be set. Similar to @ref{HyperBallBoundary},
 * new points are projected by dividing the straight line between the
 * old two points and adjusting the radius in the @p{yz}-plane.
 *
 * This class was developed to be used in conjunction with the
 * @p{cylinder} function of @ref{GridGenerator}. It should be used for
 * the hull of the cylinder only (boundary indicator 0).
 *
 *   This class is derived from @ref{StraightBoundary} rather than from
 *   @ref{Boundary}, which would seem natural, since this way we can use the
 *   @ref{StraightBoundary}@p{<dim>::in_between(neighbors)} function.
 *
 *   @author Guido Kanschat, 2001
 */
template <int dim>
class CylinderBoundary : public StraightBoundary<dim>
{
  public:
				     /**
				      * Constructor
				      */
    CylinderBoundary (const double     radius = 1.0);

				     /**
				      * Refer to the general documentation of
				      * this class and the documentation of the
				      * base class.
				      */
    virtual Point<dim>
    get_new_point_on_line (const typename Triangulation<dim>::line_iterator &line) const;

				     /**
				      * Refer to the general documentation of
				      * this class and the documentation of the
				      * base class.
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
				      * @p{get_intermediate_points_between_points}.
				      */
    virtual void
    get_intermediate_points_on_line (const typename Triangulation<dim>::line_iterator &line,
				     typename std::vector<Point<dim> > &points) const;

				     /**
				      * Refer to the general
				      * documentation of this class
				      * and the documentation of the
				      * base class.
				      *
				      * Only implemented for @p{dim=3}
				      * and for @p{points.size()==1}.
				      */
    virtual void
    get_intermediate_points_on_quad (const typename Triangulation<dim>::quad_iterator &quad,
				     typename std::vector<Point<dim> > &points) const;

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

				     /**
				      * Return the radius of the cylinder.
				      */
    double get_radius () const;

				     /**
				      * Exception. Thrown by the
				      * @p{get_radius} if the
				      * @p{compute_radius_automatically},
				      * see below, flag is set true.
				      */
    DeclException0 (ExcRadiusNotSet);
    
    
  protected:
				     /**
				      * Radius of the cylinder.
				      */
    const double radius;

  private:

				     /**
				      * Called by
				      * @p{get_intermediate_points_on_line}
				      * and by
				      * @p{get_intermediate_points_on_quad}.
				      *
				      * Refer to the general
				      * documentation of
				      * @p{get_intermediate_points_on_line}
				      * in the documentation of the
				      * base class.
				      */
    void get_intermediate_points_between_points (const Point<dim> &p0, const Point<dim> &p1,
						 typename std::vector<Point<dim> > &points) const;    
};



/**
 *   Specialisation of @ref{Boundary}<dim>, which places the new point on
 *   the boundary of a ball in arbitrary dimension. It works by projecting
 *   the point in the middle of the old points onto the ball. The middle is
 *   defined as the arithmetic mean of the points. 
 *
 *   The center of the ball and its radius may be given upon construction of
 *   an object of this type. They default to the origin and a radius of 1.0.
 *
 *   This class is derived from @ref{StraightBoundary} rather than from
 *   @ref{Boundary}, which would seem natural, since this way we can use the
 *   @ref{StraightBoundary}@p{<dim>::in_between(neighbors)} function.
 *
 *   @author Wolfgang Bangerth, 1998, Ralf Hartmann, 2001
 */
template <int dim>
class HyperBallBoundary : public StraightBoundary<dim>
{
  public:
				     /**
				      * Constructor
				      */
    HyperBallBoundary (const Point<dim> p      = Point<dim>(),
		       const double     radius = 1.0);

				     /**
				      * Refer to the general documentation of
				      * this class and the documentation of the
				      * base class.
				      */
    virtual Point<dim>
    get_new_point_on_line (const typename Triangulation<dim>::line_iterator &line) const;

				     /**
				      * Refer to the general documentation of
				      * this class and the documentation of the
				      * base class.
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
				      * @p{get_intermediate_points_between_points}.
				      */
    virtual void
    get_intermediate_points_on_line (const typename Triangulation<dim>::line_iterator &line,
				     typename std::vector<Point<dim> > &points) const;

				     /**
				      * Refer to the general
				      * documentation of this class
				      * and the documentation of the
				      * base class.
				      *
				      * Only implemented for @p{dim=3}
				      * and for @p{points.size()==1}.
				      */
    virtual void
    get_intermediate_points_on_quad (const typename Triangulation<dim>::quad_iterator &quad,
				     typename std::vector<Point<dim> > &points) const;

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

				     /**
				      * Return the center of the ball.
				      */
    Point<dim> get_center () const;

				     /**
				      * Return the radius of the ball.
				      */
    double get_radius () const;

				     /**
				      * Exception. Thrown by the
				      * @p{get_radius} if the
				      * @p{compute_radius_automatically},
				      * see below, flag is set true.
				      */
    DeclException0 (ExcRadiusNotSet);
    
    
  protected:
    
				     /**
				      * Center point of the hyperball.
				      */
    const Point<dim> center;

				     /**
				      * Radius of the hyperball.
				      */
    const double radius;

				     /**
				      * This flag is @p{false} for
				      * this class and for all derived
				      * classes that set the radius by
				      * the constructor. For example
				      * this flag is @p{false} for the
				      * @ref{HalfHyperBallBoundary}
				      * class but it is @p{true} for
				      * the @ref{HyperShellBoundary}
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
				      * @p{get_intermediate_points_on_line}
				      * and by
				      * @p{get_intermediate_points_on_quad}.
				      *
				      * Refer to the general
				      * documentation of
				      * @p{get_intermediate_points_on_line}
				      * in the documentation of the
				      * base class.
				      */
    void get_intermediate_points_between_points (const Point<dim> &p0, const Point<dim> &p1,
						 typename std::vector<Point<dim> > &points) const;    
};



/**
 * Variant of @ref{HyperBallBoundary} which denotes a half hyper ball
 * where the first coordinate is restricted to the range $x>=0$ (or
 * $x>=center(0)$). In two dimensions, this equals the right half
 * circle, in three space dimensions it is a half ball. This class
 * might be useful for computations with rotational symmetry, where
 * one dimension is the radius from the axis of rotation.
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
				      * Check if on the line @p{x==0},
				      * otherwise pass to the base
				      * class.
				      */
    virtual Point<dim>
    get_new_point_on_line (const typename Triangulation<dim>::line_iterator &line) const;

				     /**
				      * Check if on the line @p{x==0},
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
				      * @p{get_intermediate_points_between_points}.
				      */
    virtual void
    get_intermediate_points_on_line (const typename Triangulation<dim>::line_iterator &line,
				     typename std::vector<Point<dim> > &points) const;

				     /**
				      * Refer to the general
				      * documentation of this class
				      * and the documentation of the
				      * base class.
				      *
				      * Only implemented for @p{dim=3}
				      * and for @p{points.size()==1}.
				      */
    virtual void
    get_intermediate_points_on_quad (const typename Triangulation<dim>::quad_iterator &quad,
				     typename std::vector<Point<dim> > &points) const;

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
				      * base @p{HyperBallBoundary}
				      * class with a dummy radius as
				      * argument. This radius will be
				      * ignored
				      */
    HyperShellBoundary (const Point<dim> &center = Point<dim>());
};



/**
 * Variant of @ref{HyperShellBoundary} which denotes a half hyper shell
 * where the first coordinate is restricted to the range $x>=0$ (or
 * $x>=center(0)$). In two dimensions, this equals the right half arc,
 * in three space dimensions it is a half shell. This class might be
 * useful for computations with rotational symmetry, where one
 * dimension is the radius from the axis of rotation.
 *
 * @author Wolfgang Bangerth, 2000
 */
template <int dim>
class HalfHyperShellBoundary : public HyperShellBoundary<dim> 
{
  public:
				     /**
				      * Constructor. The center of the
				      * spheres defaults to the
				      * origin.
				      */
    HalfHyperShellBoundary (const Point<dim> &center = Point<dim>());
    
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
				      * @p{get_intermediate_points_between_points}.
				      */
    virtual void
    get_intermediate_points_on_line (const typename Triangulation<dim>::line_iterator &line,
				     typename std::vector<Point<dim> > &points) const;

				     /**
				      * Refer to the general
				      * documentation of this class
				      * and the documentation of the
				      * base class.
				      *
				      * Only implemented for @p{dim=3}
				      * and for @p{points.size()==1}.
				      */
    virtual void
    get_intermediate_points_on_quad (const typename Triangulation<dim>::quad_iterator &quad,
				     typename std::vector<Point<dim> > &points) const;

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
				      * Store the center of the spheres.
				      */
    const Point<dim> center;
};



/* -------------- declaration of explicit specializations ------------- */

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
get_normals_at_vertices (const Triangulation<1>::face_iterator &,
			 Boundary<1>::FaceVertexNormals &) const;
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
get_normals_at_vertices (const Triangulation<1>::face_iterator &,
			 Boundary<1>::FaceVertexNormals &) const;
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
get_normals_at_vertices (const Triangulation<1>::face_iterator &,
			 Boundary<1>::FaceVertexNormals &) const;


#endif
