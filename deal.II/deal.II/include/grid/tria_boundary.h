//----------------------------  tria_boundary.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  tria_boundary.h  ---------------------------
#ifndef __deal2__tria_boundary_h
#define __deal2__tria_boundary_h


/*----------------------------   boundary-function.h     ---------------------------*/

#include <base/subscriptor.h>
#include <base/point.h>
#include <grid/tria.h>


template <int dim> class Triangulation;



/**
 *   This class is used to represent a boundary to a triangulation.
 *   When a triangulation creates a new vertex on the boundary of the
 *   domain, it determines the new vertex' coordinates through the
 *   following code (here in two dimensions):
 *   @begin{verbatim}
 *     ...
 *     Point<2> new_vertex = boundary.get_new_point_on_line (line);
 *     ...
 *   @end{verbatim}
 *   @p{line} denotes the line at the boundary that shall be refined
 *   and for which we seek the common point of the two child lines.
 *
 *   In 3D, a new vertex may be placed on the middle of a line or on
 *   the middle of a side. Respectively, the library calls
 *   @begin{verbatim}
 *     ...
 *     Point<3> new_line_vertices[4]
 *       = { boundary.get_new_point_on_line (face->line(0)),
 *           boundary.get_new_point_on_line (face->line(1)),
 *           boundary.get_new_point_on_line (face->line(2)),
 *           boundary.get_new_point_on_line (face->line(3))  };
 *     ...
 *   @end{verbatim}
 *   to get the four midpoints of the lines bounding the quad at the
 *   boundary, and after that
 *   @begin{verbatim}
 *     ...
 *     Point<3> new_quad_vertex = boundary.get_new_point_on_quad (face);
 *     ...
 *   @end{verbatim}
 *   to get the midpoint of the face. It is guaranteed that this order
 *   (first lines, then faces) holds, so you can use information from
 *   the children of the four lines of a face, since these already exist
 *   at the time the midpoint of the face is to be computed.
 *   
 *   Since iterators are passed to the functions, you may use information
 *   about boundary indicators and the like, as well as all other information
 *   provided by these objects.
 *
 *   There are specialisations, @ref{StraightBoundary}, which places
 *   the new point right into the middle of the given points, and
 *   @ref{HyperBallBoundary} creating a hyperball with given radius
 *   around a given center point.
 *
 *   @author Wolfgang Bangerth, 1999, Ralf Hartmann, 2001
 */
template <int dim>
class Boundary : public Subscriptor
{
  public:
				     /**
				      * Destructor. Does nothing here, but
				      * needs to be declared to make it
				      * virtual.
				      */
    virtual ~Boundary ();

				     /**
				      * Return the point which shall become
				      * the new middle vertex of the two
				      * children of a regular line. In 2D,
				      * this line is a line at the boundary,
				      * while in 3d, it is bounding a face
				      * at the boundary (the lines therefore
				      * is also on the boundary).
				      */
    virtual Point<dim>
    get_new_point_on_line (const typename Triangulation<dim>::line_iterator &line) const = 0;

				     /**
				      * Return the point which shall become
				      * the common point of the four children
				      * of a quad at the boundary in three
				      * or more spatial dimensions. This
				      * function therefore is only useful in
				      * at least three dimensions and should
				      * not be called for lower dimensions.
				      *
				      * This function is called after the
				      * four lines bounding the given @p{quad}
				      * are refined, so you may want to use
				      * the information provided by
				      * @p{quad->line(i)->child(j)}, @p{i=0...3},
				      * @p{j=0,1}.
				      *
				      * Because in 2D, this function is not
				      * needed, it is not made pure virtual,
				      * to avoid the need to overload it.
				      * The default implementation throws
				      * an error in any case, however.
				      */
    virtual Point<dim>
    get_new_point_on_quad (const typename Triangulation<dim>::quad_iterator &quad) const;

				     /**
				      * Return equally spaced
				      * intermediate points on a line.
                                      *
				      * The number of points requested
				      * is given by the size of the
				      * vector @p{points}. It is the
				      * task of the derived classes to
				      * arrange the points in
				      * approximately equal distances.
				      *
				      * This function is called by the
				      * @p{MappingQ} class. This
				      * happens on each face line of a
				      * cells that has got at least
				      * one boundary line.
				      *
				      * As this function is not needed
				      * for @p{MappingQ1}, it is not
				      * made pure virtual, to avoid
				      * the need to overload it.  The
				      * default implementation throws
				      * an error in any case, however.
				      */
    virtual void
    get_intermediate_points_on_line (const typename Triangulation<dim>::line_iterator &line,
				     typename std::vector<Point<dim> > &points) const;
    
				     /**
				      * Return equally spaced
				      * intermediate points on a
				      * boundary quad.
				      *
				      * The number of points requested
				      * is given by the size of the
				      * vector @p{points}. It is
				      * required that this number is a
				      * square of another integer,
				      * i.e. @p{n=points.size()=m*m}. It
				      * is the task of the derived
				      * classes to arrange the points
				      * such they split the quad into
				      * @p{(m+1)(m+1)} approximately
				      * equal-sized subquads.
				      *
				      * This function is called by the
				      * @p{MappingQ<3>} class. This
				      * happens each face quad of
				      * cells in 3d that has got at
				      * least one boundary face quad.
				      *
				      * As this function is not needed
				      * for @p{MappingQ1}, it is not
				      * made pure virtual, to avoid
				      * the need to overload it.  The
				      * default implementation throws
				      * an error in any case, however.
				      */
    virtual void
    get_intermediate_points_on_quad (const typename Triangulation<dim>::quad_iterator &quad,
				     typename std::vector<Point<dim> > &points) const;

				     /**
				      * Exception.
				      */
    DeclException0 (ExcPureVirtualFunctionCalled);

				     /**
				      * Exception.
				      */
    DeclException1 (ExcFunctionNotUseful,
		    int,
		    << "The function called is not useful for dim=" << arg1 << ".");
};



/**
 *   Specialisation of @ref{Boundary}<dim>, which places the new point right
 *   into the middle of the given points. The middle is defined as the
 *   arithmetic mean of the points.
 *
 *   This class does not really describe a boundary in the usual sense. By
 *   placing new points in the middle of old ones, it rather assumes that the
 *   boundary of the domain is given by the polygon/polyhedron defined by the
 *   boundary of the initial coarse triangulation.
 *
 *   @author Wolfgang Bangerth, 1998, Ralf Hartmann, 2001
 */
template <int dim>
class StraightBoundary : public Boundary<dim> {
  public:
				     /**
				      * Let the new point be the arithmetic
				      * mean of the two vertices of the line.
				      *
				      * Refer to the general documentation of
				      * this class and the documentation of the
				      * base class for more information.
				      */
    virtual Point<dim>
    get_new_point_on_line (const typename Triangulation<dim>::line_iterator &line) const;

				     /**
				      * Let the new point be the arithmetic mean
				      * of the four vertices of this quad and
				      * the four midpoints of the lines, which
				      * are already created at the time of calling
				      * this function.
				      *
				      * Refer to the general documentation of
				      * this class and the documentation of the
				      * base class for more information.
				      */
    virtual Point<dim>
    get_new_point_on_quad (const typename Triangulation<dim>::quad_iterator &quad) const;

				     /**
				      * Gives @p{n=points.size()}
				      * points that splits the
				      * p{StraightBoundary} line into
				      * p{n+1} partitions of equal
				      * lengths.
				      *
				      * Refer to the general
				      * documentation of this class
				      * and the documentation of the
				      * base class.
				      */
    virtual void
    get_intermediate_points_on_line (const typename Triangulation<dim>::line_iterator &line,
				     typename std::vector<Point<dim> > &points) const;

				     /**
				      * Gives @p{n=points.size()=m*m}
				      * points that splits the
				      * p{StraightBoundary} quad into
				      * @p{(m+1)(m+1)} subquads of equal
				      * size.
				      *
				      * Refer to the general
				      * documentation of this class
				      * and the documentation of the
				      * base class.
				      */
    virtual void
    get_intermediate_points_on_quad (const typename Triangulation<dim>::quad_iterator &quad,
				     typename std::vector<Point<dim> > &points) const;

};



/*----------------------------   boundary-function.h     ---------------------------*/

#endif
/*----------------------------   boundary-function.h     ---------------------------*/
