/*----------------------------   boundary-function.h     ---------------------------*/
/*      $Id$                 */
/*      Copyright W. Bangerth, University of Heidelberg, 1998 */
#ifndef __tria_boundary_H
#define __tria_boundary_H
/*----------------------------   boundary-function.h     ---------------------------*/

#include <base/point.h>
#include <base/subscriptor.h>
#include <grid/geometry_info.h>
#include <basic/forward-declarations.h>



/**
 *   This class is used to represent a boundary to a triangulation.
 *   When a triangulation creates a new vertex on the boundary of the
 *   domain, it determines the new vertex' coordinates through the
 *   following code (here in two dimensions):
 *   \begin{verbatim}
 *     ...
 *     Point<2> new_vertex = boundary.get_new_point_on_line (line);
 *     ...
 *   \end{verbatim}
 *   #line# denotes the line at the boundary that shall be refined
 *   and for which we seek the common point of the two child lines.
 *
 *   In 3D, a new vertex may be placed on the middle of a line or on
 *   the middle of a side. Respectively, the library calls
 *   \begin{verbatim}
 *     ...
 *     Point<3> new_line_vertices[4]
 *       = { boundary.get_new_point_on_line (face->line(0)),
 *           boundary.get_new_point_on_line (face->line(1)),
 *           boundary.get_new_point_on_line (face->line(2)),
 *           boundary.get_new_point_on_line (face->line(3))  };
 *     ...
 *   \end{verbatim}
 *   to get the four midpoints of the lines bounding the quad at the
 *   boundary, and after that
 *   \begin{verbatim}
 *     ...
 *     Point<3> new_quad_vertex = boundary.get_new_point_on_quad (face);
 *     ...
 *   \end{verbatim}
 *   to get the midpoint of the face. It is guaranteed that this order
 *   (first lines, then faces) holds, so you can use information from
 *   the children of the four lines of a face, since these already exist
 *   at the time the midpoint of the face is to be computed.
 *   
 *   Since iterators are passed to the functions, you may use information
 *   about boundary indicators and the like, as well as all other information
 *   provided by these objects.
 *
 *   There are specialisations, #StraightBoundary<dim>#, which places
 *   the new point right into the middle of the given points, and
 *   #HyperBallBoundary<dim># creating a hyperball with given radius
 *   around a given center point.
 *
 *   @author Wolfgang Bangerth, 1999
 */
template <int dim>
class Boundary : public Subscriptor {
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
				      * four lines bounding the given #quad#
				      * are refined, so you may want to use
				      * the information provided by
				      * #quad->line(i)->child(j)#, #i=0...3#,
				      * #j=0,1#.
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
				      * Exception.
				      */
    DeclException0 (ExcPureVirtualFunctionCalled);
};





/**
 *   Specialisation of \Ref{Boundary}<dim>, which places the new point right
 *   into the middle of the given points. The middle is defined as the
 *   arithmetic mean of the points.
 *
 *   This class does not really describe a boundary in the usual sense. By
 *   placing new points in the middle of old ones, it rather assumes that the
 *   boundary of the domain is given by the polygon/polyhedron defined by the
 *   boundary of the initial coarse triangulation.
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
};





/*----------------------------   boundary-function.h     ---------------------------*/
/* end of #ifndef __tria_boundary_H */
#endif
/*----------------------------   boundary-function.h     ---------------------------*/
