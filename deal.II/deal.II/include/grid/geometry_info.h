//----------------------------  geometry_info.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998-2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  geometry_info.h  ---------------------------
#ifndef __deal2__geometry_info_h
#define __deal2__geometry_info_h


#include <base/config.h>
#include <base/exceptions.h>
#include <base/point.h>


template <int dim> class GeometryInfo;


/**
 * Topological description of zero dimensional cells,
 * i.e. points. This class might not look too useful but often is if
 * in a certain dimension we would like to enquire information about
 * objects with dimension one lower than the present, e.g. about
 * faces.
 *
 * This class contains as static members information on vertices and
 * faces of a @p dim-dimensional grid cell. The interface is the same
 * for all dimensions. If a value is of no use in a low dimensional
 * cell, it is (correctly) set to zero, e.g. @p subfaces_per_cell in
 * 1d.
 *
 * This information should always replace hard-coded numbers of
 * vertices, neighbors and so on, since it can be used dimension
 * independently.
 *
 * @author Wolfgang Bangerth, 1998
 */
template <>
struct GeometryInfo<0> 
{

				     /**
				      * Number of children a cell has.
				      */
    static const unsigned int children_per_cell = 1;

				     /**
				      * Number of faces a cell has.
				      */
    static const unsigned int faces_per_cell    = 0;

				     /**
				      * Number of children each face has
				      * when the adjacent cell is refined.
				      */
    static const unsigned int subfaces_per_face = 0;

				     /**
				      * Number of vertices a cell has.
				      */
    static const unsigned int vertices_per_cell = 1;

				     /**
				      * Number of vertices each face has.
				      * Since this is not useful in one
				      * dimension, we provide a useless
				      * number (in the hope that a compiler
				      * may warn when it sees constructs like
				      * <tt>for (i=0; i<vertices_per_face; ++i)</tt>,
				      * at least if @p i is an <tt>unsigned int</tt>.
				      */
    static const unsigned int vertices_per_face = 0;

				     /**
				      * Number of lines each face has.
				      */
    static const unsigned int lines_per_face    = 0;
    
				     /**
				      * Number of quads on each face.
				      */
    static const unsigned int quads_per_face    = 0;

				     /**
				      * Number of lines of a cell.
				      */
    static const unsigned int lines_per_cell    = 0;

				     /**
				      * Number of quadrilaterals of a
				      * cell.
				      */
    static const unsigned int quads_per_cell    = 0;

				     /**
				      * Number of hexahedra of a
				      * cell.
				      */
    static const unsigned int hexes_per_cell    = 0;
};



/**
 * Topological description of four dimensional cells. This class is
 * required in some exotic cases where we compute information in a
 * one-larger dimension than the present, and do so also in 3d (for
 * example, stacking the solutions of a d-dimensional time dependent
 * equation timestep after timestep in a (d+1)-dimensional space).
 *
 * This class contains as static members information on vertices and
 * faces of a @p dim-dimensional grid cell. The interface is the same
 * for all dimensions. If a value is of no use in a low dimensional
 * cell, it is (correctly) set to zero, e.g. @p subfaces_per_cell in
 * 1d.
 *
 * This information should always replace hard-coded numbers of
 * vertices, neighbors and so on, since it can be used dimension
 * independently.
 *
 * @author Wolfgang Bangerth, 1998
 */
template <>
struct GeometryInfo<4>
{

				     /**
				      * Number of children a cell has.
				      */
    static const unsigned int children_per_cell = 16;

				     /**
				      * Number of faces a cell has.
				      */
    static const unsigned int faces_per_cell    = 8;

				     /**
				      * Number of children each face has
				      * when the adjacent cell is refined.
				      */
    static const unsigned int subfaces_per_face = 8;

				     /**
				      * Number of vertices a cell has.
				      */
    static const unsigned int vertices_per_cell = 16;

				     /**
				      * Number of vertices each face
				      * has.
				      */
    static const unsigned int vertices_per_face = 8;

				     /**
				      * Number of lines each face has.
				      */
    static const unsigned int lines_per_face    = 12;
    
				     /**
				      * Number of quads on each face.
				      */
    static const unsigned int quads_per_face    = 6;

				     /**
				      * Number of lines of a cell.
				      */
    static const unsigned int lines_per_cell    = 32;

				     /**
				      * Number of quadrilaterals of a
				      * cell.
				      */
    static const unsigned int quads_per_cell    = 24;

				     /**
				      * Number of hexahedra of a
				      * cell.
				      */
    static const unsigned int hexes_per_cell    = 8;
};



/**
 * This template specifies the interface to all topological structure
 * of the mesh cells. This template class can be instantiated for
 * dim=1,2 and 3.
 *
 * @author Wolfgang Bangerth, 1998, Ralf Hartmann, 2005
 */
template <int dim>
struct GeometryInfo
{
    
				     /**
				      * Number of children of a refined cell.
				      */
    static const unsigned int children_per_cell = 1 << dim;

				     /**
				      * Number of faces of a cell.
				      */
    static const unsigned int faces_per_cell = 2 * dim;

				     /**
				      * Number of children each face has
				      * when the adjacent cell is refined.
				      */
    static const unsigned int subfaces_per_face = GeometryInfo<dim-1>::children_per_cell;

				     /**
				      * Number of vertices of a cell.
				      */
    static const unsigned int vertices_per_cell = 1 << dim;

				     /**
				      * Number of vertices on each
				      * face.
				      */
    static const unsigned int vertices_per_face = GeometryInfo<dim-1>::vertices_per_cell;

				     /**
				      * Number of lines on each face.
				      */
    static const unsigned int lines_per_face = GeometryInfo<dim-1>::lines_per_cell;
    
				     /**
				      * Number of quads on each face.
				      */
    static const unsigned int quads_per_face = (dim == 3) ? 1 : 0;

				     /**
				      * Number of lines of a cell.
				      */
    static const unsigned int lines_per_cell = (dim == 3) ? 12 : ((dim == 2) ? 4 : 1);

				     /**
				      * Number of quadrilaterals of a
				      * cell.
				      */
    static const unsigned int quads_per_cell = (dim == 3) ? 6 : ((dim == 2) ? 1 : 0);

				     /**
				      * Number of hexahedra of a
				      * cell.
				      */
    static const unsigned int hexes_per_cell = (dim == 3) ? 1 : 0;

				     /**
				      * List of numbers which
				      * denotes which face is opposite
				      * to a given face. In 1d, this
				      * list is <tt>{1,0</tt>}, in 2d <tt>{2, 3, 0, 1</tt>},
				      * in 3d <tt>{1, 0, 4, 5, 2, 3</tt>}.
				      */
    static const unsigned int opposite_face[faces_per_cell];
    

				     /**
				      * Rearrange vertices for OpenDX
				      * output.  For a cell being
				      * written in OpenDX format, each
				      * entry in this field contains
				      * the number of a vertex in
				      * <tt>deal.II</tt> that corresponds
				      * to the DX numbering at this
				      * location.
				      *
				      * Typical example: write a cell
				      * and arrange the vertices, such
				      * that OpenDX understands them.
				      *
				      * \begin{verbatim}
				      * for (i=0; i< n_vertices; ++i)
				      *   out << cell->vertex(dx_to_deal[i]);
				      * \end{verbatim}
				      */
    static const unsigned int dx_to_deal[vertices_per_cell];
    
				     /**
				      * For each face of the reference
				      * cell, this field stores the
				      * coordinate direction in which
				      * its normal vector points.
				      *
				      * Remark that this is only the
				      * coordinate number. The acual
				      * direction of the normal vector
				      * is obtained by multiplying the
				      * unit vector in this direction
				      * with #unit_normal_orientation.
				      */
    static const unsigned int unit_normal_direction[faces_per_cell];

				     /**
				      * Orientation of the unit normal
				      * vector of a face of the
				      * reference cell.
				      *
				      * Each value is either
				      * <tt>1</tt> or <tt>-1</tt>,
				      * corresponding to a normal
				      * vector pointing in the
				      * positive or negative
				      * coordinate direction,
				      * respectively.
				      */
    static const int unit_normal_orientation[faces_per_cell];
    
				     /**
				      * This field stores which child
				      * cells are adjacent to a
				      * certain face of the mother
				      * cell.
				      *
				      * For example, in 2D the layout of
				      * a cell is as follows:
				      * @verbatim
				      * .      2
				      * .   3-->--2
				      * .   |     |
				      * . 3 ^     ^ 1
				      * .   |     |
				      * .   0-->--1
				      * .      0
				      * @endverbatim
				      * Vertices and faces are indicated
				      * with their numbers, faces also with
				      * their directions.
				      *
				      * Now, when refined, the layout is
				      * like this:
				      * @verbatim
				      * *--*--*
				      * | 3|2 |
				      * *--*--*
				      * | 0|1 |
				      * *--*--*
				      * @endverbatim
				      *
				      * Thus, the child cells on face zero
				      * are (ordered in the direction of the
				      * face) 0 and 1, on face 2 they are
				      * 3 and 2, etc.
				      *
				      * For three spatial dimensions,
				      * the exact order of the
				      * children is laid down in the
				      * documentation of the
				      * Triangulation class.
				      * However, it must be noted that
				      * this class and function only
				      * deals with faces in standard
				      * orientation. In 3d, faces can
				      * exist in two orientations,
				      * though, and if a face is in
				      * the wrong orientation, then
				      * this function may not give you
				      * what you want. You can inquire
				      * about the face orientation
				      * using the
				      * <tt>cell->face_orientation</tt>
				      * function, and the function to
				      * ask for a neighbor's cell
				      * behind a given face and
				      * subface is
				      * <tt>cell->neighbor_child_on_subface</tt>.
				      * The latter function, in
				      * contrast to the present one,
				      * also takes into account the
				      * actual orientation of the
				      * faces of a cell and will
				      * return the correct result in
				      * all cases.
				      */
    static unsigned int child_cell_on_face (const unsigned int face,
					    const unsigned int subface);
    
				     /**
				      * Map line vertex number to cell
				      * vertex number, i.e. give the
				      * cell vertex number of the
				      * <tt>vertex</tt>th vertex of
				      * line <tt>line</tt>, e.g.
				      * <tt>GeometryInfo<2>::line_to_cell_vertices(2,0)=3</tt>.
				      *
				      * The order of the lines, as
				      * well as their direction (which
				      * in turn determines which is
				      * the first and which the second
				      * vertex on a line) is the
				      * canonical one in deal.II, as
				      * described in the documentation
				      * of the Triangulation
				      * class.
				      *
				      * For <tt>dim=2</tt> this call
				      * is simply passed down to the
				      * face_to_cell_vertices()
				      * function.
				      */
    static unsigned int line_to_cell_vertices (const unsigned int line,
					       const unsigned int vertex);

				     /**
				      * Map face vertex number to cell
				      * vertex number, i.e. give the
				      * cell vertex number of the
				      * <tt>vertex</tt>th vertex of
				      * face <tt>face</tt>, e.g.
				      * <tt>GeometryInfo<2>::face_to_cell_vertices(2,0)=3</tt>.
				      *
				      * As the children of a cell are
				      * ordered according to the
				      * vertices of the cell, this
				      * call is passed down to the
				      * child_cell_on_face() function.
				      * Hence this function is simply
				      * a wrapper of
				      * child_cell_on_face() giving it
				      * a suggestive name.
				      */
    static unsigned int face_to_cell_vertices (const unsigned int face,
					       const unsigned int vertex);

				     /**
				      * Map face line number to cell
				      * line number, i.e. give the
				      * cell line number of the
				      * <tt>line</tt>th line of face
				      * <tt>face</tt>, e.g.
				      * <tt>GeometryInfo<3>::face_to_cell_lines(3,1)=5</tt>.
				      *
				      * This function is useful and
				      * implemented for
				      * <tt>dim=3</tt>, only.
				      */
    static unsigned int face_to_cell_lines (const unsigned int face,
					    const unsigned int line);
    
				     /**
				      * Return the position of the
				      * @p ith vertex on the unit
				      * cell. The order of vertices is
				      * the canonical one in deal.II,
				      * as described in the
				      * documentation of the
				      * Triangulation class.
				      */
    static Point<dim> unit_cell_vertex (const unsigned int vertex);

				     /**
				      * Given a point @p p in unit
				      * coordinates, return the number
				      * of the child cell in which it
				      * would lie in. If the point
				      * lies on the interface of two
				      * children, return any one of
				      * their indices. The result is
				      * always less than
				      * <tt>GeometryInfo<dimension>::children_per_cell</tt>.
				      *
				      * The order of child cells is
				      * described the documentation of
				      * the Triangulation class.
				      */
    static unsigned int child_cell_from_point (const Point<dim> &p);

				     /**
				      * Given coordinates @p p on the
				      * unit cell, return the values
				      * of the coordinates of this
				      * point in the coordinate system
				      * of the given child. Neither
				      * original nor returned
				      * coordinates need actually be
				      * inside the cell, we simply
				      * perform a scale-and-shift
				      * operation with a shift that
				      * depends on the number of the
				      * child.
				      */
    static Point<dim> cell_to_child_coordinates (const Point<dim>    &p,
						 const unsigned int child_index);

				     /**
				      * The reverse function to the
				      * one above: take a point in the
				      * coordinate system of the
				      * child, and transform it to the
				      * coordinate system of the
				      * mother cell.
				      */
    static Point<dim> child_to_cell_coordinates (const Point<dim>    &p,
						 const unsigned int child_index);

				     /**
				      * Return true if the given point
				      * is inside the unit cell of the
				      * present space dimension.
				      */
    static bool is_inside_unit_cell (const Point<dim> &p);
    
				     /**
				      * Exception
				      */
    DeclException1 (ExcInvalidCoordinate,
		    double,
		    << "The coordinates must satisfy 0 <= x_i <= 1, "
		    << "but here we have x_i=" << arg1);
};


/* -------------- declaration of explicit specializations ------------- */

/// @if NoDoc



template <>
inline
Point<1>
GeometryInfo<1>::unit_cell_vertex (const unsigned int vertex)
{
  Assert (vertex < vertices_per_cell,
	  ExcIndexRange (vertex, 0, vertices_per_cell));

  const Point<1> vertices[vertices_per_cell] =
    { Point<1>(0.), Point<1>(1.) };
  return vertices[vertex];
}



template <>
inline
Point<2>
GeometryInfo<2>::unit_cell_vertex (const unsigned int vertex)
{
  Assert (vertex < vertices_per_cell,
	  ExcIndexRange (vertex, 0, vertices_per_cell));

  const Point<2> vertices[vertices_per_cell] =
    { Point<2>(0., 0.), Point<2>(1., 0.),
      Point<2>(1.,1.), Point<2>(0.,1.) };
  return vertices[vertex];
}



template <>
inline
Point<3>
GeometryInfo<3>::unit_cell_vertex (const unsigned int vertex)
{
  Assert (vertex < vertices_per_cell,
	  ExcIndexRange (vertex, 0, vertices_per_cell));

  const Point<3> vertices[vertices_per_cell] =
    { Point<3>(0., 0., 0.), Point<3>(1., 0., 0.),
      Point<3>(1., 0., 1.), Point<3>(0., 0., 1.),
      Point<3>(0., 1., 0.), Point<3>(1., 1., 0.),
      Point<3>(1., 1., 1.), Point<3>(0., 1., 1.) };
  return vertices[vertex];
}


template <int dim>
inline
Point<dim>
GeometryInfo<dim>::unit_cell_vertex (const unsigned int)
{
  Assert(false, ExcNotImplemented());

  return Point<dim> ();  
}



template <>
inline
unsigned int
GeometryInfo<1>::child_cell_from_point (const Point<1> &p)
{
  Assert ((p[0] >= 0) && (p[0] <= 1), ExcInvalidCoordinate(p[0]));
  
  return (p[0] <= 0.5 ? 0 : 1);
}



template <>
inline
unsigned int
GeometryInfo<2>::child_cell_from_point (const Point<2> &p)
{
  Assert ((p[0] >= 0) && (p[0] <= 1), ExcInvalidCoordinate(p[0]));
  Assert ((p[1] >= 0) && (p[1] <= 1), ExcInvalidCoordinate(p[1]));
  
  return (p[0] <= 0.5 ?
	  (p[1] <= 0.5 ? 0 : 3) :
	  (p[1] <= 0.5 ? 1 : 2));
}



template <>
inline
unsigned int
GeometryInfo<3>::child_cell_from_point (const Point<3> &p)
{
  Assert ((p[0] >= 0) && (p[0] <= 1), ExcInvalidCoordinate(p[0]));
  Assert ((p[1] >= 0) && (p[1] <= 1), ExcInvalidCoordinate(p[1]));
  Assert ((p[2] >= 0) && (p[2] <= 1), ExcInvalidCoordinate(p[2]));
  
  return (p[0] <= 0.5 ?
	  (p[1] <= 0.5 ?
	   (p[2] <= 0.5 ? 0 : 3) :
	   (p[2] <= 0.5 ? 4 : 7)) :
	  (p[1] <= 0.5 ?
	   (p[2] <= 0.5 ? 1 : 2) :
	   (p[2] <= 0.5 ? 5 : 6)));
}


template <int dim>
inline
unsigned int
GeometryInfo<dim>::child_cell_from_point (const Point<dim> &)
{
  Assert(false, ExcNotImplemented());

  return 0;
}



template <int dim>
inline
Point<dim>
GeometryInfo<dim>::cell_to_child_coordinates (const Point<dim>    &p,
					      const unsigned int child_index)
{
  Assert (child_index < GeometryInfo<dim>::children_per_cell,
	  ExcIndexRange (child_index, 0, GeometryInfo<dim>::children_per_cell));

  return 2*p - unit_cell_vertex(child_index);
}



template <int dim>
inline
Point<dim>
GeometryInfo<dim>::child_to_cell_coordinates (const Point<dim>    &p,
					      const unsigned int child_index)
{
  Assert (child_index < GeometryInfo<dim>::children_per_cell,
	  ExcIndexRange (child_index, 0, GeometryInfo<dim>::children_per_cell));

  return (p + unit_cell_vertex(child_index))/2;
}


template <>
inline
bool
GeometryInfo<1>::is_inside_unit_cell (const Point<1> &p)
{
  return (p[0] >= 0.) && (p[0] <= 1.);
}



template <>
inline
bool
GeometryInfo<2>::is_inside_unit_cell (const Point<2> &p)
{
  return (p[0] >= 0.) && (p[0] <= 1.) &&
	 (p[1] >= 0.) && (p[1] <= 1.);
}



template <>
inline
bool
GeometryInfo<3>::is_inside_unit_cell (const Point<3> &p)
{
  return (p[0] >= 0.) && (p[0] <= 1.) &&
	 (p[1] >= 0.) && (p[1] <= 1.) &&
	 (p[2] >= 0.) && (p[2] <= 1.);
}


/// @endif

#endif
