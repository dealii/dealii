/*----------------------------   geometry_info.h     ---------------------------*/
/*      $Id$                 */
/*      Copyright W. Bangerth, University of Heidelberg, 1998 */
#ifndef __geometry_info_H
#define __geometry_info_H
/*----------------------------   geometry_info.h     ---------------------------*/


#include <base/exceptions.h>

template <int _dim> struct GeometryInfo;

/**
 * Pseudo-class for recursive functions in #GeometryInfo<dim>#.
 */
struct GeometryInfo<0>
{
    static const unsigned int vertices_per_cell = 1;
    static const unsigned int lines_per_cell = 0;
    static const unsigned int quads_per_cell = 0;
    static const unsigned int hexes_per_cell = 0;
    static const unsigned int children_per_cell = 0;
};


/**
 * Topological description of cells.
 *
 * This class contains as static members information on vertices and
 * faces of a #dim#-dimensinal grid cell. The interface is the same
 * for all dimensions. If a value is of no use in a low dimensional
 * cell, it is (correctly) set to zero, e.g. #sub_faces_per_cell# in
 * 1d.
 *
 * This information should always replace hard-coded numbers of
 * vertices, neighbors and so on, since it can be used dimension
 * independent.
 *
 * @author Wolfgang Bangerth, Guido Kanschat, 1998, 1999
 */
template <int _dim>
struct GeometryInfo
{
				     /**
				      * Present dimension. Does not look useful, but might be.
				      */
    static const unsigned int dim               = _dim;

				     /**
				      * Number of children a cell has.
				      */
    static const unsigned int children_per_cell = (1<<_dim);

				     /**
				      * Number of faces a cell has.
				      */
    static const unsigned int faces_per_cell    = 2*_dim;

				     /**
				      * Number of children each face has
				      * when the adjacent cell is refined.
				      */
    static const unsigned int subfaces_per_face = GeometryInfo<_dim-1>::children_per_cell;

				     /**
				      * Number of vertices a cell has.
				      */
    static const unsigned int vertices_per_cell = (1<<_dim);

				     /**
				      * Number of vertices each face has.
				      * Since this is not useful in one
				      * dimension, we provide a useless
				      * number (in the hope that a compiler
				      * may warn when it sees constructs like
				      * #for (i=0; i<vertices_per_face; ++i)#,
				      * at least if #i# is an #unsigned int#.
				      */
    static const unsigned int vertices_per_face = GeometryInfo<_dim-1>::vertices_per_cell;

				     /**
				      * Number of lines each face has.
				      */
    static const unsigned int lines_per_face = GeometryInfo<_dim-1>::lines_per_cell;
    
				     /**
				      * Number of quads on each face.
				      */
    static const unsigned int quads_per_face = GeometryInfo<_dim-1>::quads_per_cell;

				     /**
				      * Number of lines of a cell.
				      *
				      * Computation of this value
				      * follows the idea, that
				      * building a hypercube of
				      * dimension #dim# from one of
				      * dimension #dim#-1 can be done
				      * in the following two steps:
				      *
				      * 1. Duplicated it in the new coordinate direction.
				      *
				      * 2. Connect all corresponding
				      * vertices of the original
				      * hypercube and the copy by
				      * lines.
				      */
    static const unsigned int lines_per_cell = (2*GeometryInfo<dim-1>::lines_per_cell
						+ GeometryInfo<dim-1>::vertices_per_cell);
    
				     /**
				      * Number of quadrilaterals of a cell.
				      *
				      * Computation is analogous to #lines_per_cell#.
				      */
    static const unsigned int quads_per_cell = (2*GeometryInfo<dim-1>::quads_per_cell
						+ GeometryInfo<dim-1>::lines_per_cell);

				     /**
				      * Number of hexahedra of a cell.
				      *
				      * Computation is analogous to
				      * #lines_per_cell#. Very
				      * important for more than three
				      * dimensions!
				      */
    static const unsigned int hexes_per_cell = (2*GeometryInfo<dim-1>::hexes_per_cell
						+ GeometryInfo<dim-1>::quads_per_cell);

				     /**
				      * List of numbers which is
				      * denote which face is opposite
				      * to a given face. In 1d, this
				      * list is #{1,0}#, in 2d #{2, 3, 0, 1}#,
				      * in 3d #{1, 0, 4, 5, 2, 3}#.
				      */
    static const unsigned int opposite_face[faces_per_cell];
    
				     /**
				      * This field store which child cells
				      * are adjacent to a certain face of
				      * the mother cell.
				      *
				      * For example, in 2D the layout of
				      * a cell is as follows:
				      * \begin{verbatim}
				      * .      2
				      * .   3-->--2
				      * .   |     |
				      * . 3 ^     ^ 1
				      * .   |     |
				      * .   0-->--1
				      * .      0
				      * \end{verbatim}
				      * Vertices and faces are indicated
				      * with their numbers, faces also with
				      * their directions.
				      *
				      * Now, when refined, the layout is
				      * like this:
				      * \begin{verbatim}
				      * *--*--*
				      * | 3|2 |
				      * *--*--*
				      * | 0|1 |
				      * *--*--*
				      * \end{verbatim}
				      *
				      * Thus, the child cells on face zero
				      * are (ordered in the direction of the
				      * face) 0 and 1, on face 2 they are
				      * 3 and 2, etc.
				      *
				      * For three spatial dimensions,
				      * the exact order of the children is
				      * laid down in the documentation of
				      * the #Triangulation# class.
				      */
    static unsigned int child_cell_on_face (unsigned int face,
					    unsigned int subface);

};



/*---------------------------- Inline functions --------------------------------*/

template<>
inline unsigned int
GeometryInfo<2>::child_cell_on_face (const unsigned int face,
                                     const unsigned int subface)
{
  Assert (face<faces_per_cell, ExcIndexRange(face, 0, faces_per_cell));
  Assert (subface<subfaces_per_face, ExcIndexRange(subface, 0, subfaces_per_face));
  
  const unsigned subcells[faces_per_cell][subfaces_per_face] = {{0,1},
								{1,2},
								{3,2},
								{0,3}};
  return subcells[face][subface];
};


template<>
inline
unsigned int GeometryInfo<3>::child_cell_on_face (const unsigned int face,
						  const unsigned int subface)
{
  Assert (face<faces_per_cell, ExcIndexRange(face, 0, faces_per_cell));
  Assert (subface<subfaces_per_face, ExcIndexRange(subface, 0, subfaces_per_face));
  
  const unsigned subcells[faces_per_cell][subfaces_per_face] = {{0, 1, 2, 3},
								{4, 5, 6, 7},
								{0, 1, 5, 4},
								{1, 5, 6, 2},
								{3, 2, 6, 7},
								{0, 4, 7, 3}};
  return subcells[face][subface];
};



/*----------------------------   geometry_info.h     ---------------------------*/
/* end of #ifndef __geometry_info_H */
#endif
/*----------------------------   geometry_info.h     ---------------------------*/
