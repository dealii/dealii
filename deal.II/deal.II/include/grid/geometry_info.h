/*----------------------------   geometry_info.h     ---------------------------*/
/*      $Id$                 */
#ifndef __geometry_info_H
#define __geometry_info_H
/*----------------------------   geometry_info.h     ---------------------------*/


#include <base/exceptions.h>


template <int dim> struct GeometryInfo;


/**
 *  Publish some information about geometrical interconnections to the
 *  outside world, for one spacial dimension in this case. These are,
 *  for example the numbers of children per cell, faces per cell, etc,
 *  but also neighborship information, It is especially useful if you
 *  want to loop over all faces in any space dimension, but don't want
 *  to think about their number in a dimension independent expression.
 *  This not only reduces thinking effort but also error possibilities.
 */
struct GeometryInfo<1> {
  public:
				     /**
				      * Present dimension.
				      */
    static const unsigned int dim               = 1;

				     /**
				      * Number of children a cell has.
				      */
    static const unsigned int children_per_cell = (1<<dim);

				     /**
				      * Number of faces a cell has.
				      */
    static const unsigned int faces_per_cell    = 2*dim;

				     /**
				      * Number of vertices a cell has.
				      */
    static const unsigned int vertices_per_cell = (1<<dim);

				     /**
				      * Number of vertices each face has.
				      * Since this is not useful in one
				      * dimension, we provide a useless
				      * number (in the hope that a compiler
				      * may warn when it sees constructs like
				      * #for (i=0; i<vertices_per_face; ++i)#,
				      * at least if #i# is an #unsigned int#.
				      */
    static const unsigned int vertices_per_face = 0;
    
};



/**
 *  Publish some information about geometrical interconnections to the
 *  outside world, for two spacial dimensions in this case. These are,
 *  for example the numbers of children per cell, faces per cell, etc,
 *  but also neighborship information, It is especially useful if you
 *  want to loop over all faces in any space dimension, but don't want
 *  to think about their number in a dimension independent expression.
 *  This not only reduces thinking effort but also error possibilities.
 */
struct GeometryInfo<2> {
  public:
				     /**
				      * Present dimension
				      */
    static const unsigned int dim               = 2;

				     /**
				      * Number of children a cell has.
				      */
    static const unsigned int children_per_cell = (1<<dim);

				     /**
				      * Number of faces a cell has.
				      */
    static const unsigned int faces_per_cell    = 2*dim;

				     /**
				      * Number of children each face has
				      * when the adjacent cell is refined.
				      */
    static const unsigned int subfaces_per_face = (1<<(dim-1));

				     /**
				      * Number of vertices a cell has.
				      */
    static const unsigned int vertices_per_cell = (1<<dim);

				     /**
				      * Number of vertices each face has.
				      */
    static const unsigned int vertices_per_face = (1<<(dim-1));

				     /**
				      * This field store which child cells
				      * are adjacent to a certain face of
				      * the mother cell.
				      *
				      * For example, in 2D the layout of
				      * a cell is as follows:
				      *      2
				      *   3-->--2
				      *   |     |
				      * 3 ^     ^ 1
				      *   |     |
				      *   0-->--1
				      *      0
				      * Vertices and faces are indicated
				      * with their numbers, faces also with
				      * their directions.
				      *
				      * Now, when refined, the layout is
				      * like this:
				      * *--*--*
				      * | 3|2 |
				      * *--*--*
				      * | 0|1 |
				      * *--*--*
				      *
				      * Thus, the child cells on face zero
				      * are (ordered in the direction of the
				      * face) 0 and 1, on face 2 they are
				      * 3 and 2, etc.
				      */
    static unsigned int child_cell_on_face (unsigned int face,
					    unsigned int subface);

				     /**
				      * Exception
				      */
    DeclException2 (ExcInvalidIndex,
		    int,
		    int,
		    << "Invalid index " << arg1
		    << ", index must be between 0 and " << arg2 << ".");
};





/*---------------------------- Inline functions --------------------------------*/

inline
unsigned int GeometryInfo<2>::child_cell_on_face (const unsigned int face,
						  const unsigned int subface) {
  Assert (face<faces_per_cell, ExcInvalidIndex(face,faces_per_cell));
  Assert (subface<subfaces_per_face, ExcInvalidIndex(subface,subfaces_per_face));
  
  static const unsigned subcells[faces_per_cell][subfaces_per_face] = {{0,1},
								       {1,2},
								       {3,2},
								       {0,3}};
  return subcells[face][subface];
};



/*----------------------------   geometry_info.h     ---------------------------*/
/* end of #ifndef __geometry_info_H */
#endif
/*----------------------------   geometry_info.h     ---------------------------*/
