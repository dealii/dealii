/*----------------------------   tria_quad.h     ---------------------------*/
/*      $Id$                 */
#ifndef __tria_quad_H
#define __tria_quad_H
/*----------------------------   tria_quad.h     ---------------------------*/


#include <base/exceptions.h>

/**
 *   #Quad#s denote the fundamental entities of triangulations in two dimensions
 *   and the boundaries of hexaeders in three dimensions. They are
 *   characterized by the (global) indices of the corner points.
 *
 *   A quad itself has one index, as far as the topological part handled in
 *   the triangulation is concerned: the index in the level
 *   it belongs to. The level index is implicitely given by the position
 *   in the #quads.quads# list attached to the information of each level
 *   of the triangulation.
 *
 *   @author Wolfgang Bangerth, 1998
 */
class Quad {
  public:

				     /**
				      *  Construct a Quad with line
				      *  indices #i0# through #i3#. By default,
				      *  indices are set to -1, i.e. an
				      *  invalid value.
				      *
				      *  By convention, the four lines must
				      *  be numbered in counterclockwise sense!
				      */
    Quad (const int i0 = -1,
	  const int i1 = -1,
	  const int i2 = -1,
	  const int i3 = -1);
    
				     /**
				      *  Return the index of line #i#=0
				      *  through 3.
				      */
    int line (const int i) const;
    
				     /**
				      *  Set the index of the #i#th line to
				      *  #index#. #i#=0 through 3.
				      */
    void set_line (const int i, const int index);
    
				     /**
				      *  Exception
				      */ 
    DeclException1 (ExcRange,
		    int,
		    << "Indices for the point number must be 0, 1, 2 or 3, "
		    << "but you gave " << arg1); 
  protected:
      int lines[4];
};




/*----------------------------- Inline Function: Quad ------------------------*/


inline
Quad::Quad (const int i0, const int i1, const int i2, const int i3) {
  lines[0] = i0;
  lines[1] = i1;
  lines[2] = i2;
  lines[3] = i3;
};



inline
int Quad::line (const int i) const {
  Assert ((i>=0) && (i<4),
	  ExcRange(i));
  return lines[i];
};



inline
void Quad::set_line (const int i, const int index) {
  Assert ((i>=0) && (i<4),
	  ExcRange(i));
  lines[i] = index;
};



/*----------------------------   tria_quad.h     ---------------------------*/
/* end of #ifndef __tria_quad_H */
#endif
/*----------------------------   tria_quad.h     ---------------------------*/
