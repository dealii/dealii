/*----------------------------   tria_hex.h     ---------------------------*/
/*      $Id$                 */
#ifndef __tria_hex_H
#define __tria_hex_H
/*----------------------------   tria_hex.h     ---------------------------*/



#include <base/exceptions.h>

/**
 *   #Hexahedron#s denote the fundamental entities of triangulations in three
 *   dimensions. They are
 *   characterized by the (global) indices of the six bounding quadrilaterals.
 *
 *   A heaxhedron itself has one index, as far as the topological part handled in
 *   the triangulation is concerned: the index in the level
 *   it belongs to. The level index is implicitely given by the position
 *   in the #hexes.hexes# list attached to the information of each level
 *   of the triangulation.
 *
 *   @author Wolfgang Bangerth, 1998
 */
class Hexahedron {
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
				      *  through 5.
				      */
    int quad (const int i) const;
    
				     /**
				      *  Set the index of the #i#th quad to
				      *  #index#. #i#=0 through 5.
				      */
    void set_quad (const int i, const int index);
    
				     /**
				      *  Exception
				      */ 
    DeclException1 (ExcRange,
		    int,
		    << "Indices for the quad number must be 0 through 5, "
		    << "but you gave " << arg1); 
  protected:
      int quads[6];
};




/*----------------------------- Inline Function: Hexahedron ------------------------*/


inline
Hexahedron::Hexahedron (const int i0, const int i1,
			const int i2, const int i3,
			const int i4, const int i5) {
  quads[0] = i0;
  quads[1] = i1;
  quads[2] = i2;
  quads[3] = i3;
  quads[4] = i4;
  quads[5] = i5;
};



inline
int Hexahedron::quad (const int i) const {
  Assert ((i>=0) && (i<6),
	  ExcRange(i));
  return quads[i];
};



inline
void Hexahedron::set_quad (const int i, const int index) {
  Assert ((i>=0) && (i<6),
	  ExcRange(i));
  quads[i] = index;
};




/*----------------------------   tria_hex.h     ---------------------------*/
/* end of #ifndef __tria_hex_H */
#endif
/*----------------------------   tria_hex.h     ---------------------------*/
