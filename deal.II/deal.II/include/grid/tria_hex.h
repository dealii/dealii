//----------------------------  tria_hex.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  tria_hex.h  ---------------------------
#ifndef __deal2__tria_hex_h
#define __deal2__tria_hex_h


#include <base/exceptions.h>

/**
 *   @p{Hexahedron}s denote the fundamental entities of triangulations in three
 *   dimensions. They are
 *   characterized by the (global) indices of the six bounding quadrilaterals.
 *
 *   A heaxhedron itself has one index, as far as the topological part handled in
 *   the triangulation is concerned: the index in the level
 *   it belongs to. The level index is implicitely given by the position
 *   in the @p{hexes.hexes} list attached to the information of each level
 *   of the triangulation.
 *
 *   For conventions on the numbering of faces, lines and vertices within a
 *   hexahedron, refer to the documentation of the @ref{Triangulation} class.
 *
 *   @author Wolfgang Bangerth, 1998
 */
class Hexahedron
{
  public:

				     /**
				      *  Construct a Hex with quad
				      *  indices @p{i0} through @p{i5}. By default,
				      *  indices are set to -1, i.e. an
				      *  invalid value.
				      *
				      *  No convention is set as of now on the
				      *  order of quads
				      */
    Hexahedron (const int i0 = -1,
		const int i1 = -1,
		const int i2 = -1,
		const int i3 = -1,
		const int i4 = -1,
		const int i5 = -1);
    
				     /**
				      *  Return the index of quad @p{i}=0
				      *  through 5.
				      */
    int quad (const int i) const;
    
				     /**
				      *  Set the index of the @p{i}th quad to
				      *  @p{index}. @p{i}=0 through 5.
				      */
    void set_quad (const int i, const int index);
    
				     /**
				      * Determine an estimate for the
				      * memory consumption (in bytes)
				      * of this object.
				      */
    static unsigned int memory_consumption ();

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
			const int i4, const int i5)
{
  quads[0] = i0;
  quads[1] = i1;
  quads[2] = i2;
  quads[3] = i3;
  quads[4] = i4;
  quads[5] = i5;
};



inline
int Hexahedron::quad (const int i) const
{
  Assert ((i>=0) && (i<6),
	  ExcRange(i));
  return quads[i];
};



inline
void Hexahedron::set_quad (const int i, const int index)
{
  Assert ((i>=0) && (i<6),
	  ExcRange(i));
  quads[i] = index;
};



inline
unsigned int
Hexahedron::memory_consumption ()
{
  return sizeof(Hexahedron);
};

#endif
