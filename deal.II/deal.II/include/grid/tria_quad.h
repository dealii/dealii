//----------------------------  tria_quad.h  ---------------------------
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
//----------------------------  tria_quad.h  ---------------------------
#ifndef __deal2__tria_quad_h
#define __deal2__tria_quad_h


#include <base/config.h>
#include <base/exceptions.h>


/**
 *   @p{Quad}s denote the fundamental entities of triangulations in two dimensions
 *   and the boundaries of hexaeders in three dimensions. They are
 *   characterized by the (global) indices of the dour bounding lines.
 *
 *   A quad itself has one index, as far as the topological part handled in
 *   the triangulation is concerned: the index in the level
 *   it belongs to. The level index is implicitly given by the position
 *   in the @p{quads.quads} list attached to the information of each level
 *   of the triangulation.
 *
 *   @author Wolfgang Bangerth, 1998
 */
class Quad
{
  public:

				     /**
				      *  Construct a Quad with line
				      *  indices @p{i0} through @p{i3}. By default,
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
				      *  Return the index of line @p{i}=0
				      *  through 3.
				      */
    int line (const int i) const;
    
				     /**
				      *  Set the index of the @p{i}th line to
				      *  @p{index}. @p{i}=0 through 3.
				      */
    void set_line (const int i, const int index);
    
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
		    << "Indices for the line number must be 0, 1, 2 or 3, "
		    << "but you gave " << arg1); 
  protected:
    int lines[4];
};


/*----------------------------- Inline Function: Quad ------------------------*/


inline
Quad::Quad (const int i0, const int i1, const int i2, const int i3)
{
  lines[0] = i0;
  lines[1] = i1;
  lines[2] = i2;
  lines[3] = i3;
};



inline
int Quad::line (const int i) const
{
  Assert ((i>=0) && (i<4),
	  ExcRange(i));
  return lines[i];
};



inline
void Quad::set_line (const int i, const int index)
{
  Assert ((i>=0) && (i<4),
	  ExcRange(i));
  lines[i] = index;
};



inline
unsigned int
Quad::memory_consumption ()
{
  return sizeof(Quad);
};

#endif
