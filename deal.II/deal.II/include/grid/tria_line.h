//----------------------------  tria_line.h  ---------------------------
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
//----------------------------  tria_line.h  ---------------------------
#ifndef __deal2__tria_line_h
#define __deal2__tria_line_h


#include <base/exceptions.h>


/**
 *   Lines denote the boundaries of quads and the edges of hexaeders. They are
 *   characterized by the (global) indices of the endpoints.
 *
 *   A line itself has one index, as far as the topological part handled in
 *   the triangulation is concerned: the index in the level
 *   it belongs to. The level index is implicitely given by the position
 *   in the @p{lines.lines} list attached to the information of each level.
 *
 *   @author Wolfgang Bangerth, 1998
 */
class Line
{
  public:
				     /**
				      *  Construct a line with end point
				      *  indices @p{i0} and @p{i1}. By default,
				      *  indices are set to -1, i.e. an
				      *  invalid value.
				      */
    Line (const int i0 = -1,
	  const int i1 = -1);
    
				     /**
				      *  Return the index of end point @p{i}=0
				      *  or 1.
				      */
    int vertex (const int i) const;
    
				     /**
				      *  Set the index of the @p{i}th vertex to
				      *  @p{index}. @p{i}=0 or 1.
				      */
    void set_vertex (const int i, const int index);
    
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
		    << "Indices for the point number must be 0 or 1, "
		    << "but you gave " << arg1);
      
  protected:
				     /**
				      *  Global indices of the two end points.
				      */
      int end_points[2];
};


/*----------------------------- Inline Function: Line ------------------------*/


inline                               // truly in-line here!
Line::Line (const int i0, const int i1)
{
  end_points[0] = i0;
  end_points[1] = i1;
};



inline
int Line::vertex (const int i) const
{
  Assert ((i==0) || (i==1),
	  ExcRange(i));
  return end_points[i];
};



inline
void Line::set_vertex (const int i, const int index)
{
  Assert ((i==0) || (i==1),
	  ExcRange(i));
  end_points[i] = index;
};



inline
unsigned int
Line::memory_consumption ()
{
  return sizeof(Line);
};


#endif
