/*----------------------------   tria_line.h     ---------------------------*/
/*      $Id$                 */
#ifndef __tria_line_H
#define __tria_line_H
/*----------------------------   tria_line.h     ---------------------------*/

#include <base/exceptions.h>


/**
 *   Lines denote the boundaries of quads and the edges of hexaeders. They are
 *   characterized by the (global) indices of the endpoints.
 *
 *   A line itself has one index, as far as the topological part handled in
 *   the triangulation is concerned: the index in the level
 *   it belongs to. The level index is implicitely given by the position
 *   in the #lines.lines# list attached to the information of each level.
 *
 *   @author Wolfgang Bangerth, 1998
 */
class Line {
  public:
				     /**
				      *  Construct a line with end point
				      *  indices #i0# and #i1#. By default,
				      *  indices are set to -1, i.e. an
				      *  invalid value.
				      */
    Line (const int i0 = -1,
	  const int i1 = -1);
    
				     /**
				      *  Return the index of end point #i#=0
				      *  or 1.
				      */
    int vertex (const int i) const;
    
				     /**
				      *  Set the index of the #i#th vertex to
				      *  #index#. #i#=0 or 1.
				      */
    void set_vertex (const int i, const int index);
    
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


inline                               // wahrlich inline hier!
Line::Line (const int i0, const int i1) {
  end_points[0] = i0;
  end_points[1] = i1;
};


inline
int Line::vertex (const int i) const {
  Assert ((i==0) || (i==1),
	  ExcRange(i));
  return end_points[i];
};


inline
void Line::set_vertex (const int i, const int index) {
  Assert ((i==0) || (i==1),
	  ExcRange(i));
  end_points[i] = index;
};



/*----------------------------   tria_line.h     ---------------------------*/
/* end of #ifndef __tria_line_H */
#endif
/*----------------------------   tria_line.h     ---------------------------*/
