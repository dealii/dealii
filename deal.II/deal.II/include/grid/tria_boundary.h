/*----------------------------   boundary-function.h     ---------------------------*/
/*      $Id$                 */
#ifndef __tria_boundary_H
#define __tria_boundary_H
/*----------------------------   boundary-function.h     ---------------------------*/

#include <grid/point.h>

/**
    This class is used to represent a boundary to a triangulation.
    When a triangulation creates a new vertex on the boundary of the
    domain, it determines the new vertex' coordinates through the
    following code (here in two dimensions):
    \begin{verbatim}
      ...
      const Point<2> *neighbors[2] = {&neighbor1, &neighbor2};
      Point<2> new_vertex = boundary.in_between (neighbors);
      ...
    \end{verbatim}
    #neighbor1# and #neighbor2# are the two vertices bounding the old
    line on the boundary, which is to be subdivided. #boundary# is an
    object of type #Boundary<dim>#.

    In 3D, a new vertex may be placed on the middle of a line or on
    the middle of a side. In the both cases, an array with four points
    has to be passed to #in_between#; in the latter case the two end
    points of the line have to be given consecutively twice, as
    elements 0 and 1, and 2 and 3, respectively.
    
    There is a specialisation, #StraightBoundary<dim>#, which places
    the new point right into the middle of the given points.
    */
template <int dim>
class Boundary {
  public:
				     /**
				      *  Typedef an array of the needed number
				      *  of old points.
				      */
    typedef const Point<dim>* PointArray[1<<(dim-1)];

				     /**
				      *  This function calculates the position
				      *  of the new vertex.
				      */
    virtual Point<dim> in_between (const PointArray &neighbors) const = 0;
};





/**
    Specialisation of \Ref{Boundary}<dim>, which places the new point right
    into the middle of the given points. The middle is defined as the
    arithmetic mean of the points.
    */
template <int dim>
class StraightBoundary : public Boundary<dim> {
  public:
				     /**
				      *  This function calculates the position
				      *  of the new vertex.
				      */
    virtual Point<dim> in_between (const PointArray &neighbors) const {
      Point<dim> p;
      for (int i=0; i<(1<<(dim-1)); ++i)
	p += *neighbors[i];
      p /= (1<<(dim-1))*1.0;
      return p;
    };
};




/*----------------------------   boundary-function.h     ---------------------------*/
/* end of #ifndef __tria_boundary_H */
#endif
/*----------------------------   boundary-function.h     ---------------------------*/
