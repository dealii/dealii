/*----------------------------   boundary-function.h     ---------------------------*/
/*      $Id$                 */
/*      Copyright W. Bangerth, University of Heidelberg, 1998 */
#ifndef __tria_boundary_H
#define __tria_boundary_H
/*----------------------------   boundary-function.h     ---------------------------*/

#include <base/point.h>
#include <grid/geometry_info.h>



/**
 * Workaround for a bug in egcs snapshot 1998/08/03.
 */
template <int dim> struct BoundaryHelper;
template <> struct BoundaryHelper<2> {
    typedef const Point<2> *PointArray[GeometryInfo<2>::vertices_per_face];
};
    


/**
 *   This class is used to represent a boundary to a triangulation.
 *   When a triangulation creates a new vertex on the boundary of the
 *   domain, it determines the new vertex' coordinates through the
 *   following code (here in two dimensions):
 *   \begin{verbatim}
 *     ...
 *     const Point<2> *neighbors[2] = {&neighbor1, &neighbor2};
 *     Point<2> new_vertex = boundary.in_between (neighbors);
 *     ...
 *   \end{verbatim}
 *   #neighbor1# and #neighbor2# are the two vertices bounding the old
 *   line on the boundary, which is to be subdivided. #boundary# is an
 *   object of type #Boundary<dim>#.
 *
 *   In 3D, a new vertex may be placed on the middle of a line or on
 *   the middle of a side. In the both cases, an array with four points
 *   has to be passed to #in_between#; in the latter case the two end
 *   points of the line have to be given consecutively twice, as
 *   elements 0 and 1, and 2 and 3, respectively.
 *   
 *   There are specialisations, #StraightBoundary<dim>#, which places
 *   the new point right into the middle of the given points, and
 *   #HyperBallBoundary<dim># creating a hyperball with given radius
 *   around a given center point.
 *
 *   @author Wolfgang Bangerth, 1998
 */
template <int dim>
class Boundary {
  public:
				     /**
				      *  Typedef an array of the needed number
				      *  of old points.
				      */
    typedef typename BoundaryHelper<dim>::PointArray PointArray;
    
				     /**
				      *  This function calculates the position
				      *  of the new vertex.
				      */
    virtual Point<dim> in_between (const PointArray &neighbors) const = 0;
};





/**
 *   Specialisation of \Ref{Boundary}<dim>, which places the new point right
 *   into the middle of the given points. The middle is defined as the
 *   arithmetic mean of the points.
 *
 *   This class does not really describe a boundary in the usual sense. By
 *   placing new points in teh middle of old ones, it rather assumes that the
 *   boundary of the domain is given by the polygon/polyhedron defined by the
 *   boundary of the initial coarse triangulation.
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




/**
 *   Specialisation of \Ref{Boundary}<dim>, which places the new point on
 *   the boundary of a ball in arbitrary dimension. It works by projecting
 *   the point in the middle of the old points onto the ball. The middle is
 *   defined as the arithmetic mean of the points. 
 *
 *   The center of the ball and its radius may be given upon construction of
 *   an object of this type. They default to the origin and a radius of 1.0.
 *
 *   This class is derived from #StraightBoundary# rather than from
 *   #Boundary#, which would seem natural, since this way we can use the
 *   #StraightBoundary<dim>::in_between(neighbors)# function.
 */
template <int dim>
class HyperBallBoundary : public StraightBoundary<dim> {
  public:
				     /**
				      * Constructor
				      */
    HyperBallBoundary (const Point<dim> p=Point<dim>(), const double radius=1.0) :
		    center(p), radius(radius) {};

				     /**
				      *  This function calculates the position
				      *  of the new vertex.
				      */
    virtual Point<dim> in_between (const PointArray &neighbors) const {
      Point<dim> middle = StraightBoundary<dim>::in_between(neighbors);

      middle -= center;
				       // project to boundary
      middle *= radius / sqrt(middle.square());

      middle += center;
      return middle;
    };


  private:
				     /**
				      * Center point of the hyperball.
				      */
    const Point<dim> center;

				     /**
				      * Radius of the hyperball.
				      */
    const double radius;
};

    


/*----------------------------   boundary-function.h     ---------------------------*/
/* end of #ifndef __tria_boundary_H */
#endif
/*----------------------------   boundary-function.h     ---------------------------*/
