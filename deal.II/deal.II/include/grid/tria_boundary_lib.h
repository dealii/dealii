/*----------------------------   tria_boundary_lib.h     ---------------------------*/
/*      $Id$                 */
#ifndef __tria_boundary_lib_H
#define __tria_boundary_lib_H
/*----------------------------   tria_boundary_lib.h     ---------------------------*/


#include <grid/tria_boundary.h>




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
 *
 * @author Wolfgang Bangerth, 1998
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
				      * Refer to the general documentation of
				      * this class and the documentation of the
				      * base class.
				      */
    virtual Point<dim>
    get_new_point_on_line (const typename Triangulation<dim>::line_iterator &line) const;

				     /**
				      * Refer to the general documentation of
				      * this class and the documentation of the
				      * base class.
				      */
    virtual Point<dim>
    get_new_point_on_quad (const typename Triangulation<dim>::quad_iterator &quad) const;


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

    


/**
 * Class describing the boundaries of a hyper shell. Only the center
 * of the two spheres needs to be given, the radii of inner and outer
 * sphere are computed automatically upon calling one of the virtual
 * functions.
 *
 * @author Wolfgang Bangerth, 1999
 */
template <int dim>
class HyperShellBoundary : public StraightBoundary<dim> 
{
  public:
				     /**
				      * Constructor. The center of the
				      * spheres defaults to the
				      * origin.
				      */
    HyperShellBoundary (const Point<dim> &center = Point<dim>());
    
				     /**
				      * Construct a new point on a line.
				      */
    virtual Point<dim>
    get_new_point_on_line (const typename Triangulation<dim>::line_iterator &line) const;  
    
				     /**
				      * Construct a new point on a quad.
				      */
    virtual Point<dim>
    get_new_point_on_quad (const typename Triangulation<dim>::quad_iterator &quad) const;
    
  private:
				     /**
				      * Store the center of the spheres.
				      */
    const Point<dim> center;
};




/*----------------------------   tria_boundary_lib.h     ---------------------------*/
/* end of #ifndef __tria_boundary_lib_H */
#endif
/*----------------------------   tria_boundary_lib.h     ---------------------------*/
