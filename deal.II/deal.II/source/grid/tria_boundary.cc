/* $Id$ */
/* Copyright W. Bangerth, University of Heidelberg, 1998 */


#include <grid/tria_boundary.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <grid/tria_accessor.h>
#include <cmath>



template <int dim>
Boundary<dim>::~Boundary ()
{};




template <int dim>
Point<dim>
Boundary<dim>::get_new_point_on_quad (const typename Triangulation<dim>::quad_iterator &) const 
{
  Assert (false, ExcPureVirtualFunctionCalled());
  return Point<dim>();
};



template <int dim>
Point<dim>
StraightBoundary<dim>::get_new_point_on_line (const typename Triangulation<dim>::line_iterator &line) const 
{
  return (line->vertex(0) + line->vertex(1)) / 2;
};



#if deal_II_dimension < 3

template <int dim>
Point<dim>
StraightBoundary<dim>::get_new_point_on_quad (const typename Triangulation<dim>::quad_iterator &) const 
{
  Assert (false, ExcPureVirtualFunctionCalled());
  return Point<dim>();
};


#else


template <int dim>
Point<dim>
StraightBoundary<dim>::get_new_point_on_quad (const typename Triangulation<dim>::quad_iterator &quad) const 
{
  return (quad->vertex(0) + quad->vertex(1) +
	  quad->vertex(2) + quad->vertex(3)) / 4;
};

#endif



template <int dim>
Point<dim>
HyperBallBoundary<dim>::get_new_point_on_line (const typename Triangulation<dim>::line_iterator &line) const
{
  Point<dim> middle = StraightBoundary<dim>::get_new_point_on_line (line);
  
  middle -= center;
				   // project to boundary
  middle *= radius / sqrt(middle.square());
  
  middle += center;
  return middle;
};



template <int dim>
Point<dim>
HyperBallBoundary<dim>::get_new_point_on_quad (const typename Triangulation<dim>::quad_iterator &quad) const
{
  Point<dim> middle = StraightBoundary<dim>::get_new_point_on_quad (quad);
  
  middle -= center;
				   // project to boundary
  middle *= radius / sqrt(middle.square());
  
  middle += center;
  return middle;
};




// explicit instantiations
template class Boundary<deal_II_dimension>;
template class StraightBoundary<deal_II_dimension>;
template class HyperBallBoundary<deal_II_dimension>;
