/* $Id$ */
/* Copyright W. Bangerth, University of Heidelberg, 1998 */

#include <grid/tria_boundary_lib.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <grid/tria_accessor.h>
#include <cmath>




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



#if deal_II_dimension == 1

template <>
Point<1>
HyperBallBoundary<1>::
get_new_point_on_quad (const Triangulation<1>::quad_iterator &) const
{
  Assert (false, ExcInternalError());
  return Point<1>();
};

#endif



template <int dim>
Point<dim>
HyperBallBoundary<dim>::
get_new_point_on_quad (const typename Triangulation<dim>::quad_iterator &quad) const
{
  Point<dim> middle = StraightBoundary<dim>::get_new_point_on_quad (quad);
  
  middle -= center;
				   // project to boundary
  middle *= radius / sqrt(middle.square());
  
  middle += center;
  return middle;
};





template <int dim>
HyperShellBoundary<dim>::HyperShellBoundary (const Point<dim> &center) :
		center (center) 
{};


template <int dim>
Point<dim>
HyperShellBoundary<dim>::
get_new_point_on_line (const typename Triangulation<dim>::line_iterator &line) const 
{
  const Point<dim> middle = StraightBoundary<dim>::get_new_point_on_line (line);
				   // compute the position of the points relative to the origin
  const Point<dim> middle_relative = middle - center,
		   vertex_relative = line->vertex(0) - center;
  
				   // take vertex(0) to gauge the
				   // radius corresponding to the line
				   // under consideration
  const double radius = sqrt(vertex_relative.square());

				   // scale and shift back to the
				   // original coordinate system
  return (middle_relative * (radius / sqrt(middle_relative.square()))) + center;
};



#if deal_II_dimension == 1

template <>
Point<1>
HyperShellBoundary<1>::
get_new_point_on_quad (const Triangulation<1>::quad_iterator &) const
{
  Assert (false, ExcInternalError());
  return Point<1>();
};

#endif



template <int dim>
Point<dim>
HyperShellBoundary<dim>::
get_new_point_on_quad (const typename Triangulation<dim>::quad_iterator &quad) const
{
  const Point<dim> middle = StraightBoundary<dim>::get_new_point_on_quad (quad);
				   // compute the position of the points relative to the origin
  const Point<dim> middle_relative = middle - center,
		   vertex_relative = quad->vertex(0) - center;
  
				   // take vertex(0) to gauge the
				   // radius corresponding to the line
				   // under consideration
  const double radius = sqrt(vertex_relative.square());

				   // scale and shift back to the
				   // original coordinate system
  return (middle_relative * (radius / sqrt(middle_relative.square()))) + center;
};





// explicit instantiations
template class HyperBallBoundary<deal_II_dimension>;
template class HyperShellBoundary<deal_II_dimension>;
