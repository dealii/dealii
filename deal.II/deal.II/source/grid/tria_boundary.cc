//----------------------------  tria_boundary.cc  ---------------------------
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
//----------------------------  tria_boundary.cc  ---------------------------


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
void
Boundary<dim>::get_intermediate_points_on_line (
  const typename Triangulation<dim>::line_iterator &,
  vector<Point<dim> > &) const
{
  Assert (false, ExcPureVirtualFunctionCalled());
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
  Assert (false, typename Boundary<dim>::ExcPureVirtualFunctionCalled());
  return Point<dim>();
};


#else


template <int dim>
Point<dim>
StraightBoundary<dim>::get_new_point_on_quad (const typename Triangulation<dim>::quad_iterator &quad) const 
{
  return (quad->vertex(0) + quad->vertex(1) +
	  quad->vertex(2) + quad->vertex(3) +
	  quad->line(0)->child(0)->vertex(1) +
	  quad->line(1)->child(0)->vertex(1) +
  	  quad->line(2)->child(0)->vertex(1) +
  	  quad->line(3)->child(0)->vertex(1)) / 8;
};

#endif


template <int dim>
void
StraightBoundary<dim>::get_intermediate_points_on_line (
  const typename Triangulation<dim>::line_iterator &line,
  vector<Point<dim> > &points) const
{
  const unsigned int n=points.size();
  Assert(n>0, ExcInternalError());
  const double part=1./(n+1);
  double position=part;
  for (unsigned int i=0; i<n; ++i, position+=part)
    points[i]=(1-position)*line->vertex(0) + position*line->vertex(1);
};


// explicit instantiations
template class Boundary<deal_II_dimension>;
template class StraightBoundary<deal_II_dimension>;

