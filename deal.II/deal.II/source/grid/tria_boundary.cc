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
  std::vector<Point<dim> > &) const
{
  Assert (false, ExcPureVirtualFunctionCalled());
};


template <int dim>
void
Boundary<dim>::get_intermediate_points_on_quad (
  const typename Triangulation<dim>::quad_iterator &,
  std::vector<Point<dim> > &) const
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
  Assert (false, typename Boundary<dim>::ExcFunctionNotUseful(dim));
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


#if deal_II_dimension < 2

template <int dim>
void
StraightBoundary<dim>::get_intermediate_points_on_line (
  const typename Triangulation<dim>::line_iterator &,
  typename std::vector<Point<dim> > &) const
{
  Assert(false, typename Boundary<dim>::ExcFunctionNotUseful(dim));
}


#else

template <int dim>
void
StraightBoundary<dim>::get_intermediate_points_on_line (
  const typename Triangulation<dim>::line_iterator &line,
  typename std::vector<Point<dim> > &points) const
{
  const unsigned int n=points.size();
  Assert(n>0, ExcInternalError());
  
  const double dx=1./(n+1);
  double x=dx;

  const Point<dim> vertices[2] = { line->vertex(0),
				   line->vertex(1) };
  
  for (unsigned int i=0; i<n; ++i, x+=dx)
    points[i] = (1-x)*vertices[0] + x*vertices[1];
};

#endif



#if deal_II_dimension < 3

template <int dim>
void
StraightBoundary<dim>::get_intermediate_points_on_quad (
  const typename Triangulation<dim>::quad_iterator &,
  typename std::vector<Point<dim> > &) const
{
  Assert(false, typename Boundary<dim>::ExcFunctionNotUseful(dim));
}

#else

template <int dim>
void
StraightBoundary<dim>::get_intermediate_points_on_quad (
  const typename Triangulation<dim>::quad_iterator &quad,
  typename std::vector<Point<dim> > &points) const
{
  const unsigned int n=points.size(),
		     m=static_cast<unsigned int>(sqrt(n));
				   // is n a square number
  Assert(m*m==n, ExcInternalError());

  const double ds=1./(m+1);
  double y=ds;

  const Point<dim> vertices[4] = { quad->vertex(0),
				   quad->vertex(1),
				   quad->vertex(2),
				   quad->vertex(3) };

  for (unsigned int i=0; i<n; ++i, y+=ds)
    {
      double x=ds;
      for (unsigned int j=0; j<n; ++j, x+=ds)
	points[i*m+j]=((1-x) * vertices[0] +
		       x     * vertices[1]) * (1-y) +
		      ((1-x) * vertices[3] +
		       x     * vertices[2]) * y;
    }
}

#endif


// explicit instantiations
template class Boundary<deal_II_dimension>;
template class StraightBoundary<deal_II_dimension>;

