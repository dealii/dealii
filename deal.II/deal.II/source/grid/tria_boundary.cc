/* $Id$ */
/* Copyright W. Bangerth, University of Heidelberg, 1998 */


#include <grid/tria_boundary.h>
#include <cmath>



template <int dim>
Boundary<dim>::~Boundary ()
{};


template <int dim>
Point<dim> StraightBoundary<dim>::in_between (const PointArray &neighbors) const 
{
  Point<dim> p;
  for (int i=0; i<(1<<(dim-1)); ++i)
    p += *neighbors[i];
  p /= (1<<(dim-1))*1.0;
  return p;
};



template <int dim>
Point<dim> HyperBallBoundary<dim>::in_between (const PointArray &neighbors) const
{
  Point<dim> middle = StraightBoundary<dim>::in_between(neighbors);
  
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
