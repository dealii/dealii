//----------------------------  support_point_map.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  support_point_map.cc  ---------------------------


/* Author: Wolfgang Bangerth, University of Heidelberg, 2001 */



#include <base/logstream.h>
#include <base/function_lib.h>
#include <lac/vector.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <grid/tria_boundary_lib.h>
#include <dofs/dof_accessor.h>
#include <grid/grid_generator.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_tools.h>
#include <fe/fe_q.h>
#include <fe/fe_dgq.h>
#include <fe/fe_system.h>
#include <fe/mapping_q.h>

#include <fstream>


template <int dim>
struct PointComp 
{
    bool operator () (const Point<dim> &, const Point<dim> &) const;
};


template <>
bool PointComp<1>::operator () (const Point<1> &p1,
				const Point<1> &p2) const
{
  return p1(0) < p2(0);
};


// have somewhat weird orderings in 2d and 3d
template <>
bool PointComp<2>::operator () (const Point<2> &p1,
				const Point<2> &p2) const
{
  return ((p1(0)+p1(1) < p2(0)+p2(1)) ||
	  ((p1(0)+p1(1) == p2(0)+p2(1)) &&
	   (p1(0)-p1(1) < p2(0)-p2(1))));
};



template <>
bool PointComp<3>::operator () (const Point<3> &p1,
				const Point<3> &p2) const
{
  return ((p1(2) < p2(2)) ||
	  (p1(2) == p2(2)) &&
	  ((p1(0)+p1(1) < p2(0)+p2(1)) ||
	   ((p1(0)+p1(1) == p2(0)+p2(1)) &&
	    (p1(0)-p1(1) < p2(0)-p2(1)))));
};



template <int dim>
void
check ()
{
  Triangulation<dim> tr;  
  if (dim==2)
    {
      GridGenerator::hyper_ball(tr);
    }
  else
    GridGenerator::hyper_cube(tr, -1./sqrt(dim),1./sqrt(dim));
  if (dim != 1)
    {
      static const HyperBallBoundary<dim> boundary;
      tr.set_boundary (0, boundary);
    };
  tr.refine_global (1);
  tr.begin_active()->set_refine_flag ();
  tr.execute_coarsening_and_refinement ();
  if (dim==1)
    tr.refine_global(2);

  MappingQ<dim> mapping(3);
  
  FESystem<dim> element (FE_Q<dim>(2), 1, FE_DGQ<dim>(1), 1);
  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(element);

				   // get the forward map
  std::vector<Point<dim> > support_points (dof.n_dofs());
  DoFTools::map_dofs_to_support_points (mapping, dof, support_points);
  for (unsigned int i=0; i<dof.n_dofs(); ++i)
    deallog << i << ": " << support_points[i] << std::endl;

				   // now get the backward map
  std::map<Point<dim>,unsigned int,PointComp<dim> > point_map;
  DoFTools::map_support_points_to_dofs (mapping, dof, point_map);
  typename std::map<Point<dim>,unsigned int,PointComp<dim> >::const_iterator
    i = point_map.begin(),
    e = point_map.end();
  for (; i!=e; ++i)
    deallog << i->first << ',' << i->second << std::endl;
}


int main ()
{
  ofstream logfile ("support_point_map.output");
  logfile.precision (2);
  logfile.setf(ios::fixed);  
  deallog.attach(logfile);
  deallog.depth_console (0);

  deallog.push ("1d");
  check<1> ();
  deallog.pop ();
  deallog.push ("2d");
  check<2> ();
  deallog.pop ();
  deallog.push ("3d");
  check<3> ();
  deallog.pop ();
}
