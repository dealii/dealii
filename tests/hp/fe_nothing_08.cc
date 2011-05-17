//----------------------------  fe_nothing_08.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2009, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  fe_nothing_08.cc  ---------------------------


// test that FE_Nothing can be called with interpolate_boundary_values
// with scalar elements


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/numerics/vectors.h>


#include <fstream>


template <int dim>
void test ()
{
  Triangulation<dim>       triangulation;
  GridGenerator :: hyper_cube (triangulation, -1, 1);
  triangulation.refine_global(1);

  hp::FECollection<dim>    fe_collection;
  fe_collection.push_back (FE_Q<dim>(1));
  fe_collection.push_back (FE_Nothing<dim>());

  hp::DoFHandler<dim>      dof_handler (triangulation);

  typename hp::DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

  for(; cell != endc; cell++)
    if (cell->center()[0] > 0)
      cell->set_active_fe_index(1);
    else
      cell->set_active_fe_index(0);

  dof_handler.distribute_dofs (fe_collection);
  deallog << dof_handler.n_dofs() << " dofs" << std::endl;

  std::map<unsigned int, double> bv;
  VectorTools::interpolate_boundary_values (dof_handler,
					    0,
					    ZeroFunction<dim>(1),
					    bv);
  for (std::map<unsigned int, double>::iterator
	 p = bv.begin(); p!=bv.end(); ++p)
    deallog << p->first << ' ' << p->second << std::endl;
}



int main ()
{
  std::ofstream logfile("fe_nothing_08/output");
  logfile.precision(2);

  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<1> ();
  test<2> ();
  test<3> ();

  deallog << "OK" << std::endl;
}
