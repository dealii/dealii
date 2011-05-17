//----------------------------  crash_04.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  crash_04.cc  ---------------------------


// trigger a problem in DerivativeApproximation with memory handling


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/numerics/derivative_approximation.h>

#include <fstream>


template <int dim>
void test ()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global (1);

  hp::FECollection<dim> fe_collection (FE_DGQ<dim> (1));

  hp::DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs (fe_collection);

  Vector<double> v(dof_handler.n_dofs());
  Vector<float> e(tria.n_active_cells());
  DerivativeApproximation::approximate_gradient (dof_handler, v, e);
}



int main ()
{
  std::ofstream logfile("crash_04/output");
  logfile.precision(2);
  
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  
  test<1> ();
  test<2> ();
  test<3> ();
  
  deallog << "OK" << std::endl;
}
