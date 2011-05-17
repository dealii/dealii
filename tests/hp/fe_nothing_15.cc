//----------------------------  fe_nothing_15.cc  ---------------------------
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
//----------------------------  fe_nothing_15.cc  ---------------------------


// make sure we can extract shape functions from an FE_Nothing.


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
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
#include <deal.II/fe/fe_system.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/lac/constraint_matrix.h>


#include <fstream>

template <int dim>
void test ()
{
  Triangulation<dim>       triangulation;
  GridGenerator :: hyper_cube (triangulation, -0.5, 0.5);

  FESystem<dim> fe (FE_Q<dim>(1), 1,
		    FE_Nothing<dim>(), 1);
  DoFHandler<dim> dof_handler (triangulation);
  dof_handler.distribute_dofs (fe);

  FEValues<dim> fe_values(fe, QGauss<dim>(2), update_values);
  FEValuesExtractors::Scalar nothing(1);
  fe_values.reinit (dof_handler.begin_active());
  for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
    for (unsigned int q=0; q<fe_values.n_quadrature_points; ++q)
      deallog << "i=" << i
	      << ", q=" << q
	      << ", value="
	      << fe_values[nothing].value(i,q)
	      << std::endl;
}



int main ()
{
  std::ofstream logfile("fe_nothing_15/output");
  logfile.precision(2);

  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  deallog << "Try dim == 1" << std::flush << std::endl;
  test<1> ();

  deallog << "Try dim == 2" << std::flush << std::endl;
  test<2> ();

  deallog << "Try dim == 3" << std::flush << std::endl;
  test<3> ();

  deallog << "OK" << std::endl;
}
