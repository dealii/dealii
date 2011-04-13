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
#include <base/logstream.h>
#include <base/quadrature_lib.h>
#include <fe/fe_nothing.h>
#include <fe/fe_q.h>
#include <hp/fe_collection.h>
#include <hp/dof_handler.h>
#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/grid_refinement.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <fe/fe_q.h>
#include <fe/fe_system.h>
#include <hp/dof_handler.h>
#include <hp/fe_values.h>
#include <lac/constraint_matrix.h>


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
