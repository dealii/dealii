//----------------------------  fe_nothing_16.cc  ---------------------------
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
//----------------------------  fe_nothing_16.cc  ---------------------------


// test a problem we used to have: FESystem would delete an internal object
// after reinitialization for the first time if it determined that it was no
// longer necessary. yet, somehow, it was still referenced. the point seems to
// have been that the base element always only had update_default for the
// values that need to be updated on each cell, which is rather uncommon (the
// base element is FE_Nothing)
//
// an extract of this bug is fe/crash_01


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

  QGauss<dim-1> q(2);
  FEFaceValues<dim> fe_values(fe, q, update_values);
  FEValuesExtractors::Scalar nothing(1);
  fe_values.reinit (dof_handler.begin_active(), 0);

				   // the following (second) call to reinit
				   // used to abort
  fe_values.reinit (dof_handler.begin_active(), 1);
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
  std::ofstream logfile("fe_nothing_16/output");
  logfile.precision(2);

  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<2> ();

  deallog << "OK" << std::endl;
}
