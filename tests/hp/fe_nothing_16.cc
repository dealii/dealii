// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2014 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



// test a problem we used to have: FESystem would delete an internal object
// after reinitialization for the first time if it determined that it was no
// longer necessary. yet, somehow, it was still referenced. the point seems to
// have been that the base element always only had update_default for the
// values that need to be updated on each cell, which is rather uncommon (the
// base element is FE_Nothing)
//
// an extract of this bug is fe/crash_01


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
  std::ofstream logfile("output");
  logfile.precision(2);

  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<2> ();

  deallog << "OK" << std::endl;
}
