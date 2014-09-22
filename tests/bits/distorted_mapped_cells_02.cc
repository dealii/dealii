// ---------------------------------------------------------------------
//
// Copyright (C) 2013 by the deal.II authors
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


// Test that an assertion is thrown when MappingQ distorts the geometry too
// much on a thin hyper shell

#include "../tests.h"
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>


void test()
{
  const int dim = 2;
  Triangulation<dim> tria(Triangulation<dim>::none,true);

  // shell is so narrow that MappingQ(2) distorts the cell
  GridGenerator::quarter_hyper_shell (tria, Point<dim>(), 0.95, 1, 1);
  static HyperShellBoundary<dim> boundary;
  tria.set_boundary (0, boundary);

  FE_Nothing<dim> dummy;
  MappingQ<dim> mapping(2);
  QGauss<dim> quad(2);
  FEValues<dim> fe_val (mapping, dummy, quad, update_JxW_values);
  double integral = 0.;
  /*typename*/ Triangulation<dim>::active_cell_iterator
  cell = tria.begin_active(), endc = tria.end();
  for ( ; cell != endc; ++cell)
    {
      try
        {
          fe_val.reinit (cell);
          for (unsigned int q=0; q<quad.size(); ++q)
            integral += fe_val.JxW(q);
        }
      catch (ExceptionBase &e)
        {
          deallog << e.get_exc_name() << std::endl;
        }
    }
  deallog << "Integral = " << integral << std::endl;
}


int
main()
{
  deal_II_exceptions::disable_abort_on_exception();

  std::ofstream logfile ("output");
  deallog << std::setprecision(4) << std::fixed;
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test();

  return 0;
}



