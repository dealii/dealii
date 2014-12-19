// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
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


// turns out that cylinder_shell_01 wasn't enough: measure just takes the
// positive value it computes. we need to ask a FEValues object for the JxW
// values to get the sign of the Jacobian


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <fstream>



int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  deallog << std::setprecision (2);

  // generate a hyperball in 3d
  Triangulation<3> tria;
  GridGenerator::cylinder_shell (tria, 1, .8, 1);

  FE_Q<3> fe(1);
  DoFHandler<3> dof_handler (tria);
  dof_handler.distribute_dofs (fe);

  QMidpoint<3> q;
  FEValues<3> fe_values (fe, q, update_JxW_values);

  // make sure that all cells have positive
  // volume
  for (DoFHandler<3>::active_cell_iterator cell=dof_handler.begin_active();
       cell!=dof_handler.end(); ++cell)
    {
      fe_values.reinit (cell);
      deallog << cell << ' ' << fe_values.JxW(0) << std::endl;

      Assert (fe_values.JxW(0) > 0, ExcInternalError());
    }
}
