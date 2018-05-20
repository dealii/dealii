// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2017 by the deal.II authors
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

// check the P1NC element's gradient on a square

#include "../tests.h"
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/fe/fe_p1nc.h>

#include <string>

template <int dim>
void
check()
{
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation, 0, 5);

  FE_P1NC         fe;
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);

  QGauss<dim>   quadrature(3);
  FEValues<dim> fe_values(fe, quadrature, update_gradients | update_q_points);
  fe_values.reinit(dof_handler.begin_active());

  for(unsigned int q = 0; q < quadrature.size(); ++q)
    {
      deallog << "index=" << q << " position=" << fe_values.quadrature_point(q)
              << " values=";
      for(unsigned int i = 0; i < fe.dofs_per_cell; ++i)
        deallog << '[' << fe_values.shape_grad(i, q) << "] ";
      deallog << std::endl;
    }
}

int
main()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(5);
  deallog << std::fixed;
  deallog.attach(logfile);
  deallog.depth_console(0);

  check<2>();
}
