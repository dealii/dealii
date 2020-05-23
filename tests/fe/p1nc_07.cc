// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------



// check the P1NC element's gradient on a parallelogram

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_p1nc.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <string>

#include "../tests.h"



template <int dim>
Point<dim>
affine(const Point<dim> &p)
{
  Point<dim> q = p;
  if (dim >= 2)
    {
      q(0) = 2. * p(0) + p(1);
      q(1) = p(0) + 3. * p(1);
    }

  return q;
}



template <int dim>
void
check()
{
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation, 0., 1.);
  GridTools::transform(&affine<dim>, triangulation);

  FE_P1NC         fe;
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);

  QGauss<dim>   quadrature(3);
  FEValues<dim> fe_values(fe, quadrature, update_gradients | update_q_points);
  fe_values.reinit(dof_handler.begin_active());

  for (unsigned int q = 0; q < quadrature.size(); ++q)
    {
      deallog << "index=" << q << " position=" << fe_values.quadrature_point(q)
              << " values=";
      for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
        deallog << '[' << fe_values.shape_grad(i, q) << "] ";
      deallog << std::endl;
    }
}

int
main()
{
  initlog();
  deallog << std::setprecision(5) << std::fixed;
  deallog.depth_console(0);

  check<2>();
}
