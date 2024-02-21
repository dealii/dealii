// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2011 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Check the documented property of the DGP element that the first shape
// function is constant in space

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <string>

#include "../tests.h"

#define PRECISION 3



template <int dim>
void
test(const unsigned int degree)
{
  FE_DGP<dim> fe(degree);
  deallog << fe.get_name() << std::endl;

  Triangulation<dim> tr;
  if (dim > 1)
    GridGenerator::hyper_ball(tr);
  else
    GridGenerator::hyper_cube(tr);

  QGauss<dim>   q(degree + 1);
  FEValues<dim> fe_values(fe, q, update_values);
  for (typename Triangulation<dim>::active_cell_iterator cell =
         tr.begin_active();
       cell != tr.end();
       ++cell)
    {
      fe_values.reinit(cell);
      for (const auto q : fe_values.quadrature_point_indices())
        AssertThrow(fe_values.shape_value(0, q) == 1, ExcInternalError());
    }

  deallog << "OK" << std::endl;
}


int
main()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(PRECISION);
  deallog << std::fixed;
  deallog.attach(logfile);

  for (unsigned int degree = 0; degree <= 4; ++degree)
    test<1>(degree);

  for (unsigned int degree = 0; degree <= 4; ++degree)
    test<2>(degree);

  for (unsigned int degree = 0; degree <= 4; ++degree)
    test<3>(degree);

  return 0;
}
