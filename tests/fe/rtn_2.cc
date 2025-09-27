// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// rt(2) had some problems with shape functions (tested in rt_7.cc). because
// it is so simple, just run the same test for the RT-Nodal element as well

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/vector.h>

#include <string>
#include <vector>

#include "../tests.h"

#define PRECISION 8



template <int dim>
void
plot_shape_functions(const unsigned int degree)
{
  FE_RaviartThomasNodal<dim> fe_rt(degree);
  Triangulation<dim>         tr;
  GridGenerator::hyper_cube(tr, 0., 1.);

  DoFHandler<dim>                         dof(tr);
  typename DoFHandler<dim>::cell_iterator c = dof.begin();
  dof.distribute_dofs(fe_rt);

  QTrapezoid<1>      q_trapez;
  const unsigned int div = 10;
  QIterated<dim>     q(q_trapez, div);
  FEValues<dim>      fe(fe_rt,
                   q,
                   update_values | update_gradients | update_quadrature_points);
  fe.reinit(c);

  for (unsigned int q_point = 0; q_point < q.size(); ++q_point)
    {
      if (q_point % QIterated<1>(q_trapez, div).size() == 0)
        deallog << std::endl;

      deallog << fe.quadrature_point(q_point) << ' ';

      for (unsigned int i = 0; i < fe_rt.dofs_per_cell; ++i)
        for (unsigned int c = 0; c < fe.get_fe().n_components(); ++c)
          deallog << ' ' << fe.shape_value_component(i, q_point, c);

      deallog << std::endl;
    }
}


int
main()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(PRECISION) << std::fixed;
  deallog.attach(logfile);

  plot_shape_functions<2>(2);

  return 0;
}
