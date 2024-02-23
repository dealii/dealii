// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2003 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Show the shape functions of the BDM element on the unit cell
// Plots are gnuplot compatible if lines with desired prefix are selected.

#include <deal.II/fe/fe_bdm.h>

#include <sstream>
#include <string>
#include <vector>

#include "../tests.h"

#define PRECISION 8



template <int dim>
inline void
plot_shape_functions(const unsigned int degree)
{
  FE_BDM<dim> fe_bdm(degree);
  deallog.push(fe_bdm.get_name());

  const unsigned int div = 2;
  for (unsigned int mz = 0; mz <= ((dim > 2) ? div : 0); ++mz)
    for (unsigned int my = 0; my <= ((dim > 1) ? div : 0); ++my)
      {
        for (unsigned int mx = 0; mx <= div; ++mx)
          {
            const Point<dim> p =
              (dim == 2 ?
                 Point<dim>(1. * mx / div, 1. * my / div) :
                 Point<dim>(1. * mx / div, 1. * my / div, 1. * mz / div));

            // Lines with function
            // values contain
            // quadrature point and one
            // vector of dim entries
            // for each shape function
            deallog << "value " << p;
            for (unsigned int i = 0; i < fe_bdm.dofs_per_cell; ++i)
              {
                for (unsigned int c = 0; c < dim; ++c)
                  deallog << ' ' << fe_bdm.shape_value_component(i, p, c);
                deallog << "  ";
              }
            deallog << std::endl << "grad " << p;
            for (unsigned int i = 0; i < fe_bdm.dofs_per_cell; ++i)
              {
                for (unsigned int c = 0; c < dim; ++c)
                  {
                    deallog << ' ';
                    for (unsigned int d = 0; d < dim; ++d)
                      deallog << ' ' << fe_bdm.shape_grad_component(i, p, c)[d];
                  }
              }
            deallog << std::endl;
          }
        deallog << "value " << std::endl;
        deallog << "grad " << std::endl;
      }

  deallog.pop();
}


int
main()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(PRECISION);
  deallog << std::fixed;
  deallog.attach(logfile);

  for (unsigned int degree = 1; degree < 4; ++degree)
    {
      plot_shape_functions<2>(degree);
      plot_shape_functions<3>(degree);
    }

  return 0;
}
