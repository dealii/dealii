// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// a simplified version of shapes_system.cc

#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q1.h>

#include <string>

#include "../tests.h"

#include "shapes.h"

#define PRECISION 8

template <int dim>
void
plot_FE_System_shape_functions()
{
  MappingQ<dim> m(1);

  FESystem<dim> p3(FE_Q<dim>(1), 1, FESystem<dim>(FE_Q<dim>(1), 2), 2);
  test_compute_functions(m, p3, "System_1");

  FESystem<dim> p4(p3, 2, FESystem<dim>(p3, 3), 1, p3, 1);
  test_compute_functions(m, p4, "System_2");
}


int
main()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(PRECISION) << std::fixed;
  deallog.attach(logfile);
  deallog << "FE_System<1>" << std::endl;
  plot_FE_System_shape_functions<1>();
  deallog << "FE_System<2>" << std::endl;
  plot_FE_System_shape_functions<2>();
  deallog << "FE_System<3>" << std::endl;
  plot_FE_System_shape_functions<3>();

  return 0;
}
