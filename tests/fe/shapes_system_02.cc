// ---------------------------------------------------------------------
//
// Copyright (C) 2013 - 2020 by the deal.II authors
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
  MappingQGeneric<dim> m(1);

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
