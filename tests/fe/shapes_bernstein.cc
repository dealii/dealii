// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/fe/fe_bernstein.h>
#include <deal.II/fe/mapping_q1.h>

#include <string>

#include "../tests.h"

#include "shapes.h"
// #include "../../include/fe_bernstein.h"

#define PRECISION 8


template <int dim>
void
plot_FE_Bernstein_shape_functions()
{
  MappingQ<dim> m(1);

  FE_Bernstein<dim> b1(1);
  plot_shape_functions(m, b1, "B1");
  plot_face_shape_functions(m, b1, "B1");
  test_compute_functions(m, b1, "B1");

  FE_Bernstein<dim> b2(2);
  plot_shape_functions(m, b2, "B2");
  plot_face_shape_functions(m, b2, "B2");
  test_compute_functions(m, b2, "B2");

  // skip the following tests to
  // reduce run-time
  if (dim < 3)
    {
      FE_Bernstein<dim> b3(3);
      plot_shape_functions(m, b3, "B3");
      plot_face_shape_functions(m, b3, "B3");
      test_compute_functions(m, b3, "B3");

      FE_Bernstein<dim> b4(4);
      plot_shape_functions(m, b4, "B4");
      plot_face_shape_functions(m, b4, "B4");
      test_compute_functions(m, b4, "B4");
    };
}


int
main()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(PRECISION) << std::fixed;
  deallog.attach(logfile);

  plot_FE_Bernstein_shape_functions<1>();
  plot_FE_Bernstein_shape_functions<2>();
  plot_FE_Bernstein_shape_functions<3>();

  return 0;
}
