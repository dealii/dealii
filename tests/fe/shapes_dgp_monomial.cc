// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/fe/fe_dgp_monomial.h>
#include <deal.II/fe/mapping_q1.h>

#include <string>

#include "../tests.h"

#include "shapes.h"

#define PRECISION 8


template <int dim>
void
plot_FE_DGPMonomial_shape_functions()
{
  MappingQ<dim> m(1);

  FE_DGPMonomial<dim> p1(1);
  plot_shape_functions(m, p1, "DGPMonomial1");
  plot_face_shape_functions(m, p1, "DGPMonomial1");
  test_compute_functions(m, p1, "DGPMonomial1");

  FE_DGPMonomial<dim> p2(2);
  plot_shape_functions(m, p2, "DGPMonomial2");
  plot_face_shape_functions(m, p2, "DGPMonomial2");
  test_compute_functions(m, p2, "DGPMonomial2");

  if (dim < 3)
    {
      FE_DGPMonomial<dim> p3(3);
      plot_shape_functions(m, p3, "DGPMonomial3");
      plot_face_shape_functions(m, p3, "DGPMonomial3");
      test_compute_functions(m, p3, "DGPMonomial3");
    }
}


int
main()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(PRECISION) << std::fixed;
  deallog.attach(logfile);

  plot_FE_DGPMonomial_shape_functions<1>();
  plot_FE_DGPMonomial_shape_functions<2>();
  plot_FE_DGPMonomial_shape_functions<3>();

  return 0;
}
