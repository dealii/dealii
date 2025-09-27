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


#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/mapping_q1.h>

#include <string>

#include "../tests.h"

#include "shapes.h"

#define PRECISION 8


template <int dim>
void
plot_FE_DGP_shape_functions()
{
  MappingQ<dim> m(1);

  FE_DGP<dim> p1(1);
  plot_shape_functions(m, p1, "DGP1");
  plot_face_shape_functions(m, p1, "DGP1");
  test_compute_functions(m, p1, "DGP1");

  FE_DGP<dim> p2(2);
  plot_shape_functions(m, p2, "DGP2");
  plot_face_shape_functions(m, p2, "DGP2");
  test_compute_functions(m, p2, "DGP2");

  FE_DGP<dim> p3(3);
  plot_shape_functions(m, p3, "DGP3");
  plot_face_shape_functions(m, p3, "DGP3");
  test_compute_functions(m, p3, "DGP3");

  //    FE_DGP<dim> p4(4);
  //    plot_shape_functions(m, p4, "DGP4");
  //    plot_face_shape_functions(m, p4, "DGP4");
  //    test_compute_functions(m, p4, "DGP4");

  //    FE_DGP<dim> p5(5);
  //    plot_shape_functions(m, p5, "DGP5");
  //    FE_DGP<dim> p6(6);
  //    plot_shape_functions(m, p6, "DGP6");
  //    FE_DGP<dim> p7(7);
  //    plot_shape_functions(m, p7, "DGP7");
  //    FE_DGP<dim> p8(8);
  //    plot_shape_functions(m, p8, "DGP8");
  //    FE_DGP<dim> p9(9);
  //    plot_shape_functions(m, p9, "DGP9");
  //    FE_DGP<dim> p10(10);
  //    plot_shape_functions(m, p10, "DGP10");
}


int
main()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(PRECISION) << std::fixed;
  deallog.attach(logfile);

  plot_FE_DGP_shape_functions<1>();
  plot_FE_DGP_shape_functions<2>();
  plot_FE_DGP_shape_functions<3>();

  return 0;
}
