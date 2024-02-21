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


#include <deal.II/fe/fe_dgp_nonparametric.h>
#include <deal.II/fe/mapping_q1.h>

#include <string>

#include "../tests.h"

#include "shapes.h"

#define PRECISION 8


template <int dim>
void
plot_FE_DGPNonparametric_shape_functions()
{
  MappingQ<dim> m(1);

  FE_DGPNonparametric<dim> p0(0);
  plot_shape_functions(m, p0, "DGPNonparametric0");
  plot_face_shape_functions(m, p0, "DGPNonparametric0");

  FE_DGPNonparametric<dim> p1(1);
  plot_shape_functions(m, p1, "DGPNonparametric1");
  plot_face_shape_functions(m, p1, "DGPNonparametric1");

  FE_DGPNonparametric<dim> p2(2);
  plot_shape_functions(m, p2, "DGPNonparametric2");
  plot_face_shape_functions(m, p2, "DGPNonparametric2");

  //    FE_DGPNonparametric<dim> p3(3);
  //    plot_shape_functions(m, p3, "DGPNonparametric3");
  //    plot_face_shape_functions(m, p3, "DGPNonparametric3");

  //    FE_DGPNonparametric<dim> p4(4);
  //    plot_shape_functions(m, p4, "DGPNonparametric4");
  //    plot_face_shape_functions(m, p4, "DGPNonparametric4");
}


int
main()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(PRECISION) << std::fixed;
  deallog.attach(logfile);

  plot_FE_DGPNonparametric_shape_functions<1>();
  plot_FE_DGPNonparametric_shape_functions<2>();
  plot_FE_DGPNonparametric_shape_functions<3>();

  return 0;
}
