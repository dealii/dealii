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


#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/mapping_q1.h>

#include <string>

#include "../tests.h"

#include "shapes.h"

#define PRECISION 8


template <int dim>
void
plot_FE_Nedelec_shape_functions()
{
  MappingQ<dim>   m(1);
  FE_Nedelec<dim> p0(0);
  //   plot_shape_functions(m, p1, "Nedelec1");
  //   plot_face_shape_functions(m, p1, "Nedelec1");
  test_compute_functions(m, p0, "Nedelec0");
  FE_Nedelec<dim> p1(1);
  test_compute_functions(m, p1, "Nedelec1");
}


int
main()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(PRECISION) << std::fixed;
  deallog.attach(logfile);
  deallog << "FE_Nedelec<2>" << std::endl;
  plot_FE_Nedelec_shape_functions<2>();
  deallog << "FE_Nedelec<3>" << std::endl;
  plot_FE_Nedelec_shape_functions<3>();

  return 0;
}
