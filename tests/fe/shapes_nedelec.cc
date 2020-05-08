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
  MappingQGeneric<dim> m(1);
  FE_Nedelec<dim>      p0(0);
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
