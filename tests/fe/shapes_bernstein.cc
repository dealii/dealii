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


#include <deal.II/fe/fe_bernstein.h>
#include <deal.II/fe/mapping_q1.h>

#include <string>

#include "../tests.h"

#include "shapes.h"
//#include "../../include/fe_bernstein.h"

#define PRECISION 8


template <int dim>
void
plot_FE_Bernstein_shape_functions()
{
  MappingQGeneric<dim> m(1);

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
