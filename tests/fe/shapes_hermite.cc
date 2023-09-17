// ---------------------------------------------------------------------
//
// Copyright (C) 2023 - 2023 by the deal.II authors
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

#include "../tests.h"

#define PRECISION 8


#include <deal.II/fe/fe_hermite.h>
#include <deal.II/fe/mapping_cartesian.h>

#include <string>

#include "shapes.h"



/*
 * Test case on various hypercube reference cells to check that the shape
 * values of <code>FE_Hermite<dim>(reg)<\code> are correct. Values of
 * derivatives at the vertices are checked in the test case
 * derivatives_hermite.cc.
 */

template <int dim>
void
plot_FE_Hermite_shape_functions()
{
  MappingCartesian<dim> m;

  FE_Hermite<dim> herm_1(1);
  plot_shape_functions(m, herm_1, "Hermite-1");
  plot_face_shape_functions(m, herm_1, "Hermite-1");
  test_compute_functions(m, herm_1, "Hermite-1");

  FE_Hermite<dim> herm_3(3);
  plot_shape_functions(m, herm_3, "Hermite-3");
  plot_face_shape_functions(m, herm_3, "Hermite-3");
  test_compute_functions(m, herm_3, "Hermite-3");

  // skip the following tests to
  // reduce run-time
  if (dim < 3)
    {
      FE_Hermite<dim> herm_5(5);
      plot_shape_functions(m, herm_5, "Hermite-5");
      plot_face_shape_functions(m, herm_5, "Hermite-5");
      test_compute_functions(m, herm_5, "Hermite-5");
    }

  if (dim == 1)
    {
      FE_Hermite<dim> herm_7(7);
      plot_shape_functions(m, herm_7, "Hermite-7");
      plot_face_shape_functions(m, herm_7, "Hermite-7");
      test_compute_functions(m, herm_7, "Hermite-7");

      FE_Hermite<dim> herm_9(9);
      plot_shape_functions(m, herm_9, "Hermite-9");
      plot_face_shape_functions(m, herm_9, "Hermite-9");
      test_compute_functions(m, herm_9, "Hermite-9");
    };
}



int
main()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(PRECISION) << std::fixed;
  deallog.attach(logfile);

  plot_FE_Hermite_shape_functions<1>();
  plot_FE_Hermite_shape_functions<2>();
  plot_FE_Hermite_shape_functions<3>();

  return 0;
}
