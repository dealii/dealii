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

#include "derivatives.h"



/**
 * Test case to check the derivative values  up to @p regularity
 * of <code>FE_Hermite<dim>(2 \times regularity + 1)<\code> at the vertices
 * of a hypercube reference element are calculated correctly.
 * This checks that the basis is assembled correctly and will not
 * accidentally enforce discontinuities rather than continuity
 * in derivatives across element boundaries.
 */


template <int dim>
void
print_hermite_endpoint_derivatives()
{
  MappingCartesian<dim> m;

  FE_Hermite<dim> herm_1(1);
  plot_function_derivatives<dim>(m, herm_1, "Hermite-1");

  FE_Hermite<dim> herm_3(3);
  plot_function_derivatives<dim>(m, herm_3, "Hermite-3");

  // Skip the following for dim 3 or greater
  if (dim < 3)
    {
      FE_Hermite<dim> herm_5(5);
      plot_function_derivatives<dim>(m, herm_5, "Hermite-5");
    }
  if (dim == 1)
    {
      FE_Hermite<dim> herm_7(7);
      plot_function_derivatives<dim>(m, herm_7, "Hermite-7");

      FE_Hermite<dim> herm_9(9);
      plot_function_derivatives<dim>(m, herm_9, "Hermite-9");
    }
}



int
main()
{
  std::ofstream logfile("output");

  deallog << std::setprecision(PRECISION) << std::fixed;
  deallog.attach(logfile);

  print_hermite_endpoint_derivatives<1>();
  print_hermite_endpoint_derivatives<2>();
  print_hermite_endpoint_derivatives<3>();

  return 0;
}
