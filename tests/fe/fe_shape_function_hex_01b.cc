// ---------------------------------------------------------------------
//
// Copyright (C) 2022 by the deal.II authors
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


// Plot shape functions for scalar and vector shape functions
// - Hexahedra
// - FE_Q_Hierarchical

#include <deal.II/fe/fe_q_hierarchical.h>

#include "../tests.h"

#include "common/fe_shape_function_plotter.h"


template <int dim>
void
plot_all(const unsigned int n_divs = 4)
{
  std::cout << "Plotting in dim " << dim << "." << std::endl;

  Assert(dim == 2 || dim == 3, ExcDimensionMismatch2(dim, 2, 3));
  const ReferenceCell &reference_cell = ReferenceCells::get_hypercube<dim>();

  plot_one<FE_Q_Hierarchical<dim>>(
    "FE_Q_Hierarchical", 1, 2, reference_cell, n_divs);
}


int
main()
{
  initlog();

  plot_all<3>();

  deallog << "OK" << std::endl;

  return 0;
}
