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
// - Pyramid
// - FE_PyramidDGP

#include <deal.II/fe/fe_pyramid_p.h>

#include "../tests.h"

#include "common/fe_shape_function_plotter.h"


template <int dim>
void
plot_all(const unsigned int n_divs = 4)
{
  std::cout << "Plotting in dim " << dim << "." << std::endl;

  Assert(dim == 3, ExcDimensionMismatch(dim, 3));
  const ReferenceCell &reference_cell = ReferenceCells::Pyramid;

  // Not yet implemented for higher order elements
  const unsigned int max_order = 1;
  plot_one<FE_PyramidDGP<dim>>(
    "FE_PyramidDGP", 1, max_order, reference_cell, n_divs);
}


int
main()
{
  initlog();

  // Not yet implemented for isotropic refinement
  const unsigned int n_divs = 0;
  plot_all<3>(n_divs);

  deallog << "OK" << std::endl;

  return 0;
}
