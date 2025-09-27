// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/bounding_box.h>

#include "../tests.h"

/*
 * Lock all possible coordinates in dim + 1. For each locked coordinate, print
 * how the coordinates in dim dimensions map to the coordinates in dim + 1
 * dimensions.
 */
template <int dim>
void
print_what_each_coordinate_maps_to()
{
  deallog << "dim = " << dim << std::endl;

  for (unsigned int locked = 0; locked < dim + 1; ++locked)
    {
      deallog << "locked coordinate = " << locked << std::endl;
      for (unsigned int i = 0; i < dim; ++i)
        {
          deallog << i << " -> "
                  << internal::coordinate_to_one_dim_higher<dim>(locked, i)
                  << std::endl;
        }
    }
}



int
main()
{
  initlog();

  print_what_each_coordinate_maps_to<1>();
  deallog << std::endl;
  print_what_each_coordinate_maps_to<2>();

  return 0;
}
