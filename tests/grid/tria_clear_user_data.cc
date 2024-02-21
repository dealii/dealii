// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// We assess the possibility of using clear user data in any dim.

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

template <int dim>
void
test()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.clear_user_data();
  deallog << "Cleared user data in " << dim << 'D' << std::endl;
}


int
main()
{
  initlog();
  deallog << std::setprecision(2);

  test<1>();
  test<2>();
  test<3>();

  return 0;
}
