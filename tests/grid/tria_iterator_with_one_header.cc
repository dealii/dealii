// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// It used to be that if you wanted to use
// Triangulation::active_cell_iterator in a context where that type
// was actually used, you also had to #include
// <deal.II/grid/tria_iterator.h> and <deal.II/grid/tria_accessor.h>.
// This was changed in r25531 in such a way that these types are now
// no longer only forward declared. test that this continues to work

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


template <int dim>
void
test()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  deallog << tria.begin_active()->center() << std::endl;
}


int
main()
{
  initlog();
  deallog << std::setprecision(4);

  test<1>();
  test<2>();
  test<3>();

  return 0;
}
