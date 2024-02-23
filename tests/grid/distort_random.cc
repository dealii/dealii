// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2001 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check GridTools::distort_random


#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


template <int dim>
void
test1(const bool keep_boundary)
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(2);

  GridTools::distort_random(0.1, tria, keep_boundary);

  std::ostream &logfile = deallog.get_file_stream();
  deallog << "dim=" << dim << ", keep_boundary=" << keep_boundary << std::endl;
  GridOut().write_gnuplot(tria, logfile);
}



int
main()
{
  initlog();
  deallog.get_file_stream() << std::setprecision(4);

  test1<1>(true);
  test1<1>(false);
  test1<2>(true);
  test1<2>(false);
  test1<3>(true);
  test1<3>(false);

  return 0;
}
