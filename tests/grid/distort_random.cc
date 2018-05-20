// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

// check GridTools::distort_random

#include "../tests.h"
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

template <int dim>
void
test1(const bool keep_boundary)
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(2);

  GridTools::distort_random(0.1, tria, keep_boundary);

  std::ostream& logfile = deallog.get_file_stream();
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
