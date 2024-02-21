// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test GridGenerator::extrude

#include <deal.II/base/tensor.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"



template <int dim>
void
test(std::ostream &out)
{
  Triangulation<2> triangulation;
  Triangulation<3> tr;
  GridGenerator::hyper_rectangle(triangulation,
                                 Point<2>(0, 0),
                                 Point<2>(1, 1),
                                 true);

  GridGenerator::extrude_triangulation(triangulation, 3, 2.0, tr);
  GridOut go;
  go.set_flags(GridOutFlags::Ucd(true));
  go.write_ucd(tr, out);
}


int
main()
{
  initlog();

  test<3>(deallog.get_file_stream());
}
