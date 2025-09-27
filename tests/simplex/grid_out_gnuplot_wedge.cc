// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Create a triangulation with a wedge and output it in gmsh format

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <fstream>

#include "../tests.h"

template <int dim>
void
check()
{
  Triangulation<dim> triangulation;
  GridGenerator::reference_cell(triangulation, ReferenceCells::Wedge);
  {
    GridOut go;
    go.set_flags(GridOutFlags::Msh(true));

    go.write_gnuplot(triangulation, deallog.get_file_stream());
  }

  deallog << "OK!" << std::endl;
}


int
main()
{
  initlog();
  check<3>();
}
