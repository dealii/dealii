// ---------------------------------------------------------------------
//
// Copyright (C) 2021 by the deal.II authors
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


// Create a triangulation with a pyramid and output it in gmsh format

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
  GridGenerator::reference_cell(triangulation, ReferenceCells::Pyramid);
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
