// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Read in a ugrid mesh and save the results to an output file. The
// mesh corresponds to the three-piece high-lift configuration of a
// wing, with a very large area around it.

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <string>

#include "../tests.h"


template <int dim, int spacedim = dim>
void
read_and_print(const std::string &filename)
{
  Triangulation<dim, spacedim> tria;
  GridIn<dim, spacedim>        gi;
  gi.attach_triangulation(tria);
  std::ifstream in(filename);
  Assert(in, ExcIO());
  gi.read_ugrid(in);

  deallog << "Number of vertices: " << tria.get_vertices().size() << std::endl;
  deallog << "Number of cells: " << tria.n_cells() << std::endl;

  // Output diagnostic information. Do so for all cells if the mesh is
  // small, or every 50th cell if the mesh is large:
  for (const auto &cell : tria.active_cell_iterators())
    if ((tria.n_active_cells() < 20) || (cell->active_cell_index() % 50 == 0))
      deallog << "cell " << cell->index()
              << " type = " << cell->reference_cell().to_string()
              << " volume = " << cell->measure() << std::endl;

  // Same with faces, for every 10th cell:
  for (const auto &cell : tria.active_cell_iterators())
    if ((tria.n_active_cells() < 20) || (cell->active_cell_index() % 10 == 0))
      for (const auto &face : cell->face_iterators())
        if (face->at_boundary())
          deallog << "boundary face " << face
                  << " type = " << face->reference_cell().to_string()
                  << " at position=" << face->center()
                  << " boundary_id = " << face->boundary_id() << std::endl;

  GridOut go;
  go.write_vtk(tria, deallog.get_file_stream());
}

int
main()
{
  initlog();
  deallog.get_file_stream() << std::setprecision(2);

  // Read the whole mesh file describing the 3-part high-lift wing:
  read_and_print<2, 2>(SOURCE_DIR "/grids/ugrid_01.txt");
  read_and_print<2, 3>(SOURCE_DIR "/grids/ugrid_01.txt");

  // Also read a very much reduced mesh with just two cells:
  read_and_print<2, 2>(SOURCE_DIR "/grids/ugrid_02.txt");
  read_and_print<2, 3>(SOURCE_DIR "/grids/ugrid_02.txt");

  deallog << "OK" << std::endl;
}
