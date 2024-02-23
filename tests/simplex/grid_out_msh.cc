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


// Create a triangulation with simplices and output it in gmsh format

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

  {
    Triangulation<dim> x;
    GridGenerator::subdivided_cylinder(x, 4, 1, 2);

    // For the cylinder geometry
    for (auto cell : x.active_cell_iterators())
      for (const unsigned int face : cell->face_indices())
        if (cell->at_boundary(face))
          {
            if (std::fabs(cell->face(face)->center()[0] - 0) <
                0.1 * cell->diameter())
              cell->face(face)->set_boundary_id(1);
            else if (std::fabs(cell->face(face)->center()[0] - 4) <
                     0.1 * cell->diameter())
              cell->face(face)->set_boundary_id(2);
          }
    GridGenerator::convert_hypercube_to_simplex_mesh(x, triangulation);
  }

  {
    GridOut go;
    go.set_flags(GridOutFlags::Msh(true));

    go.write_msh(triangulation, deallog.get_file_stream());
  }

  deallog << "OK!" << std::endl;
}


int
main()
{
  initlog();
  check<3>();
}
