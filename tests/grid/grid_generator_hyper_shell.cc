// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test if colorizing a thin quarter_hyper_shell yields correct
// boundary indicators. The function used to not colorize
// boundaries correctly if the hyper shell was too thin.

#include <deal.II/base/tensor.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"



template <int dim>
void
test(const Point<dim> origin = Point<dim>())
{
  SphericalManifold<dim> manifold(origin);
  {
    deallog << "quarter_hyper_shell with origin: " << origin << std::endl;
    Triangulation<dim> tr;
    GridGenerator::quarter_hyper_shell(tr, origin, 1460800, 1560800, 0, true);
    tr.set_manifold(0, manifold);

    bool                     found_top_boundary = false;
    const types::boundary_id top_id             = 1;
    deallog << std::endl << "Boundary ids of boundary faces:";

    for (const auto &cell : tr.active_cell_iterators())
      for (const unsigned int face_no : cell->face_indices())
        {
          const typename Triangulation<dim>::face_iterator face =
            cell->face(face_no);

          if (face->at_boundary())
            {
              const types::boundary_id boundary_id = face->boundary_id();
              deallog << " " << boundary_id;

              if (boundary_id == top_id)
                {
                  found_top_boundary = true;
                }
            }
        }

    if (dim > 2)
      {
        deallog << std::endl << std::endl << "Boundary ids of boundary lines:";

        for (const auto &cell : tr.active_cell_iterators())
          for (const unsigned int face_no : cell->face_indices())
            {
              const typename Triangulation<dim>::face_iterator face =
                cell->face(face_no);

              if (face->at_boundary())
                {
                  for (const unsigned l : face->line_indices())
                    deallog << " " << face->line(l)->boundary_id();
                }
            }
      }

    deallog << std::endl;

    if (found_top_boundary)
      deallog << "Ok. Found top boundary." << std::endl;
    else
      deallog << "Error. Failed to find top boundary." << std::endl;
  }
}


int
main()
{
  initlog();

  deallog.get_file_stream() << std::setprecision(9);

  deallog.push("2d");
  test<2>();
  test<2>(Point<2>(1.0, 0.0));
  deallog.pop();

  deallog.push("3d");
  test<3>();
  test<3>(Point<3>(1.0, 0.0, 0.0));

  deallog.pop();
}
