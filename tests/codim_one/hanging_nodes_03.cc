// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2011 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// a further extract of the _02 test. the results here are correct and
// show the relationship between the various cells of the mesh

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


int
main()
{
  initlog();

  const unsigned int spacedim = 3;
  const unsigned int dim      = spacedim - 1;

  Triangulation<dim, spacedim> boundary_mesh;

  /*****************************************************************/
  // Create Surface Mesh:  Boundary of hypercube without one face
  {
    Triangulation<spacedim> volume_mesh;
    GridGenerator::hyper_cube(volume_mesh);
    Triangulation<spacedim>::active_cell_iterator cell =
      volume_mesh.begin_active();

    cell->face(0)->set_all_boundary_ids(1);
    const std::set<types::boundary_id> boundary_ids = {0};
    GridGenerator::extract_boundary_mesh(volume_mesh,
                                         boundary_mesh,
                                         boundary_ids);
  }

  Triangulation<dim, spacedim>::active_cell_iterator cell =
    boundary_mesh.begin_active();
  for (; cell != boundary_mesh.end(); ++cell)
    {
      deallog << "Cell = " << cell << std::endl;
      deallog << "  direction_flag = " << cell->direction_flag() << std::endl;

      for (const unsigned int face : GeometryInfo<dim>::face_indices())
        {
          deallog << "  face = " << face
                  << "  (neighbor = " << cell->neighbor(face) << ')'
                  << std::endl;

          if (cell->face(face)->has_children())
            for (unsigned int c = 0; c < cell->face(face)->n_children(); ++c)
              {
                deallog << "    subface = " << c << std::endl;
              }
        }
    }

  return 0;
}
