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



// an extract of the _01 test that shows the essence of what's going
// wrong. compared to the _03 testcase, upon refinement of the middle
// cell (0.0), cell 0.1 forgets who some of its neighbors are. this is
// clearly not good
//
// the underlying reason was that when we update neighborship
// information in the triangulation class, we have to take into
// account the direction_flag of cells

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
  boundary_mesh.begin_active()->set_refine_flag();
  boundary_mesh.execute_coarsening_and_refinement();

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
                deallog << "              "
                        << cell->neighbor_child_on_subface(face, c)
                        << std::endl;
              }
        }
    }

  return 0;
}
