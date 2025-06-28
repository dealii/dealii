/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2020 - 2025 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Part of the source code is dual licensed under Apache-2.0 WITH
 * LLVM-exception OR LGPL-2.1-or-later. Detailed license information
 * governing the source code and code contributions can be found in
 * LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
 *
 * ------------------------------------------------------------------------
 *
 * Test function GridGenerator::alfeld_split_of_simplex_mesh() in 2D and 3D.
 */

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <iostream>

#include "../tests.h"


template <int dim, int spacedim>
void
create_triangulation(Triangulation<dim, spacedim> &triangulation)
{
  GridGenerator::subdivided_hyper_cube(triangulation, 4);
}

template <int dim>
void
create_triangulation(Triangulation<dim, dim> &triangulation)
{
  GridGenerator::quarter_hyper_ball(triangulation);
}

template <int dim, int spacedim>
void
check()
{
  Triangulation<dim, spacedim> in_tria, out_tria;
  GridGenerator::subdivided_hyper_cube_with_simplices(
    in_tria, 2, 0.0, 1.0, true);

  // make each cell a different material id
  unsigned int m_id = 0;
  for (const auto &cell : in_tria)
    {
      cell.set_material_id(m_id++);
    }

  GridGenerator::alfeld_split_of_simplex_mesh(in_tria, out_tria);

  // write 2 outputs (total mesh and only surface mesh)
  const auto grid_out = [](const auto &tria, const bool surface_mesh_only) {
    GridOutFlags::Vtk flags;

    if (surface_mesh_only)
      {
        flags.output_cells         = false;
        flags.output_faces         = true;
        flags.output_edges         = false;
        flags.output_only_relevant = false;
      }

    GridOut grid_out;
    grid_out.set_flags(flags);

    grid_out.write_vtk(tria, deallog.get_file_stream());
  };

  grid_out(out_tria, false); // total mesh
  grid_out(out_tria, true);  // only surface mesh

  deallog << "OK!" << std::endl;
}


int
main()
{
  initlog();

  deallog.push("2d");
  check<2, 2>();
  deallog.pop();

  deallog.push("3d");
  check<3, 3>();
  deallog.pop();
}
