// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check whether we can read in with the gmsh format and obtain the same results
// for the GMSH-2 format and the GMSH-4 format

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"



template <int dim>
void
gmsh_grid(const char *name_v2, const char *name_v4)
{
  Triangulation<dim> tria_v2;
  {
    GridIn<dim> grid_in;
    grid_in.attach_triangulation(tria_v2);
    std::ifstream input_file(name_v2);
    grid_in.read_msh(input_file);
  }

  Triangulation<dim> tria_v4;
  {
    GridIn<dim> grid_in;
    grid_in.attach_triangulation(tria_v4);
    std::ifstream input_file(name_v4);
    grid_in.read_msh(input_file);
  }

  // The vertex index is different between the meshes parsed from the GMSH-2 and
  // the GMSH-4 fil but all other information should match.
  AssertThrow(tria_v2.n_active_cells() == tria_v4.n_active_cells(),
              ExcInternalError());
  deallog << "  " << tria_v2.n_active_cells() << " active cells" << std::endl;

  auto       cell_v2 = tria_v2.begin_active();
  auto       cell_v4 = tria_v4.begin_active();
  const auto end_v2  = tria_v2.end();
  for (; cell_v2 != end_v2; ++cell_v2, ++cell_v4)
    {
      AssertThrow(cell_v2->material_id() == cell_v4->material_id(),
                  ExcInternalError());
      for (const unsigned int i : cell_v2->vertex_indices())
        {
          AssertThrow((cell_v2->vertex(i) - cell_v4->vertex(i)).norm() < 1.e-10,
                      ExcInternalError());
        }
      for (const unsigned int i : GeometryInfo<dim>::face_indices())
        {
          AssertThrow(cell_v2->face(i)->boundary_id() ==
                        cell_v4->face(i)->boundary_id(),
                      ExcInternalError());
        }
      for (unsigned int i = 0; i < GeometryInfo<dim>::lines_per_cell; ++i)
        {
          AssertThrow(cell_v2->line(i)->boundary_id() ==
                        cell_v4->line(i)->boundary_id(),
                      ExcInternalError());
        }
    }
  deallog << "  OK" << std::endl;
}


int
main()
{
  initlog();

  try
    {
      deallog << "/grid_in_msh_01.2d.v4.msh" << std::endl;
      gmsh_grid<2>(SOURCE_DIR "/grids/grid_in_msh_01.2d.msh",
                   SOURCE_DIR "/grids/grid_in_msh_01.2d.v4.msh");
      deallog << "/grid_in_msh_01.2da.v4.msh" << std::endl;
      gmsh_grid<2>(SOURCE_DIR "/grids/grid_in_msh_01.2da.msh",
                   SOURCE_DIR "/grids/grid_in_msh_01.2da.v4.msh");
      deallog << "/grid_in_msh_01.3d.v4.msh" << std::endl;
      gmsh_grid<3>(SOURCE_DIR "/grids/grid_in_msh_01.3d.msh",
                   SOURCE_DIR "/grids/grid_in_msh_01.3d.v4.msh");
      deallog << "/grid_in_msh_01.3da.v4.msh" << std::endl;
      gmsh_grid<3>(SOURCE_DIR "/grids/grid_in_msh_01.3da.msh",
                   SOURCE_DIR "/grids/grid_in_msh_01.3da.v4.msh");
      deallog << "/grid_in_msh_01.3d_neg.v4.msh" << std::endl;
      gmsh_grid<3>(SOURCE_DIR "/grids/grid_in_msh_01.3d_neg.msh",
                   SOURCE_DIR "/grids/grid_in_msh_01.3d_neg.v4.msh");
    }
  catch (const std::exception &exc)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    }
  catch (...)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    };

  return 0;
}
