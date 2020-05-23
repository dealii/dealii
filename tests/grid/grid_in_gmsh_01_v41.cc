// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2020 by the deal.II authors
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
gmsh_grid(const char *name_v2, const char *name_v41)
{
  Triangulation<dim> tria_v2;
  {
    GridIn<dim> grid_in;
    grid_in.attach_triangulation(tria_v2);
    std::ifstream input_file(name_v2);
    grid_in.read_msh(input_file);
  }

  Triangulation<dim> tria_v41;
  {
    GridIn<dim> grid_in;
    grid_in.attach_triangulation(tria_v41);
    std::ifstream input_file(name_v41);
    grid_in.read_msh(input_file);
  }

  // The vertex index is different between the meshes parsed from the GMSH-2 and
  // the GMSH-4 fil but all other information should match.
  AssertThrow(tria_v2.n_active_cells() == tria_v41.n_active_cells(),
              ExcInternalError());
  deallog << "  " << tria_v2.n_active_cells() << " active cells" << std::endl;

  auto       cell_v2  = tria_v2.begin_active();
  auto       cell_v41 = tria_v41.begin_active();
  const auto end_v2   = tria_v2.end();
  for (; cell_v2 != end_v2; ++cell_v2, ++cell_v41)
    {
      AssertThrow(cell_v2->material_id() == cell_v41->material_id(),
                  ExcInternalError());
      for (const unsigned int i : GeometryInfo<dim>::vertex_indices())
        {
          AssertThrow((cell_v2->vertex(i) - cell_v41->vertex(i)).norm() <
                        1.e-10,
                      ExcInternalError());
        }
      for (const unsigned int i : GeometryInfo<dim>::face_indices())
        {
          AssertThrow(cell_v2->face(i)->boundary_id() ==
                        cell_v41->face(i)->boundary_id(),
                      ExcInternalError());
        }
      for (unsigned int i = 0; i < GeometryInfo<dim>::lines_per_cell; ++i)
        {
          AssertThrow(cell_v2->line(i)->boundary_id() ==
                        cell_v41->line(i)->boundary_id(),
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
      deallog << "/grid_in_msh_01.2d.v41.msh" << std::endl;
      gmsh_grid<2>(SOURCE_DIR "/grids/grid_in_msh_01.2d.msh",
                   SOURCE_DIR "/grids/grid_in_msh_01.2d.v41.msh");
      deallog << "/grid_in_msh_01.2da.v41.msh" << std::endl;
      gmsh_grid<2>(SOURCE_DIR "/grids/grid_in_msh_01.2da.msh",
                   SOURCE_DIR "/grids/grid_in_msh_01.2da.v41.msh");
      deallog << "/grid_in_msh_01.3d.v41.msh" << std::endl;
      gmsh_grid<3>(SOURCE_DIR "/grids/grid_in_msh_01.3d.msh",
                   SOURCE_DIR "/grids/grid_in_msh_01.3d.v41.msh");
      deallog << "/grid_in_msh_01.3da.v41.msh" << std::endl;
      gmsh_grid<3>(SOURCE_DIR "/grids/grid_in_msh_01.3da.msh",
                   SOURCE_DIR "/grids/grid_in_msh_01.3da.v41.msh");
      deallog << "/grid_in_msh_01.3d_neg.v41.msh" << std::endl;
      gmsh_grid<3>(SOURCE_DIR "/grids/grid_in_msh_01.3d_neg.msh",
                   SOURCE_DIR "/grids/grid_in_msh_01.3d_neg.v41.msh");
    }
  catch (std::exception &exc)
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
