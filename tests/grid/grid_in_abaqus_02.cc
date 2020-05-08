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



// Check whether we can read in with the abaqus format
// Specifically, this tests for the codimension 1 case; previously the
// third dimension was being ignored, effectively projecting the mesh to
// the XY plane.

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"



template <int dim, int spacedim>
void
abaqus_grid(const std::string path_and_name,
            const bool        output_grid  = false,
            const std::string gridout_name = "")
{
  Triangulation<dim, spacedim> tria;
  GridIn<dim, spacedim>        grid_in;
  grid_in.attach_triangulation(tria);
  std::ifstream input_file(path_and_name);
  grid_in.read_abaqus(input_file);

  deallog << "  " << tria.n_active_cells() << " active cells" << std::endl;

  int hash  = 0;
  int index = 0;
  for (typename Triangulation<dim, spacedim>::active_cell_iterator c =
         tria.begin_active();
       c != tria.end();
       ++c, ++index)
    for (const unsigned int i : GeometryInfo<dim>::vertex_indices())
      hash += (index * i * c->vertex_index(i)) % (tria.n_active_cells() + 1);
  deallog << "  hash=" << hash << std::endl;

  if (output_grid)
    {
      // Write in uncompressed VTK so that its possible to verify by eye
      // that the last coordinate of all nodes is non-zero.
      std::ofstream     fout_vtk(gridout_name + ".vtk");
      GridOutFlags::Vtk flags_vtk;
      flags_vtk.compression_level = DataOutBase::VtkFlags::no_compression;
      GridOut gridout;
      gridout.set_flags(flags_vtk);
      gridout.write_vtk(tria, fout_vtk);
    }
}


int
main()
{
  initlog();

  try
    {
      deallog << "codim_1-tube" << std::endl;
      abaqus_grid<2, 3>(SOURCE_DIR "/grids/abaqus/3d/codim_1-tube.inp",
                        false,
                        "codim_1-tube");
      deallog << "codim_1-jagged_surface" << std::endl;
      abaqus_grid<2, 3>(SOURCE_DIR
                        "/grids/abaqus/3d/codim_1-jagged_surface.inp",
                        false,
                        "codim_1-jagged_surface");
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
