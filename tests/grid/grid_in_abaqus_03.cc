// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2004 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



// Like the _02 test, but check that we can read files that have DOS
// line endings. To this end, we read the files from the _02 test, and
// convert them to ones that have DOS line endings. Then pass those on
// to the reader. This used to cause total chaos with error messages
// that were really not useful.

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <sstream>

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

  std::ostringstream file_with_dos_line_endings;
  while (input_file)
    {
      std::string line;
      std::getline(input_file, line);

      // Put the line back into a stream, but with a DOS line ending:
      file_with_dos_line_endings << line << '\r' << '\n';
    }
  file_with_dos_line_endings.flush();

  {
    std::ofstream x("x.abc");
    x << file_with_dos_line_endings.str();
  }

  std::ifstream input("x.abc");
  grid_in.read_abaqus(input);

  deallog << "  " << tria.n_active_cells() << " active cells" << std::endl;

  int hash  = 0;
  int index = 0;
  for (typename Triangulation<dim, spacedim>::active_cell_iterator c =
         tria.begin_active();
       c != tria.end();
       ++c, ++index)
    for (const unsigned int i : c->vertex_indices())
      hash += (index * i * c->vertex_index(i)) % (tria.n_active_cells() + 1);
  deallog << "  hash=" << hash << std::endl;

  if (output_grid)
    {
      // Write in uncompressed VTK so that its possible to verify by eye
      // that the last coordinate of all nodes is non-zero.
      std::ofstream     fout_vtk(gridout_name + ".vtk");
      GridOutFlags::Vtk flags_vtk;
      flags_vtk.compression_level =
        DataOutBase::CompressionLevel::no_compression;
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
      abaqus_grid<3, 3>(SOURCE_DIR "/grids/abaqus/3d/ellipsoid_cavity.inp",
                        false);
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
