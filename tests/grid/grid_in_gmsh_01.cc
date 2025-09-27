// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2004 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check whether we can read in with the gmsh format

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"



template <int dim>
void
gmsh_grid(const char *name)
{
  Triangulation<dim> tria;
  GridIn<dim>        grid_in;
  grid_in.attach_triangulation(tria);
  std::ifstream input_file(name);
  grid_in.read_msh(input_file);

  deallog << "  " << tria.n_active_cells() << " active cells" << std::endl;

  int hash  = 0;
  int index = 0;
  for (typename Triangulation<dim>::active_cell_iterator c =
         tria.begin_active();
       c != tria.end();
       ++c, ++index)
    for (const unsigned int i : c->vertex_indices())
      hash += (index * i * c->vertex_index(i)) % (tria.n_active_cells() + 1);
  deallog << "  hash=" << hash << std::endl;
}


int
main()
{
  initlog();

  try
    {
      gmsh_grid<2>(SOURCE_DIR "/grids/grid_in_msh_01.2d.msh");
      gmsh_grid<2>(SOURCE_DIR "/grids/grid_in_msh_01.2da.msh");
      gmsh_grid<3>(SOURCE_DIR "/grids/grid_in_msh_01.3d.msh");
      gmsh_grid<3>(SOURCE_DIR "/grids/grid_in_msh_01.3da.msh");
      gmsh_grid<3>(SOURCE_DIR "/grids/grid_in_msh_01.3d_neg.msh");
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
