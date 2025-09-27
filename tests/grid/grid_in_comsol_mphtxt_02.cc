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



// Check whether we can read in meshes in COMSOL's .mphtxt format.
//
// Compared to the _01 test, we here read in a whole bunch of other
// meshes but output less detailed information about them.

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"



template <int dim>
void
comsol_grid(const char *name)
{
  Triangulation<dim> tria;
  GridIn<dim>        grid_in;
  grid_in.attach_triangulation(tria);
  std::ifstream input_file(name);
  grid_in.read_comsol_mphtxt(input_file);

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

#if 0
  {
    static int    output_file = 0;
    GridOut       go;
    std::ofstream o("x-" + std::to_string(output_file) + ".gnuplot");
    go.write_gnuplot(tria, o);

    ++output_file;
  }

#endif
}


int
main()
{
  initlog();

  try
    {
      comsol_grid<3>(SOURCE_DIR "/grids/comsol/mesh_example_intro.mphtxt");
      comsol_grid<3>(
        SOURCE_DIR "/grids/comsol/mesh_example_with_complete_geom_info.mphtxt");
      comsol_grid<3>(
        SOURCE_DIR
        "/grids/comsol/mesh_example_with_domain_geom_info_1_7.mphtxt");
      comsol_grid<3>(SOURCE_DIR
                     "/grids/comsol/mesh_example_with_domain_geom_info.mphtxt");
      comsol_grid<3>(
        SOURCE_DIR
        "/grids/comsol/mesh_example_without_domain_geom_info.mphtxt");
      comsol_grid<3>(SOURCE_DIR "/grids/comsol/busbar_with_selections.mphtxt");
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
