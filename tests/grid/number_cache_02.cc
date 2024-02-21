// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2008 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// test that we correctly compute the number of lines, quads, hexes, etc and
// store them in the number cache of the Triangulation class. use the meshes
// from grid_in (2d) and grid_in_3d (3d) for this

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"


template <int dim>
void
output(const Triangulation<dim> &tria)
{
  deallog << "  " << tria.n_active_cells() << std::endl;
  deallog << "  " << tria.n_cells() << std::endl;
  deallog << "  " << tria.n_active_lines() << std::endl;
  deallog << "  " << tria.n_lines() << std::endl;
  deallog << "  " << tria.n_active_quads() << std::endl;
  deallog << "  " << tria.n_quads() << std::endl;
  deallog << "  " << tria.n_active_hexs() << std::endl;
  deallog << "  " << tria.n_hexs() << std::endl;

  for (unsigned int i = 0; i < tria.n_levels(); ++i)
    {
      deallog << "  " << tria.n_active_cells(i) << std::endl;
      deallog << "  " << tria.n_cells(i) << std::endl;
      if (dim == 1)
        {
          deallog << "  " << tria.n_active_lines(i) << std::endl;
          deallog << "  " << tria.n_lines(i) << std::endl;
        }

      if (dim == 2)
        {
          deallog << "  " << tria.n_active_quads(i) << std::endl;
          deallog << "  " << tria.n_quads(i) << std::endl;
        }

      if (dim == 3)
        {
          deallog << "  " << tria.n_active_hexs(i) << std::endl;
          deallog << "  " << tria.n_hexs(i) << std::endl;
        }
    }
}


template <int dim>
void
test(const char *filename)
{
  Triangulation<dim> tria(Triangulation<dim>::none, true);
  GridIn<dim>        gi;
  gi.attach_triangulation(tria);
  std::ifstream in(filename);

  try
    {
      gi.read_xda(in);
    }
  catch (typename Triangulation<dim>::DistortedCellList &dcv)
    {
      // ignore the exception that we
      // get because the mesh has
      // distorted cells
      deallog << dcv.distorted_cells.size() << " cells are distorted."
              << std::endl;
    }
  catch (const std::exception &exc)
    {
      deallog << "  caught exception:" << std::endl << exc.what() << std::endl;
      return;
    }

  output(tria);

  // now refine a few cells and output again
  deallog << "  Refining..." << std::endl;
  typename Triangulation<dim>::active_cell_iterator cell = tria.begin_active();
  for (unsigned int i = 0; i <= std::min(tria.n_active_cells() / 3, 10U);
       ++i, ++cell)
    cell->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  output(tria);
}


int
main()
{
  initlog();
  deallog << std::setprecision(2);

  test<2>(SOURCE_DIR "/grid_in/2d.xda");

  test<3>(SOURCE_DIR "/grid_in_3d/1.in");
  test<3>(SOURCE_DIR "/grid_in_3d/2.in");
  test<3>(SOURCE_DIR "/grid_in_3d/3.in");
  test<3>(SOURCE_DIR "/grid_in_3d/4.in");

  test<3>(SOURCE_DIR "/grid_in_3d/evil_0.in");
  test<3>(SOURCE_DIR "/grid_in_3d/evil_1.in");
  test<3>(SOURCE_DIR "/grid_in_3d/evil_2.in");
  test<3>(SOURCE_DIR "/grid_in_3d/evil_3.in");
  test<3>(SOURCE_DIR "/grid_in_3d/evil_4.in");
}
