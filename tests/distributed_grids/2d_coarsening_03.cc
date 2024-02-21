// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Like coarsening_02, but with a complex grid

#include <deal.II/base/tensor.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/intergrid_map.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"

#include "coarse_grid_common.h"



template <int dim>
void
test(std::ostream & /*out*/)
{
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);
  Triangulation<dim>                        tr2(
    Triangulation<dim>::limit_level_difference_at_vertices);

  {
    GridIn<dim> gi;
    gi.attach_triangulation(tr);
    std::ifstream in(SOURCE_DIR "/../grid/grid_in_02/2d.xda");
    try
      {
        gi.read_xda(in);
      }
    catch (
      const typename Triangulation<dim>::DistortedCellList &distorted_cells)
      {
        // ignore distorted cells
        deallog << distorted_cells.distorted_cells.size()
                << " distorted cells after creating mesh." << std::endl;
      }
  }

  {
    GridIn<dim> gi;
    gi.attach_triangulation(tr2);
    std::ifstream in(SOURCE_DIR "/../grid/grid_in_02/2d.xda");
    try
      {
        gi.read_xda(in);
      }
    catch (
      const typename Triangulation<dim>::DistortedCellList &distorted_cells)
      {
        // ignore distorted cells
        deallog << distorted_cells.distorted_cells.size()
                << " distorted cells after creating mesh." << std::endl;
      }
  }

  deallog << tr.n_active_cells() << ' ' << tr2.n_active_cells() << std::endl;
  Assert(tr.n_active_cells() == tr2.n_active_cells(), ExcInternalError());

  try
    {
      tr.refine_global(1);
    }
  catch (const typename Triangulation<dim>::DistortedCellList &distorted_cells)
    {
      // ignore distorted cells
      deallog << distorted_cells.distorted_cells.size()
              << " distorted cells after refining mesh." << std::endl;
    }

  try
    {
      tr2.refine_global(1);
    }
  catch (const typename Triangulation<dim>::DistortedCellList &distorted_cells)
    {
      // ignore distorted cells
      deallog << distorted_cells.distorted_cells.size()
              << " distorted cells after refining mesh." << std::endl;
    }

  deallog << tr.n_active_cells() << ' ' << tr2.n_active_cells() << std::endl;
  assert_tria_equal(tr, tr2);

  for (unsigned int i = 0; i < 1; ++i)
    {
      // refine one-fifth of cells randomly
      std::vector<bool> flags(tr.n_active_cells(), false);
      for (unsigned int k = 0; k < flags.size() / 5 + 1; ++k)
        flags[Testing::rand() % flags.size()] = true;
      // make sure there's at least one that
      // will be refined
      flags[0] = true;

      InterGridMap<Triangulation<dim>> intergrid_map;
      intergrid_map.make_mapping(tr, tr2);

      // refine tr and tr2
      unsigned int index = 0;
      for (typename Triangulation<dim>::active_cell_iterator cell =
             tr.begin_active();
           cell != tr.end();
           ++cell, ++index)
        if (flags[index])
          {
            cell->set_refine_flag();
            intergrid_map[cell]->set_refine_flag();
          }
      Assert(index == tr.n_active_cells(), ExcInternalError());

      // flag all other cells for coarsening
      // (this should ensure that at least
      // some of them will actually be
      // coarsened)
      index = 0;
      for (typename Triangulation<dim>::active_cell_iterator cell =
             tr.begin_active();
           cell != tr.end();
           ++cell, ++index)
        if (!flags[index])
          {
            cell->set_coarsen_flag();
            intergrid_map[cell]->set_coarsen_flag();
          }

      try
        {
          tr.execute_coarsening_and_refinement();
        }
      catch (
        const typename Triangulation<dim>::DistortedCellList &distorted_cells)
        {
          // ignore distorted cells
          deallog << distorted_cells.distorted_cells.size()
                  << " distorted cells after adaptively refining mesh."
                  << std::endl;
        }
      try
        {
          tr2.execute_coarsening_and_refinement();
        }
      catch (
        const typename Triangulation<dim>::DistortedCellList &distorted_cells)
        {
          // ignore distorted cells
          deallog << distorted_cells.distorted_cells.size()
                  << " distorted cells after adaptively refining mesh."
                  << std::endl;
        }

      deallog << std::endl;

      deallog << i << " Number of cells: " << tr.n_active_cells() << ' '
              << tr2.n_active_cells() << std::endl;
      assert_tria_equal(tr, tr2);
    }
}


int
main(int argc, char *argv[])
{
  initlog();
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  deallog.push("2d");
  test<2>(deallog.get_file_stream());
  deallog.pop();
}
