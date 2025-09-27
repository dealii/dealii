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



// like _05 but do the test with a complex mesh

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
    std::ifstream in(SOURCE_DIR "/../grid/grid_in_3d_02/747.ucd");
    gi.read(in);
  }

  {
    GridIn<dim> gi;
    gi.attach_triangulation(tr2);
    std::ifstream in(SOURCE_DIR "/../grid/grid_in_3d_02/747.ucd");
    gi.read(in);
  }

  while (tr.n_active_cells() < 70000)
    {
      std::vector<bool> flags(tr.n_active_cells(), false);

      // refine one 1/50 of all cells each time (but at least one)
      deallog << "Refining cells ... " << std::endl;
      for (unsigned int i = 0; i < tr.n_active_cells() / 50 + 1; ++i)
        {
          const unsigned int x = Testing::rand() % flags.size();
          flags[x]             = true;
        }

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
      tr.execute_coarsening_and_refinement();
      tr2.execute_coarsening_and_refinement();

      deallog << " Number of cells: " << tr.n_active_cells() << ' '
              << tr2.n_active_cells() << std::endl;

      assert_tria_equal(tr, tr2);
    }
}


int
main(int argc, char *argv[])
{
  initlog();
#ifdef DEAL_II_WITH_MPI
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
#else
  (void)argc;
  (void)argv;
#endif

  test<3>(deallog.get_file_stream());
}
