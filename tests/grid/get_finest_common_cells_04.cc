// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/mpi.h>

#include <deal.II/distributed/shared_tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>

#include "../tests.h"



template <int dim>
void
test()
{
  // create 2 triangulations with the
  // same coarse grid, and refine
  // them differently
  parallel::shared::Triangulation<dim> tria_0(MPI_COMM_WORLD);
  parallel::shared::Triangulation<dim> tria_1(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria_0);
  GridGenerator::hyper_cube(tria_1);

  tria_0.refine_global(2);
  tria_1.refine_global(2);

  tria_0.begin_active()->set_refine_flag();
  tria_0.execute_coarsening_and_refinement();

  tria_1.last_active()->set_refine_flag();
  tria_1.execute_coarsening_and_refinement();

  tria_1.last_active()->set_refine_flag();
  tria_1.execute_coarsening_and_refinement();

#if 0 // turn this block back on to get plots of grids and locally owned cells
  if (dim == 2)
    {
      std::string tria_0_out =
        "tria-0-" + std::to_string(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)) + ".svg";
      std::ofstream tria_0_out_stream(tria_0_out);
      std::string tria_1_out =
        "tria-1-" + std::to_string(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)) + ".svg";
      std::ofstream tria_1_out_stream(tria_1_out);
      GridOut go;
      go.set_flags (GridOutFlags::Svg(2, 4, true, GridOutFlags::Svg::white, 0, 0, GridOutFlags::Svg::level_number, false, true, true));
      go.write_svg(tria_0, tria_0_out_stream);
      go.write_svg(tria_1, tria_1_out_stream);
    }

  deallog << "tria 0 cells:" << std::endl;
  for (const auto &cell : tria_0.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        deallog << cell << std::endl;
    }

  deallog << "tria 1 cells:" << std::endl;
  for (const auto &cell : tria_1.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        deallog << cell << std::endl;
    }
#endif

  using CellList =
    std::list<std::pair<typename Triangulation<dim>::cell_iterator,
                        typename Triangulation<dim>::cell_iterator>>;

  deallog << "number of locally owned cells in tria 0 and tria 1: "
          << tria_0.n_active_cells() << ' ' << tria_1.n_active_cells()
          << std::endl;
  const CellList cell_list = GridTools::get_finest_common_cells(tria_0, tria_1);
  for (typename CellList::const_iterator cell_pair = cell_list.begin();
       cell_pair != cell_list.end();
       ++cell_pair)
    deallog << cell_pair->first << ' ' << cell_pair->second << std::endl;
}


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  MPILogInitAll all;

  test<2>();
  test<3>();
}
