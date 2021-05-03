// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2018 by the deal.II authors
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


// The purpose of this test is to ensure that the p::s::T is
// refined/coarsened in the same way on all participating processors.

#include <deal.II/distributed/shared_tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>

#include "../tests.h"



template <int dim>
void
check()
{
  parallel::shared::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(3);

  Vector<float> estimated_error_per_cell;
  estimated_error_per_cell.reinit(tria.n_active_cells());
  for (unsigned int i = 0; i < estimated_error_per_cell.size(); ++i)
    estimated_error_per_cell(i) = i + 1;

  GridRefinement::refine_and_coarsen_fixed_number(tria,
                                                  estimated_error_per_cell,
                                                  0.2,
                                                  0.1);

  std::vector<bool> refined_cells(tria.n_active_cells() * dim);
  tria.save_refine_flags(refined_cells);
  int n_refined_cells =
    std::count(refined_cells.begin(), refined_cells.end(), true);

  std::vector<bool> coarsened_cells(tria.n_active_cells() * dim);
  tria.save_coarsen_flags(coarsened_cells);
  int n_coarsened_cells =
    std::count(coarsened_cells.begin(), coarsened_cells.end(), true);

  tria.execute_coarsening_and_refinement();
  int n_cells = tria.n_active_cells();

  if (Utilities::MPI::max(n_refined_cells, MPI_COMM_WORLD) ==
      Utilities::MPI::min(n_refined_cells, MPI_COMM_WORLD))
    deallog << "OK" << std::endl;

  if (Utilities::MPI::max(n_coarsened_cells, MPI_COMM_WORLD) ==
      Utilities::MPI::min(n_coarsened_cells, MPI_COMM_WORLD))
    deallog << "OK" << std::endl;

  if (Utilities::MPI::max(n_cells, MPI_COMM_WORLD) ==
      Utilities::MPI::min(n_cells, MPI_COMM_WORLD))
    deallog << "OK" << std::endl;
}


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  mpi_initlog();

  deallog.push("2d");
  check<2>();
  deallog.pop();

  deallog.push("3d");
  check<3>();
  deallog.pop();

  return 0;
}
