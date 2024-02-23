// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Verify fixed fraction algorithm with l1-norm and l2-norm
// Equidistant indicators

#include <deal.II/base/utilities.h>

#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/vector.h>

#include "../tests.h"


template <int dim>
void
count_flags(const parallel::distributed::Triangulation<dim> &tria)
{
  unsigned int n_refine = 0, n_coarsen = 0;
  for (const auto &cell : tria.active_cell_iterators())
    if (cell->is_locally_owned())
      {
        if (cell->refine_flag_set())
          ++n_refine;
        else if (cell->coarsen_flag_set())
          ++n_coarsen;
      }

  const unsigned int n_refine_global =
                       Utilities::MPI::sum(n_refine, MPI_COMM_WORLD),
                     n_coarsen_global =
                       Utilities::MPI::sum(n_coarsen, MPI_COMM_WORLD);

  deallog << "n_refine_flags: " << n_refine_global
          << ", n_coarsen_flags: " << n_coarsen_global;
}



template <int dim>
void
test()
{
  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(2);

  Vector<float> indicator(tria.n_active_cells());
  // assign each cell a globally unique cellid
  for (const auto &cell :
       tria.active_cell_iterators() | IteratorFilters::LocallyOwnedCell())
    {
      const std::string  cellid = cell->id().to_string();
      const unsigned int fine_cellid =
        std::stoul(cellid.substr(cellid.find(':') + 1, std::string::npos));

      indicator[cell->active_cell_index()] = fine_cellid + 1;
    }

  deallog << "l1-norm: ";
  parallel::distributed::GridRefinement::refine_and_coarsen_fixed_fraction(
    tria, indicator, 0.3, 0.3, VectorTools::L1_norm);
  count_flags(tria);
  deallog << std::endl;

  // reset refinement flags
  for (const auto &cell :
       tria.active_cell_iterators() | IteratorFilters::LocallyOwnedCell())
    {
      cell->clear_refine_flag();
      cell->clear_coarsen_flag();
    }

  deallog << "l2-norm: ";
  parallel::distributed::GridRefinement::refine_and_coarsen_fixed_fraction(
    tria, indicator, 0.3, 0.3, VectorTools::L2_norm);
  count_flags(tria);
  deallog << std::endl;
}



int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  const unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  if (myid == 0)
    {
      initlog();

      test<2>();
    }
  else
    test<2>();
}
