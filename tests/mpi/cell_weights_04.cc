// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// just create a 16x16 coarse mesh, refine it once, and partition it
// like _03, but with a larger spread of weights

#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>

#include "../tests.h"



template <int dim>
unsigned int
cell_weight(
  const typename parallel::distributed::Triangulation<dim>::cell_iterator &cell,
  const CellStatus status)
{
  const unsigned int cell_weight =
    (cell->center()[0] < 0.5 || cell->center()[1] < 0.5 ? 1 : 1000);

  return cell_weight;
}

template <int dim>
void
test()
{
  unsigned int myid    = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  parallel::distributed::Triangulation<dim> tr(
    MPI_COMM_WORLD, dealii::Triangulation<dim>::none);

  GridGenerator::subdivided_hyper_cube(tr, 16);

  tr.signals.weight.connect(&cell_weight<dim>);

  tr.refine_global(1);

  const auto n_locally_owned_active_cells_per_processor =
    Utilities::MPI::all_gather(tr.get_mpi_communicator(),
                               tr.n_locally_owned_active_cells());
  if (myid == 0)
    for (unsigned int p = 0; p < numproc; ++p)
      deallog << "processor " << p << ": "
              << n_locally_owned_active_cells_per_processor[p]
              << " locally owned active cells" << std::endl;

  // let each processor sum up its weights
  std::vector<unsigned int> integrated_weights(numproc, 0);
  for (const auto &cell :
       tr.active_cell_iterators() | IteratorFilters::LocallyOwnedCell())
    integrated_weights[myid] +=
      cell_weight<dim>(cell, CellStatus::cell_will_persist);

  Utilities::MPI::sum(integrated_weights, MPI_COMM_WORLD, integrated_weights);
  if (myid == 0)
    for (unsigned int p = 0; p < numproc; ++p)
      deallog << "processor " << p << ": " << integrated_weights[p] << " weight"
              << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  if (myid == 0)
    {
      initlog();

      deallog.push("2d");
      test<2>();
      deallog.pop();
      deallog.push("3d");
      test<3>();
      deallog.pop();
    }
  else
    {
      test<2>();
      test<3>();
    }
}
