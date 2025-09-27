// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2024 by the deal.II authors
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
//
// like _01, but first partition it with one sets of weights, and then
// partition it again with all equal weights. this should yield the
// same mesh as if there had been no weights at all to begin with

#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>

#include "../tests.h"



unsigned int current_cell_weight;

template <int dim>
unsigned int
cell_weight_1(
  const typename parallel::distributed::Triangulation<dim>::cell_iterator &cell,
  const CellStatus status)
{
  return current_cell_weight++;
}

template <int dim>
unsigned int
cell_weight_2(
  const typename parallel::distributed::Triangulation<dim>::cell_iterator &cell,
  const CellStatus status)
{
  return 1;
}


template <int dim>
void
test()
{
  unsigned int myid    = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  parallel::distributed::Triangulation<dim> tr(
    MPI_COMM_WORLD,
    dealii::Triangulation<dim>::none,
    parallel::distributed::Triangulation<dim>::no_automatic_repartitioning);

  GridGenerator::subdivided_hyper_cube(tr, 16);
  tr.refine_global(1);

  current_cell_weight = 1;

  // repartition the mesh as described above, first in some arbitrary
  // way, and then with all equal weights
  tr.signals.weight.connect(&cell_weight_1<dim>);
  tr.repartition();

  tr.signals.weight.disconnect_all_slots();

  tr.signals.weight.connect(&cell_weight_2<dim>);
  tr.repartition();

  const auto n_locally_owned_active_cells_per_processor =
    Utilities::MPI::all_gather(tr.get_mpi_communicator(),
                               tr.n_locally_owned_active_cells());
  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    for (unsigned int p = 0; p < numproc; ++p)
      deallog << "processor " << p << ": "
              << n_locally_owned_active_cells_per_processor[p]
              << " locally owned active cells" << std::endl;
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
