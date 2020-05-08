// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2020 by the deal.II authors
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



// just create a 16x16 coarse mesh, and partition it
//
// like _03, but with a spread of weights that ensures that one
// processor only should carry one cell. we do this by giving all
// cells a weight of one, except for one cell that has a weight equal
// to the total number of cells in the mesh (i.e., a weight one larger
// than that of all other cells together).
//
// p4est does not seem to partition this correctly right now. it
// should have one processor that stores exactly one cell, and the
// remaining ones should be split equally among the other processors;
// at present, however, the remaining processors do not split the
// remaining cells even close to equally

#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>

#include "../tests.h"



template <int dim>
unsigned int
cell_weight(
  const typename parallel::distributed::Triangulation<dim>::cell_iterator &cell,
  const typename parallel::distributed::Triangulation<dim>::CellStatus status)
{
  return (
    // bottom left corner
    (cell->center()[0] < 1) && (cell->center()[1] < 1) &&
        (dim == 3 ? (cell->center()[2] < 1) : true) ?
      // one cell has more weight than all others together
      Utilities::pow(16, dim) * 1000 :
      0);
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

  // create a 16x16 or 16x16x16 mesh where each cell has size 1x1 ir 1x1x1
  GridGenerator::subdivided_hyper_cube(tr, 16, 0, 16);


  // repartition the mesh; attach different weights to all cells
  tr.signals.cell_weight.connect(
    std::bind(&cell_weight<dim>, std::placeholders::_1, std::placeholders::_2));

  tr.repartition();

  const auto n_locally_owned_active_cells_per_processor =
    Utilities::MPI::all_gather(tr.get_communicator(),
                               tr.n_locally_owned_active_cells());
  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    for (unsigned int p = 0; p < numproc; ++p)
      deallog << "processor " << p << ": "
              << n_locally_owned_active_cells_per_processor[p]
              << " locally owned active cells" << std::endl;

  // let each processor sum up its weights
  std::vector<double> integrated_weights(numproc, 0.0);
  for (typename Triangulation<dim>::active_cell_iterator cell =
         tr.begin_active();
       cell != tr.end();
       ++cell)
    if (cell->is_locally_owned())
      integrated_weights[myid] +=
        1000 + cell_weight<dim>(
                 cell, parallel::distributed::Triangulation<dim>::CELL_PERSIST);
  Utilities::MPI::sum(integrated_weights, MPI_COMM_WORLD, integrated_weights);
  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
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
