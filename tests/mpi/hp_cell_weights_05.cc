// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Same as test _01, but with precomputed weights.


#include <deal.II/distributed/cell_weights.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include "../tests.h"


template <int dim>
void
test()
{
  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);

  GridGenerator::hyper_cube(tria);
  tria.refine_global(2);

  // Apply ndof cell weights.
  hp::FECollection<dim> fes;
  fes.push_back(FE_Q<dim>(1));
  fes.push_back(FE_Q<dim>(5));

  DoFHandler<dim> dh(tria);

  // default: active_fe_index = 0
  for (auto &cell : dh.active_cell_iterators())
    if (cell->is_locally_owned())
      if (cell->id().to_string() == "0_2:00")
        cell->set_active_fe_index(1);

  dh.distribute_dofs(fes);

  deallog << "Number of cells before repartitioning: "
          << tria.n_locally_owned_active_cells() << std::endl;
  {
    unsigned int dof_counter = 0;
    for (auto &cell : dh.active_cell_iterators())
      if (cell->is_locally_owned())
        dof_counter += cell->get_fe().dofs_per_cell;
    deallog << "  Cumulative dofs per cell: " << dof_counter << std::endl;
  }


  const auto weighting_function =
    parallel::CellWeights<dim>::ndofs_weighting({1, 1});
  const auto precomputed_weights =
    parallel::CellWeights<dim>::precompute_weights(fes, weighting_function);
  const parallel::CellWeights<dim> cell_weights(dh, precomputed_weights);

  tria.repartition();


  deallog << "Number of cells after repartitioning: "
          << tria.n_locally_owned_active_cells() << std::endl;
  {
    unsigned int dof_counter = 0;
    for (auto &cell : dh.active_cell_iterators())
      if (cell->is_locally_owned())
        dof_counter += cell->get_fe().dofs_per_cell;
    deallog << "  Cumulative dofs per cell: " << dof_counter << std::endl;
  }

  // make sure no processor is hanging
  MPI_Barrier(MPI_COMM_WORLD);

  deallog << "OK" << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  deallog.push("2d");
  test<2>();
  deallog.pop();
  deallog.push("3d");
  test<3>();
  deallog.pop();
}
