// ---------------------------------------------------------------------
//
// Copyright (C) 2021 - 2022 by the deal.II authors
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


// Test RepartitioningPolicyTools::DefaultPolicy's tightening functionality.

#include <deal.II/base/mpi_consensus_algorithms.h>

#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/distributed/repartitioning_policy_tools.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria_description.h>

#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/numerics/data_out.h>

#include "../grid/tests.h"

using namespace dealii;

template <int dim>
void
test(const MPI_Comm comm)
{
  parallel::shared::Triangulation<dim> tria(
    comm,
    Triangulation<dim>::none,
    true,
    parallel::shared::Triangulation<dim>::partition_custom_signal);

  tria.signals.create.connect([&tria]() {
    for (const auto &cell : tria.active_cell_iterators())
      {
        switch (cell->active_cell_index())
          {
            case 0:
              cell->set_subdomain_id(0);
              break;
            case 1:
            case 2:
              cell->set_subdomain_id(2);
              break;
            case 3:
              cell->set_subdomain_id(3);
              break;
          }
      }
  });

  GridGenerator::subdivided_hyper_cube(tria, 2);

  const auto construction_data =
    TriangulationDescription::Utilities::create_description_from_triangulation(
      tria,
      RepartitioningPolicyTools::DefaultPolicy<dim>(true).partition(tria));

  parallel::fullydistributed::Triangulation<dim> tria_pft(comm);
  tria_pft.create_triangulation(construction_data);

  deallog << tria.n_locally_owned_active_cells() << ' '
          << tria_pft.n_locally_owned_active_cells() << std::endl;
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);
  MPILogInitAll                    all;

  MPI_Comm comm = MPI_COMM_WORLD;

  deallog.push("all");
  test<2>(comm);
  deallog.pop();
}
