// ---------------------------------------------------------------------
//
// Copyright (C) 2021 by the deal.II authors
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


// Test
// TriangulationDescription::Utilities::create_description_from_triangulation()
// with repartitioning capabilities.

#include <deal.II/base/mpi_consensus_algorithms.h>

#include <deal.II/distributed/fully_distributed_tria.h>
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


template <int dim, int spacedim>
LinearAlgebra::distributed::Vector<double>
partition_distributed_triangulation(const Triangulation<dim, spacedim> &tria_in,
                                    const unsigned int n_partitions)
{
  const auto tria =
    dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(&tria_in);

  Assert(tria, ExcNotImplemented());

  LinearAlgebra::distributed::Vector<double> partition(
    tria->global_active_cell_index_partitioner().lock());

  for (const auto &cell :
       tria_in.active_cell_iterators() | IteratorFilters::LocallyOwnedCell())
    partition[cell->global_active_cell_index()] =
      std::floor(cell->center()[0] * n_partitions);

  partition.update_ghost_values();

  return partition;
}


template <int dim>
void
test(const MPI_Comm comm, const unsigned int n_partitions)
{
  parallel::distributed::Triangulation<dim> tria(
    comm,
    Triangulation<dim>::none,
    parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy);
  GridGenerator::subdivided_hyper_cube(tria, 4);
  tria.refine_global(3);

  const auto partition_new =
    partition_distributed_triangulation(tria, n_partitions);

  // repartition triangulation so that it has strided partitioning
  const auto construction_data =
    TriangulationDescription::Utilities::create_description_from_triangulation(
      tria,
      partition_new,
      TriangulationDescription::Settings::construct_multigrid_hierarchy);

  parallel::fullydistributed::Triangulation<dim> tria_pft(comm);
  tria_pft.create_triangulation(construction_data);

  FE_Q<dim>       fe(2);
  DoFHandler<dim> dof_handler(tria_pft);
  dof_handler.distribute_dofs(fe);
  dof_handler.distribute_mg_dofs();

  // print statistics
  print_statistics(tria_pft);
  print_statistics(dof_handler);
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);
  MPILogInitAll                    all;

  MPI_Comm comm = MPI_COMM_WORLD;

  deallog.push("all");
  test<2>(comm, Utilities::MPI::n_mpi_processes(comm));
  deallog.pop();

  // test that we can eliminate processes
  deallog.push("reduced");
  test<2>(comm,
          std::max<unsigned int>(1, Utilities::MPI::n_mpi_processes(comm) / 2));
  deallog.pop();
}
