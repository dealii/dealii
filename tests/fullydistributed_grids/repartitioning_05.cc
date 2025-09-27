// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test
// TriangulationDescription::Utilities::create_description_from_triangulation()
// so that it also works for p:d:T set up on a subcommunicator.

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


MPI_Comm
create_sub_comm(const MPI_Comm comm, const unsigned int size)
{
  const auto rank = Utilities::MPI::this_mpi_process(comm);

  int color = rank < size;

  MPI_Comm sub_comm;
  int      ierr = MPI_Comm_split(comm, color, rank, &sub_comm);
  AssertThrowMPI(ierr);

  if (rank < size)
    return sub_comm;
  else
    {
      ierr = MPI_Comm_free(&sub_comm);
      AssertThrowMPI(ierr);
      return MPI_COMM_SELF;
    }
}



template <int dim, int spacedim>
LinearAlgebra::distributed::Vector<double>
partition_distributed_triangulation(const Triangulation<dim, spacedim> &tria_in,
                                    const MPI_Comm                      comm)
{
  const auto comm_tria = tria_in.get_mpi_communicator();

  const auto n_global_active_cells = Utilities::MPI::max(
    comm_tria == MPI_COMM_SELF ? 0 : tria_in.n_global_active_cells(), comm);

  if (comm_tria == MPI_COMM_SELF)
    {
      LinearAlgebra::distributed::Vector<double> partition{
        IndexSet(n_global_active_cells), IndexSet(n_global_active_cells), comm};
      partition.update_ghost_values();
      return partition;
    }

  const auto tria =
    dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(&tria_in);

  Assert(tria, ExcNotImplemented());

  LinearAlgebra::distributed::Vector<double> partition(
    tria->global_active_cell_index_partitioner().lock()->locally_owned_range(),
    tria->global_active_cell_index_partitioner().lock()->ghost_indices(),
    comm);

  const unsigned int n_partitions = Utilities::MPI::n_mpi_processes(comm);

  for (const auto &cell :
       tria_in.active_cell_iterators() | IteratorFilters::LocallyOwnedCell())
    partition[cell->global_active_cell_index()] =
      cell->global_active_cell_index() * n_partitions /
      tria_in.n_global_active_cells();

  partition.update_ghost_values();

  return partition;
}


template <int dim>
void
test(const MPI_Comm comm, const unsigned int n_partitions)
{
  auto sub_comm = create_sub_comm(comm, n_partitions);

  parallel::distributed::Triangulation<dim> tria(sub_comm);

  if (sub_comm != MPI_COMM_SELF)
    {
      GridGenerator::subdivided_hyper_cube(tria, 4);
      tria.refine_global(3);
    }

  const auto partition_new = partition_distributed_triangulation(tria, comm);
  partition_new.update_ghost_values();

  // repartition triangulation
  const auto construction_data =
    TriangulationDescription::Utilities::create_description_from_triangulation(
      tria, partition_new);

  parallel::fullydistributed::Triangulation<dim> tria_pft(comm);
  tria_pft.create_triangulation(construction_data);

  FE_Q<dim>       fe(2);
  DoFHandler<dim> dof_handler(tria_pft);
  dof_handler.distribute_dofs(fe);

  // print statistics
  print_statistics(tria_pft);
  print_statistics(dof_handler);

  if (sub_comm != MPI_COMM_SELF)
    {
      const int ierr = MPI_Comm_free(&sub_comm);
      AssertThrowMPI(ierr);
    }
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
