// ---------------------------------------------------------------------
//
// Copyright (C) 2021 - 2023 by the deal.II authors
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


template <int dim, int spacedim = dim>
class MyPolicy : public RepartitioningPolicyTools::Base<dim, spacedim>
{
public:
  MyPolicy(const MPI_Comm comm, const unsigned int direction)
    : comm(comm)
    , direction(direction)
  {}

  virtual LinearAlgebra::distributed::Vector<double>
  partition(const Triangulation<dim, spacedim> &tria_in) const override
  {
    const auto tria =
      dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
        &tria_in);

    Assert(tria, ExcNotImplemented());

    LinearAlgebra::distributed::Vector<double> partition(
      tria->global_active_cell_index_partitioner().lock());

    const unsigned int n_partitions = Utilities::MPI::n_mpi_processes(comm);

    for (const auto &cell :
         tria_in.active_cell_iterators() | IteratorFilters::LocallyOwnedCell())
      partition[cell->global_active_cell_index()] =
        std::floor(cell->center()[direction] * n_partitions);

    partition.update_ghost_values();

    return partition;
  }

private:
  const MPI_Comm     comm;
  const unsigned int direction;
};


template <int dim>
void
test(const MPI_Comm comm)
{
  parallel::distributed::Triangulation<dim> tria(
    comm,
    Triangulation<dim>::none,
    parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy);
  GridGenerator::subdivided_hyper_cube(tria, 4);
  tria.refine_global(2);

  MyPolicy<dim> policy_0(comm, 0);
  MyPolicy<dim> policy_1(comm, 1);

  const auto partition_0 = policy_0.partition(tria);

  const auto settings =
    TriangulationDescription::Settings::construct_multigrid_hierarchy;

  // repartition triangulation so that it has strided partitioning
  const auto construction_data =
    TriangulationDescription::Utilities::create_description_from_triangulation(
      tria, partition_0, settings);

  parallel::fullydistributed::Triangulation<dim> tria_pft(comm);
  tria_pft.create_triangulation(construction_data);

  {
    FE_Q<dim>       fe(2);
    DoFHandler<dim> dof_handler(tria_pft);
    dof_handler.distribute_dofs(fe);
    dof_handler.distribute_mg_dofs();

    // print statistics
    print_statistics(tria_pft);
    print_statistics(dof_handler);
  }

  GridOut go;
  go.write_mesh_per_processor_as_vtu(tria_pft, "mesh_old");

  tria_pft.set_partitioner(policy_1, settings);

  tria_pft.repartition();

  go.write_mesh_per_processor_as_vtu(tria_pft, "mesh_new");

  {
    FE_Q<dim>       fe(2);
    DoFHandler<dim> dof_handler(tria_pft);
    dof_handler.distribute_dofs(fe);
    dof_handler.distribute_mg_dofs();

    // print statistics
    print_statistics(tria_pft);
    print_statistics(dof_handler);
  }
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);
  MPILogInitAll                    all;

  MPI_Comm comm = MPI_COMM_WORLD;

  test<2>(comm);
  deallog.pop();
}
