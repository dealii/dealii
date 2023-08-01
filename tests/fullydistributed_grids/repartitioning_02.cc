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


// Test different partitioning policies from RepartitioningPolicyTools.


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

template <int dim, int spacedim>
void
create_triangulation(Triangulation<dim, spacedim> &tria)
{
  const unsigned int n_ref_global = 3;
  const unsigned int n_ref_local  = 3;
  GridGenerator::hyper_cube(tria, -1.0, +1.0);
  tria.refine_global(n_ref_global);

  for (unsigned int i = 0; i < n_ref_local; ++i)
    {
      for (auto cell : tria.active_cell_iterators())
        if (cell->is_locally_owned())
          {
            bool flag = true;
            for (unsigned int d = 0; d < dim; ++d)
              if (cell->center()[d] > 0.0)
                flag = false;
            if (flag)
              cell->set_refine_flag();
          }
      tria.execute_coarsening_and_refinement();
    }
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);
  MPILogInitAll                    all;

  const unsigned int dim  = 2;
  const MPI_Comm     comm = MPI_COMM_WORLD;

  parallel::distributed::Triangulation<dim> tria(comm);
  create_triangulation(tria);

  const auto test = [&](const auto &policy, const std::string &label) {
    const auto partition_new = policy.partition(tria);
    partition_new.update_ghost_values();

    // repartition triangulation so that it has strided partitioning
    const auto construction_data =
      partition_new.size() == 0 ?
        TriangulationDescription::Utilities::
          create_description_from_triangulation(tria, comm) :
        TriangulationDescription::Utilities::
          create_description_from_triangulation(tria, partition_new);

    parallel::fullydistributed::Triangulation<dim> tria_pft(comm);
    tria_pft.create_triangulation(construction_data);

    FE_Q<dim>       fe(2);
    DoFHandler<dim> dof_handler(tria_pft);
    dof_handler.distribute_dofs(fe);

    // print statistics
    deallog.push(label);
    print_statistics(tria_pft);
    print_statistics(dof_handler);
    deallog.pop();
  };

  // default policy (simply copy p:d:T)
  test(RepartitioningPolicyTools::DefaultPolicy<dim>(), "grid_policy_default");

  // first-child policy
  test(RepartitioningPolicyTools::FirstChildPolicy<dim>(tria),
       "grid_policy_first");

  // first-child policy
  test(RepartitioningPolicyTools::MinimalGranularityPolicy<dim>(4),
       "grid_policy_minimal_4");
  test(RepartitioningPolicyTools::MinimalGranularityPolicy<dim>(100),
       "grid_policy_minimal_100");
  test(RepartitioningPolicyTools::MinimalGranularityPolicy<dim>(500),
       "grid_policy_minimal_500");
}
