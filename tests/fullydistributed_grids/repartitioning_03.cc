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


// Test MGTransferGlobalCoarseningTools::create_geometric_coarsening_sequence
// with RepartitioningPolicyTools.


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

#include <deal.II/multigrid/mg_transfer_global_coarsening.h>

#include <deal.II/numerics/data_out.h>

#include "tests.h"

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

template <int dim, int spacedim>
void
output_grid(
  const std::vector<std::shared_ptr<const Triangulation<dim, spacedim>>> &trias,
  const std::string &                                                     label)
{
  deallog.push(label);
  const auto comm    = trias.front()->get_communicator();
  const auto my_rank = Utilities::MPI::this_mpi_process(comm);

  for (unsigned int i = 0; i < trias.size(); ++i)
    {
      const auto &tria = *trias[i];

      FE_Q<dim>       fe(2);
      DoFHandler<dim> dof_handler(tria);
      dof_handler.distribute_dofs(fe);

      // print statistics
      print_statistics(tria);
      print_statistics(dof_handler);

      if (false)
        {
          Vector<double> ranks(tria.n_active_cells());
          ranks = my_rank;

          DataOut<dim> data_out;
          data_out.attach_triangulation(tria);
          data_out.add_data_vector(ranks, "ranks");
          data_out.build_patches();

          data_out.write_vtu_with_pvtu_record("./", label, i, comm, 2, 0);
        }
    }
  deallog.pop();
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
    const auto trias =
      MGTransferGlobalCoarseningTools::create_geometric_coarsening_sequence(
        tria, policy);
    output_grid(trias, label);
  };

  // default policy (simply copy p:d:T)
  test(RepartitioningPolicyTools::DefaultPolicy<dim>(), "grid_policy_default");

  // first-child policy
  test(RepartitioningPolicyTools::FirstChildPolicy<dim>(tria),
       "grid_policy_first");

  // first-child policy
  test(RepartitioningPolicyTools::MinimalGranularityPolicy<dim>(4),
       "grid_policy_minimal");

  // cell-weight policy
  test(RepartitioningPolicyTools::CellWeightPolicy<dim>(
         [](const auto &cell, const auto) {
           return (cell->global_active_cell_index() <
                   cell->get_triangulation().n_global_active_cells() / 2) ?
                    1.0 :
                    2.0;
         }),
       "grid_policy_weight");
}
