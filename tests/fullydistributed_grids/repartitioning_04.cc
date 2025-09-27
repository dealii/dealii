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

#include "../grid/tests.h"


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
output_grid(const Triangulation<dim, spacedim> &tria, const std::string &label)
{
  deallog.push(label);

  FE_Q<dim>       fe(2);
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  // print statistics
  print_statistics(tria);
  print_statistics(dof_handler);

  deallog.pop();
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);
  MPILogInitAll                    all;

  const unsigned int dim  = 2;
  const MPI_Comm     comm = MPI_COMM_WORLD;

  const auto test = [&](const bool keep_fine_triangulation,
                        const bool repartition_fine_triangulation) {
    parallel::distributed::Triangulation<dim> tria(comm);
    create_triangulation(tria);

    const RepartitioningPolicyTools::DefaultPolicy<dim> policy;

    std::string label =
      (keep_fine_triangulation ? std::string("true") : std::string("false")) +
      "_" +
      (repartition_fine_triangulation ? std::string("true") :
                                        std::string("false"));

    const auto trias =
      MGTransferGlobalCoarseningTools::create_geometric_coarsening_sequence(
        tria, policy, keep_fine_triangulation, repartition_fine_triangulation);
    output_grid(tria, label);
    output_grid(*trias.back(), label);
  };

  test(true, false);
  test(true, true);
  test(false, false);
  test(false, true);
}
