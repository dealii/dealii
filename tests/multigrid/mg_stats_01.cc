// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test geometric metrics from MGTools.

#include <deal.II/base/utilities.h>

#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/repartitioning_policy_tools.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_transfer_global_coarsening.h>
#include <deal.II/multigrid/multigrid.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

template <int dim>
void
create_quadrant(Triangulation<dim> &tria, const unsigned int n_refinements)
{
  GridGenerator::subdivided_hyper_cube(tria, 2, -1.0, +1.0);

  if (n_refinements == 0)
    return;

  for (unsigned int i = 0; i < n_refinements; i++)
    {
      for (auto cell : tria.active_cell_iterators())
        if (cell->is_locally_owned())
          {
            bool flag = true;
            for (int d = 0; d < dim; d++)
              if (cell->center()[d] > 0.0)
                flag = false;
            if (flag)
              cell->set_refine_flag();
          }
      tria.execute_coarsening_and_refinement();
    }
}


template <typename T>
void
print_stats(const T &tria)
{
  deallog << "Workload efficiency:               "
          << 1.0 / MGTools::workload_imbalance(tria) << std::endl;

  const auto local_workload = MGTools::local_workload(tria);
  const auto local_vertical_communication_cost =
    MGTools::local_vertical_communication_cost(tria);

  deallog << "Local workload:                    ";
  for (unsigned i = 0; i < local_workload.size(); ++i)
    deallog << local_workload[i] << " ";
  deallog << std::endl;

  deallog << "Vertical communication efficiency: "
          << MGTools::vertical_communication_efficiency(tria) << std::endl;

  deallog << "Vertical communication costs:      ";
  for (unsigned i = 0; i < local_workload.size(); ++i)
    deallog << "(" << local_vertical_communication_cost[i].first << ", "
            << local_vertical_communication_cost[i].second << ") ";
  deallog << std::endl;
}


template <int dim>
void
test()
{
  parallel::distributed::Triangulation<dim> tria(
    MPI_COMM_WORLD,
    Triangulation<dim>::limit_level_difference_at_vertices,
    parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy);

  create_quadrant(tria, 2);

  print_stats(tria);

  const RepartitioningPolicyTools::FirstChildPolicy<dim> policy(tria);
  const auto                                             trias =
    MGTransferGlobalCoarseningTools::create_geometric_coarsening_sequence(
      tria, policy);

  print_stats(trias);
}



int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  deallog.push("2d");
  test<2>();
  deallog.pop();
}
