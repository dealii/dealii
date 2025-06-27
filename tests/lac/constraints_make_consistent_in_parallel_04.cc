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



// test AffineConstraints<double>::make_consistent_in_parallel for case
// where constraints need to be combined between different ranks

#include <deal.II/base/index_set.h>

#include <deal.II/lac/affine_constraints.h>

#include "../tests.h"

using namespace dealii;

void
test()
{
  IndexSet locally_owned_dofs(3);
  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    locally_owned_dofs.add_range(0, 2);
  else if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 1)
    locally_owned_dofs.add_range(2, 3);

  IndexSet                  locally_relevant_dofs = complete_index_set(3);
  AffineConstraints<double> constraints(locally_owned_dofs,
                                        locally_relevant_dofs);

  const auto show_constraints_1_2 = [&]() {
    deallog << "What process "
            << Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)
            << " believes about DoF 1:" << std::endl;
    for (const auto &c : constraints.get_lines())
      if (c.index == 1)
        for (const auto &entry : c.entries)
          deallog << "    constrained against " << entry.first
                  << " with weight " << entry.second << std::endl;
    deallog << "What process "
            << Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)
            << " believes about DoF 2:" << std::endl;
    for (const auto &c : constraints.get_lines())
      if (c.index == 2)
        for (const auto &entry : c.entries)
          deallog << "    constrained against " << entry.first
                  << " with weight " << entry.second << std::endl;
  };


  deallog << "------------- make_hanging_node_constraints():" << std::endl;
  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    {
      constraints.add_constraint(1, {{0, 0.5}});
      constraints.add_constraint(2, {{0, 0.25}});
    }
  else if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 1)
    {
      constraints.add_constraint(2, {{1, 0.5}});
    }
  show_constraints_1_2();

  deallog << "------------- make_consistent_in_parallel():" << std::endl;
  constraints.make_consistent_in_parallel(locally_owned_dofs,
                                          locally_relevant_dofs,
                                          MPI_COMM_WORLD);
  show_constraints_1_2();
}



int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);
  MPILogInitAll                    all;

  test();
}
