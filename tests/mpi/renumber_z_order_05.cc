// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test that DofRenumbering::hierarchical() gives the same DoF numbers
// independent of their previous ordering. This was broken before where
// the hierarchical numbering on each processor re-used the indices
// the processor previously owned. This worked just fine with the
// original ordering used by default by DoFHandler, but not if there
// was another reordering step in between.
//
// This test is like the _04 test, but for dim!=spacedim


#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/trilinos_vector.h>

#include <boost/concept_check.hpp>

#include "../tests.h"


void
test()
{
  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  parallel::distributed::Triangulation<2, 3> tr(MPI_COMM_WORLD);
  GridGenerator::subdivided_hyper_rectangle(tr,
                                            std::vector<unsigned int>{{2, 2}},
                                            Point<2>(),
                                            Point<2>(2, 2));

  const FE_Q<2, 3> fe(1);

  for (unsigned int test = 0; test < 2; ++test)
    {
      DoFHandler<2, 3> dof_handler(tr);
      dof_handler.distribute_dofs(fe);

      // in the second test run, revert the global order of DoF
      // indices. the point is to make sure that processors do not
      // have strictly increasing, contiguous groups of DoF indices.
      if (test == 1)
        {
          IndexSet locally_owned_dofs = dof_handler.locally_owned_dofs();
          std::vector<types::global_dof_index> new_numbers(
            locally_owned_dofs.n_elements());
          for (auto i : locally_owned_dofs)
            new_numbers[locally_owned_dofs.index_within_set(i)] =
              dof_handler.n_dofs() - i - 1;
          dof_handler.renumber_dofs(new_numbers);
        }

      // next, let the hierarchical numbering work on it. the
      // important aspect is that the result of hierarchical
      // renumbering should not depend on the *previous* numbering. so
      // we should be getting the same output as in the first run
      DoFRenumbering::hierarchical(dof_handler);

      // output DoF indices
      deallog << (test == 0 ? "Without " : "With ")
              << "prior reordering:" << std::endl;
      const unsigned int                     dofs_per_cell = fe.dofs_per_cell;
      std::vector<types::global_dof_index>   local_dof_indices(dofs_per_cell);
      DoFHandler<2, 3>::active_cell_iterator cell = dof_handler.begin_active(),
                                             endc = dof_handler.end();
      for (; cell != endc; ++cell)
        if (cell->subdomain_id() == tr.locally_owned_subdomain())
          {
            deallog << "Cell=" << cell << std::endl;
            cell->get_dof_indices(local_dof_indices);
            for (auto i : local_dof_indices)
              deallog << i << ' ';
            deallog << std::endl;
          }
    }
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  test();
}
