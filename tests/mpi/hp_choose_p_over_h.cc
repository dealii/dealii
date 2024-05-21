// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// verify hp::Refinement::choose_p_over_h() on strongly distributed meshes


#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/refinement.h>

#include "../tests.h"



void
test()
{
  constexpr unsigned int dim = 2;

  // Setup distributed triangulation.
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tr, -1, 1);
  tr.refine_global(1);

  tr.begin_active()->set_refine_flag();
  tr.execute_coarsening_and_refinement();

  // For two dimensions in this particular scenario, each quadrant in the
  // coordinate system will be assigned to single MPI process.
  //
  // See also the corresponding output file and this ASCII sketch on how p4est
  // decided to partition this mesh with four MPI processes.
  //
  //  +---+---+
  //  |   |   |
  //  | 3 | 3 |
  //  |   |   |
  //  +-+-+---+
  //  |1|1|   |
  //  +-+-+ 2 |
  //  |1|1|   |
  //  +-+-+---+

  // Set h-coarsen and p-refinement flags on the coarser cells.
  DoFHandler<dim> dh(tr);
  for (const auto &cell : dh.active_cell_iterators_on_level(1) |
                            IteratorFilters::LocallyOwnedCell())
    {
      cell->set_coarsen_flag();
      cell->set_future_fe_index(1);
    }

  // Setup FEs.
  hp::FECollection<dim> fes;
  for (unsigned int d = 1; d <= 2; ++d)
    fes.push_back(FE_Q<dim>(d));
  dh.distribute_dofs(fes);

  // Check whether choose_p_over_h() finishes.
  // (It triggered an assertion in the past.)
  hp::Refinement::choose_p_over_h(dh);

  // Check the new flags.
  for (const auto &cell :
       dh.active_cell_iterators() | IteratorFilters::LocallyOwnedCell())
    deallog << cell->id().to_string()
            << ": future_fe=" << cell->future_fe_index()
            << " coarsen=" << cell->coarsen_flag_set() << std::endl;

  deallog << "OK" << std::endl;
}



int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  test();
}
