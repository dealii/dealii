// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// hp::Refinement::choose_p_over_h called future_fe_index_set() on cells that
// are not locally owned and triggered an assertion at some point.


#include <deal.II/distributed/shared_tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/filtered_iterator.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/refinement.h>

#include "../tests.h"

#include "../test_grids.h"


template <int dim>
void
test()
{
  // setup
  parallel::shared::Triangulation<dim> tr(
    MPI_COMM_WORLD,
    ::Triangulation<dim>::none,
    false,
    parallel::shared::Triangulation<dim>::partition_custom_signal);
  tr.signals.post_refinement.connect([&tr]() {
    // partition the triangulation by hand
    for (const auto &cell : tr.active_cell_iterators())
      cell->set_subdomain_id(cell->active_cell_index() %
                             Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD));
  });

  TestGrids::hyper_line(tr, 2);
  tr.refine_global(1);

  hp::FECollection<dim> fes;
  for (unsigned int d = 1; d <= 2; ++d)
    fes.push_back(FE_Q<dim>(d));

  DoFHandler<dim> dh(tr);

  // set flags
  for (auto cell = dh.begin(0); cell != dh.end(0); ++cell)
    {
      if (cell->id().to_string() == "0_0:")
        {
          // all children will be flagged for both h- and p-adaptation.
          // choose_p_over_h() will decide in favor of p-adaptation.
          for (unsigned int i = 0; i < cell->n_children(); ++i)
            {
              const auto &child = cell->child(i);

              child->set_coarsen_flag();

              if (child->is_locally_owned())
                child->set_future_fe_index(1);
            }
        }
      else if (cell->id().to_string() == "1_0:")
        {
          // all children will be flagged for both h-adaptation
          // and only one of them for p-adaptation.
          // choose_p_over_h() will decide in favor of h-adaptation.
          for (unsigned int i = 0; i < cell->n_children(); ++i)
            {
              const auto &child = cell->child(i);

              child->set_coarsen_flag();

              if (child->is_locally_owned())
                if (i == 0)
                  child->set_future_fe_index(1);
            }
        }
    }

  dh.distribute_dofs(fes);

  // decide between p and h flags
  hp::Refinement::choose_p_over_h(dh);

  // verify
  unsigned int h_flagged_cells = 0, p_flagged_cells = 0;
  for (const auto &cell :
       dh.active_cell_iterators() | IteratorFilters::LocallyOwnedCell())
    {
      if (cell->coarsen_flag_set())
        ++h_flagged_cells;
      if (cell->future_fe_index_set())
        ++p_flagged_cells;
    }
  const unsigned int global_h_flagged_cells =
                       Utilities::MPI::sum(h_flagged_cells, MPI_COMM_WORLD),
                     global_p_flagged_cells =
                       Utilities::MPI::sum(p_flagged_cells, MPI_COMM_WORLD);

  deallog << "h-flags:" << global_h_flagged_cells << std::endl
          << "p-flags:" << global_p_flagged_cells << std::endl
          << "OK" << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    {
      initlog();

      deallog.push("2d");
      test<2>();
      deallog.pop();
      deallog.push("3d");
      test<3>();
      deallog.pop();
    }
  else
    {
      test<2>();
      test<3>();
    }
}
