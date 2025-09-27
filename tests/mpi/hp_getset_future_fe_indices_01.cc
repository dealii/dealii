// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Transfer future FE indices between DoFHandler objects.


#include <deal.II/base/mpi.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <algorithm>

#include "../tests.h"

#include "../test_grids.h"


template <int dim, int spacedim>
void
test()
{
  MPI_Comm mpi_comm = MPI_COMM_WORLD;

  const unsigned int n_procs = Utilities::MPI::n_mpi_processes(mpi_comm);

  parallel::distributed::Triangulation<dim, spacedim> tria(mpi_comm);
  TestGrids::hyper_line(tria, n_procs);
  tria.refine_global(1);

  DoFHandler<dim, spacedim> dofh_src(tria);
  DoFHandler<dim, spacedim> dofh_dst(tria);

  // set future fe index on first cell of every process
  for (const auto &cell : dofh_src.active_cell_iterators())
    if (cell->is_locally_owned())
      {
        cell->set_future_fe_index(1);
        break;
      }

  // transfer future FE indices
  dofh_dst.set_future_fe_indices(dofh_src.get_future_fe_indices());

  // count copied future FE indices
  const unsigned int n_indices_copied =
    std::count_if(dofh_dst.active_cell_iterators().begin(),
                  dofh_dst.active_cell_iterators().end(),
                  [](const auto &cell) {
                    return (cell->is_locally_owned() &&
                            cell->future_fe_index_set());
                  });
  deallog << "indices copied: " << n_indices_copied << std::endl;

  // verify that future FE indices match
  for (const auto &cell : tria.active_cell_iterators())
    if (cell->is_locally_owned())
      {
        const auto cell_src = cell->as_dof_handler_iterator(dofh_src);
        const auto cell_dst = cell->as_dof_handler_iterator(dofh_dst);

        AssertThrow(cell_src->future_fe_index_set() ==
                      cell_dst->future_fe_index_set(),
                    ExcInternalError());
        AssertThrow(cell_src->future_fe_index() == cell_dst->future_fe_index(),
                    ExcInternalError());
      }

  deallog << "OK" << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  test<2, 2>();
}
