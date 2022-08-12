// ---------------------------------------------------------------------
//
// Copyright (C) 2022 by the deal.II authors
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



// Check functionality of hp::Refinement::copy_future_fe_indices.


#include <deal.II/base/mpi.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/hp/refinement.h>

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


  deallog << "indices copied: "
          << hp::Refinement::copy_future_fe_indices(dofh_src, dofh_dst)
          << std::endl;


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
