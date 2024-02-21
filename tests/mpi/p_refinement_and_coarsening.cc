// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Check execution of p-refinement and p-coarsening via
// parallel::distributed::Triangulation::execute_coarsening_and_refinement()


#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/hp/fe_collection.h>

#include "../tests.h"



template <int dim>
void
test()
{
  // we pick a subdivided hyper cube so that the domain is
  // distributed on more than one processor.
  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::subdivided_hyper_cube(tria, 2);
  tria.refine_global(1);

  hp::FECollection<dim> fe;
  for (unsigned int i = 0; i < std::pow(2, dim); ++i)
    fe.push_back(FE_Q<dim>(1));

  DoFHandler<dim> dh(tria);

  // set future_fe_indices
  unsigned int future_feidx = 0;
  for (const auto &cell :
       dh.active_cell_iterators() | IteratorFilters::LocallyOwnedCell())
    {
      // check if cell is initialized correctly
      Assert(cell->active_fe_index() == 0, ExcInternalError());

      cell->set_future_fe_index(future_feidx);
      future_feidx = ((future_feidx + 1) < fe.size()) ? future_feidx + 1 : 0;
    }

  dh.distribute_dofs(fe);
  tria.execute_coarsening_and_refinement();

  // check if all flags were cleared and verify fe_indices
  for (const auto &cell :
       dh.active_cell_iterators() | IteratorFilters::LocallyOwnedCell())
    {
      Assert(cell->future_fe_index_set() == false, ExcInternalError());

      deallog << "cell:" << cell->id().to_string()
              << ", fe_index:" << cell->active_fe_index() << std::endl;
    }
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  deallog.push("2d");
  test<2>();
  deallog.pop();

  deallog.push("3d");
  test<3>();
  deallog.pop();

  deallog << "OK" << std::endl;
}
