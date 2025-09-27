// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// We check if active_fe_indices are set correctly
// before and after coarsening and refinement.
//
// We use a parallel::share::Triangulation that allows artificial cells.
//
// This test is inspired by 'mpi/hp_active_fe_indices_transfer_01.cc'.

#include <deal.II/distributed/shared_tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>

#include <deal.II/hp/fe_collection.h>

#include "../tests.h"


template <int dim>
void
test()
{
  // ------ setup ------
  parallel::shared::Triangulation<dim> tria(
    MPI_COMM_WORLD,
    ::Triangulation<dim>::none,
    true,
    parallel::shared::Triangulation<dim>::partition_zoltan);

  GridGenerator::subdivided_hyper_cube(tria, 2);
  tria.refine_global(1);
  deallog << "cells before: " << tria.n_global_active_cells() << std::endl;

  DoFHandler<dim>       dh(tria);
  hp::FECollection<dim> fe_collection;

  // prepare FECollection with arbitrary number of entries
  const unsigned int max_degree = 1 + Utilities::pow(2, dim);
  for (unsigned int i = 0; i < max_degree; ++i)
    fe_collection.push_back(FE_Q<dim>(max_degree - i));

  typename DoFHandler<dim, dim>::active_cell_iterator cell;
  unsigned int                                        i = 0;

  for (auto &cell : dh.active_cell_iterators())
    {
      // set refinement/coarsening flags
      if (cell->id().to_string() == "0_1:0")
        cell->set_refine_flag();
      else if (cell->parent()->id().to_string() ==
               ((dim == 2) ? "3_0:" : "7_0:"))
        cell->set_coarsen_flag();

      if (cell->is_locally_owned())
        {
          // set active FE index
          if (i >= fe_collection.size())
            i = 0;
          cell->set_active_fe_index(i++);

          deallog << " cellid=" << cell->id()
                  << " fe_index=" << cell->active_fe_index()
                  << " feq_degree=" << max_degree - cell->active_fe_index();
          if (cell->coarsen_flag_set())
            deallog << " coarsening";
          else if (cell->refine_flag_set())
            deallog << " refining";
          deallog << std::endl;
        }
    }

  dh.distribute_dofs(fe_collection);

  // ----- refine -----
  tria.execute_coarsening_and_refinement();
  deallog << "cells after: " << tria.n_global_active_cells() << std::endl;

  // for further calculations, distribute dofs after refinement, i.e.
  // dh.distribute_dofs(fe_collection);

  // ------ verify ------
  // check if all children adopted the correct id
  for (auto &cell : dh.active_cell_iterators())
    if (cell->is_locally_owned())
      deallog << " cellid=" << cell->id()
              << " fe_index=" << cell->active_fe_index()
              << " feq_degree=" << max_degree - cell->active_fe_index()
              << std::endl;

  // make sure no processor is hanging
  MPI_Barrier(MPI_COMM_WORLD);

  deallog << "OK" << std::endl;
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
}
