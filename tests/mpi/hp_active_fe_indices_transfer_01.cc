// ---------------------------------------------------------------------
//
// Copyright (C) 2018 by the deal.II authors
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



// ActiveFEIndicesTransfer Test


#include <deal.II/distributed/active_fe_indices_transfer.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/hp/dof_handler.h>

#include "../tests.h"


template <int dim>
void
test()
{
  const unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  // ------ setup ------
  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);

  GridGenerator::subdivided_hyper_cube(tria, 2);
  tria.refine_global(1);
  deallog << "cells before: " << tria.n_global_active_cells() << std::endl;

  hp::DoFHandler<dim>   dh(tria);
  hp::FECollection<dim> fe_collection;

  // prepare FECollection with arbitrary number of entries
  const unsigned int max_degree = 1 + std::pow(2, dim);
  for (unsigned int i = 0; i < max_degree; ++i)
    fe_collection.push_back(FE_Q<dim>(max_degree - i));

  typename hp::DoFHandler<dim, dim>::active_cell_iterator cell;
  unsigned int                                            i = 0;

  for (cell = dh.begin_active(); cell != dh.end(); ++cell)
    {
      if (cell->is_locally_owned())
        {
          // set active fe index
          if (!(cell->is_artificial()))
            {
              if (i >= fe_collection.size())
                i = 0;
              cell->set_active_fe_index(i++);
            }

          // set refinement/coarsening flags
          if (cell->id().to_string() == "0_1:0")
            cell->set_refine_flag();
          else if (cell->parent()->id().to_string() ==
                   ((dim == 2) ? "3_0:" : "7_0:"))
            cell->set_coarsen_flag();

          deallog << "myid=" << myid << " cellid=" << cell->id()
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

  // ----- transfer -----
  parallel::distributed::ActiveFEIndicesTransfer<dim, dim> feidx_transfer(dh);

  feidx_transfer.prepare_for_transfer();
  tria.execute_coarsening_and_refinement();
  deallog << "cells after: " << tria.n_global_active_cells() << std::endl;

  feidx_transfer.unpack();

  // for further calculations, distribute dofs after unpacking, i.e.
  // dh.distribute_dofs(fe_collection);

  // ------ verify ------
  // check if all children adopted the correct id
  for (cell = dh.begin_active(); cell != dh.end(); ++cell)
    {
      if (cell->is_locally_owned())
        {
          deallog << "myid=" << myid << " cellid=" << cell->id()
                  << " fe_index=" << cell->active_fe_index()
                  << " feq_degree=" << max_degree - cell->active_fe_index()
                  << std::endl;
        }
    }

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
