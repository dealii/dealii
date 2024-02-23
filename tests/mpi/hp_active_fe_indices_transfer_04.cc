// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// active FE indices serialization with a different number of cpus


#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include "../tests.h"


template <int dim>
void
test()
{
  unsigned int myid    = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  MPI_Comm     com_all = MPI_COMM_WORLD;
  MPI_Comm     com_small;

  // split the communicator in proc 0,1,2 and 3,4
  MPI_Comm_split(com_all, (myid < 3) ? 0 : 1, myid, &com_small);

  // prepare FECollection with arbitrary number of entries
  hp::FECollection<dim> fe_collection;
  const unsigned int    max_degree = 1 + Utilities::pow(2, dim);
  for (unsigned int i = 0; i < max_degree; ++i)
    fe_collection.push_back(FE_Q<dim>(max_degree - i));

  // write with small com
  if (myid < 3)
    {
      deallog << "writing with " << Utilities::MPI::n_mpi_processes(com_small)
              << std::endl;

      // ------ setup ------
      parallel::distributed::Triangulation<dim> tria(com_small);
      GridGenerator::subdivided_hyper_cube(tria, 2);
      tria.refine_global(1);

      DoFHandler<dim> dh(tria);

      unsigned int i = 0;
      for (auto &cell : dh.active_cell_iterators())
        if (cell->is_locally_owned())
          {
            // set active FE index
            if (i >= fe_collection.size())
              i = 0;
            cell->set_active_fe_index(i++);

            deallog << "cellid=" << cell->id()
                    << " fe_index=" << cell->active_fe_index() << std::endl;
          }

      dh.distribute_dofs(fe_collection);

      // ----- transfer -----
      dh.prepare_for_serialization_of_active_fe_indices();
      tria.save("file");
    }

  // make sure no processor is hanging
  MPI_Barrier(MPI_COMM_WORLD);

  {
    deallog << "reading with " << Utilities::MPI::n_mpi_processes(com_all)
            << std::endl;

    // ------ setup ------
    parallel::distributed::Triangulation<dim> tria(com_all);
    GridGenerator::subdivided_hyper_cube(tria, 2);
    // triangulation has to be initialized with correct coarse cells

    // we need to introduce dof_handler to its fe_collection first
    DoFHandler<dim> dh(tria);

    // ----- transfer -----
    tria.load("file");
    dh.deserialize_active_fe_indices();

    // ------ verify ------
    // check if all children adopted the correct id
    for (auto &cell : dh.active_cell_iterators())
      if (!cell->is_artificial())
        {
          deallog << "cellid=" << cell->id()
                  << " fe_index=" << cell->active_fe_index();
          if (cell->is_ghost())
            deallog << " ghost";
          deallog << std::endl;
        }

    // distribute dofs again for further calculations, i.e.
    dh.distribute_dofs(fe_collection);
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
