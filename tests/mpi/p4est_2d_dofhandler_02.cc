// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2020 by the deal.II authors
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



// create a parallel DoFHandler on a 2d adaptively refined mesh

#include <deal.II/base/tensor.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"



template <int dim>
void
test()
{
  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    deallog << "hyper_cube" << std::endl;

  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

  GridGenerator::hyper_cube(tr);
  tr.refine_global(1);

  while (tr.n_locally_owned_active_cells() < 2000)
    {
      if (tr.n_locally_owned_active_cells())
        {
          std::vector<bool> flags(tr.n_locally_owned_active_cells(), false);
          for (unsigned int i = 0;
               i < tr.n_locally_owned_active_cells() / 5 + 1;
               ++i)
            {
              const unsigned int x = Testing::rand() % flags.size();
              flags[x]             = true;
            }

          unsigned int index = 0;
          for (typename Triangulation<dim>::active_cell_iterator cell =
                 tr.begin_active();
               cell != tr.end();
               ++cell)
            if (cell->subdomain_id() == myid)
              {
                if (flags[index])
                  {
                    cell->set_refine_flag();
                  }
                ++index;
              }
        }

      tr.execute_coarsening_and_refinement();
      if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
        deallog << "#local cells:" << tr.n_locally_owned_active_cells()
                << std::endl;

      DoFHandler<dim> dofh(tr);

      static const FE_Q<dim> fe(2);
      dofh.distribute_dofs(fe);

      std::vector<types::global_dof_index> n_locally_owned_dofs_per_processor =
        Utilities::MPI::all_gather(MPI_COMM_WORLD, dofh.n_locally_owned_dofs());
      if (myid == 0)
        {
          deallog << "dofh.n_dofs() " << n_locally_owned_dofs_per_processor
                  << std::endl;
          deallog << "dofh.n_locally_owned_dofs() "
                  << dofh.n_locally_owned_dofs() << std::endl;
        }
    }
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);


  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      initlog();

      deallog.push("2d");
      test<2>();
      deallog.pop();
    }
  else
    test<2>();
}
