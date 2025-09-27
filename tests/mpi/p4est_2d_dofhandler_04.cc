// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// create a parallel DoFHandler on a 2d adaptively refined mesh
// check DoFTools::extract_locally_active_dofs

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

      const IndexSet dof_set = DoFTools::extract_locally_active_dofs<>(dofh);
      std::set<unsigned int> control_dof_set;

      typename DoFHandler<dim>::active_cell_iterator cell = dofh.begin_active();

      const unsigned int dofs_per_cell = dofh.get_fe().dofs_per_cell;
      std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

      for (; cell != dofh.end(); ++cell)
        if (!cell->is_artificial())
          {
            cell->get_dof_indices(local_dof_indices);
            for (std::vector<types::global_dof_index>::iterator it =
                   local_dof_indices.begin();
                 it != local_dof_indices.end();
                 ++it)
              {
                types::global_dof_index invalid_dofindex =
                  numbers::invalid_dof_index;
                Assert(*it != invalid_dofindex, ExcInternalError());
                control_dof_set.insert(*it);
              }
          }

      if (myid == 0)
        {
          deallog << "# in locally_active_dofs:" << dof_set.n_elements()
                  << std::endl;
          dof_set.print(deallog);
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
