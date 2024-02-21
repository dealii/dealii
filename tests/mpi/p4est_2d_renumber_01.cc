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



// create a parallel DoFHandler on a 2d mesh and check componentwise
// renumbering


#include <deal.II/base/tensor.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

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
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

  GridGenerator::hyper_cube(tr);

  const FE_Q<dim> fe_q(1);
  FESystem<dim>   fe(fe_q, 2);
  {
    tr.refine_global(2);
    tr.begin_active()->set_refine_flag();
    tr.execute_coarsening_and_refinement();

    DoFHandler<dim> dofh(tr);
    dofh.distribute_dofs(fe);

    {
      const IndexSet dof_set = DoFTools::extract_locally_active_dofs(dofh);
      if (myid == 0)
        dof_set.print(deallog);
    }

    unsigned int n_dofs = dofh.n_dofs();
    if (myid == 0)
      deallog << "**** n_dofs = " << n_dofs << std::endl;

    DoFRenumbering::component_wise(dofh);
    // check if n_dofs() is still
    // correct. This was a bug at some point.
    Assert(n_dofs == dofh.n_dofs(), ExcInternalError());
    {
      const IndexSet dof_set = DoFTools::extract_locally_active_dofs(dofh);

      if (myid == 0)
        dof_set.print(deallog);

      if (myid == 0)
        {
          std::vector<types::global_dof_index>           local_dof_indices;
          typename DoFHandler<dim>::active_cell_iterator cell,
            endc = dofh.end();

          if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
            for (cell = dofh.begin_active(); cell != endc; ++cell)
              if (!cell->is_artificial() && !cell->is_ghost())
                {
                  local_dof_indices.resize(cell->get_fe().dofs_per_cell);
                  cell->get_dof_indices(local_dof_indices);
                  for (unsigned int i = 0; i < cell->get_fe().dofs_per_cell;
                       ++i)
                    deallog << local_dof_indices[i] << ' ';
                  deallog << std::endl;
                }
        }
    }
  }

  if (0)
    {
      tr.refine_global(1);
      DoFHandler<dim> dofh(tr);
      dofh.distribute_dofs(fe);
    }

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    deallog << "OK" << std::endl;
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

      deallog.push("3d");
      test<3>();
      deallog.pop();
    }
  else
    {
      deallog.push("2d");
      test<2>();
      deallog.pop();

      deallog.push("3d");
      test<3>();
      deallog.pop();
    }
}
