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



// Test DoFTools::extract_locally_active_dofs and ensure that it returns the
// same result as DoFTools::extract_dofs_with_subdomain_association()

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
  parallel::distributed::Triangulation<dim> tr(
    MPI_COMM_WORLD,
    Triangulation<dim>::limit_level_difference_at_vertices,
    parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy);

  GridGenerator::hyper_ball(tr);
  tr.refine_global(2);

  const FE_Q<dim> fe(2);
  DoFHandler<dim> dofh(tr);
  dofh.distribute_dofs(fe);
  dofh.distribute_mg_dofs();

  for (unsigned int level = 0; level < tr.n_levels(); ++level)
    {
      const IndexSet locally_level_active =
        DoFTools::extract_locally_active_level_dofs(dofh, level);

      if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
        {
          deallog << "Level: " << level << ", " << locally_level_active.size()
                  << ", " << locally_level_active.n_elements() << std::endl;

          for (unsigned int i = 0; i < locally_level_active.size(); ++i)
            if (locally_level_active.is_element(i))
              deallog << i << ' ';
          deallog << "OK" << std::endl;
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
