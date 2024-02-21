/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2016 - 2023 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Part of the source code is dual licensed under Apache-2.0 WITH
 * LLVM-exception OR LGPL-2.1-or-later. Detailed license information
 * governing the source code and code contributions can be found in
 * LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
 *
 * ------------------------------------------------------------------------
 */

#include <deal.II/distributed/shared_tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <iostream>

using namespace dealii;

static const unsigned int dim = 2;

int
main(int argc, char **argv)
{
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      {
        parallel::shared::Triangulation<dim> triangulation(
          MPI_COMM_WORLD,
          ::Triangulation<dim>::none,
          false,
          parallel::shared::Triangulation<dim>::partition_metis);
        FE_Q<dim>       fe(1);
        DoFHandler<dim> dof_handler(triangulation);

        GridGenerator::hyper_cube(triangulation, -1, 1);
        triangulation.refine_global(2);
        dof_handler.distribute_dofs(fe);
        IndexSet locally_owned_dofs = dof_handler.locally_owned_dofs();
        deallog << locally_owned_dofs.n_elements() << std::endl;
        dof_handler.clear();
        deallog << "OK" << std::endl;
      }
    }

  catch (const std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
