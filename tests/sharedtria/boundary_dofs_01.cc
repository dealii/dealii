// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Check DoFTools::extract_boundary_dofs() with shared triangulations.

#include <deal.II/distributed/shared_tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include "../tests.h"



template <int dim>
void
test()
{
  const MPI_Comm                       mpi_communicator = MPI_COMM_WORLD;
  parallel::shared::Triangulation<dim> triangulation(mpi_communicator);

  GridGenerator::hyper_cube(triangulation, 0, 1);
  triangulation.refine_global(4);

  FE_Q<dim>       fe(1);
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);

  deallog << "On process " << Utilities::MPI::this_mpi_process(mpi_communicator)
          << std::endl;
  deallog << "  centers of locally owned cells=";
  for (const auto &cell : triangulation.active_cell_iterators())
    if (cell->is_locally_owned())
      deallog << cell->center() << ' ';
  deallog << "\n  locally_owned_dofs=";
  for (const auto i : dof_handler.locally_owned_dofs())
    deallog << i << ' ';
  deallog << "\n  boundary_dofs=";
  for (const auto i : DoFTools::extract_boundary_dofs(dof_handler))
    deallog << i << ' ';
  deallog << std::endl;
  deallog << std::endl;
}


int
main(int argc, char **argv)
{
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

      MPILogInitAll all;

      test<1>();
    }
  catch (std::exception &exc)
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
