/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2013 - 2015 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------
 */

#include <deal.II/base/logstream.h>
#include <deal.II/distributed/shared_tria.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>

#include <fstream>
#include <iostream>

using namespace dealii;

static const unsigned int dim = 2;

int main (int argc, char **argv)
{
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, 1);
      {
        parallel::shared::Triangulation<dim> triangulation(MPI_COMM_WORLD);
        FE_Q<dim> fe(1);
        DoFHandler<dim> dof_handler (triangulation);

        GridGenerator::hyper_cube (triangulation, -1, 1);
        triangulation.refine_global (2);
        dof_handler.distribute_dofs (fe);
        IndexSet locally_owned_dofs = dof_handler.locally_owned_dofs();
        deallog << locally_owned_dofs.n_elements() << std::endl;
        dof_handler.clear ();
        deallog << "OK" << std::endl;
      }
    }

  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
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
      std::cerr << std::endl << std::endl
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
