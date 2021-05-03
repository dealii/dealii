// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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


// Test that DataOut::write_vtu_with_pvtu_record() prints the path name
// if a file cannot be opened.

#include <deal.II/base/mpi.h>

#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/distributed/shared_tria.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_dgq.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_description.h>

#include <deal.II/numerics/data_out.h>

#include <sstream>

#include "../tests.h"

using namespace dealii;


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  const unsigned int dim = 2;

  const MPI_Comm comm = MPI_COMM_WORLD;

  parallel::distributed::Triangulation<dim> tria(comm);

  GridGenerator::subdivided_hyper_cube(tria, 2);

  FE_DGQ<dim>     fe(2);
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  Vector<double> vec(dof_handler.n_dofs());

  DataOut<dim> data_out;

  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(vec, "solution");

  const auto print_line = [](const std::string &str) {
    std::stringstream ss(str);
    std::string       token;

    for (unsigned int i = 0; i < 8; ++i)
      std::getline(ss, token, '\n');

    deallog << token << std::endl;
  };

  try
    {
      data_out.write_vtu_with_pvtu_record("./folder/", "file", 0, comm, 1, 0);
    }
  catch (const std::exception &exc)
    {
      print_line(exc.what());
    }

  try
    {
      data_out.write_vtu_with_pvtu_record("./folder/", "file", 0, comm, 1, 1);
    }
  catch (const std::exception &exc)
    {
      print_line(exc.what());
    }

  try
    {
      data_out.write_vtu_with_pvtu_record(
        "./folder/", "file", 0, comm, 1, Utilities::MPI::n_mpi_processes(comm));
    }
  catch (const std::exception &exc)
    {
      print_line(exc.what());
    }
}
