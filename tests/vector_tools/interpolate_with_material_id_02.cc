// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Tests if proper checks are run when interpolating on distributed vectors

#include <deal.II/base/function.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/mpi.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_cache.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


template <int dim>
void
test_interpolating_distributed(unsigned int ref_cube)
{
  MPI_Comm     mpi_communicator = MPI_COMM_WORLD;
  unsigned int n_procs = Utilities::MPI::n_mpi_processes(mpi_communicator);

  FE_Q<dim> fe(1);
  // Implementing the constant function
  Functions::ConstantFunction<dim> f_constant(1.0);

  // Creating a distributed triangulation

  parallel::distributed::Triangulation<dim> cube_d(mpi_communicator);
  GridGenerator::hyper_cube(cube_d);
  cube_d.refine_global(ref_cube);
  DoFHandler<dim> dof_handler_d(cube_d);
  dof_handler_d.distribute_dofs(fe);

  LinearAlgebra::distributed::Vector<double> data_d(dof_handler_d.n_dofs());
  VectorTools::interpolate(dof_handler_d, f_constant, data_d);
  data_d.print(deallog.get_file_stream());
  deallog << "Test passed" << std::endl;
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;
  test_interpolating_distributed<2>(3);
}
