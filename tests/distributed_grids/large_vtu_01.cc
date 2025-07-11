/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2022 - 2023 by the deal.II authors
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

// Test writing of large vtu files with a total file size above 4 GB

// This test is running a tiny version by default. Set the following
// flag to true to run a test that generates a file larger than 4GB
// when running with 2 MPI ranks. The total file size is about
// 7 GB. Warning, you need quite a bit of RAM and patience to run this.
const bool run_big = false;

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

#include "../tests.h"

template <int dim>
void
run()
{
  MPI_Comm                                  mpi_communicator = MPI_COMM_WORLD;
  parallel::distributed::Triangulation<dim> triangulation(mpi_communicator);

  GridGenerator::subdivided_hyper_cube(triangulation, 1);
  unsigned int n_cycles_global = (run_big) ? 4 : 1;
  triangulation.refine_global(n_cycles_global);

  if (dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    deallog << "n_global_active_cells: "
            << triangulation.n_global_active_cells()
            << " n_global_levels: " << triangulation.n_global_levels()
            << std::endl;

  const unsigned int n_vectors = 1;

  FE_DGQ<dim>     fe(2);
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);
  // Make FE vector
  const IndexSet locally_relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(dof_handler);

  using VectorType = LinearAlgebra::distributed::Vector<double>;
  VectorType global_dof_vector;
  global_dof_vector.reinit(dof_handler.locally_owned_dofs(),
                           locally_relevant_dofs,
                           MPI_COMM_WORLD);

  VectorTools::interpolate(dof_handler,
                           Functions::SquareFunction<dim>(),
                           global_dof_vector);

  {
    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    for (unsigned int i = 0; i < n_vectors; ++i)
      data_out.add_data_vector(
        dof_handler,
        global_dof_vector,
        std::vector<std::string>(1, "data" + Utilities::to_string(i)));

    data_out.build_patches((run_big) ? 42 : 2);
    data_out.write_vtu_in_parallel("base.vtu", MPI_COMM_WORLD);
  }

  if (dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    {
      // Let's print 10 lines from vtu (exclude the header and the data):
      std::ifstream in("base.vtu");
      for (int i = 0; i < 10; ++i)
        {
          std::string line;
          std::getline(in, line);
          if (line[0] == '<')
            deallog << line << '\n';
        }

      deallog << "OK" << std::endl;
    }
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  run<3>();
}
