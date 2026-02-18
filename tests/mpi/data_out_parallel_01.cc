// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2003 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

// tests DataOut with more MPI ranks than cells / patches and writing
// a TrilinosWrappers::Vector

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/trilinos_parallel_block_vector.h>

#include <deal.II/numerics/data_out.h>

#include "../tests.h"



template <int dim>
void
test()
{
  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria, 0., 1.);

  FESystem<dim>   fe(FE_Q<dim>(2), 1);
  DoFHandler<dim> dof_handler(tria);

  dof_handler.distribute_dofs(fe);

  std::vector<unsigned int> stokes_sub_blocks(1, 0);
  DoFRenumbering::component_wise(dof_handler, stokes_sub_blocks);

  deallog << "Total number of DoFs: " << dof_handler.n_dofs() << std::endl;

  const std::vector<types::global_dof_index> stokes_dofs_per_block =
    DoFTools::count_dofs_per_fe_block(dof_handler, stokes_sub_blocks);

  const types::global_dof_index n_u = stokes_dofs_per_block[0];

  const IndexSet &stokes_locally_owned_index_set =
    dof_handler.locally_owned_dofs();
  const IndexSet stokes_locally_relevant_set =
    DoFTools::extract_locally_relevant_dofs(dof_handler);

  std::vector<IndexSet> stokes_partitioning;
  stokes_partitioning.push_back(
    stokes_locally_owned_index_set.get_view(0, n_u));

  std::vector<IndexSet> stokes_relevant_partitioning;
  stokes_relevant_partitioning.push_back(
    stokes_locally_relevant_set.get_view(0, n_u));

  TrilinosWrappers::MPI::Vector v(stokes_partitioning[0],
                                  stokes_relevant_partitioning[0],
                                  tria.get_mpi_communicator());

  DataOut<dim> data_out;
  data_out.add_data_vector(dof_handler, v, "linear");
  data_out.build_patches();

  std::stringstream output;
  data_out.write(output, DataOutBase::OutputFormat::gnuplot);

  deallog << "Output: " << output.str() << std::endl << std::endl;
  deallog << "Ok" << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  try
    {
      test<2>();

      return 0;
    }
  catch (const std::exception &exc)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    }
  catch (...)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    };
}
