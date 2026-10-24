// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

// Test DataOutInterface::write_vtu_with_pvtu_record() with n_groups = 2
// and ParallelWriteMode::send_to_rank_0. The mesh consists of a single
// cell, so all but one rank have no patches and one of the two groups has
// no patches at all. Also check that all ranks throw if the file cannot be
// opened.

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/numerics/data_out.h>

#include "../tests.h"


template <int dim>
void
test()
{
  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria);

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(FE_Q<dim>(2));

  LinearAlgebra::distributed::Vector<double> x(
    dof_handler.locally_owned_dofs(),
    DoFTools::extract_locally_relevant_dofs(dof_handler),
    MPI_COMM_WORLD);
  x = 2.0;
  x.update_ghost_values();

  DataOutBase::VtkFlags vtk_flags;
  vtk_flags.compression_level = DataOutBase::CompressionLevel::best_compression;

  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(x, "x");
  data_out.build_patches();
  data_out.set_flags(vtk_flags);

  deallog << "n_patches: " << data_out.get_patches().size() << std::endl;

  const std::string pvtu_filename = data_out.write_vtu_with_pvtu_record(
    "",
    "output",
    0,
    MPI_COMM_WORLD,
    2,
    /*n_groups =*/2,
    DataOutBase::ParallelWriteMode::send_to_rank_0);
  MPI_Barrier(MPI_COMM_WORLD);

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    {
      cat_file("output_00.0.vtu");
      cat_file("output_00.1.vtu");
      cat_file(pvtu_filename.c_str());
    }

  try
    {
      data_out.write_vtu_in_parallel(
        "nonexistent_directory/output.vtu",
        MPI_COMM_WORLD,
        DataOutBase::ParallelWriteMode::send_to_rank_0);
    }
  catch (const ExcFileNotOpen &)
    {
      deallog << "caught ExcFileNotOpen" << std::endl;
    }

  deallog << "OK" << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  test<2>();
}
