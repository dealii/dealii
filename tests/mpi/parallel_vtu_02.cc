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

// Test DataOutInterface::write_vtu_in_parallel(): Check that writing a
// vtu file with and without MPI IO results in identical files.

#include <deal.II/base/function_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <sstream>

#include "../tests.h"


std::string
read_file(const std::string &filename)
{
  std::ifstream     in(filename);
  std::stringstream ss;
  ss << in.rdbuf();
  return ss.str();
}


template <int dim>
void
test()
{
  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(2);

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(FE_Q<dim>(1));

  LinearAlgebra::distributed::Vector<double> x(
    dof_handler.locally_owned_dofs(),
    DoFTools::extract_locally_relevant_dofs(dof_handler),
    MPI_COMM_WORLD);
  VectorTools::interpolate(dof_handler, Functions::SquareFunction<dim>(), x);
  x.update_ghost_values();

  DataOutBase::VtkFlags vtk_flags;
  vtk_flags.compression_level = DataOutBase::CompressionLevel::best_compression;

  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(x, "x");
  data_out.build_patches();
  data_out.set_flags(vtk_flags);

  data_out.write_vtu_in_parallel("output_mpi_io.vtu",
                                 MPI_COMM_WORLD,
                                 /*use_mpi_io =*/true);
  data_out.write_vtu_in_parallel("output.vtu",
                                 MPI_COMM_WORLD,
                                 /*use_mpi_io =*/false);
  MPI_Barrier(MPI_COMM_WORLD);

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    {
      deallog << "files identical: " << std::boolalpha
              << (read_file("output.vtu") == read_file("output_mpi_io.vtu"))
              << std::endl;
      cat_file("output.vtu");
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
