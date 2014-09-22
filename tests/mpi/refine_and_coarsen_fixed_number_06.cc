// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



// document bug in parallel::distributed::GridRefinement
// ::refine_and_coarsen_fixed_number() and fixed_fraction() with one CPU with
// 0 cells:
//#0  0x00007ffff67297e2 in dealii::(anonymous namespace)::min_element<float> (
//    criteria=...)
//    at /w/heister/deal-trunk/deal.II/source/distributed/grid_refinement.cc:57
//#1  0x00007ffff6727672 in dealii::(anonymous namespace)::compute_global_min_and_max_at_root<float> (criteria=..., mpi_communicator=0x6fc860)
//    at /w/heister/deal-trunk/deal.II/source/distributed/grid_refinement.cc:82
//#2  0x00007ffff672b4ef in dealii::parallel::distributed::GridRefinement::refine_and_coarsen_fixed_number<2, dealii::Vector<float>, 2>

#include "../tests.h"

#include <deal.II/lac/compressed_simple_sparsity_pattern.h>
#include <deal.II/lac/vector.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_generator.h>
#include <fstream>
//#include <mpi.h>

template <int dim>
void test()
{
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  unsigned int numprocs = Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);

  parallel::distributed::Triangulation<dim> triangulation(MPI_COMM_WORLD);
  GridGenerator::hyper_cube (triangulation);

  Vector<float> estimated_error_per_cell (triangulation.n_active_cells());
  parallel::distributed::GridRefinement::
  refine_and_coarsen_fixed_number (triangulation,
                                   estimated_error_per_cell,
                                   0.3, 0.03);
  parallel::distributed::GridRefinement::
  refine_and_coarsen_fixed_fraction (triangulation,
                                     estimated_error_per_cell,
                                     0.3, 0.03);
  triangulation.execute_coarsening_and_refinement ();

  if (myid==0)
    deallog << "n_global_active_cells="
            << triangulation.n_global_active_cells()
            << std::endl;

  if (myid==0)
    deallog << "OK" << std::endl;
}


int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    {
      std::ofstream logfile("output");
      deallog.attach(logfile);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      test<2>();
    }
  else
    test<2>();

}
