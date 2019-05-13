// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2018 by the deal.II authors
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


// Test DoFTools::locally_owned_dofs_per_subdomain (for a standard
// Triangulation)

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


template <int dim>
void
test()
{
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation, -1.0, 1.0);
  triangulation.refine_global(1);

  const unsigned int nproc = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  GridTools::partition_triangulation(nproc, triangulation);

  FE_Q<dim>       fe(1);
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);

  std::vector<IndexSet> locally_owned_dofs_per_proc =
    DoFTools::locally_owned_dofs_per_subdomain(dof_handler);

  for (unsigned int p = 0; p < nproc; ++p)
    {
      deallog << "proc " << p << ": " << std::endl;
      locally_owned_dofs_per_proc[p].print(deallog.get_file_stream());
    }
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  MPILogInitAll all;

  test<2>();
}
