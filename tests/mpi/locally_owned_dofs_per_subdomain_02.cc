// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


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
