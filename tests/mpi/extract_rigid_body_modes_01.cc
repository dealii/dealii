// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// test DoFTools::rigid_body_modes

#include <deal.II/base/tensor.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"



template <int dim>
void
test()
{
  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);
  GridGenerator::subdivided_hyper_cube(tr, 2);

  MappingQ1<dim> mapping;

  FESystem<dim>   fe(FE_Q<dim>(1), dim);
  DoFHandler<dim> dofh(tr);
  dofh.distribute_dofs(fe);

  deallog << "Total dofs=" << dofh.n_dofs() << std::endl;

  // extract constant modes and print

  std::vector<std::vector<double>> constant_modes =
    DoFTools::extract_rigid_body_modes(mapping, dofh);

  for (unsigned int i = 0; i < constant_modes.size(); ++i)
    {
      for (unsigned int j = 0; j < constant_modes[i].size(); ++j)
        deallog << constant_modes[i][j] << ' ';
      deallog << std::endl;
    }
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  deallog.push("2d");
  test<2>();
  deallog.pop();

  deallog.push("3d");
  test<3>();
  deallog.pop();
}
