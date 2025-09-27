// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check AffineConstraints<double> for a distributed mesh,
// also compare with/without sparse line_cache via IndexSet.
// Mesh: 3d Random refinement.

#include <deal.II/base/tensor.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

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
  unsigned int myid    = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    deallog << "hyper_cube" << std::endl;

  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

  GridGenerator::hyper_cube(tr);
  tr.refine_global(2);
  for (unsigned int step = 0; step < 8; ++step)
    {
      typename Triangulation<dim>::active_cell_iterator cell =
                                                          tr.begin_active(),
                                                        endc = tr.end();

      for (; cell != endc; ++cell)
        if (Testing::rand() % 42 == 1)
          cell->set_refine_flag();

      tr.execute_coarsening_and_refinement();
    }

  DoFHandler<dim> dofh(tr);

  static const FE_Q<dim> fe(1);
  dofh.distribute_dofs(fe);

  const IndexSet dof_set = DoFTools::extract_locally_relevant_dofs(dofh);

  AffineConstraints<double> cm;
  DoFTools::make_hanging_node_constraints(dofh, cm);
  AffineConstraints<double> cm2(dofh.locally_owned_dofs(), dof_set);
  DoFTools::make_hanging_node_constraints(dofh, cm2);

  {
    std::ofstream file(
      (std::string("dat.") + Utilities::int_to_string(myid)).c_str());
    file << "**** proc " << myid << std::endl;
    cm.print(file);
    file << "****" << std::endl;
    cm2.print(file);
  }

  MPI_Barrier(MPI_COMM_WORLD);



  if (myid == 0)
    {
      for (unsigned int i = 0; i < numproc; ++i)
        {
          cat_file((std::string("dat.") + Utilities::int_to_string(i)).c_str());
        }
    }
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);


  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      initlog();

      deallog.push("3d");
      test<3>();
      deallog.pop();
    }
  else
    test<3>();
}
