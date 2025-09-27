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



// check AffineConstraints<double>.distribute() for a distributed mesh
// with Trilinos

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

#include <deal.II/lac/trilinos_vector.h>

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
  for (unsigned int step = 0; step < 15; ++step)
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

  const IndexSet &owned_set    = dofh.locally_owned_dofs();
  const IndexSet  dof_set      = DoFTools::extract_locally_active_dofs(dofh);
  const IndexSet  relevant_set = DoFTools::extract_locally_relevant_dofs(dofh);

  TrilinosWrappers::MPI::Vector x;
  x.reinit(owned_set, MPI_COMM_WORLD);
  x = 2.0;

  TrilinosWrappers::MPI::Vector x_rel;
  x_rel.reinit(relevant_set, MPI_COMM_WORLD);

  AffineConstraints<double> cm(owned_set, relevant_set);
  DoFTools::make_hanging_node_constraints(dofh, cm);
  cm.close();

  cm.distribute(x);
  x_rel = x;

  TrilinosWrappers::MPI::Vector x_dub;
  x_dub.reinit(complete_index_set(dof_set.size()));
  x_dub.reinit(x_rel, false, true);

  {
    std::stringstream out;
    out << "**** proc " << myid << std::endl;
    x_dub.print(out);

    if (myid == 0)
      deallog << out.str() << std::endl;
    else
      MPI_Send((void *)out.str().c_str(),
               out.str().size() + 1,
               MPI_CHAR,
               0,
               1,
               MPI_COMM_WORLD);
  }

  if (myid == 0)
    {
      for (unsigned int i = 1; i < numproc; ++i)
        {
          MPI_Status status;
          int        msglen;
          MPI_Probe(i, 1, MPI_COMM_WORLD, &status);
          MPI_Get_count(&status, MPI_CHAR, &msglen);
          std::vector<char> buf(msglen);
          MPI_Recv(&buf[0],
                   msglen,
                   MPI_CHAR,
                   status.MPI_SOURCE,
                   1,
                   MPI_COMM_WORLD,
                   &status);
          deallog << &buf[0] << std::endl;
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

      deallog.push("2d");
      test<2>();
      deallog.pop();
    }
  else
    test<2>();
}
