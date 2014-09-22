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



// check ConstraintMatrix.distribute() for a distributed mesh
// with Trilinos

#include "../tests.h"
#include "coarse_grid_common.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor.h>
#include <deal.II/grid/tria.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/lac/trilinos_vector.h>

#include <deal.II/fe/fe_q.h>

#include <fstream>
#include <sstream>

template<int dim>
void test()
{
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);

  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    deallog << "hyper_cube" << std::endl;

  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

  GridGenerator::hyper_cube(tr);
  tr.refine_global (2);
  for (unsigned int step=0; step<15; ++step)
    {
      typename Triangulation<dim>::active_cell_iterator
      cell = tr.begin_active(),
      endc = tr.end();

      for (; cell!=endc; ++cell)
        if (Testing::rand()%42==1)
          cell->set_refine_flag ();

      tr.execute_coarsening_and_refinement ();
    }

  DoFHandler<dim> dofh(tr);

  static const FE_Q<dim> fe(1);
  dofh.distribute_dofs (fe);

  IndexSet owned_set = dofh.locally_owned_dofs();

  IndexSet dof_set;
  DoFTools::extract_locally_active_dofs (dofh, dof_set);

  IndexSet relevant_set;
  DoFTools::extract_locally_relevant_dofs (dofh, relevant_set);

  TrilinosWrappers::MPI::Vector x;
  x.reinit(owned_set, MPI_COMM_WORLD);
  x=2.0;

  TrilinosWrappers::MPI::Vector x_rel;
  x_rel.reinit(relevant_set, MPI_COMM_WORLD);

  ConstraintMatrix cm(relevant_set);
  DoFTools::make_hanging_node_constraints (dofh, cm);
  cm.close ();

  cm.distribute(x);
  x_rel = x;

  //x.print(std::cout);

//  x_rel.print(std::cout);

  TrilinosWrappers::Vector x_dub;
  x_dub.reinit(dof_set.size());

  x_dub = x_rel;

  {
    std::stringstream out;
    out << "**** proc " << myid << std::endl;
    x_dub.print (out);

    if (myid==0)
      deallog << out.str() << std::endl;
    else
      MPI_Send((void *)out.str().c_str(),out.str().size()+1, MPI_CHAR, 0, 1, MPI_COMM_WORLD);
  }

  if (myid == 0)
    {
      for (unsigned int i = 1; i < numproc; ++i)
        {
          MPI_Status status;
          int msglen;
          MPI_Probe(i, 1, MPI_COMM_WORLD, &status);
          MPI_Get_count(&status, MPI_CHAR, &msglen);
          std::vector<char> buf(msglen);
          MPI_Recv(&buf[0], msglen, MPI_CHAR, status.MPI_SOURCE, 1,
                   MPI_COMM_WORLD, &status);
          deallog << &buf[0] << std::endl;
        }
    }

}


int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);


  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      std::ofstream logfile("output");
      deallog.attach(logfile);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      deallog.push("2d");
      test<2>();
      deallog.pop();
    }
  else
    test<2>();

}
