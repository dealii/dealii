// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2013 by the deal.II authors
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


// check FETools::interpolate on parallel vector

#include "../tests.h"
#include <deal.II/base/utilities.h>
#include <deal.II/base/index_set.h>
#include <deal.II/lac/parallel_vector.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/base/function.h>
#include <fstream>
#include <iostream>
#include <vector>


void test ()
{
  const unsigned int dim = 2;

  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(3);

  FE_Q<dim> fe1(1), fe2(2);
  DoFHandler<dim> dof1(tria), dof2(tria);
  dof1.distribute_dofs(fe1);
  dof2.distribute_dofs(fe2);

  IndexSet locally_relevant_dofs1, locally_relevant_dofs2;
  DoFTools::extract_locally_relevant_dofs (dof1,
                                           locally_relevant_dofs1);
  DoFTools::extract_locally_relevant_dofs (dof2,
                                           locally_relevant_dofs2);

  parallel::distributed::Vector<double> v1(dof1.locally_owned_dofs(),
                                           locally_relevant_dofs1,
                                           MPI_COMM_WORLD),
                                           v2(dof2.locally_owned_dofs(), locally_relevant_dofs2, MPI_COMM_WORLD);

  // set first vector to 1
  VectorTools::interpolate(dof1, ConstantFunction<dim>(1.), v1);
  for (unsigned int i=0; i<v1.local_size(); ++i)
    Assert(v1.local_element(i) == 1., ExcInternalError());

  v1.update_ghost_values();
  FETools::interpolate(dof1, v1, dof2, v2);
  for (unsigned int i=0; i<v2.local_size(); ++i)
    Assert(v2.local_element(i) == 1., ExcInternalError());

  v2.update_ghost_values();
  for (unsigned int i=0; i<v2.local_size()+v2.n_ghost_entries(); ++i)
    Assert(v2.local_element(i) == 1., ExcInternalError());
}



int main (int argc, char **argv)
{
  Utilities::System::MPI_InitFinalize mpi_initialization(argc, argv);

  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      std::ofstream logfile("output");
      deallog.attach(logfile);
      deallog << std::setprecision(4);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      test();
    }
  else
    test();

}
