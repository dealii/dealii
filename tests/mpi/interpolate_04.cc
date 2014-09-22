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



// test FETools::interpolation_difference

#include "../tests.h"
#include "coarse_grid_common.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/tensor.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/lac/petsc_parallel_vector.h>

#include <fstream>




template<int dim>
void test()
{
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

  GridGenerator::hyper_cube(tr);

//tr.refine_global (2);

  const FE_Q<dim> fe(2);
  DoFHandler<dim> dofh1(tr);
  DoFHandler<dim> dofh2(tr);
  dofh1.distribute_dofs (fe);
  dofh2.distribute_dofs (fe);

  ConstraintMatrix cm1;
  cm1.close();
  ConstraintMatrix cm2;
  cm2.close();

  IndexSet  dof1_locally_owned_dofs = dofh1.locally_owned_dofs();
  IndexSet  dof2_locally_owned_dofs = dofh2.locally_owned_dofs();
  IndexSet  dof1_locally_relevant_dofs;
  IndexSet  dof2_locally_relevant_dofs;
  DoFTools::extract_locally_relevant_dofs (dofh2,
                                           dof2_locally_relevant_dofs);
  DoFTools::extract_locally_relevant_dofs (dofh1,
                                           dof1_locally_relevant_dofs);

  PETScWrappers::MPI::Vector u1(dof1_locally_owned_dofs, dof1_locally_relevant_dofs, MPI_COMM_WORLD);

  PETScWrappers::MPI::Vector out(dof1_locally_owned_dofs, MPI_COMM_WORLD);

  FETools::interpolation_difference
  (dofh1, cm1, u1, dofh2, cm2, out);

  double norm = out.l2_norm();

  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    deallog << "norm = " << norm
            << std::endl;
}


int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);

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

      deallog.push("3d");
      //    test<3>();
      deallog.pop();
    }
  else
    {
      deallog.push("2d");
      test<2>();
      deallog.pop();

      deallog.push("3d");
//      test<3>();
      deallog.pop();
    }
}
