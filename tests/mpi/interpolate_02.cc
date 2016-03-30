// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2015 by the deal.II authors
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



// Like _01 but with PETSc vectors. this program hangs at the time of writing
// in much the same way solution_transfer_01 does because one processor has no
// cells and all the others decide that PETScWrappers::Vector::compress is a
// no-op when of course the one who does have cells doesn't think so


#include "../tests.h"
#include "coarse_grid_common.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/function.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/grid/filtered_iterator.h>

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
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/lac/petsc_parallel_vector.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <fstream>
#include <sstream>


template<int dim>
void test()
{
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);

  // create a mesh so that all but one
  // processor are empty
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);
  GridGenerator::hyper_cube (tr);
  FE_Q<dim> fe (1);
  DoFHandler<dim> dofh(tr);
  dofh.distribute_dofs (fe);

  IndexSet owned_set = dofh.locally_owned_dofs();
  PETScWrappers::MPI::Vector x (MPI_COMM_WORLD,dofh.n_dofs(),owned_set.n_elements());

  VectorTools::interpolate(dofh,
                           ConstantFunction<dim,PetscScalar>(1),
                           x);
  const double norm = x.l2_norm();
  if (myid == 0)
    deallog << dofh.n_locally_owned_dofs() << ' ' << dofh.n_dofs()
            << std::endl
            << norm
            << std::endl;
}


int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, 1);

  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      std::ofstream logfile("output");
      deallog.attach(logfile);
      deallog.threshold_double(1.e-10);

      deallog.push("2d");
      test<2>();
      deallog.pop();
      deallog.push("3d");
      test<3>();
      deallog.pop();
    }
  else
    {
      test<2>();
      test<3>();
    }


}
