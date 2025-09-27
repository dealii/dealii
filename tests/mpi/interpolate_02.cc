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



// Like _01 but with PETSc vectors. this program hangs at the time of writing
// in much the same way solution_transfer_01 does because one processor has no
// cells and all the others decide that PETScWrappers::MPI::Vector::compress
// is a no-op when of course the one who does have cells doesn't think so


#include <deal.II/base/function.h>
#include <deal.II/base/tensor.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/petsc_vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


template <int dim>
void
test()
{
  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  // create a mesh so that all but one
  // processor are empty
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tr);
  FE_Q<dim>       fe(1);
  DoFHandler<dim> dofh(tr);
  dofh.distribute_dofs(fe);

  IndexSet                   owned_set = dofh.locally_owned_dofs();
  PETScWrappers::MPI::Vector x(MPI_COMM_WORLD,
                               dofh.n_dofs(),
                               owned_set.n_elements());

  VectorTools::interpolate(dofh,
                           Functions::ConstantFunction<dim, PetscScalar>(1),
                           x);
  const double norm = x.l2_norm();
  if (myid == 0)
    deallog << dofh.n_locally_owned_dofs() << ' ' << dofh.n_dofs() << std::endl
            << norm << std::endl;
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
