// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check that LinearAlgebra::distributed::Vector::operator= does not carry
// over any state from a duplicated MPI communicator in conjunction with a
// distributed mesh


#include <deal.II/base/mpi.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgq.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/vector_memory.h>

#include "../tests.h"

template <typename VectorType>
void
do_test(MPI_Comm communicator)
{
  const int                                 dim = 2;
  parallel::distributed::Triangulation<dim> tria(communicator);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(2);
  FE_DGQ<dim>     fe(0);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);

  VectorType v1;
  v1.reinit(dof.locally_owned_dofs(), communicator);
  if (dof.locally_owned_dofs().n_elements() > 0)
    v1.local_element(0) = 1;

  GrowingVectorMemory<VectorType>            memory;
  typename VectorMemory<VectorType>::Pointer v3(memory);
  *v3 = v1;

  deallog << v1.l2_norm() << ' ' << v3->l2_norm() << std::endl;
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  mpi_initlog();

  do_test<LinearAlgebra::distributed::Vector<double>>(MPI_COMM_WORLD);

  {
    MPI_Comm communicator =
      Utilities::MPI::duplicate_communicator(MPI_COMM_WORLD);
    do_test<LinearAlgebra::distributed::Vector<double>>(communicator);
    MPI_Comm_free(&communicator);
  }
  {
    MPI_Comm communicator =
      Utilities::MPI::duplicate_communicator(MPI_COMM_WORLD);
    do_test<LinearAlgebra::distributed::Vector<double>>(communicator);
    MPI_Comm_free(&communicator);
  }

  do_test<LinearAlgebra::distributed::Vector<double>>(MPI_COMM_WORLD);
}
