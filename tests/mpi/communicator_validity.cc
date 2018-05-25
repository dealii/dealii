// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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

  deallog << v1.l2_norm() << " " << v3->l2_norm() << std::endl;
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
