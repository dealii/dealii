// ---------------------------------------------------------------------
//
// Copyright (C) 2015 by the deal.II authors
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

// Test that it is possible to instantiate a LinearOperator object for all
// different kinds of Trilinos matrices and vectors
// TODO: A bit more tests...

#include "../tests.h"


#include <deal.II/lac/linear_operator.h>

#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_sparse_matrix.h>

#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_parallel_sparse_matrix.h>

using namespace dealii;

int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, 1);

  initlog();
  deallog << std::setprecision(10);

  {
    PETScWrappers::SparseMatrix a;
    auto op_a  = linear_operator<PETScWrappers::Vector>(a);
  }

  {
    PETScWrappers::MPI::SparseMatrix a;
    auto op_a  = linear_operator<PETScWrappers::MPI::Vector>(a);
  }


  deallog << "OK" << std::endl;

  return 0;
}


