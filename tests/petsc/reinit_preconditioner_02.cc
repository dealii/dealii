// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2017 by the deal.II authors
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



// check re-initializing a preconditioner (serial version)

#include "../tests.h"
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/base/index_set.h>
#include <iostream>
#include <vector>

template <class PRE>
void
test ()
{
  DynamicSparsityPattern csp (5, 5);

  for (unsigned int i=0; i<5; ++i)
    csp.add(i,i);

  csp.add(0,1);
  csp.add(1,0);

  PETScWrappers::SparseMatrix mat;
  mat.reinit (csp);

  for (unsigned int i=0; i<5; ++i)
    mat.set(i,i, 1.0+i*2.0);
  mat.set(0,1, 0.1);
  mat.set(1,0, 0.1);

  mat.compress(VectorOperation::insert);

  {
    IndexSet indices(5);
    indices.add_range(0, 5);
    PETScWrappers::MPI::Vector src, dst;
    src.reinit(indices, MPI_COMM_WORLD);
    dst.reinit(indices, MPI_COMM_WORLD);
    src(0) = 1.0;
    src(1) = 2.0;
    src.compress(VectorOperation::insert);

    PRE pre;
    pre.initialize(mat);
    pre.vmult(dst, src);
    dst.print(deallog.get_file_stream());

    mat.add(0,0,1.0);
    mat.compress(VectorOperation::add);

    pre.initialize(mat);
    pre.vmult(dst, src);
    dst.print(deallog.get_file_stream());
  }

  deallog << "OK" << std::endl;
}



int
main (int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, 1);
  MPILogInitAll log;

  test<PETScWrappers::PreconditionJacobi> ();
  test<PETScWrappers::PreconditionBlockJacobi> ();
  test<PETScWrappers::PreconditionSOR> ();
  test<PETScWrappers::PreconditionSSOR> ();
  // todo: this crashes test<PETScWrappers::PreconditionEisenstat> ();
  test<PETScWrappers::PreconditionICC> ();
  test<PETScWrappers::PreconditionILU> ();
  test<PETScWrappers::PreconditionLU> ();
  test<PETScWrappers::PreconditionBoomerAMG> ();
  test<PETScWrappers::PreconditionParaSails> ();
  test<PETScWrappers::PreconditionNone> ();
}
