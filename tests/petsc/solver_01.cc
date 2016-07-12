// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2016 by the deal.II authors
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

// test the PETSc Richardson solver

#include "../tests.h"
#include "../lac/testmatrix.h"

#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_precondition.h>

int main(int argc, char **argv)
{
  initlog();

  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, 1);
  {
    const unsigned int size = 32;
    unsigned int dim = (size-1)*(size-1);

    deallog << "Size " << size << " Unknowns " << dim << std::endl;

    // Make matrix
    FDMatrix testproblem(size, size);
    PETScWrappers::SparseMatrix  A(dim, dim, 5);
    testproblem.five_point(A);

    PETScWrappers::Vector  f(dim);
    PETScWrappers::Vector  u(dim);

    f = 1.;
    u = 0.;

    A.compress (VectorOperation::insert);
    f.compress (VectorOperation::insert);
    u.compress (VectorOperation::insert);

    // Richardson is a tricky smoother for the kind of FD matrix we use in
    // this test. So, simply test that we're able to reduce the residual to
    // a reasonably small value of 1.e-4.
    SolverControl control(2500, 1.e-4);

    PETScWrappers::SolverRichardson solver(control);
    PETScWrappers::PreconditionJacobi preconditioner(A);

    check_solver_within_range(
      solver.solve(A,u,f, preconditioner),
      control.last_step(), 2295, 2300);
  }
}
