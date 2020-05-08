// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

// test the PETSc Chebychev solver

#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>

#include "../tests.h"

#include "../testmatrix.h"

int
main(int argc, char **argv)
{
  initlog();

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  {
    const unsigned int size = 32;
    unsigned int       dim  = (size - 1) * (size - 1);

    deallog << "Size " << size << " Unknowns " << dim << std::endl;

    // Make matrix
    FDMatrix                    testproblem(size, size);
    PETScWrappers::SparseMatrix A(dim, dim, 5);
    testproblem.five_point(A);

    IndexSet indices(dim);
    indices.add_range(0, dim);
    PETScWrappers::MPI::Vector f(indices, MPI_COMM_WORLD);
    PETScWrappers::MPI::Vector u(indices, MPI_COMM_WORLD);

    f = 1.;
    u = 0.;

    A.compress(VectorOperation::insert);
    f.compress(VectorOperation::insert);
    u.compress(VectorOperation::insert);

    // Chebychev is a tricky smoother for the kind of FD matrix we use in
    // this test. So, simply test that we're able to reduce the residual to
    // a reasonably small value of 1.e-3.
    SolverControl control(2500, 1.5e-3);

    PETScWrappers::SolverChebychev    solver(control);
    PETScWrappers::PreconditionJacobi preconditioner(A);

    check_solver_within_range(solver.solve(A, u, f, preconditioner),
                              control.last_step(),
                              1050,
                              1150);
  }
}
