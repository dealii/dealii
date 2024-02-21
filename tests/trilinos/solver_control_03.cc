// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2004 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test IterationNumberControl with Trilinos solver
// This test is adapted from tests/trilinos/solver_03.cc


#include <deal.II/base/utilities.h>

#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/vector_memory.h>

#include <iostream>
#include <typeinfo>

#include "../tests.h"

#include "../testmatrix.h"


int
main(int argc, char **argv)
{
  initlog();
  deallog << std::setprecision(4);

  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());


  {
    const unsigned int size = 32;
    unsigned int       dim  = (size - 1) * (size - 1);

    deallog << "Size " << size << " Unknowns " << dim << std::endl;

    // Make matrix
    FDMatrix               testproblem(size, size);
    DynamicSparsityPattern csp(dim, dim);
    testproblem.five_point_structure(csp);
    TrilinosWrappers::SparseMatrix A;
    A.reinit(csp);
    testproblem.five_point(A);

    TrilinosWrappers::MPI::Vector f;
    f.reinit(complete_index_set(dim), MPI_COMM_WORLD);
    TrilinosWrappers::MPI::Vector u;
    u.reinit(complete_index_set(dim), MPI_COMM_WORLD);
    f = 1.;
    A.compress(VectorOperation::insert);
    f.compress(VectorOperation::insert);
    u.compress(VectorOperation::insert);

    TrilinosWrappers::PreconditionJacobi preconditioner;
    preconditioner.initialize(A);

    deallog.push("Abs tol");
    {
      // Expects success
      IterationNumberControl     control(2000, 1.e-3);
      TrilinosWrappers::SolverCG solver(control);
      check_solver_within_range(solver.solve(A, u, f, preconditioner),
                                control.last_step(),
                                42,
                                44);
    }
    deallog.pop();

    u = 0.;
    deallog.push("Iterations");
    {
      // Expects success
      IterationNumberControl     control(20, 1.e-3);
      TrilinosWrappers::SolverCG solver(control);
      check_solver_within_range(solver.solve(A, u, f, preconditioner),
                                control.last_step(),
                                19,
                                21);
    }
    deallog.pop();
  }
}
