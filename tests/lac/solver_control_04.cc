// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test SolverControl's get_history_data function
// This test is adapted from tests/trilinos/solver_03.cc


#include <deal.II/base/utilities.h>

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_relaxation.h>
#include <deal.II/lac/solver_richardson.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_memory.h>

#include <iostream>

#include "../tests.h"

#include "../testmatrix.h"

template <typename MatrixType, typename VectorType, class PRECONDITION>
void
check_solve(SolverControl      &solver_control,
            const MatrixType   &A,
            VectorType         &u,
            VectorType         &f,
            const PRECONDITION &P,
            const bool          expected_result)
{
  SolverCG<VectorType> solver(solver_control);

  u            = 0.;
  f            = 1.;
  bool success = false;
  try
    {
      solver.solve(A, u, f, P);
      deallog << "Success. " << std::endl;
      success = true;
    }
  catch (const std::exception &e)
    {
      deallog << "Failure. " << std::endl;
    }

  deallog << "Solver history data: " << solver_control.get_history_data()
          << std::endl;
  Assert(success == expected_result, ExcMessage("Incorrect result."));
}


int
main(int argc, char **argv)
{
  initlog();
  deallog.get_file_stream() << std::setprecision(4);

  {
    const unsigned int size = 32;
    unsigned int       dim  = (size - 1) * (size - 1);

    deallog << "Size " << size << " Unknowns " << dim << std::endl;

    // Make matrix
    FDMatrix        testproblem(size, size);
    SparsityPattern structure(dim, dim, 5);
    testproblem.five_point_structure(structure);
    structure.compress();
    SparseMatrix<double> A(structure);
    testproblem.five_point(A);

    Vector<double> f(dim);
    Vector<double> u(dim);
    f = 1.;

    PreconditionJacobi<> preconditioner;
    preconditioner.initialize(A);

    deallog.push("Abs tol");
    {
      // Expects success
      SolverControl solver_control(2000, 1.e-3);
      solver_control.enable_history_data();

      check_solve(solver_control, A, u, f, preconditioner, true);
    }
    deallog.pop();
    deallog.push("Iterations");
    {
      // Expects failure
      SolverControl solver_control(20, 1.e-3);
      solver_control.enable_history_data();

      check_solve(solver_control, A, u, f, preconditioner, false);
    }
    deallog.pop();
    deallog.push("Reuse");
    {
      // Expects success
      SolverControl solver_control(200, 1.e-1);
      solver_control.enable_history_data();

      check_solve(solver_control, A, u, f, preconditioner, true);
      check_solve(solver_control, A, u, f, preconditioner, true);
    }
    deallog.pop();
  }
}
