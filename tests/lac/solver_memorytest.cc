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



// Test all solver with a start value that is already a solution (i.e. 0
// iterations). This caused a memory leak in FGMRES. In the same situation
// test GMRES without using the default residual, which caused a double
// memory freeing (see issue #886).


#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_bicgstab.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/solver_minres.h>
#include <deal.II/lac/solver_qmrs.h>
#include <deal.II/lac/solver_richardson.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_memory.h>

#include "../tests.h"

#include "../testmatrix.h"

template <typename SolverType, typename MatrixType, typename VectorType>
void
check_solve(const MatrixType                          &A,
            VectorType                                &u,
            VectorType                                &f,
            const typename SolverType::AdditionalData &additional_data =
              typename SolverType::AdditionalData())
{
  GrowingVectorMemory<> mem;
  SolverControl         control(100, 1.e-3);
  SolverType            solver(control, mem, additional_data);
  PreconditionIdentity  prec_no;
  u = 0.;
  f = 0.;

  try
    {
      solver.solve(A, u, f, prec_no);
    }
  catch (const std::exception &e)
    {
      deallog << e.what() << std::endl;
    }
}

int
main()
{
  std::ofstream logfile("output");
  //  logfile.setf(std::ios::fixed);
  deallog << std::setprecision(4);
  deallog.attach(logfile);

  for (unsigned int size = 4; size <= 30; size *= 3)
    {
      unsigned int dim = (size - 1) * (size - 1);

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
      Vector<double> res(dim);

      deallog.push("alreadydone");
      check_solve<SolverCG<>>(A, u, f);
      check_solve<SolverGMRES<>>(A, u, f);
      //      check_solve<SolverFGMRES<> >(A,u,f);
      check_solve<SolverBicgstab<>>(A, u, f);
      check_solve<SolverQMRS<>>(A, u, f);

      // test use_default_residual=false case
      check_solve<SolverGMRES<>>(
        A, u, f, SolverGMRES<>::AdditionalData(30, false, false));

      deallog.pop();
    }
}
