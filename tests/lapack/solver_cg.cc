// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Check eigenvalue capabilities of SolverCG

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_matrix.templates.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_memory.h>

#include <iostream>

#include "../tests.h"

#include "../testmatrix.h"

void
output_double_number(double input, const std::string &text)
{
  deallog << text << input << std::endl;
}

template <class NUMBER>
void
output_eigenvalues(const std::vector<NUMBER> &eigenvalues,
                   const std::string         &text)
{
  deallog << text;
  for (unsigned int j = 0; j < eigenvalues.size(); ++j)
    {
      deallog << ' ' << eigenvalues.at(j);
    }
  deallog << std::endl;
}


template <typename SolverType,
          typename MatrixType,
          typename VectorType,
          class PRECONDITION>
void
check_solve(SolverType         &solver,
            const MatrixType   &A,
            VectorType         &u,
            VectorType         &f,
            const PRECONDITION &P)
{
  u = 0.;
  f = 1.;
  try
    {
      solver.solve(A, u, f, P);
    }
  catch (dealii::SolverControl::NoConvergence &e)
    {
      deallog << "Failure step " << e.last_step << " value " << e.last_residual
              << std::endl;
      deallog << "Exception: " << e.get_exc_name() << std::endl;
    }
}

template <typename SolverType,
          typename MatrixType,
          typename VectorType,
          class PRECONDITION>
void
check_Tsolve(SolverType         &solver,
             const MatrixType   &A,
             VectorType         &u,
             VectorType         &f,
             const PRECONDITION &P)
{
  u = 0.;
  f = 1.;
  try
    {
      solver.Tsolve(A, u, f, P);
    }
  catch (dealii::SolverControl::NoConvergence &e)
    {
      deallog << "Exception: " << e.get_exc_name() << std::endl;
    }
}

int
main()
{
  initlog();
  deallog << std::setprecision(4);

  GrowingVectorMemory<> mem;
  SolverControl         control(100, 1.e-3);
  SolverControl         verbose_control(100, 1.e-3, true);
  SolverCG<>            cg(control, mem);
  cg.connect_condition_number_slot(std::bind(output_double_number,
                                             std::placeholders::_1,
                                             "Condition number estimate: "),
                                   true);
  cg.connect_eigenvalues_slot(std::bind(output_eigenvalues<double>,
                                        std::placeholders::_1,
                                        "Final Eigenvalues: "));


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

      PreconditionIdentity   prec_no;
      PreconditionRichardson prec_richardson;
      prec_richardson.initialize(0.6);
      PreconditionSSOR<> prec_ssor;
      prec_ssor.initialize(A, 1.2);

      std::vector<unsigned int> permutation(dim);
      std::vector<unsigned int> inverse_permutation(dim);

      Vector<double> f(dim);
      Vector<double> u(dim);
      Vector<double> res(dim);

      f = 1.;
      u = 1.;

      try
        {
          deallog.push("no-fail");
          control.set_max_steps(10);
          check_solve(cg, A, u, f, prec_no);
          control.set_max_steps(100);
          deallog.pop();

          deallog.push("no");
          check_solve(cg, A, u, f, prec_no);
          deallog.pop();

          deallog.push("rich");
          check_solve(cg, A, u, f, prec_richardson);
          deallog.pop();

          deallog.push("ssor");
          check_solve(cg, A, u, f, prec_ssor);
          deallog.pop();
        }
      catch (const std::exception &e)
        {
          std::cerr << "Exception: " << e.what() << std::endl;
        }
    }
}
