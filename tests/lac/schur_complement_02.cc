// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2018 by the deal.II authors
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

// Test internal preconditioner and solver options

#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/linear_operator.h>
#include <deal.II/lac/packaged_operation.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/schur_complement.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include "../tests.h"

#define PRINTME(name, var) deallog << "RHS vector: " << name << ": " << var;



int
main()
{
  initlog();
  deallog.depth_console(0);
  deallog << std::setprecision(10);

  // deal.II SparseMatrix
  {
    deallog << "Schur complement" << std::endl;
    deallog.push("SC_SparseMatrix");

    {
      deallog << "SparseMatrix 1" << std::endl;

      /* MATLAB / Gnu Octave code

        clear all;
        printf("SparseMatrix 1")
        A = [1];
        B = [2];
        C = [3];
        D = [4];
        y = [6];

        S = D - C*inv(A)*B

        % vmult
        g1 = S*y
        % Tvmult
        g2 = S'*y

        g = [2];

        % vmult_add
        g3 = S*y + g
        % Tvmult_add
        g4 = S'*y + g

       */

      const unsigned int rc = 1;
      SparsityPattern    sparsity_pattern(rc, rc, 0);
      sparsity_pattern.compress();

      SparseMatrix<double> A(sparsity_pattern);
      SparseMatrix<double> B(sparsity_pattern);
      SparseMatrix<double> C(sparsity_pattern);
      SparseMatrix<double> D(sparsity_pattern);
      Vector<double>       y(rc);
      Vector<double>       g(rc);
      for (unsigned int i = 0; i < rc; ++i)
        {
          A.diag_element(i) = 1.0 * (i + 1);
          B.diag_element(i) = 2.0 * (i + 1);
          C.diag_element(i) = 3.0 * (i + 1);
          D.diag_element(i) = 4.0 * (i + 1);
          y(i)              = 6.0 * (i + 1);
          g(i)              = 2.0 * (i + 1);
        }

      const auto lo_A = linear_operator(A);
      const auto lo_B = linear_operator(B);
      const auto lo_C = linear_operator(C);
      const auto lo_D = linear_operator(D);

      SolverControl                            solver_control_A(100, 1.0e-10);
      SolverCG<Vector<double>>                 solver_A(solver_control_A);
      PreconditionJacobi<SparseMatrix<double>> preconditioner_A;
      preconditioner_A.initialize(A);
      const auto lo_A_inv = inverse_operator(lo_A, solver_A, preconditioner_A);

      const auto lo_S   = schur_complement(lo_A_inv, lo_B, lo_C, lo_D);
      const auto lo_S_t = transpose_operator(lo_S);

      const Vector<double> g1 = lo_S * y;
      const Vector<double> g2 = lo_S_t * y;
      const Vector<double> g3 = lo_S * y + g;
      const Vector<double> g4 = lo_S_t * y + g;

      PRINTME("g1", g1);
      PRINTME("g2", g2);
      PRINTME("g3", g3);
      PRINTME("g4", g4);
    }

    deallog << "SparseMatrix OK" << std::endl;
  }
}
