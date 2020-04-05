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

#define PRINTME(name, var) \
  deallog << "Solution vector: " << name << ": " << var;



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
        A = [1,2;3,4]
        b = [5;6]
        x = A\b

       */

      const unsigned int rc = 1;
      SparsityPattern    sparsity_pattern(rc, rc, 0);
      sparsity_pattern.compress();

      SparseMatrix<double> A(sparsity_pattern);
      SparseMatrix<double> B(sparsity_pattern);
      SparseMatrix<double> C(sparsity_pattern);
      SparseMatrix<double> D(sparsity_pattern);
      Vector<double>       x(rc);
      Vector<double>       y(rc);
      Vector<double>       f(rc);
      Vector<double>       g(rc);
      for (unsigned int i = 0; i < rc; ++i)
        {
          A.diag_element(i) = 1.0 * (i + 1);
          B.diag_element(i) = 2.0 * (i + 1);
          C.diag_element(i) = 3.0 * (i + 1);
          D.diag_element(i) = 4.0 * (i + 1);
          f(i)              = 5.0 * (i + 1);
          g(i)              = 6.0 * (i + 1);
        }

      const auto lo_A = linear_operator(A);
      const auto lo_B = linear_operator(B);
      const auto lo_C = linear_operator(C);
      const auto lo_D = linear_operator(D);

      SolverControl            solver_control_A(1, 1.0e-10, false, false);
      SolverCG<Vector<double>> solver_A(solver_control_A);
      PreconditionJacobi<SparseMatrix<double>> preconditioner_A;
      preconditioner_A.initialize(A);
      const auto lo_A_inv = inverse_operator(lo_A, solver_A, preconditioner_A);

      const auto lo_S = schur_complement(lo_A_inv, lo_B, lo_C, lo_D);

      SolverControl            solver_control_S(1, 1.0e-10, false, false);
      SolverCG<Vector<double>> solver_S(solver_control_S);
      PreconditionJacobi<SparseMatrix<double>> preconditioner_S;
      preconditioner_S.initialize(D); // Same space as S
      const auto lo_S_inv = inverse_operator(lo_S, solver_S, preconditioner_S);

      auto rhs = condense_schur_rhs(lo_A_inv, lo_C, f, g);
      check_solver_within_range(y = lo_S_inv * rhs, // Solve for y
                                solver_control_S.last_step(),
                                1,
                                1);

      x = postprocess_schur_solution(lo_A_inv, lo_B, y, f);

      PRINTME("x", x);
      PRINTME("y", y);
    }

    deallog << "SparseMatrix OK" << std::endl;
  }

  // deal.II BlockSparseMatrix
  {
    deallog.push("SC_BlockSparseMatrix");

    {
      deallog << "BlockSparseMatrix 1" << std::endl;

      /* MATLAB / Gnu Octave code

         clear all;
         printf("BlockSparseMatrix 1")
         blks=2;
         rc=10;
         for (i=0:rc-1)
           for (bi=0:blks-1)
             b(bi*rc+i+1,1) = bi*rc + i;
             for (bj=0:blks-1)
               el_i = 1 + i + bi*rc;
               el_j = 1 + i + bj*rc;
               A(el_i,el_j) = 2.0*bi + 1.5*bj + (i+1);
             endfor
           endfor
         endfor
         A
         b
         x = A\b

       */

      const unsigned int   blks = 2;
      const unsigned int   rc   = 10;
      BlockSparsityPattern sparsity_pattern;
      {
        BlockDynamicSparsityPattern csp(blks, blks);
        for (unsigned int bi = 0; bi < blks; ++bi)
          for (unsigned int bj = 0; bj < blks; ++bj)
            csp.block(bi, bj).reinit(rc, rc);

        csp.collect_sizes();
        sparsity_pattern.copy_from(csp);
      }

      BlockSparseMatrix<double> A(sparsity_pattern);
      BlockVector<double>       b(blks, rc);
      for (unsigned int i = 0; i < rc; ++i)
        {
          for (unsigned int bi = 0; bi < blks; ++bi)
            {
              b.block(bi)(i) = bi * rc + i;
              for (unsigned int bj = 0; bj < blks; ++bj)
                A.block(bi, bj).diag_element(i) = 2.0 * bi + 1.5 * bj + (i + 1);
            }
        }

      const auto lo_A = linear_operator(A.block(1, 1));
      const auto lo_B = linear_operator(A.block(1, 0));
      const auto lo_C = linear_operator(A.block(0, 1));
      const auto lo_D = linear_operator(A.block(0, 0));

      Vector<double> &f = b.block(1);
      Vector<double> &g = b.block(0);

      BlockVector<double> s(blks, rc);
      Vector<double> &    x = s.block(1);
      Vector<double> &    y = s.block(0);

      SolverControl            solver_control_A(1, 1.0e-10, false, false);
      SolverCG<Vector<double>> solver_A(solver_control_A);
      PreconditionJacobi<SparseMatrix<double>> preconditioner_A;
      preconditioner_A.initialize(A.block(1, 1));
      const auto lo_A_inv = inverse_operator(lo_A, solver_A, preconditioner_A);

      const auto lo_S = schur_complement(lo_A_inv, lo_B, lo_C, lo_D);

      // Preconditinoed by D
      {
        SolverControl            solver_control_S(11, 1.0e-10, false, false);
        SolverCG<Vector<double>> solver_S(solver_control_S);
        PreconditionJacobi<SparseMatrix<double>> preconditioner_S;
        preconditioner_S.initialize(A.block(0, 0)); // Same space as S
        const auto lo_S_inv =
          inverse_operator(lo_S, solver_S, preconditioner_S);

        auto rhs = condense_schur_rhs(lo_A_inv, lo_C, f, g);
        check_solver_within_range(y = lo_S_inv * rhs, // Solve for y
                                  solver_control_S.last_step(),
                                  11,
                                  11);

        x = postprocess_schur_solution(lo_A_inv, lo_B, y, f);

        PRINTME("x = s.block(1)", x);
        PRINTME("y = s.block(0)", y);
      }

      // Preconditinoed by S_approx_inv
      {
        const auto lo_A_inv_approx = linear_operator(preconditioner_A);
        const auto lo_S_approx =
          schur_complement(lo_A_inv_approx, lo_B, lo_C, lo_D);

        // Setup inner solver: Approximation of inverse of Schur complement
        IterationNumberControl solver_control_S_approx(
          1, 1.0e-10, false, false); // Perform only a limited number of sweeps
        SolverCG<Vector<double>> solver_S_approx(solver_control_S_approx);
        PreconditionJacobi<SparseMatrix<double>> preconditioner_S_approx;
        preconditioner_S_approx.initialize(A.block(0, 0)); // Same space as S
        const auto lo_S_inv_approx = inverse_operator(lo_S_approx,
                                                      solver_S_approx,
                                                      preconditioner_S_approx);

        // Setup outer solver: Exact inverse of Schur complement
        SolverControl            solver_control_S(11, 1.0e-10, false, false);
        SolverCG<Vector<double>> solver_S(solver_control_S);
        const auto lo_S_inv = inverse_operator(lo_S, solver_S, lo_S_inv_approx);

        auto rhs = condense_schur_rhs(lo_A_inv, lo_C, f, g);
        check_solver_within_range(y = lo_S_inv * rhs, // Solve for y
                                  solver_control_S.last_step(),
                                  11,
                                  11);

        x = postprocess_schur_solution(lo_A_inv, lo_B, y, f);

        PRINTME("x = s.block(1)", x);
        PRINTME("y = s.block(0)", y);
      }

      //      A.print(std::cout);
      //      b.print(std::cout);
      //      s.print(std::cout);
    }

    deallog << "BlockSparseMatrix OK" << std::endl;
  }
}
