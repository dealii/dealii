// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2020 by the deal.II authors
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

// Test for composite operations

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/linear_operator.h>
#include <deal.II/lac/packaged_operation.h>
#include <deal.II/lac/precondition.h>
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
  deallog << std::setprecision(10);

  // SparseMatrix:
  {
    const unsigned int rc = 2;
    SparsityPattern    sparsity_pattern(rc, rc, 0);
    sparsity_pattern.compress();

    SparseMatrix<double> A(sparsity_pattern);
    Vector<double>       b(rc);
    for (unsigned int i = 0; i < rc; ++i)
      {
        A.diag_element(i) = 5.0;
        b(i)              = 1.0;
      }

    const auto lo_A   = linear_operator(A);
    const auto lo_A_t = transpose_operator(lo_A);

    // build transpose of inverse
    SolverControl                            solver_control_A_1(100, 1.0e-10);
    SolverCG<Vector<double>>                 solver_A_1(solver_control_A_1);
    PreconditionJacobi<SparseMatrix<double>> preconditioner_A_1;
    preconditioner_A_1.initialize(A);
    const auto lo_A_inv =
      inverse_operator(lo_A, solver_A_1, preconditioner_A_1);
    const auto lo_A_inv_t = transpose_operator(lo_A_inv);

    // build inverse of transpose
    SolverControl                            solver_control_A_2(100, 1.0e-10);
    SolverCG<Vector<double>>                 solver_A_2(solver_control_A_2);
    PreconditionJacobi<SparseMatrix<double>> preconditioner_A_2;
    preconditioner_A_2.initialize(A);
    const auto lo_A_t_inv =
      inverse_operator(lo_A_t, solver_A_2, preconditioner_A_2);



    deallog << "Normal and inverse multiplication operations" << std::endl;

    const Vector<double> x1  = lo_A * b;
    const Vector<double> x2  = lo_A_t * b;
    const Vector<double> x3  = lo_A_inv * b;
    const Vector<double> x4a = lo_A_inv_t * b;
    const Vector<double> x4b = lo_A_t_inv * b;

    //   PRINTME("x1", x1);
    //   PRINTME("x2", x2);
    //   PRINTME("x3", x3);
    //   PRINTME("x4a", x4a);
    //   PRINTME("x4b", x4b);

    deallog << "x4a==x4b : " << (x4a == x4b) << std::endl;

    // Schur-type composition
    SparseMatrix<double> B(sparsity_pattern);
    SparseMatrix<double> C(sparsity_pattern);
    SparseMatrix<double> D(sparsity_pattern);

    for (unsigned int i = 0; i < rc; ++i)
      {
        B.diag_element(i) = 4.0;
        C.diag_element(i) = 4.0; // K = [A,B ; C,D] is symmetric
        D.diag_element(i) = 3.0;
      }

    const auto lo_B   = linear_operator(B);
    const auto lo_C   = linear_operator(C);
    const auto lo_D   = linear_operator(D);
    const auto lo_B_t = transpose_operator(lo_B);
    const auto lo_C_t = transpose_operator(lo_C);
    const auto lo_D_t = transpose_operator(lo_D);

    deallog << "Single packaged operation" << std::endl;
    {
      const auto S = lo_D - lo_C * lo_A_inv * lo_B;
      const auto S_t_1 =
        lo_D_t - lo_B_t * lo_A_inv_t * lo_C_t; // using transpose of inverse
      const auto S_t_2 =
        lo_D_t - lo_B_t * lo_A_t_inv * lo_C_t; // using inverse of transpose

      const Vector<double> x5  = S * b;
      const Vector<double> x6a = S_t_1 * b;
      const Vector<double> x6b = S_t_2 * b;

      deallog << "x5==x6a : " << (x5 == x6a)
              << std::endl; // using transpose of inverse
      deallog << "x5==x6b : " << (x5 == x6b)
              << std::endl; // using inverse of transpose
                            //      PRINTME("x5", x5);
                            //      PRINTME("x6a", x6a);
                            //      PRINTME("x6b", x6b);
    }

    deallog << "Manual operations" << std::endl;
    {
      const Vector<double> x5a = lo_B * b;
      const Vector<double> x5b = lo_A_inv * x5a;
      const Vector<double> x5c = lo_C * x5b;
      const Vector<double> x5d = lo_D * b - x5c;
      const Vector<double> x6a = lo_C_t * b;
      const Vector<double> x6b_1 =
        lo_A_inv_t * x6a; // using transpose of inverse
      const Vector<double> x6c_1 = lo_B_t * x6b_1;
      const Vector<double> x6d_1 = lo_D_t * b - x6c_1;
      const Vector<double> x6b_2 =
        lo_A_t_inv * x6a; // using inverse of transpose
      const Vector<double> x6c_2 = lo_B_t * x6b_2;
      const Vector<double> x6d_2 = lo_D_t * b - x6c_2;

      deallog << "x5a==x6a : " << (x5a == x6a) << std::endl;
      deallog << "x5b==x6b_1 : " << (x5b == x6b_1)
              << std::endl; // using transpose of inverse
      deallog << "x5c==x6c_1 : " << (x5c == x6c_1) << std::endl;
      deallog << "x5d==x6d_1 : " << (x5d == x6d_1) << std::endl;
      deallog << "x5b==x6b_2 : " << (x5b == x6b_2)
              << std::endl; // using inverse of transpose
      deallog << "x5c==x6c_2 : " << (x5c == x6c_2) << std::endl;
      deallog << "x5d==x6d_2 : " << (x5d == x6d_2) << std::endl;
      //     PRINTME("x5a", x5a);
      //     PRINTME("x6a", x6a);
      //     PRINTME("x5b", x5b);
      //     PRINTME("x6b_1", x6b_1);
      //     PRINTME("x6b_2", x6b_2);
      //     PRINTME("x5c", x5c);
      //     PRINTME("x6c_1", x6c_1);
      //     PRINTME("x6c_2", x6c_2);
      //     PRINTME("x5d", x5d);
      //     PRINTME("x6d_1", x6d_1);
      //     PRINTME("x6d_2", x6d_2);
    }

    deallog << "OK" << std::endl;
  }
}
