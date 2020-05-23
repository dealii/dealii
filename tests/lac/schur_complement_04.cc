// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2020 by the deal.II authors
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

// Test schur complement for Trilinos vectors
// Should give the equivalent output as schur_complement_02

#include <deal.II/lac/linear_operator.h>
#include <deal.II/lac/packaged_operation.h>
#include <deal.II/lac/schur_complement.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_linear_operator.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>
#include <deal.II/lac/trilinos_vector.h>

#include <cctype>
#include <sstream>

#include "../tests.h"

// Have to remove the control commands returned from
// Vector::print in order for the deallog to print properly
// http://www.cplusplus.com/forum/general/76900/
#define PRINTME(name, var)                                         \
  {                                                                \
    std::ostringstream stream;                                     \
    var.print(stream);                                             \
    std::string str = stream.str();                                \
    str.resize(remove_if(str.begin(),                              \
                         str.end(),                                \
                         [](char x) { return std::iscntrl(x); }) - \
               str.begin());                                       \
    deallog << "RHS vector: " << name << ": " << str << std::endl; \
  }



int
main(int argc, char **argv)
{
  initlog();
  deallog.depth_console(0);
  deallog << std::setprecision(10);

  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

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

      const unsigned int                    rc = 1;
      typedef TrilinosWrappers::MPI::Vector VectorType;

      TrilinosWrappers::SparseMatrix A(rc, rc, rc);
      TrilinosWrappers::SparseMatrix B(rc, rc, rc);
      TrilinosWrappers::SparseMatrix C(rc, rc, rc);
      TrilinosWrappers::SparseMatrix D(rc, rc, rc);

      VectorType y;
      y.reinit(complete_index_set(rc));
      VectorType g;
      g.reinit(complete_index_set(rc));
      for (unsigned int i = 0; i < rc; ++i)
        {
          A.set(i, i, 1.0 * (i + 1));
          B.set(i, i, 2.0 * (i + 1));
          C.set(i, i, 3.0 * (i + 1));
          D.set(i, i, 4.0 * (i + 1));
          y(i) = 6.0 * (i + 1);
          g(i) = 2.0 * (i + 1);
        }
      A.compress(VectorOperation::insert);
      B.compress(VectorOperation::insert);
      C.compress(VectorOperation::insert);
      D.compress(VectorOperation::insert);

      // Note that the vector type has to be specified
      // Also note that this call to linear_operator is
      // actually to TrilinosWrappers::linear_operator
      const auto lo_A = linear_operator<VectorType>(A);
      const auto lo_B = linear_operator<VectorType>(B);
      const auto lo_C = linear_operator<VectorType>(C);
      const auto lo_D = linear_operator<VectorType>(D);

      SolverControl                        solver_control_A(100, 1.0e-10);
      TrilinosWrappers::SolverCG           solver_A(solver_control_A);
      TrilinosWrappers::PreconditionJacobi preconditioner_A;
      preconditioner_A.initialize(A);
      const auto lo_A_inv = inverse_operator(lo_A, solver_A, preconditioner_A);

      const auto lo_S   = schur_complement(lo_A_inv, lo_B, lo_C, lo_D);
      const auto lo_S_t = transpose_operator(lo_S);

      const VectorType g1 = lo_S * y;
      const VectorType g2 = lo_S_t * y;
      const VectorType g3 = lo_S * y + g;
      const VectorType g4 = lo_S_t * y + g;

      PRINTME("g1", g1);
      PRINTME("g2", g2);
      PRINTME("g3", g3);
      PRINTME("g4", g4);
    }

    deallog << "SparseMatrix OK" << std::endl;
  }
}
