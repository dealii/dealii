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

#include <deal.II/lac/packaged_operation.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_linear_operator.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>
#include <deal.II/lac/trilinos_vector.h>

#include "../tests.h"



template <typename VECTOR>
void
print(const VECTOR &vec)
{
  for (types::global_dof_index i = 0; i < vec.size(); ++i)
    {
      deallog << vec(i) << " ";
    }
  deallog << std::endl;
}


template <class PRECONDITIONER,
          class MATRIX,
          class VECTOR,
          class ADDITIONAL_DATA = typename PRECONDITIONER::AdditionalData>
void
test_preconditioner(const MATRIX &         A,
                    const VECTOR &         b,
                    const ADDITIONAL_DATA &data = ADDITIONAL_DATA())
{
  const auto lo_A = linear_operator<VECTOR>(A);
  // Note: The above should be equivalent to the following:
  //
  //  typedef
  //  dealii::TrilinosWrappers::internal::LinearOperatorImplementation::TrilinosPayload
  //  PAYLOAD; const auto lo_A = linear_operator<VECTOR,VECTOR,PAYLOAD>(A);

  PRECONDITIONER preconditioner;
  preconditioner.initialize(A, data);

  typedef TrilinosWrappers::SolverCG SOLVER;
  SolverControl                      solver_control(100, 1.0e-10);
  SOLVER                             solver(solver_control);

  // Exact inverse
  const auto lo_A_inv = inverse_operator(lo_A, solver, preconditioner);
  // Note: The above should be equivalent to the following:
  //
  //  const auto lo_A_inv = inverse_operator<PAYLOAD,SOLVER,
  //             PRECONDITIONER,
  //             VECTOR,VECTOR>(lo_A,
  //                            solver,
  //                            preconditioner);

  // Singular operation
  {
    deallog.push("S_Op");
    const VECTOR x = lo_A_inv * b;
    print(x);
    deallog.pop();
  }

  // Composite operation
  {
    deallog.push("C_Op");
    const VECTOR x = (lo_A_inv * lo_A * lo_A_inv) * b;
    print(x);
    deallog.pop();
  }

  // Approximate inverse
  deallog.push("Approx");
  {
    // Using exemplar matrix
    deallog.push("Exemp");
    const auto lo_A_inv_approx =
      linear_operator<VECTOR, VECTOR>(A, preconditioner);
    // Note: The above should be equivalent to the following:
    //
    //    typedef
    //    dealii::TrilinosWrappers::internal::LinearOperatorImplementation::TrilinosPayload
    //    PAYLOAD; const auto lo_A_inv_approx =
    //    linear_operator<VECTOR,VECTOR,PAYLOAD>(A, preconditioner);

    // Singular operation
    {
      deallog.push("S_Op");
      const VECTOR x_approx = lo_A_inv_approx * b;
      print(x_approx);
      deallog.pop();
    }

    // Composite operation
    {
      deallog.push("C_Op");
      const VECTOR x_approx = (lo_A_inv_approx * lo_A * lo_A_inv_approx) * b;
      print(x_approx);
      deallog.pop();
    }

    deallog.pop();
  }
  {
    // Stand-alone
    deallog.push("S.A.");
    typedef dealii::TrilinosWrappers::internal::LinearOperatorImplementation::
      TrilinosPayload PAYLOAD;
    const auto        lo_A_inv_approx =
      linear_operator<VECTOR, VECTOR, PAYLOAD>(preconditioner);

    // Singular operation
    {
      deallog.push("S_Op");
      const VECTOR x_approx = lo_A_inv_approx * b;
      print(x_approx);
      deallog.pop();
    }

    // Composite operation
    {
      deallog.push("C_Op");
      const VECTOR x_approx = (lo_A_inv_approx * lo_A * lo_A_inv_approx) * b;
      print(x_approx);
      deallog.pop();
    }

    deallog.pop();
  }
  deallog.pop();
}

template <class SOLVER, class MATRIX, class VECTOR>
void
test_solver(const MATRIX &A, const VECTOR &b)
{
  const auto lo_A = linear_operator<VECTOR>(A);
  // Note: The above should be equivalent to the following:
  //
  //  typedef
  //  dealii::TrilinosWrappers::internal::LinearOperatorImplementation::TrilinosPayload
  //  PAYLOAD; const auto lo_A = linear_operator<VECTOR,VECTOR,PAYLOAD>(A);

  SolverControl solver_control(100, 1.0e-10);
  SOLVER        solver(solver_control);

  typedef TrilinosWrappers::PreconditionJacobi PRECONDITIONER;
  PRECONDITIONER                               preconditioner;
  preconditioner.initialize(A);

  {
    const auto lo_A_inv = inverse_operator(lo_A, solver, preconditioner);
    // Note: The above should be equivalent to the following:
    //
    //  const auto lo_A_inv = inverse_operator<PAYLOAD,SOLVER,
    //             PRECONDITIONER,
    //             VECTOR,VECTOR>(lo_A,
    //                            solver,
    //                            preconditioner);

    // Singular operation
    {
      deallog.push("S_Op");
      const VECTOR x_approx = lo_A_inv * b;
      print(x_approx);
      deallog.pop();
    }

    // Composite operation
    {
      deallog.push("C_Op");
      const VECTOR x_approx = (lo_A_inv * lo_A * lo_A_inv) * b;
      print(x_approx);
      deallog.pop();
    }
  }

  // Composite operation 2
  {
    deallog.push("C_Op2");
    SolverControl solver_control_1(100, 1.0e-10);
    SOLVER        solver_1(solver_control_1);
    const auto    lo_A_inv_1 = inverse_operator(lo_A, solver_1, preconditioner);
    SolverControl solver_control_2(100, 1.0e-10);
    SOLVER        solver_2(solver_control_2);
    const auto    lo_A_inv_2 = inverse_operator(lo_A, solver_2, preconditioner);
    const VECTOR  x_approx   = (lo_A_inv_2 * lo_A * lo_A_inv_1) * b;
    print(x_approx);
    deallog.pop();
  }
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  initlog();

  deallog.depth_console(0);
  deallog << std::setprecision(10);

  // TrilinosWrappers::SparseMatrix
  {
    const unsigned int                rc = 10;
    TrilinosWrappers::SparsityPattern sparsity_pattern(
      rc, rc, /*n_entries_per_row =*/1);
    for (unsigned int i = 0; i < rc; ++i)
      {
        sparsity_pattern.add(i, i);
      }
    sparsity_pattern.compress();

    TrilinosWrappers::SparseMatrix A(sparsity_pattern);
    TrilinosWrappers::MPI::Vector  b;
    b.reinit(A.locally_owned_domain_indices());
    TrilinosWrappers::MPI::Vector c;
    c.reinit(A.locally_owned_domain_indices());
    for (unsigned int i = 0; i < rc; ++i)
      {
        A.set(i, i, 2.0);
        b(i) = i;
      }

    // === PRECONDITIONERS ===
    deallog << "PRECONDITIONERS" << std::endl;
    deallog.push("Preconditioners");

    {
      deallog.push("PreconditionAMG");
      typedef TrilinosWrappers::PreconditionAMG PREC;
      test_preconditioner<PREC>(A, b);
      test_preconditioner<PREC>(A, c);
      deallog.pop();
    }

#ifdef DEAL_II_TRILINOS_WITH_MUELU
    {
      deallog.push("PreconditionAMGMueLu");
      typedef TrilinosWrappers::PreconditionAMGMueLu PREC;
      test_preconditioner<PREC>(A, b);
      test_preconditioner<PREC>(A, c);
      deallog.pop();
    }
#endif

    {
      deallog.push("PreconditionChebyshev");
      typedef TrilinosWrappers::PreconditionChebyshev PREC;
      test_preconditioner<PREC>(A, b);
      test_preconditioner<PREC>(A, c);
      deallog.pop();
    }

    {
      deallog.push("PreconditionIC");
      typedef TrilinosWrappers::PreconditionIC PREC;
      test_preconditioner<PREC>(A, b);
      test_preconditioner<PREC>(A, c);
      deallog.pop();
    }

    {
      deallog.push("PreconditionIdentity");
      typedef TrilinosWrappers::PreconditionIdentity PREC;
      test_preconditioner<PREC>(A, b);
      test_preconditioner<PREC>(A, c);
      deallog.pop();
    }

    {
      deallog.push("PreconditionILU");
      typedef TrilinosWrappers::PreconditionILU PREC;
      test_preconditioner<PREC>(A, b);
      test_preconditioner<PREC>(A, c);
      deallog.pop();
    }

    {
      deallog.push("PreconditionILUT");
      typedef TrilinosWrappers::PreconditionILUT PREC;
      test_preconditioner<PREC>(A, b);
      test_preconditioner<PREC>(A, c);
      deallog.pop();
    }

    {
      deallog.push("PreconditionJacobi");
      typedef TrilinosWrappers::PreconditionJacobi PREC;
      test_preconditioner<PREC>(A, b);
      test_preconditioner<PREC>(A, c);
      deallog.pop();
    }

    {
      deallog.push("PreconditionSOR");
      typedef TrilinosWrappers::PreconditionSOR PREC;
      test_preconditioner<PREC>(A, b);
      test_preconditioner<PREC>(A, c);
      deallog.pop();
    }

    {
      deallog.push("PreconditionSSOR");
      typedef TrilinosWrappers::PreconditionSSOR PREC;
      test_preconditioner<PREC>(A, b);
      test_preconditioner<PREC>(A, c);
      deallog.pop();
    }

    deallog.pop();

    // === SOLVERS ===
    deallog << std::endl;
    deallog << "SOLVERS" << std::endl;
    deallog.push("Solvers");

    {
      try
        {
          // TODO: AztecOO::Iterate error code -3: loss of precision
          // Note: Even though this test fails, we would still like to check
          //       that we can build LinearOperators with this solver
          deallog.push("SolverBicgstab");
          typedef TrilinosWrappers::SolverBicgstab SLVR;
          test_solver<SLVR>(A, b);
          deallog.pop();
        }
      catch (...)
        {
          deallog.pop();
          deallog << "Known AztecOO error in SolverBicgstab" << std::endl;
          deallog.pop();
        }
    }

    {
      deallog.push("SolverCG");
      typedef TrilinosWrappers::SolverCG SLVR;
      test_solver<SLVR>(A, b);
      deallog.pop();
    }

    {
      deallog.push("SolverCGS");
      typedef TrilinosWrappers::SolverCGS SLVR;
      test_solver<SLVR>(A, b);
      deallog.pop();
    }

    // The TrilinosWrappers::SolverDirect class is not (yet) compatible with
    // the LinearOperator class.
    // {
    //   deallog.push("SolverDirect");
    //   typedef TrilinosWrappers::SolverDirect SLVR;
    //   // test_solver<SLVR> (A, b);
    //
    //   typedef dealii::TrilinosWrappers::MPI::Vector VectorType;
    //   typedef
    //   dealii::TrilinosWrappers::internal::LinearOperatorImplementation::TrilinosPayload
    //   PayloadType;
    //   // TODO: Full template expansion required for composite operator. Can
    //   one prevent this?
    //   //       i.e. is 'const auto lo_A = linear_operator<VectorType>(A);'
    //   possible? const auto lo_A =
    //   linear_operator<VectorType,VectorType,PayloadType>(A);
    //
    //   SolverControl solver_control (100, 1.0e-10);
    //   SLVR solver (solver_control);
    //   solver.initialize(A);
    //
    //   // TODO: Full template expansion required for composite operator. Can
    //   one prevent this?
    //   //       i.e. is 'const auto lo_A_inv =
    //   linear_operator<VectorType>(solver);' possible? const auto lo_A_inv =
    //   linear_operator<VectorType,VectorType,PayloadType>(solver);
    //
    //   // Singular operation
    //   {
    //     deallog.push("S_Op");
    //     const TrilinosWrappers::MPI::Vector x_approx = lo_A_inv*b;
    //     print(x_approx);
    //     deallog.pop();
    //   }
    //
    //   // Composite operation
    //   {
    //     deallog.push("C_Op");
    //     const TrilinosWrappers::MPI::Vector x_approx =
    //     (lo_A_inv*lo_A*lo_A_inv)*b; print(x_approx); deallog.pop();
    //   }
    //
    //   deallog.pop();
    // }

    {
      deallog.push("SolverGMRES");
      typedef TrilinosWrappers::SolverGMRES SLVR;
      test_solver<SLVR>(A, b);
      deallog.pop();
    }

    {
      try
        {
          // TODO: AztecOO::Iterate error code -2: numerical breakdown
          // Note: Even though this test fails, we would still like to check
          //       that we can build LinearOperators with this solver
          deallog.push("SolverTFQMR");
          typedef TrilinosWrappers::SolverTFQMR SLVR;
          test_solver<SLVR>(A, b);
          deallog.pop();
        }
      catch (...)
        {
          deallog.pop();
          deallog << "Known AztecOO error in SolverTFQMR" << std::endl;
          deallog.pop();
        }
    }

    deallog.pop();
    deallog << "TrilinosWrappers::SparseMatrix OK" << std::endl;
  } // TrilinosWrappers::SparseMatrix

  // {
  //   // TODO: Implement these checks
  //   TrilinosWrappers::PreconditionBlockJacobi
  //   TrilinosWrappers::PreconditionBlockSOR
  //   TrilinosWrappers::PreconditionBlockSSOR
  //   TrilinosWrappers::PreconditionBlockwiseDirect
  // }
}
