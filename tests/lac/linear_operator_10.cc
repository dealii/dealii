// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

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
      deallog << vec(i) << ' ';
    }
  deallog << std::endl;
}


template <class PRECONDITIONER,
          class MATRIX,
          class VECTOR,
          class ADDITIONAL_DATA = typename PRECONDITIONER::AdditionalData>
void
test_preconditioner(const MATRIX          &A,
                    const VECTOR          &b,
                    const ADDITIONAL_DATA &data = ADDITIONAL_DATA())
{
  // This test might trigger spurious floating point exceptions in Trilinos
  // despite functioning properly. Simply disable floating point exceptions
  // again (after they had been enabled int tests.h).
#if defined(DEBUG) && defined(DEAL_II_HAVE_FP_EXCEPTIONS)
  {
    const int current_fe_except = fegetexcept();
    fedisableexcept(current_fe_except);
  }
#endif

  const auto lo_A = linear_operator<VECTOR>(A);
  // Note: The above should be equivalent to the following:
  //
  //  typedef
  //  dealii::TrilinosWrappers::internal::LinearOperatorImplementation::TrilinosPayload
  //  PAYLOAD; const auto lo_A = linear_operator<VECTOR,VECTOR,PAYLOAD>(A);

  PRECONDITIONER preconditioner;
  preconditioner.initialize(A, data);

  using SOLVER = TrilinosWrappers::SolverCG;
  SolverControl solver_control(100, 1.0e-10);
  SOLVER        solver(solver_control);

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
    using PAYLOAD = dealii::TrilinosWrappers::internal::
      LinearOperatorImplementation::TrilinosPayload;
    const auto lo_A_inv_approx =
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

  using PRECONDITIONER = TrilinosWrappers::PreconditionJacobi;
  PRECONDITIONER preconditioner;
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
      using PREC = TrilinosWrappers::PreconditionAMG;
      test_preconditioner<PREC>(A, b);
      test_preconditioner<PREC>(A, c);
      deallog.pop();
    }

#ifdef DEAL_II_TRILINOS_WITH_MUELU
    {
      deallog.push("PreconditionAMGMueLu");
      using PREC = TrilinosWrappers::PreconditionAMGMueLu;
      test_preconditioner<PREC>(A, b);
      test_preconditioner<PREC>(A, c);
      deallog.pop();
    }
#endif

    {
      deallog.push("PreconditionChebyshev");
      using PREC = TrilinosWrappers::PreconditionChebyshev;
      test_preconditioner<PREC>(A, b);
      test_preconditioner<PREC>(A, c);
      deallog.pop();
    }

    {
      deallog.push("PreconditionIC");
      using PREC = TrilinosWrappers::PreconditionIC;
      test_preconditioner<PREC>(A, b);
      test_preconditioner<PREC>(A, c);
      deallog.pop();
    }

    {
      deallog.push("PreconditionIdentity");
      using PREC = TrilinosWrappers::PreconditionIdentity;
      test_preconditioner<PREC>(A, b);
      test_preconditioner<PREC>(A, c);
      deallog.pop();
    }

    {
      deallog.push("PreconditionILU");
      using PREC = TrilinosWrappers::PreconditionILU;
      test_preconditioner<PREC>(A, b);
      test_preconditioner<PREC>(A, c);
      deallog.pop();
    }

    {
      deallog.push("PreconditionILUT");
      using PREC = TrilinosWrappers::PreconditionILUT;
      test_preconditioner<PREC>(A, b);
      test_preconditioner<PREC>(A, c);
      deallog.pop();
    }

    {
      deallog.push("PreconditionJacobi");
      using PREC = TrilinosWrappers::PreconditionJacobi;
      test_preconditioner<PREC>(A, b);
      test_preconditioner<PREC>(A, c);
      deallog.pop();
    }

    {
      deallog.push("PreconditionSOR");
      using PREC = TrilinosWrappers::PreconditionSOR;
      test_preconditioner<PREC>(A, b);
      test_preconditioner<PREC>(A, c);
      deallog.pop();
    }

    {
      deallog.push("PreconditionSSOR");
      using PREC = TrilinosWrappers::PreconditionSSOR;
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
          using SLVR = TrilinosWrappers::SolverBicgstab;
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
      using SLVR = TrilinosWrappers::SolverCG;
      test_solver<SLVR>(A, b);
      deallog.pop();
    }

    {
      deallog.push("SolverCGS");
      using SLVR = TrilinosWrappers::SolverCGS;
      test_solver<SLVR>(A, b);
      deallog.pop();
    }

    // The TrilinosWrappers::SolverDirect class is not (yet) compatible with
    // the LinearOperator class.
    // {
    //   deallog.push("SolverDirect");
    //   using SLVR = TrilinosWrappers::SolverDirect;
    //   // test_solver<SLVR> (A, b);
    //
    //   using VectorType = dealii::TrilinosWrappers::MPI::Vector;
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
      using SLVR = TrilinosWrappers::SolverGMRES;
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
          using SLVR = TrilinosWrappers::SolverTFQMR;
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
