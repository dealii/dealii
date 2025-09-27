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

// Test internal preconditioner and solver options

#include <deal.II/lac/packaged_operation.h>
#include <deal.II/lac/solver_bicgstab.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/solver_minres.h>
#include <deal.II/lac/solver_qmrs.h>
#include <deal.II/lac/solver_richardson.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_linear_operator.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_precondition.h>
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
  const auto lo_A = linear_operator<VECTOR>(A);
  // Note: The above should be equivalent to the following:
  //
  //  typedef
  //  dealii::TrilinosWrappers::internal::LinearOperatorImplementation::TrilinosPayload
  //  PAYLOAD; const auto lo_A = linear_operator<VECTOR,VECTOR,PAYLOAD>(A);

  PRECONDITIONER preconditioner;
  preconditioner.initialize(A, data);

  using SOLVER = SolverCG<VECTOR>;
  SolverControl solver_control(100, 1.0e-10, false, false);
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

  SolverControl solver_control(100, 1.0e-10, false, false);
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
    SolverControl solver_control_1(100, 1.0e-10, false, false);
    SOLVER        solver_1(solver_control_1);
    const auto    lo_A_inv_1 = inverse_operator(lo_A, solver_1, preconditioner);
    SolverControl solver_control_2(100, 1.0e-10, false, false);
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
      deallog.push("SolverBicgstab");
      using SLVR = SolverBicgstab<TrilinosWrappers::MPI::Vector>;
      test_solver<SLVR>(A, b);
      deallog.pop();
    }

    {
      deallog.push("SolverCG");
      using SLVR = SolverCG<TrilinosWrappers::MPI::Vector>;
      test_solver<SLVR>(A, b);
      deallog.pop();
    }

    {
      deallog.push("SolverGMRES");
      using SLVR = SolverGMRES<TrilinosWrappers::MPI::Vector>;
      test_solver<SLVR>(A, b);
      deallog.pop();
    }

    {
      deallog.push("SolverFGMRES");
      using SLVR = SolverFGMRES<TrilinosWrappers::MPI::Vector>;
      test_solver<SLVR>(A, b);
      deallog.pop();
    }

    {
      deallog.push("SolverMinRes");
      using SLVR = SolverMinRes<TrilinosWrappers::MPI::Vector>;
      test_solver<SLVR>(A, b);
      deallog.pop();
    }

    {
      deallog.push("SolverQMRS");
      using SLVR = SolverQMRS<TrilinosWrappers::MPI::Vector>;
      test_solver<SLVR>(A, b);
      deallog.pop();
    }

    {
      deallog.push("SolverRichardson");
      using SLVR = SolverRichardson<TrilinosWrappers::MPI::Vector>;
      test_solver<SLVR>(A, b);
      deallog.pop();
    }

    deallog.pop();
    deallog << "TrilinosWrappers::SparseMatrix OK" << std::endl;
  } // TrilinosWrappers::SparseMatrix
}
