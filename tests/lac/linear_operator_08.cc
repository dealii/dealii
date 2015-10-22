// ---------------------------------------------------------------------
//
// Copyright (C) 2015 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

// Test internal preconditioner and solver options

#include "../tests.h"

#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/linear_operator.h>
#include <deal.II/lac/packaged_operation.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/precondition_block.h>
#include <deal.II/lac/precondition_selector.h>
#include <deal.II/lac/sparse_ilu.h>
#include <deal.II/lac/sparse_decomposition.h>
#include <deal.II/lac/sparse_mic.h>
#include <deal.II/lac/sparse_vanka.h>

#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_bicgstab.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/solver_minres.h>
#include <deal.II/lac/solver_qmrs.h>
#include <deal.II/lac/solver_relaxation.h>
#include <deal.II/lac/solver_richardson.h>
#include <deal.II/lac/solver_selector.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/iterative_inverse.h>


using namespace dealii;


template<class PRECONDITIONER, class MATRIX, class VECTOR,
         class ADDITIONAL_DATA = typename PRECONDITIONER::AdditionalData>
void
test_preconditioner (const MATRIX &A,
                     const VECTOR &b,
                     const ADDITIONAL_DATA &data = ADDITIONAL_DATA())
{
  const auto lo_A = linear_operator<VECTOR>(A);

  PRECONDITIONER preconditioner;
  preconditioner.initialize(A, data);

  SolverControl solver_control (100, 1.0e-10);
  SolverCG<VECTOR> solver (solver_control);

  // Exact inverse
  const auto lo_A_inv = inverse_operator(lo_A,
                                         solver,
                                         preconditioner);

  const VECTOR x = lo_A_inv*b;

  // Approximate inverse
  {
    // Using exemplar matrix
    const auto lo_A_inv_approx = linear_operator<VECTOR>(A, preconditioner);
    const VECTOR x_approx = lo_A_inv_approx*b;
  }
  {
    // Stand-alone
    const auto lo_A_inv_approx = linear_operator<VECTOR>(preconditioner);
    const VECTOR x_approx = lo_A_inv_approx*b;
  }
}

// For Vector <double>
// Cannot use more generic function as Vector <double>
// does not define vector_type
template<class PRECONDITIONER>
void
test_preconditioner (const SparseMatrix<double> &A,
                     const Vector<double> &b,
                     const typename PRECONDITIONER::AdditionalData &data
                     = typename PRECONDITIONER::AdditionalData ())
{
  const auto lo_A = linear_operator(A);

  PRECONDITIONER preconditioner;
  preconditioner.initialize(A, data);

  // Exact inverse
  {
    deallog.push("Exact inverse");
    SolverControl solver_control (100, 1.0e-10);
    SolverCG< Vector<double> > solver (solver_control);
    const auto lo_A_inv = inverse_operator(lo_A,
                                           solver,
                                           preconditioner);

    const Vector<double> x = lo_A_inv*b;
    deallog.pop();
  }

  // Approximate inverses:
  // Using an exemplar matrix
  {
    deallog.push("Exemplar matrix");
    const auto lo_A_inv_approx = linear_operator(A, preconditioner);
    const Vector<double> x_approx = lo_A_inv_approx*b;
    deallog.pop();
  }
  // Stand-alone
  {
    deallog.push("Stand-alone");
    const auto lo_A_inv_approx = linear_operator(preconditioner);
    const Vector<double> x_approx = lo_A_inv_approx*b;
    deallog.pop();
  }
}

template<class SOLVER>
void
test_solver (const SparseMatrix<double> &A,
             const Vector<double> &b)
{
  // Standard solver
  {
    deallog.push("Standard solver");
    SolverControl solver_control (100, 1.0e-10);
    SOLVER solver (solver_control);

    PreconditionJacobi< SparseMatrix<double> > preconditioner;
    preconditioner.initialize(A);

    Vector<double> x;
    x.reinit(b);

    solver.solve(A,x,b,preconditioner);
    deallog.pop();
  }

  // Linear operator
  {
    deallog.push("Linear operator");
    const auto lo_A = linear_operator(A);

    SolverControl solver_control (100, 1.0e-10);
    SOLVER solver (solver_control);

    PreconditionJacobi< SparseMatrix<double> > preconditioner;
    preconditioner.initialize(A);

    const auto lo_A_inv = inverse_operator(lo_A,
                                           solver,
                                           preconditioner);
    const Vector<double> x = lo_A_inv*b;
    deallog.pop();
  }
}

template <typename MatrixType>
class PreconditionBlockIdentity : public BlockMatrixBase<MatrixType>
{
public:

  struct AdditionalData
  {
    AdditionalData () {}
  };

  virtual ~PreconditionBlockIdentity ()
  {}

  void
  initialize (const BlockMatrixBase<MatrixType> &matrix,
              const AdditionalData &additional_data = AdditionalData())
  {
    this->row_block_indices    = matrix.get_row_indices();
    this->column_block_indices = matrix.get_column_indices();

    this->sub_objects.reinit (matrix.n_block_rows(),
                              matrix.n_block_cols());
  }

  template<typename VECTOR>
  void
  vmult(VECTOR &dst, const VECTOR &src) const
  {
    dst = src;
  }

  template<class VECTOR>
  void
  Tvmult(VECTOR &dst, const VECTOR &src) const
  {
    dst = src;
  }
};

// Not tested:
// PreconditionLU
// PreconditionMG
// PreconditionUseMatrix
// The following don't work as expected: vmult acts on Vectors, not BlockVectors
// PreconditionBlockJacobi
// PreconditionBlockSOR
// PreconditionBlockSSOR
// SparseBlockVanka

int main()
{
  initlog();
  deallog.depth_console(0);
  deallog << std::setprecision(10);

  // deal.II SparseMatrix
  {
    const unsigned int rc=10;
    SparsityPattern sparsity_pattern (rc, rc, 0);
    sparsity_pattern.compress();

    SparseMatrix<double> A (sparsity_pattern);
    Vector<double> b (rc);
    for (unsigned int i=0; i < rc; ++i)
      {
        A.diag_element(i) = 2.0;
        b(i) = i;
      }

    // === PRECONDITIONERS ===
    deallog << "Preconditioners" << std::endl;
    deallog.push("Preconditioners");

    {
      deallog << "PreconditionChebyshev" << std::endl;
      typedef PreconditionChebyshev< SparseMatrix<double> > PREC;
      test_preconditioner<PREC>(A, b);
    }
    {
      deallog << "PreconditionIdentity" << std::endl;
      typedef PreconditionIdentity PREC;
      test_preconditioner<PREC>(A, b);
    }
    {
      deallog << "PreconditionJacobi" << std::endl;
      typedef PreconditionJacobi< SparseMatrix<double> > PREC;
      test_preconditioner<PREC>(A, b);
    }
    {
      deallog << "PreconditionPSOR" << std::endl;
      typedef PreconditionPSOR< SparseMatrix<double> > PREC;
      std::vector<unsigned int> permutation(b.size());
      std::vector<unsigned int> inverse_permutation(b.size());
      test_preconditioner<PREC>(A, b,
                                typename PREC::AdditionalData(permutation,
                                                              inverse_permutation));
    }
    {
      deallog << "PreconditionRichardson" << std::endl;
      typedef PreconditionRichardson PREC;
      test_preconditioner<PREC>(A, b);
    }
    {
      deallog << "PreconditionSelector" << std::endl;
      const auto lo_A = linear_operator(A);

      PreconditionSelector< SparseMatrix<double>, Vector<double> >
      preconditioner ("jacobi");
      preconditioner.use_matrix(A);

      SolverControl solver_control (100, 1.0e-10);
      SolverCG< Vector<double> > solver (solver_control);

      // Exact inverse
      const auto lo_A_inv = inverse_operator(lo_A,
                                             solver,
                                             preconditioner);

      const Vector<double> x = lo_A_inv*b;

      // Approximate inverse
      const auto lo_A_inv_approx = linear_operator(preconditioner);
      const Vector<double> x_approx = lo_A_inv_approx*b;
    }
    {
      deallog << "PreconditionSOR" << std::endl;
      typedef PreconditionSOR< SparseMatrix<double> > PREC;
      test_preconditioner<PREC>(A, b);
    }
    {
      deallog << "PreconditionSSOR" << std::endl;
      typedef PreconditionSSOR< SparseMatrix<double> > PREC;
      test_preconditioner<PREC>(A, b);
    }
    {
      deallog << "SparseILU" << std::endl;
      typedef SparseILU<double> PREC;
      test_preconditioner<PREC>(A, b);
    }
    {
      deallog << "SparseMIC" << std::endl;
      typedef SparseMIC<double> PREC;
      test_preconditioner<PREC>(A, b);
    }
    {
      deallog << "SparseVanka" << std::endl;
      typedef SparseVanka<double> PREC;
      typedef typename PREC::AdditionalData PREC_AD;
      test_preconditioner<PREC>(A, b,
                                PREC_AD(std::vector<bool>(rc, true)));
    }
    deallog.pop();

    // === SOLVERS ===
    deallog << std::endl;
    deallog << "Solvers" << std::endl;
    deallog.push("Solvers");

    {
      deallog << "SolverBicgstab" << std::endl;
      typedef SolverBicgstab< Vector<double> > SLVR;
      test_solver<SLVR> (A, b);
    }
    {
      deallog << "SolverCG" << std::endl;
      typedef SolverCG< Vector<double> > SLVR;
      test_solver<SLVR> (A, b);
    }
    {
      deallog << "SolverFGMRES" << std::endl;
      typedef SolverFGMRES< Vector<double> > SLVR;
      test_solver<SLVR> (A, b);
    }
    {
      deallog << "SolverGMRES" << std::endl;
      typedef SolverGMRES< Vector<double> > SLVR;
      test_solver<SLVR> (A, b);
    }
    {
      deallog << "SolverMinRes" << std::endl;
      typedef SolverMinRes< Vector<double> > SLVR;
      test_solver<SLVR> (A, b);
    }
    {
      deallog << "SolverQMRS" << std::endl;
      typedef SolverQMRS< Vector<double> > SLVR;
      test_solver<SLVR> (A, b);
    }
    {
      deallog << "SolverRelaxation" << std::endl;
      typedef SolverRelaxation< Vector<double> > SLVR;
      test_solver<SLVR> (A, b);
    }
    {
      deallog << "SolverRichardson" << std::endl;
      typedef SolverRichardson< Vector<double> > SLVR;
      test_solver<SLVR> (A, b);
    }
    {
      deallog << "SolverSelector" << std::endl;
      const auto lo_A = linear_operator(A);

      ReductionControl solver_control (10, 1.e-30, 1.e-2);
      SolverSelector< Vector<double> > solver;
      solver.select("cg");
      solver.set_control(solver_control);

      PreconditionJacobi< SparseMatrix<double> > preconditioner;
      preconditioner.initialize(A);

      const auto lo_A_inv = inverse_operator(lo_A,
                                             solver,
                                             preconditioner);

      const Vector<double> b (rc);
      const Vector<double> x = lo_A_inv*b;
    }
    {
      deallog << "SparseDirectUMFPACK" << std::endl;
      const auto lo_A = linear_operator(A);

      SparseDirectUMFPACK solver;
      solver.initialize(A);

      const auto lo_A_inv = linear_operator(solver);
      const Vector<double> x = lo_A_inv*b;
    }
    {
      deallog << "IterativeInverse" << std::endl;

      PreconditionJacobi< SparseMatrix<double> > preconditioner;
      preconditioner.initialize(A);

      ReductionControl solver_control (10, 1.e-30, 1.e-2);
      IterativeInverse< Vector<double> > A_inv;
      A_inv.initialize(A,preconditioner);
      A_inv.solver.select("cg");
      A_inv.solver.set_control(solver_control);

      const auto lo_A_inv = linear_operator(A_inv);
      const Vector<double> x = lo_A_inv*b;
    }
    deallog.pop();


    deallog << "SparseMatrix OK" << std::endl;
  }

  // deal.II BlockSparseMatrix
  {
    const unsigned int blks=2;
    const unsigned int rc=10;
    BlockSparsityPattern sparsity_pattern;
    {
      BlockCompressedSimpleSparsityPattern csp(blks, blks);
      for (unsigned int bi=0; bi<blks; ++bi)
        for (unsigned int bj=0; bj<blks; ++bj)
          csp.block(bi,bj).reinit(rc,rc);

      csp.collect_sizes();
      sparsity_pattern.copy_from(csp);
    }

    BlockSparseMatrix<double> A (sparsity_pattern);
    BlockVector<double> b (blks,rc);
    for (unsigned int bi=0; bi<blks; ++bi)
      for (unsigned int i=0; i<rc; ++i)
        {
          A.block(bi,bi).diag_element(i) = 2.0;
          b.block(bi)(i) = bi*rc + i;
        }

    // === PRECONDITIONERS ===
    {
      deallog << "PreconditionBlockIdentity" << std::endl;
      typedef PreconditionBlockIdentity< SparseMatrix<double> > PREC;
      test_preconditioner<PREC>(A, b);
    }

    deallog << "BlockSparseMatrix OK" << std::endl;
  }

}
