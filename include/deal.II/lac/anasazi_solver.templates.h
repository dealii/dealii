// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2009 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


#ifndef dealii_trilinos_anasazi_trilinos_solver_templates_h
#define dealii_trilinos_anasazi_trilinos_solver_templates_h

#include <deal.II/base/config.h>
#ifdef DEAL_II_TRILINOS_WITH_TPETRA
#  include <deal.II/lac/trilinos_tpetra_sparse_matrix.h>
#  include <deal.II/lac/trilinos_tpetra_vector.h>
// including this header fixed some issues i encountered
#  ifdef DEAL_II_TRILINOS_WITH_ANASAZI
#    include <deal.II/lac/anasazi_solver.h>

#    include <AnasaziBasicEigenproblem.hpp>
#    include <AnasaziBlockDavidsonSolMgr.hpp>
#    include <AnasaziBlockKrylovSchurSolMgr.hpp>
#    include <AnasaziGeneralizedDavidsonSolMgr.hpp>
#    include <AnasaziLOBPCGSolMgr.hpp>
#    include <AnasaziSolverManager.hpp>
#    include <AnasaziTpetraAdapter.hpp>
#  endif
#endif

DEAL_II_NAMESPACE_OPEN
#ifdef DEAL_II_TRILINOS_WITH_TPETRA
#  ifdef DEAL_II_TRILINOS_WITH_ANASAZI

namespace LinearAlgebra
{
  namespace TpetraWrappers
  {
    template <typename Number, typename MemorySpace>
    AnasaziSolverBase<Number, MemorySpace>::AnasaziSolverBase(
      SolverControl                              &cn,
      const Teuchos::RCP<Teuchos::ParameterList> &parameters_)
      : solver_control(cn)
      , parameters(parameters_)
    {
      // set the eigenvalue problem specifically
      // the BasicEigenproblem should be enough
      this->eigenproblem = Teuchos::rcp(
        new Anasazi::BasicEigenproblem<Number,
                                       Tpetra::MultiVector<Number>,
                                       Tpetra::Operator<Number>>());

      parameters_->set("Convergence Tolerance", solver_control.tolerance());
      parameters_->set("Maximum Iterations", solver_control.max_steps());
    }



    template <typename Number, typename MemorySpace>
    void
    AnasaziSolverBase<Number, MemorySpace>::solve(
      const LinearAlgebra::TpetraWrappers::SparseMatrix<Number, MemorySpace> &A,
      std::vector<Number> &eigenvalues,
      std::vector<LinearAlgebra::TpetraWrappers::Vector<Number, MemorySpace>>
                        &eigenvectors,
      const unsigned int n_eigenpairs)
    {
      // Panic if the number of eigenpairs wanted is out of bounds.
      AssertThrow((n_eigenpairs > 0) && (n_eigenpairs <= A.m()),
                  ExcAnasaziWrappersUsageError());

      // set the matrix of the problem
      set_matrices(A);

      // and solve
      unsigned int n_converged = 0;
      solve(n_eigenpairs, n_converged);

      if (n_converged > n_eigenpairs)
        {
          n_converged = n_eigenpairs;
        }

      AssertThrow(n_converged == n_eigenpairs,
                  ExcAnasaziEigenvectorConvergenceMismatchError(n_converged,
                                                                n_eigenpairs));

      AssertThrow(eigenvectors.size() != 0, ExcAnasaziWrappersUsageError());
      eigenvectors.resize(n_converged, eigenvectors.front());
      eigenvalues.resize(n_converged);

      for (unsigned int index = 0; index < n_converged; ++index)
        {
          get_eigenpair(index, eigenvalues[index], eigenvectors[index]);
        }
    }



    template <typename Number, typename MemorySpace>
    void
    AnasaziSolverBase<Number, MemorySpace>::solve(
      const LinearAlgebra::TpetraWrappers::SparseMatrix<Number, MemorySpace> &A,
      const LinearAlgebra::TpetraWrappers::SparseMatrix<Number, MemorySpace> &B,
      std::vector<Number> &eigenvalues,
      std::vector<LinearAlgebra::TpetraWrappers::Vector<Number, MemorySpace>>
                        &eigenvectors,
      const unsigned int n_eigenpairs)
    {
      // Guard against incompatible matrix sizes:
      AssertThrow(A.m() == B.m(), ExcDimensionMismatch(A.m(), B.m()));
      AssertThrow(A.n() == B.n(), ExcDimensionMismatch(A.n(), B.n()));

      // Panic if the number of eigenpairs wanted is out of bounds.
      AssertThrow((n_eigenpairs > 0) && (n_eigenpairs <= A.m()),
                  ExcAnasaziWrappersUsageError());

      // set the matrix of the problem
      set_matrices(A, B);

      // and solve
      unsigned int n_converged = 0;
      solve(n_eigenpairs, n_converged);

      if (n_converged > n_eigenpairs)
        n_converged = n_eigenpairs;
      AssertThrow(n_converged == n_eigenpairs,
                  ExcAnasaziEigenvectorConvergenceMismatchError(n_converged,
                                                                n_eigenpairs));

      AssertThrow(eigenvectors.size() != 0, ExcAnasaziWrappersUsageError());
      eigenvectors.resize(n_converged, eigenvectors.front());
      eigenvalues.resize(n_converged);

      for (unsigned int index = 0; index < n_converged; ++index)
        {
          get_eigenpair(index, eigenvalues[index], eigenvectors[index]);
        }
    }



    template <typename Number, typename MemorySpace>
    void
    AnasaziSolverBase<Number, MemorySpace>::set_initial_space(
      const std::vector<
        LinearAlgebra::TpetraWrappers::Vector<Number, MemorySpace>>
        &initial_space)
    {
      // this can probably also be simplified
      using MultiVector =
        Tpetra::MultiVector<Number,
                            int,
                            dealii::types::signed_global_dof_index,
                            TpetraTypes::NodeType<MemorySpace>>;
      Teuchos::RCP<const Tpetra::Map<int,
                                     types::signed_global_dof_index,
                                     TpetraTypes::NodeType<MemorySpace>>>
        map = initial_space[0].trilinos_vector().getMap();
      Teuchos::RCP<MultiVector> multivector =
        Teuchos::rcp(new MultiVector(map, initial_space.size()));

      for (std::size_t i = 0; i < initial_space.size(); ++i)
        {
          Teuchos::RCP<Tpetra::Vector<Number,
                                      int,
                                      types::signed_global_dof_index,
                                      TpetraTypes::NodeType<MemorySpace>>>
            column = multivector->getVectorNonConst(i);
          column->assign(initial_space[i].trilinos_vector());
        }

      this->eigenproblem->setInitVec(multivector);
    }



    template <typename Number, typename MemorySpace>
    void
    AnasaziSolverBase<Number, MemorySpace>::set_preconditioner(
      LinearAlgebra::TpetraWrappers::PreconditionBase<Number, MemorySpace>
        &Preconditioner)
    {
      this->eigenproblem->setPrec(Preconditioner.trilinos_rcp());
    }



    template <typename Number, typename MemorySpace>
    void
    AnasaziSolverBase<Number, MemorySpace>::is_hermitian(bool hermitian)
    {
      this->eigenproblem->setHermitian(hermitian);
    }



    template <typename Number, typename MemorySpace>
    void
    AnasaziSolverBase<Number, MemorySpace>::solve(
      const unsigned int n_eigenpairs,
      unsigned int      &n_converged)
    {
      // set the number of eigenvalues to be computed
      this->eigenproblem->setNEV(static_cast<int>(n_eigenpairs));

      if (!this->eigenproblem->getInitVec())
        {
          // determine what we do
        }

      // basically we tell Aasazi that we are ready
      bool problem_set = this->eigenproblem->setProblem();
      if (!problem_set)
        {
          AssertThrow(
            false,
            "The Eigenvalue problem could not be set inside of Trilinos. This can happen because no initial space or matrix was defined.");
        }

      // we set the solver the user needs.
      set_solver();

      // now we actually do the solve
      this->eigensolver->solve();

      // grab the solution from the eigenvalue problem
      eigensolution = eigenproblem->getSolution();

      //
      n_converged = static_cast<unsigned int>(eigensolution.numVecs);
    }



    template <typename Number, typename MemorySpace>
    void
    AnasaziSolverBase<Number, MemorySpace>::get_eigenpair(
      const unsigned int                                          index,
      Number                                                     &eigenvalues,
      LinearAlgebra::TpetraWrappers::Vector<Number, MemorySpace> &eigenvectors)
    {
      if constexpr (std::is_same_v<Number, std::complex<double>> ||
                    std::is_same_v<Number, std::complex<float>>)
        {
          eigenvalues = Number(eigensolution.Evals[index].realpart,
                               eigensolution.Evals[index].imagpart);
        }
      else
        {
          eigenvalues = eigensolution.Evals[index].realpart;
        }
      // eigenvectors.trilinos_vector() = *((*eigensolution.Evecs)(index));

      // not sure if this is correct
      auto evec = eigensolution.Evecs->getVector(index);
      eigenvectors.trilinos_vector().update(1.0, *evec, 0.0);
    }



    template <typename Number, typename MemorySpace>
    void
    AnasaziSolverBase<Number, MemorySpace>::set_matrices(
      const LinearAlgebra::TpetraWrappers::SparseMatrix<Number, MemorySpace> &A)
    {
      // in case we have a basic eigenvalue problem
      // it is easiest to just set the operator.
      this->eigenproblem->setOperator(A.trilinos_rcp());
    }



    template <typename Number, typename MemorySpace>
    void
    AnasaziSolverBase<Number, MemorySpace>::set_matrices(
      const LinearAlgebra::TpetraWrappers::SparseMatrix<Number, MemorySpace> &A,
      const LinearAlgebra::TpetraWrappers::SparseMatrix<Number, MemorySpace> &B)
    {
      // in the case where we have two matrices
      // it is a bit cleaner to set the specific operators
      this->eigenproblem->setA(A.trilinos_rcp());
      this->eigenproblem->setM(B.trilinos_rcp());
    }



    template <typename Number, typename MemorySpace>
    AnasaziSolverBlockKrylovSchur<Number, MemorySpace>::
      AnasaziSolverBlockKrylovSchur(
        SolverControl                              &cn,
        const Teuchos::RCP<Teuchos::ParameterList> &parameters_)
      : AnasaziSolverBase<Number, MemorySpace>(cn, parameters_)
    {}



    template <typename Number, typename MemorySpace>
    void
    AnasaziSolverBlockKrylovSchur<Number, MemorySpace>::set_solver()
    {
      this->eigensolver = Teuchos::rcp(
        new Anasazi::BlockKrylovSchurSolMgr<Number,
                                            Tpetra::MultiVector<Number>,
                                            Tpetra::Operator<Number>>(
          this->eigenproblem, *(this->parameters)));
    }



    template <typename Number, typename MemorySpace>
    AnasaziSolverBlockDavidson<Number, MemorySpace>::AnasaziSolverBlockDavidson(
      SolverControl                              &cn,
      const Teuchos::RCP<Teuchos::ParameterList> &parameters_)
      : AnasaziSolverBase<Number, MemorySpace>(cn, parameters_)
    {}



    template <typename Number, typename MemorySpace>
    void
    AnasaziSolverBlockDavidson<Number, MemorySpace>::set_solver()
    {
      this->eigensolver = Teuchos::rcp(
        new Anasazi::BlockDavidsonSolMgr<Number,
                                         Tpetra::MultiVector<Number>,
                                         Tpetra::Operator<Number>>(
          this->eigenproblem, *(this->parameters)));
    }



    template <typename Number, typename MemorySpace>
    AnasaziSolverGeneralizedDavidson<Number, MemorySpace>::
      AnasaziSolverGeneralizedDavidson(
        SolverControl                              &cn,
        const Teuchos::RCP<Teuchos::ParameterList> &parameters_)
      : AnasaziSolverBase<Number, MemorySpace>(cn, parameters_)
    {}



    template <typename Number, typename MemorySpace>
    void
    AnasaziSolverGeneralizedDavidson<Number, MemorySpace>::set_solver()
    {
      this->eigensolver = Teuchos::rcp(
        new Anasazi::GeneralizedDavidsonSolMgr<Number,
                                               Tpetra::MultiVector<Number>,
                                               Tpetra::Operator<Number>>(
          this->eigenproblem, *(this->parameters)));
    }



    template <typename Number, typename MemorySpace>
    AnasaziSolverLOBPCG<Number, MemorySpace>::AnasaziSolverLOBPCG(
      SolverControl                              &cn,
      const Teuchos::RCP<Teuchos::ParameterList> &parameters_)
      : AnasaziSolverBase<Number, MemorySpace>(cn, parameters_)
    {}



    template <typename Number, typename MemorySpace>
    void
    AnasaziSolverLOBPCG<Number, MemorySpace>::set_solver()
    {
      this->eigensolver =
        Teuchos::rcp(new Anasazi::LOBPCGSolMgr<Number,
                                               Tpetra::MultiVector<Number>,
                                               Tpetra::Operator<Number>>(
          this->eigenproblem, *(this->parameters)));
    }


    /* ---------------------- SolverControls ----------------------- */
    template <typename Number, typename MemorySpace>
    SolverControl &
    AnasaziSolverBase<Number, MemorySpace>::control() const
    {
      return solver_control;
    }
  } // namespace TpetraWrappers
} // namespace LinearAlgebra
#  endif
#endif
DEAL_II_NAMESPACE_CLOSE

#endif
