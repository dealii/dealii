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


#ifndef dealii_trilinos_anasazi_trilinos_solver
#define dealii_trilinos_anasazi_trilinos_solver


#include <deal.II/base/config.h>

#include <deal.II/base/types.h>

#include <deal.II/lac/trilinos_tpetra_types.h>

#ifdef DEAL_II_TRILINOS_WITH_TPETRA
#  include <deal.II/lac/trilinos_tpetra_precondition.h>
#  include <deal.II/lac/trilinos_tpetra_sparse_matrix.h>
#  include <deal.II/lac/trilinos_tpetra_vector.h>
#  ifdef DEAL_II_TRILINOS_WITH_ANASAZI
#    include <deal.II/base/types.h>

#    include <deal.II/lac/solver_control.h>

#    include <AnasaziBasicEigenproblem.hpp>
#    include <AnasaziSolverManager.hpp>
#    include <AnasaziTpetraAdapter.hpp>
#    include <AnasaziTypes.hpp>
#  endif
#endif

DEAL_II_NAMESPACE_OPEN
#ifdef DEAL_II_TRILINOS_WITH_TPETRA
#  ifdef DEAL_II_TRILINOS_WITH_ANASAZI

/**
 * To do:
 * 1. See if we can also add the other SLEPc functions
 * 2. Write tests by transforming the SLEPc tests into Anasazi tests and writing extra tests
 * 3. Generalize tests for the other solvers
 *  Useful resource: https://github.com/t-sakashita/rokko/blob/ebd49e1198c4ec9e7612ad4a9806d16a4ff0bdc9/rokko/anasazi/solver.hpp
 */
namespace LinearAlgebra
{
  namespace TpetraWrappers
  {
    template <typename Number, typename MemorySpace = MemorySpace::Host>
    class AnasaziSolverBase
    {
    public:
      /**
       * Constructor. Takes in the SolverControl
       */
      AnasaziSolverBase(
        SolverControl                              &cn,
        const Teuchos::RCP<Teuchos::ParameterList> &parameters_ =
          Teuchos::rcp(new Teuchos::ParameterList)); // implemented

      void
      solve(
        const LinearAlgebra::TpetraWrappers::SparseMatrix<Number, MemorySpace>
                            &A,
        std::vector<Number> &eigenvalues,
        std::vector<LinearAlgebra::TpetraWrappers::Vector<Number, MemorySpace>>
                          &eigenvectors,
        const unsigned int n_eigenpairs = 1); // implemented

      void
      solve(
        const LinearAlgebra::TpetraWrappers::SparseMatrix<Number, MemorySpace>
          &A,
        const LinearAlgebra::TpetraWrappers::SparseMatrix<Number, MemorySpace>
                            &B,
        std::vector<Number> &eigenvalues,
        std::vector<LinearAlgebra::TpetraWrappers::Vector<Number, MemorySpace>>
                          &eigenvectors,
        const unsigned int n_eigenpairs = 1); // implemented

      void
      set_initial_space(
        const std::vector<
          LinearAlgebra::TpetraWrappers::Vector<Number, MemorySpace>>
          &initial_space); // done

      void
      set_preconditioner(
        LinearAlgebra::TpetraWrappers::PreconditionBase<Number, MemorySpace>
          &Preconditioner);
      
      void
      is_hermitian(bool hermitian);

      /**
       * Exception. Standard exception.
       */
      DeclException0(ExcAnasaziWrappersUsageError);

      /**
       * Exception. Convergence failure on the number of eigenvectors.
       */
      DeclException2(ExcAnasaziEigenvectorConvergenceMismatchError,
                     int,
                     int,
                     << "    The number of converged eigenvectors is " << arg1
                     << " but " << arg2 << " were requested. ");

      /**
       * Access to the object that controls convergence.
       */
      SolverControl &
      control() const; // implemented

    protected:
      /**
       * Reference to the object that controls convergence of the iterative
       * solver.
       */
      SolverControl &solver_control;

      /**
       * Solve the linear system for <code>n_eigenpairs</code> eigenstates.
       * Parameter <code>n_converged</code> contains the actual number of
       * eigenstates that have  converged; this can be both fewer or more than
       * n_eigenpairs, depending on the Anasazi eigensolver used.
       * Note: This implementation does not match the SLEPc wrappes
       */
      void
      solve(const unsigned int n_eigenpairs,
            unsigned int      &n_converged); // implemented

      void
      get_eigenpair(const unsigned int index,
                    Number            &eigenvalues,
                    LinearAlgebra::TpetraWrappers::Vector<Number, MemorySpace>
                      &eigenvectors); // implemented

      virtual void
      set_solver() = 0;

      void
      set_matrices(
        const LinearAlgebra::TpetraWrappers::SparseMatrix<Number, MemorySpace>
          &A); // implemented

      void
      set_matrices(
        const LinearAlgebra::TpetraWrappers::SparseMatrix<Number, MemorySpace>
          &A,
        const LinearAlgebra::TpetraWrappers::SparseMatrix<Number, MemorySpace>
          &B); // implemented

      

    protected:
      Teuchos::RCP<Anasazi::BasicEigenproblem<Number,
                                              Tpetra::MultiVector<Number>,
                                              Tpetra::Operator<Number>>>
        eigenproblem; // exists and used
      Anasazi::Eigensolution<Number, Tpetra::MultiVector<Number>>
        eigensolution; // exists and used

      const Teuchos::RCP<Teuchos::ParameterList> parameters; // exists and used

      Teuchos::RCP<Anasazi::SolverManager<Number,
                                          Tpetra::MultiVector<Number>,
                                          Tpetra::Operator<Number>>>
        eigensolver; // exists and used
    };



    template <typename Number, typename MemorySpace = MemorySpace::Host>
    class AnasaziSolverBlockKrylovSchur
      : public AnasaziSolverBase<Number, MemorySpace>
    {
    public:
      explicit AnasaziSolverBlockKrylovSchur(
        SolverControl                              &cn,
        const Teuchos::RCP<Teuchos::ParameterList> &parameters_ =
          Teuchos::rcp(new Teuchos::ParameterList));

    protected:
      void
      set_solver() override;
    };



    template <typename Number, typename MemorySpace = MemorySpace::Host>
    class AnasaziSolverBlockDavidson
      : public AnasaziSolverBase<Number, MemorySpace>
    {
    public:
      explicit AnasaziSolverBlockDavidson(
        SolverControl                              &cn,
        const Teuchos::RCP<Teuchos::ParameterList> &parameters_ =
          Teuchos::rcp(new Teuchos::ParameterList));

    protected:
      void
      set_solver() override;
    };

    template <typename Number, typename MemorySpace = MemorySpace::Host>
    class AnasaziSolverGeneralizedDavidson
      : public AnasaziSolverBase<Number, MemorySpace>
    {
    public:
      explicit AnasaziSolverGeneralizedDavidson(
        SolverControl                              &cn,
        const Teuchos::RCP<Teuchos::ParameterList> &parameters_ =
          Teuchos::rcp(new Teuchos::ParameterList));

    protected:
      void
      set_solver() override;
    };



    template <typename Number, typename MemorySpace = MemorySpace::Host>
    class AnasaziSolverLOBPCG
      : public AnasaziSolverBase<Number, MemorySpace>
    {
    public:
      explicit AnasaziSolverLOBPCG(
        SolverControl                              &cn,
        const Teuchos::RCP<Teuchos::ParameterList> &parameters_ =
          Teuchos::rcp(new Teuchos::ParameterList));

    protected:
      void
      set_solver() override;
    };

  } // namespace TpetraWrappers
} // namespace LinearAlgebra

#  endif
#endif

DEAL_II_NAMESPACE_CLOSE

#endif
