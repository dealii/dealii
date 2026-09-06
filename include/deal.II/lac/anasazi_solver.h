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

#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>
#ifdef DEAL_II_TRILINOS_WITH_ANASAZI
#  include <deal.II/base/types.h>

#  include <deal.II/lac/solver_control.h>

#  include <AnasaziBasicEigenproblem.hpp>
#  include <AnasaziBlockKrylovSchurSolMgr.hpp>
#  include <AnasaziEpetraAdapter.hpp>
#  include <AnasaziSolverManager.hpp>
#endif

DEAL_II_NAMESPACE_OPEN
#ifdef DEAL_II_TRILINOS_WITH_ANASAZI

/**
 * To do:
 * 1. Finish the remaining functions that are not implemented yet
 * 2. See if we can also add the other SLEPc functions
 * 3. Make sure it works with Tpetra and Epetra->
 * 4. Write tests by transforming the SLEPc tests into Anasazi tests
 * 5. Add more Solvers
 * 6. Add tests for those solvers
 */
namespace TrilinosWrappers
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
    solve(const TrilinosWrappers::SparseMatrix       &A,
          std::vector<Number>                        &eigenvalues,
          std::vector<TrilinosWrappers::MPI::Vector> &eigenvectors,
          const unsigned int n_eigenpairs = 1); // implemented

    void
    solve(const TrilinosWrappers::SparseMatrix       &A,
          const TrilinosWrappers::SparseMatrix       &B,
          std::vector<Number>                        &eigenvalues,
          std::vector<TrilinosWrappers::MPI::Vector> &eigenvectors,
          const unsigned int n_eigenpairs = 1); // implemented

    void
    solve(const TrilinosWrappers::SparseMatrix       &A,
          const TrilinosWrappers::SparseMatrix       &B,
          std::vector<double>                        &real_eigenvalues,
          std::vector<double>                        &imag_eigenvalues,
          std::vector<TrilinosWrappers::MPI::Vector> &real_eigenvectors,
          std::vector<TrilinosWrappers::MPI::Vector> &imag_eigenvectors,
          const unsigned int n_eigenpairs = 1); // implemented

    void
    set_initial_space(const std::vector<TrilinosWrappers::MPI::Vector>
                        &initial_space); // started

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
          unsigned int       n_converged); // implemented

    void
    get_eigenpair(const unsigned int             index,
                  Number                        &eigenvalues,
                  TrilinosWrappers::MPI::Vector &eigenvectors); // implemented
    void
    get_eigenpair(
      const unsigned int             index,
      double                        &real_eigenvalues,
      double                        &imag_eigenvalues,
      TrilinosWrappers::MPI::Vector &real_eigenvector,
      TrilinosWrappers::MPI::Vector &
        imag_eigenvector); // started, but not completed. Extraction is annoying

    virtual void
    set_solver() = 0;

    void
    set_matrices(const TrilinosWrappers::SparseMatrix &A); // implemented

    void
    set_matrices(const TrilinosWrappers::SparseMatrix &A,
                 const TrilinosWrappers::SparseMatrix &B); // implemented

  protected:
    Teuchos::RCP<
      Anasazi::BasicEigenproblem<Number, Epetra_MultiVector, Epetra_Operator>>
      eigenproblem; // exists and used
    Anasazi::Eigensolution<Number, Epetra_MultiVector>
      eigensolution; // exists and used

    const Teuchos::RCP<Teuchos::ParameterList> parameters; // exists and used

    Teuchos::RCP<
      Anasazi::SolverManager<Number, Epetra_MultiVector, Epetra_Operator>>
      eigensolver; // exists and used
  };

  template <typename Number, typename MemorySpace>
  class SolverBlockKrylovSchur : public AnasaziSolverBase<Number, MemorySpace>
  {
    explicit SolverBlockKrylovSchur(
      SolverControl                              &cn,
      const Teuchos::RCP<Teuchos::ParameterList> &parameters_ =
        Teuchos::rcp(new Teuchos::ParameterList));

  protected:
    void
    set_solver() override;
  };

} // namespace TrilinosWrappers

#endif

DEAL_II_NAMESPACE_CLOSE

#endif
