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
 * 2. Write tests by transforming the SLEPc tests into Anasazi tests and writing
 * extra tests
 * 3. Generalize tests for the other solvers
 *  Useful resource:
 * https://github.com/t-sakashita/rokko/blob/ebd49e1198c4ec9e7612ad4a9806d16a4ff0bdc9/rokko/anasazi/solver.hpp
 */
namespace LinearAlgebra
{
  namespace TpetraWrappers
  {

    /**
     * Base class for all the solver classes using the Anasazi eigensolver
     * from Trilinos. Since all the different solvers are organized
     * in a similar way by specifying most information via the eigenproblem
     * and leaving the solver only to solve the problem, this class does most of
     * the heavy lifting.
     *
     * This class is designed very similar to the SLEPcWrappers::Solverbase
     * class, so we refer to its documentation for explanation for most of the
     * functionality. The main differences will be outlined below.
     *
     * With the exception of the desired tolerance and number of iterations,
     * most information for the solver, i.e. for which eigenvalues to solve or
     * the number of restarts, are passed via a
     * @code
     * Teuchos::RCP<Teuchos::ParameterList>.
     * @endcode
     *
     * A short demonstration with some of the default values.
     *
     * @code
     * SolverControl solver_control(1000, 1e-9);
     *
     * // create a Teuchos::ParameterList
     * Teuchos::RCP<Teuchos::ParameterList> parameters =  Teuchos::rcp(new
     * Teuchos::ParameterList);
     *
     * // Specify for which eigenvalues to solve
     * parameters->set("Which", "LM"); // set to the eigenvalues with the
     * largest magnitude
     *
     * // Sets the maximum number of restarts
     * parameters->set( "Maximum Restarts", 0);
     *
     * // employ locking of converged eigenvalues
     * parameters->set( "Use Locking", false);
     *
     * // specify the solver
     * LinearAlgebra::TpetraWrappers::AnasaziSolverGeneralizedDavidson(solver_control,
     * parameters);
     * @endcode
     *
     * For more details on which parameters to specify please see the Ansazi
     * documentation https://trilinos.github.io/docs/anasazi/index.html
     *
     */
    template <typename Number, typename MemorySpace = MemorySpace::Host>
    class AnasaziSolverBase
    {
    public:
      /**
       * Constructor. Takes in the SolverControl and a
       * Teuchos::RCP<Teuchos::ParameterList>.
       */
      AnasaziSolverBase(
        SolverControl                              &cn,
        const Teuchos::RCP<Teuchos::ParameterList> &parameters_ =
          Teuchos::rcp(new Teuchos::ParameterList));

      void



      /**
       * Composite method that solves the eigensystem $Ax=\lambda x$. The
       * eigenvector sent in has to have at least one element that we can use as
       * a template when resizing, since we do not know the parameters of the
       * specific vector class used (i.e. local_dofs for MPI vectors). However,
       * while copying eigenvectors, at least twice the memory size of
       * <tt>eigenvectors</tt> is being used (and can be more). To avoid doing
       * this, the fairly standard calling sequence executed here is used: Set
       * up matrices for solving; Actually solve the system; Gather the
       * solution(s).
       *
       * @note Note that the number of converged eigenvectors can be larger than
       * the number of eigenvectors requested; this is due to a round off error
       * (success) of the eigenproblem solver context. If this is found to be
       * the case we simply do not bother with more eigenpairs than requested,
       * but handle that it may be more than specified by ignoring any extras.
       * By default one eigenvector/eigenvalue pair is computed.
       */
      solve(
        const LinearAlgebra::TpetraWrappers::SparseMatrix<Number, MemorySpace>
                            &A,
        std::vector<Number> &eigenvalues,
        std::vector<LinearAlgebra::TpetraWrappers::Vector<Number, MemorySpace>>
                          &eigenvectors,
        const unsigned int n_eigenpairs = 1);



      /**
       * Same as above, but here a composite method for solving the system $A
       * x=\lambda B x$ with real matrices $A, B$ and imaginary eigenpairs $x,
       * \lambda$.
       */
      void
      solve(
        const LinearAlgebra::TpetraWrappers::SparseMatrix<Number, MemorySpace>
          &A,
        const LinearAlgebra::TpetraWrappers::SparseMatrix<Number, MemorySpace>
                            &B,
        std::vector<Number> &eigenvalues,
        std::vector<LinearAlgebra::TpetraWrappers::Vector<Number, MemorySpace>>
                          &eigenvectors,
        const unsigned int n_eigenpairs = 1);



      /**
       * Set the initial vector space for the solver.
       *
       * Anasazi does not create an initial space by default,
       * so one has to either manually specify an initial space
       * or generate a random initial space. See below for the random initial
       * space.
       */
      void
      set_initial_space(
        const std::vector<
          LinearAlgebra::TpetraWrappers::Vector<Number, MemorySpace>>
          &initial_space); // done



      /**
       * Generates a random initial space for the solver.
       *
       * Anasazi does not generate a default initial space.
       * This function generates a random initial space for ease of use.
       */
      void
      generate_random_initial_space(
        const LinearAlgebra::TpetraWrappers::SparseMatrix<Number, MemorySpace>
          &A);



      /**
       * Set the preconditioner of the problem.
       * It should be noted that only certain eigensolvers allow for the use of
       * preconditioners.
       *
       * The eigensolvers that allow preconditioners are:
       * @ref AnasaziSolverBlockDavidson "AnasaziSolverBlockDavidson",
       * @ref AnasaziSolverGeneralizedDavidson "AnasaziSolverGeneralizedDavidson" and
       * @ref AnasaziSolverLOBPCG "AnasaziSolverLOBPCG".
       *
       */
      void
      set_preconditioner(
        LinearAlgebra::TpetraWrappers::PreconditionBase<Number, MemorySpace>
          &Preconditioner);



      /**
       * Tells the solver that the eigenproblem is hermitian
       * or in the real case symmetric.
       */
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
       * It should be noted that as of right now the SolverControl object is
       * only used for setting the tolerance and the maximum number of
       * iterations of the simulation.
       */
      SolverControl &
      control() const;

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
       */
      void
      solve(const unsigned int n_eigenpairs, unsigned int &n_converged);



      /**
       * Access the eigenvectors and eigenvalues of the solved eigenvalue
       * problem.
       * This function differs from SLEPcWrappers::Solverbase as it will not
       * return the real part, but both the imaginary part and the real part if
       * the number type was specified to be  @code std::complex<> @endcode.
       * Otherwise it will simply return the real part.
       */
      void
      get_eigenpair(const unsigned int index,
                    Number            &eigenvalues,
                    LinearAlgebra::TpetraWrappers::Vector<Number, MemorySpace>
                      &eigenvectors);



      /**
       * Virtual function that is used to create the solvermanager.
       */
      virtual void
      set_solver() = 0;

      /**
       * Initialize solver for the linear system $Ax=\lambda x$. (Note: this is
       * required before calling solve ())
       */
      void
      set_matrices(
        const LinearAlgebra::TpetraWrappers::SparseMatrix<Number, MemorySpace>
          &A);



      /**
       * Same as above, but here initialize solver for the linear system $A
       * x=\lambda B x$.
       */
      void
      set_matrices(
        const LinearAlgebra::TpetraWrappers::SparseMatrix<Number, MemorySpace>
          &A,
        const LinearAlgebra::TpetraWrappers::SparseMatrix<Number, MemorySpace>
          &B);



    protected:
      /**
       * The eigenproblem object used to specify the eigenproblem and
       * its properties.
       */
      Teuchos::RCP<Anasazi::BasicEigenproblem<Number,
                                              Tpetra::MultiVector<Number>,
                                              Tpetra::Operator<Number>>>
        eigenproblem;


      /**
       * The solution object of the eigenproblem.
       */
      Anasazi::Eigensolution<Number, Tpetra::MultiVector<Number>> eigensolution;



      /**
       * Additional data with advanced settings.
       */
      const Teuchos::RCP<Teuchos::ParameterList> parameters;


      /**
       * The solver manager used to solve the eigenproblem.
       */
      Teuchos::RCP<Anasazi::SolverManager<Number,
                                          Tpetra::MultiVector<Number>,
                                          Tpetra::Operator<Number>>>
        eigensolver;
    };


    /**
     * The Block-Krylov solver class.
     * This solver can handle the most general types of matrices, including
     * non-hermitian matrices. This solver cannot be used together with a
     * preconditioner.
     */
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



    /**
     * The BlockDavidson solver class.
     * This solver is designed for hermitian (or in the real case symmetric)
     * problems, but is not designed for non-hermitian problems. Similar to CG
     * not being a solver designed for non-symmetric, non-positive problems.
     * This solver can be used together with a preconditioner.
     */
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



    /**
     * The GeneralizedDavidson solver class.
     * This solver is designed to handle non-hermitian problems, but
     * but is not designed for complex hermitian problems.
     * This solver can be used together with a preconditioner.
     */
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



    /**
     * The LOBPCG(Locally Optimal Block Preconditioned Conjugate Gradient)
     * solver class. This solver is designed to handle complex hermitian
     * problems, but is not designed for non-hermitian problems. As the name
     * suggests, this solver can be used together with a preconditioner.
     */
    template <typename Number, typename MemorySpace = MemorySpace::Host>
    class AnasaziSolverLOBPCG : public AnasaziSolverBase<Number, MemorySpace>
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
