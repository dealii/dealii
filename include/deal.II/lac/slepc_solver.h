// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2009 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#ifndef dealii_slepc_solver_h
#  define dealii_slepc_solver_h

#  include <deal.II/base/config.h>

#  ifdef DEAL_II_WITH_SLEPC

#    include <deal.II/lac/petsc_matrix_base.h>
#    include <deal.II/lac/slepc_spectral_transformation.h>
#    include <deal.II/lac/solver_control.h>

#    include <petscconf.h>
#    include <petscksp.h>

#    include <slepceps.h>

#    include <memory>

DEAL_II_NAMESPACE_OPEN

/**
 * Base namespace for solver classes using the SLEPc solvers which are
 * selected based on flags passed to the eigenvalue problem solver context.
 * Derived classes set the right flags to set the right solver.
 *
 * The SLEPc solvers are intended to be used for solving the generalized
 * eigenspectrum problem $(A-\lambda B)x=0$, for $x\neq0$; where $A$ is a
 * system matrix, $B$ is a @ref GlossMassMatrix "mass matrix", and $\lambda, x$ are a set of
 * eigenvalues and eigenvectors respectively. The emphasis is on methods and
 * techniques appropriate for problems in which the associated matrices are
 * sparse. Most of the methods offered by the SLEPc library are projection
 * methods or other methods with similar properties; and wrappers are provided
 * to interface to SLEPc solvers that handle both of these problem sets.
 *
 * SLEPcWrappers can be implemented in application codes in the following way:
 * @code
 * SolverControl solver_control (1000, 1e-9);
 * SolverArnoldi system (solver_control, mpi_communicator);
 * system.solve (A, B, lambda, x, size_of_spectrum);
 * @endcode
 * for the generalized eigenvalue problem $Ax=B\lambda x$, where the variable
 * <code>const unsigned int size_of_spectrum</code> tells SLEPc the number of
 * eigenvector/eigenvalue pairs to solve for. Additional options and solver
 * parameters can be passed to the SLEPc solvers before calling
 * <code>solve()</code>. For example, if the matrices of the general
 * eigenspectrum problem are not hermitian and the lower eigenvalues are
 * wanted only, the following code can be implemented before calling
 * <code>solve()</code>:
 * @code
 * system.set_problem_type (EPS_NHEP);
 * system.set_which_eigenpairs (EPS_SMALLEST_REAL);
 * @endcode
 * These options can also be set at the command line.
 *
 * See also <code>step-36</code> for a hands-on example.
 *
 * For cases when spectral transformations are used in conjunction with
 * Krylov-type solvers or Davidson-type eigensolvers are employed one can
 * additionally specify which linear solver and preconditioner to use. This
 * can be achieved as follows
 * @code
 * PETScWrappers::PreconditionBoomerAMG::AdditionalData data;
 * data.symmetric_operator = true;
 * PETScWrappers::PreconditionBoomerAMG preconditioner(mpi_communicator, data);
 * SolverControl linear_solver_control (dof_handler.n_dofs(),
 *                                      1e-12, false, false);
 * PETScWrappers::SolverCG linear_solver(linear_solver_control,
 *                                       mpi_communicator);
 * linear_solver.initialize(preconditioner);
 * SolverControl solver_control (100, 1e-9,false,false);
 * SLEPcWrappers::SolverKrylovSchur eigensolver(solver_control,
 *                                              mpi_communicator);
 * SLEPcWrappers::TransformationShift spectral_transformation(mpi_communicator);
 * spectral_transformation.set_solver(linear_solver);
 * eigensolver.set_transformation(spectral_transformation);
 * eigensolver.solve(stiffness_matrix, mass_matrix,
 *                   eigenvalues, eigenfunctions, eigenfunctions.size());
 * @endcode
 *
 * In order to support this usage case, different from PETSc wrappers, the
 * classes in this namespace are written in such a way that the underlying
 * SLEPc objects are initialized in constructors. By doing so one also avoid
 * caching of different settings (such as target eigenvalue or type of the
 * problem); instead those are applied straight away when the corresponding
 * functions of the wrapper classes are called.
 *
 * An alternative implementation to the one above is to use the API internals
 * directly within the application code. In this way the calling sequence
 * requires calling several of SolverBase functions rather than just one. This
 * freedom is intended for use of the SLEPcWrappers that require a greater
 * handle on the eigenvalue problem solver context. See also the API of, for
 * example:
 * @code
 * template <typename OutputVector>
 * void
 * SolverBase::solve (const PETScWrappers::MatrixBase &A,
 *                    const PETScWrappers::MatrixBase &B,
 *                    std::vector<PetscScalar>        &eigenvalues,
 *                    std::vector<OutputVector>       &eigenvectors,
 *                    const unsigned int               n_eigenpairs)
 * {
 *   ...
 * }
 * @endcode
 * as an example on how to do this.
 *
 * For further information and explanations on handling the
 * @ref SLEPcWrappers "SLEPcWrappers",
 * see also the
 * @ref PETScWrappers "PETScWrappers",
 * on which they depend.
 *
 * @ingroup SLEPcWrappers
 */
namespace SLEPcWrappers
{
  /**
   * Base class for solver classes using the SLEPc solvers. Since solvers in
   * SLEPc are selected based on flags passed to a generic solver object,
   * basically all the actual solver calls happen in this class, and derived
   * classes simply set the right flags to select one solver or another, or to
   * set certain parameters for individual solvers.
   *
   * For examples of how this and its derived classes can be used, including
   * how to provide preconditioners to the matrix of which eigenvalues are
   * to be computed, see the documentation of the SolverBase class as well
   * as the extensive discussions in the documentation of the SLEPcWrappers
   * namespace.
   */
  class SolverBase
  {
  public:
    /**
     * Constructor. Takes the MPI communicator over which parallel
     * computations are to happen.
     */
    SolverBase(SolverControl &cn, const MPI_Comm mpi_communicator);

    /**
     * Destructor.
     */
    virtual ~SolverBase();

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
    template <typename OutputVector>
    void
    solve(const PETScWrappers::MatrixBase &A,
          std::vector<PetscScalar>        &eigenvalues,
          std::vector<OutputVector>       &eigenvectors,
          const unsigned int               n_eigenpairs = 1);

    /**
     * Same as above, but here a composite method for solving the system $A
     * x=\lambda B x$, for real matrices, vectors, and values $A, B, x,
     * \lambda$.
     */
    template <typename OutputVector>
    void
    solve(const PETScWrappers::MatrixBase &A,
          const PETScWrappers::MatrixBase &B,
          std::vector<PetscScalar>        &eigenvalues,
          std::vector<OutputVector>       &eigenvectors,
          const unsigned int               n_eigenpairs = 1);

    /**
     * Same as above, but here a composite method for solving the system $A
     * x=\lambda B x$ with real matrices $A, B$ and imaginary eigenpairs $x,
     * \lambda$.
     */
    template <typename OutputVector>
    void
    solve(const PETScWrappers::MatrixBase &A,
          const PETScWrappers::MatrixBase &B,
          std::vector<double>             &real_eigenvalues,
          std::vector<double>             &imag_eigenvalues,
          std::vector<OutputVector>       &real_eigenvectors,
          std::vector<OutputVector>       &imag_eigenvectors,
          const unsigned int               n_eigenpairs = 1);

    /**
     * Set the initial vector space for the solver.
     *
     * By default, SLEPc initializes the starting vector or the initial
     * subspace randomly.
     */
    template <typename Vector>
    void
    set_initial_space(const std::vector<Vector> &initial_space);

    /**
     * Set the spectral transformation to be used.
     */
    void
    set_transformation(SLEPcWrappers::TransformationBase &this_transformation);

    /**
     * Set target eigenvalues in the spectrum to be computed. By default, no
     * target is set.
     */
    void
    set_target_eigenvalue(const PetscScalar &this_target);

    /**
     * Indicate which part of the spectrum is to be computed. By default
     * largest magnitude eigenvalues are computed.
     *
     * @note For other allowed values see the SLEPc documentation.
     */
    void
    set_which_eigenpairs(EPSWhich set_which);

    /**
     * Specify the type of the eigenspectrum problem. This can be used to
     * exploit known symmetries of the matrices that make up the
     * standard/generalized eigenspectrum problem.  By default a non-Hermitian
     * problem is assumed.
     *
     * @note For other allowed values see the SLEPc documentation.
     */
    void
    set_problem_type(EPSProblemType set_problem);

    /**
     * Take the information provided from SLEPc and checks it against
     * deal.II's own SolverControl objects to see if convergence has been
     * reached.
     */
    void
    get_solver_state(const SolverControl::State state);

    /**
     * Exception. Standard exception.
     */
    DeclException0(ExcSLEPcWrappersUsageError);

    /**
     * Exception. SLEPc error with error number.
     */
    DeclException1(ExcSLEPcError,
                   int,
                   << "    An error with error number " << arg1
                   << " occurred while calling a SLEPc function");

    /**
     * Exception. Convergence failure on the number of eigenvectors.
     */
    DeclException2(ExcSLEPcEigenvectorConvergenceMismatchError,
                   int,
                   int,
                   << "    The number of converged eigenvectors is " << arg1
                   << " but " << arg2 << " were requested. ");

    /**
     * Access to the object that controls convergence.
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
     * Copy of the MPI communicator object to be used for the solver.
     */
    const MPI_Comm mpi_communicator;

    /**
     * Solve the linear system for <code>n_eigenpairs</code> eigenstates.
     * Parameter <code>n_converged</code> contains the actual number of
     * eigenstates that have  converged; this can be both fewer or more than
     * n_eigenpairs, depending on the SLEPc eigensolver used.
     */
    void
    solve(const unsigned int n_eigenpairs, unsigned int *n_converged);

    /**
     * Access the real parts of solutions for a solved eigenvector problem,
     * pair index solutions, $\text{index}\,\in\,0\dots
     * \mathrm{n\_converged}-1$.
     */
    void
    get_eigenpair(const unsigned int         index,
                  PetscScalar               &eigenvalues,
                  PETScWrappers::VectorBase &eigenvectors);

    /**
     * Access the real and imaginary parts of solutions for a solved
     * eigenvector problem, pair index solutions, $\text{index}\,\in\,0\dots
     * \mathrm{n\_converged}-1$.
     */
    void
    get_eigenpair(const unsigned int         index,
                  double                    &real_eigenvalues,
                  double                    &imag_eigenvalues,
                  PETScWrappers::VectorBase &real_eigenvectors,
                  PETScWrappers::VectorBase &imag_eigenvectors);

    /**
     * Initialize solver for the linear system $Ax=\lambda x$. (Note: this is
     * required before calling solve ())
     */
    void
    set_matrices(const PETScWrappers::MatrixBase &A);

    /**
     * Same as above, but here initialize solver for the linear system $A
     * x=\lambda B x$.
     */
    void
    set_matrices(const PETScWrappers::MatrixBase &A,
                 const PETScWrappers::MatrixBase &B);

  protected:
    /**
     * Objects for Eigenvalue Problem Solver.
     */
    EPS eps;

  private:
    /**
     * Convergence reason.
     */
    EPSConvergedReason reason;


    /**
     * A function that can be used in SLEPc as a callback to check on
     * convergence.
     *
     * @note This function is not used currently.
     */
    static int
    convergence_test(EPS         eps,
                     PetscScalar real_eigenvalue,
                     PetscScalar imag_eigenvalue,
                     PetscReal   residual_norm,
                     PetscReal  *estimated_error,
                     void       *solver_control);
  };



  /**
   * An implementation of the solver interface using the SLEPc Krylov-Schur
   * solver. Usage: All spectrum, all problem types, complex.
   *
   * For examples of how this and its sibling classes can be used, including
   * how to provide preconditioners to the matrix of which eigenvalues are
   * to be computed, see the documentation of the SolverBase class as well
   * as the extensive discussions in the documentation of the SLEPcWrappers
   * namespace.
   *
   * @ingroup SLEPcWrappers
   */
  class SolverKrylovSchur : public SolverBase
  {
  public:
    /**
     * Standardized data struct to pipe additional data to the solver, should
     * it be needed.
     */
    struct AdditionalData
    {};

    /**
     * SLEPc solvers will want to have an MPI communicator context over which
     * computations are parallelized. By default, this carries the same
     * behavior as the PETScWrappers, but you can change that.
     */
    explicit SolverKrylovSchur(
      SolverControl        &cn,
      const MPI_Comm        mpi_communicator = PETSC_COMM_SELF,
      const AdditionalData &data             = AdditionalData());

  protected:
    /**
     * Store a copy of the flags for this particular solver.
     */
    const AdditionalData additional_data;
  };



  /**
   * An implementation of the solver interface using the SLEPc Arnoldi solver.
   * Usage: All spectrum, all problem types, complex.
   *
   * For examples of how this and its sibling classes can be used, including
   * how to provide preconditioners to the matrix of which eigenvalues are
   * to be computed, see the documentation of the SolverBase class as well
   * as the extensive discussions in the documentation of the SLEPcWrappers
   * namespace.
   *
   * @ingroup SLEPcWrappers
   */
  class SolverArnoldi : public SolverBase
  {
  public:
    /**
     * Standardized data struct to pipe additional data to the solver, should
     * it be needed.
     */
    struct AdditionalData
    {
      /**
       * Constructor. By default, set the option of delayed
       * reorthogonalization to false, i.e. don't do it.
       */
      explicit AdditionalData(const bool delayed_reorthogonalization = false);

      /**
       * Flag for delayed reorthogonalization.
       */
      bool delayed_reorthogonalization;
    };

    /**
     * SLEPc solvers will want to have an MPI communicator context over which
     * computations are parallelized. By default, this carries the same
     * behavior as the PETScWrappers, but you can change that.
     */
    explicit SolverArnoldi(SolverControl &cn,
                           const MPI_Comm mpi_communicator = PETSC_COMM_SELF,
                           const AdditionalData &data      = AdditionalData());

  protected:
    /**
     * Store a copy of the flags for this particular solver.
     */
    const AdditionalData additional_data;
  };



  /**
   * An implementation of the solver interface using the SLEPc Lanczos solver.
   * Usage: All spectrum, all problem types, complex.
   *
   * For examples of how this and its sibling classes can be used, including
   * how to provide preconditioners to the matrix of which eigenvalues are
   * to be computed, see the documentation of the SolverBase class as well
   * as the extensive discussions in the documentation of the SLEPcWrappers
   * namespace.
   *
   * @ingroup SLEPcWrappers
   */
  class SolverLanczos : public SolverBase
  {
  public:
    /**
     * Standardized data struct to pipe additional data to the solver, should
     * it be needed.
     */
    struct AdditionalData
    {
      /**
       * The type of reorthogonalization used during the Lanczos iteration.
       */
      EPSLanczosReorthogType reorthog;

      /**
       * Constructor. By default sets the type of reorthogonalization used
       * during the Lanczos iteration to full.
       */
      explicit AdditionalData(
        const EPSLanczosReorthogType r = EPS_LANCZOS_REORTHOG_FULL);
    };

    /**
     * SLEPc solvers will want to have an MPI communicator context over which
     * computations are parallelized. By default, this carries the same
     * behavior as the PETScWrappers, but you can change that.
     */
    explicit SolverLanczos(SolverControl &cn,
                           const MPI_Comm mpi_communicator = PETSC_COMM_SELF,
                           const AdditionalData &data      = AdditionalData());

  protected:
    /**
     * Store a copy of the flags for this particular solver.
     */
    const AdditionalData additional_data;
  };



  /**
   * An implementation of the solver interface using the SLEPc Power solver.
   * Usage: Largest values of spectrum only, all problem types, complex.
   *
   * For examples of how this and its sibling classes can be used, including
   * how to provide preconditioners to the matrix of which eigenvalues are
   * to be computed, see the documentation of the SolverBase class as well
   * as the extensive discussions in the documentation of the SLEPcWrappers
   * namespace.
   *
   * @ingroup SLEPcWrappers
   */
  class SolverPower : public SolverBase
  {
  public:
    /**
     * Standardized data struct to pipe additional data to the solver, should
     * it be needed.
     */
    struct AdditionalData
    {};

    /**
     * SLEPc solvers will want to have an MPI communicator context over which
     * computations are parallelized. By default, this carries the same
     * behavior as the PETScWrappers, but you can change that.
     */
    explicit SolverPower(SolverControl &cn,
                         const MPI_Comm mpi_communicator = PETSC_COMM_SELF,
                         const AdditionalData &data      = AdditionalData());

  protected:
    /**
     * Store a copy of the flags for this particular solver.
     */
    const AdditionalData additional_data;
  };



  /**
   * An implementation of the solver interface using the SLEPc Davidson
   * solver. Usage: All problem types.
   *
   * For examples of how this and its sibling classes can be used, including
   * how to provide preconditioners to the matrix of which eigenvalues are
   * to be computed, see the documentation of the SolverBase class as well
   * as the extensive discussions in the documentation of the SLEPcWrappers
   * namespace.
   *
   * @ingroup SLEPcWrappers
   */
  class SolverGeneralizedDavidson : public SolverBase
  {
  public:
    /**
     * Standardized data struct to pipe additional data to the solver, should
     * it be needed.
     */
    struct AdditionalData
    {
      /**
       * Use double expansion in search subspace.
       */
      bool double_expansion;

      /**
       * Constructor. By default set double_expansion to false.
       */
      explicit AdditionalData(bool double_expansion = false);
    };

    /**
     * SLEPc solvers will want to have an MPI communicator context over which
     * computations are parallelized. By default, this carries the same
     * behavior as the PETScWrappers, but you can change that.
     */
    explicit SolverGeneralizedDavidson(
      SolverControl        &cn,
      const MPI_Comm        mpi_communicator = PETSC_COMM_SELF,
      const AdditionalData &data             = AdditionalData());

  protected:
    /**
     * Store a copy of the flags for this particular solver.
     */
    const AdditionalData additional_data;
  };



  /**
   * An implementation of the solver interface using the SLEPc Jacobi-Davidson
   * solver. Usage: All problem types.
   *
   * For examples of how this and its sibling classes can be used, including
   * how to provide preconditioners to the matrix of which eigenvalues are
   * to be computed, see the documentation of the SolverBase class as well
   * as the extensive discussions in the documentation of the SLEPcWrappers
   * namespace.
   *
   * @ingroup SLEPcWrappers
   */
  class SolverJacobiDavidson : public SolverBase
  {
  public:
    /**
     * Standardized data struct to pipe additional data to the solver, should
     * it be needed.
     */
    struct AdditionalData
    {};

    /**
     * SLEPc solvers will want to have an MPI communicator context over which
     * computations are parallelized. By default, this carries the same
     * behavior as the PETScWrappers, but you can change that.
     */
    explicit SolverJacobiDavidson(
      SolverControl        &cn,
      const MPI_Comm        mpi_communicator = PETSC_COMM_SELF,
      const AdditionalData &data             = AdditionalData());

  protected:
    /**
     * Store a copy of the flags for this particular solver.
     */
    const AdditionalData additional_data;
  };



  /**
   * An implementation of the solver interface using the SLEPc LAPACK direct
   * solver.
   *
   * For examples of how this and its sibling classes can be used, including
   * how to provide preconditioners to the matrix of which eigenvalues are
   * to be computed, see the documentation of the SolverBase class as well
   * as the extensive discussions in the documentation of the SLEPcWrappers
   * namespace.
   *
   * @ingroup SLEPcWrappers
   */
  class SolverLAPACK : public SolverBase
  {
  public:
    /**
     * Standardized data struct to pipe additional data to the solver, should
     * it be needed.
     */
    struct AdditionalData
    {};

    /**
     * SLEPc solvers will want to have an MPI communicator context over which
     * computations are parallelized. By default, this carries the same
     * behavior as the PETScWrappers, but you can change that.
     */
    explicit SolverLAPACK(SolverControl &cn,
                          const MPI_Comm mpi_communicator = PETSC_COMM_SELF,
                          const AdditionalData &data      = AdditionalData());

  protected:
    /**
     * Store a copy of the flags for this particular solver.
     */
    const AdditionalData additional_data;
  };



  // --------------------------- inline and template functions -----------
  /**
   * This is declared here to make it possible to take a std::vector of
   * different PETScWrappers vector types
   */
  // todo: The logic of these functions can be simplified without breaking
  // backward compatibility...

  template <typename OutputVector>
  void
  SolverBase::solve(const PETScWrappers::MatrixBase &A,
                    std::vector<PetscScalar>        &eigenvalues,
                    std::vector<OutputVector>       &eigenvectors,
                    const unsigned int               n_eigenpairs)
  {
    // Panic if the number of eigenpairs wanted is out of bounds.
    AssertThrow((n_eigenpairs > 0) && (n_eigenpairs <= A.m()),
                ExcSLEPcWrappersUsageError());

    // Set the matrices of the problem
    set_matrices(A);

    // and solve
    unsigned int n_converged = 0;
    solve(n_eigenpairs, &n_converged);

    if (n_converged > n_eigenpairs)
      n_converged = n_eigenpairs;
    AssertThrow(n_converged == n_eigenpairs,
                ExcSLEPcEigenvectorConvergenceMismatchError(n_converged,
                                                            n_eigenpairs));

    AssertThrow(eigenvectors.size() != 0, ExcSLEPcWrappersUsageError());
    eigenvectors.resize(n_converged, eigenvectors.front());
    eigenvalues.resize(n_converged);

    for (unsigned int index = 0; index < n_converged; ++index)
      get_eigenpair(index, eigenvalues[index], eigenvectors[index]);
  }

  template <typename OutputVector>
  void
  SolverBase::solve(const PETScWrappers::MatrixBase &A,
                    const PETScWrappers::MatrixBase &B,
                    std::vector<PetscScalar>        &eigenvalues,
                    std::vector<OutputVector>       &eigenvectors,
                    const unsigned int               n_eigenpairs)
  {
    // Guard against incompatible matrix sizes:
    AssertThrow(A.m() == B.m(), ExcDimensionMismatch(A.m(), B.m()));
    AssertThrow(A.n() == B.n(), ExcDimensionMismatch(A.n(), B.n()));

    // Panic if the number of eigenpairs wanted is out of bounds.
    AssertThrow((n_eigenpairs > 0) && (n_eigenpairs <= A.m()),
                ExcSLEPcWrappersUsageError());

    // Set the matrices of the problem
    set_matrices(A, B);

    // and solve
    unsigned int n_converged = 0;
    solve(n_eigenpairs, &n_converged);

    if (n_converged >= n_eigenpairs)
      n_converged = n_eigenpairs;

    AssertThrow(n_converged == n_eigenpairs,
                ExcSLEPcEigenvectorConvergenceMismatchError(n_converged,
                                                            n_eigenpairs));
    AssertThrow(eigenvectors.size() != 0, ExcSLEPcWrappersUsageError());

    eigenvectors.resize(n_converged, eigenvectors.front());
    eigenvalues.resize(n_converged);

    for (unsigned int index = 0; index < n_converged; ++index)
      get_eigenpair(index, eigenvalues[index], eigenvectors[index]);
  }

  template <typename OutputVector>
  void
  SolverBase::solve(const PETScWrappers::MatrixBase &A,
                    const PETScWrappers::MatrixBase &B,
                    std::vector<double>             &real_eigenvalues,
                    std::vector<double>             &imag_eigenvalues,
                    std::vector<OutputVector>       &real_eigenvectors,
                    std::vector<OutputVector>       &imag_eigenvectors,
                    const unsigned int               n_eigenpairs)
  {
    // Guard against incompatible matrix sizes:
    AssertThrow(A.m() == B.m(), ExcDimensionMismatch(A.m(), B.m()));
    AssertThrow(A.n() == B.n(), ExcDimensionMismatch(A.n(), B.n()));

    // and incompatible eigenvalue/eigenvector sizes
    AssertThrow(real_eigenvalues.size() == imag_eigenvalues.size(),
                ExcDimensionMismatch(real_eigenvalues.size(),
                                     imag_eigenvalues.size()));
    AssertThrow(real_eigenvectors.size() == imag_eigenvectors.size(),
                ExcDimensionMismatch(real_eigenvectors.size(),
                                     imag_eigenvectors.size()));

    // Panic if the number of eigenpairs wanted is out of bounds.
    AssertThrow((n_eigenpairs > 0) && (n_eigenpairs <= A.m()),
                ExcSLEPcWrappersUsageError());

    // Set the matrices of the problem
    set_matrices(A, B);

    // and solve
    unsigned int n_converged = 0;
    solve(n_eigenpairs, &n_converged);

    if (n_converged >= n_eigenpairs)
      n_converged = n_eigenpairs;

    AssertThrow(n_converged == n_eigenpairs,
                ExcSLEPcEigenvectorConvergenceMismatchError(n_converged,
                                                            n_eigenpairs));
    AssertThrow((real_eigenvectors.size() != 0) &&
                  (imag_eigenvectors.size() != 0),
                ExcSLEPcWrappersUsageError());

    real_eigenvectors.resize(n_converged, real_eigenvectors.front());
    imag_eigenvectors.resize(n_converged, imag_eigenvectors.front());
    real_eigenvalues.resize(n_converged);
    imag_eigenvalues.resize(n_converged);

    for (unsigned int index = 0; index < n_converged; ++index)
      get_eigenpair(index,
                    real_eigenvalues[index],
                    imag_eigenvalues[index],
                    real_eigenvectors[index],
                    imag_eigenvectors[index]);
  }

  template <typename VectorType>
  void
  SolverBase::set_initial_space(
    const std::vector<VectorType> &this_initial_space)
  {
    std::vector<Vec> vecs(this_initial_space.size());

    for (unsigned int i = 0; i < this_initial_space.size(); ++i)
      {
        Assert(this_initial_space[i].l2_norm() > 0.0,
               ExcMessage("Initial vectors should be nonzero."));
        vecs[i] = this_initial_space[i];
      }

    // if the eigensolver supports only a single initial vector, but several
    // guesses are provided, then all except the first one will be discarded.
    // One could still build a vector that is rich in the directions of all
    // guesses, by taking a linear combination of them. (TODO: make function
    // virtual?)

    const PetscErrorCode ierr =
      EPSSetInitialSpace(eps, vecs.size(), vecs.data());
    AssertThrow(ierr == 0, ExcSLEPcError(ierr));
  }

} // namespace SLEPcWrappers

DEAL_II_NAMESPACE_CLOSE

#  else

// Make sure the scripts that create the C++20 module input files have
// something to latch on if the preprocessor #ifdef above would
// otherwise lead to an empty content of the file.
DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE

#  endif // DEAL_II_WITH_SLEPC

/*----------------------------   slepc_solver.h  ---------------------------*/

#endif

/*----------------------------   slepc_solver.h  ---------------------------*/
