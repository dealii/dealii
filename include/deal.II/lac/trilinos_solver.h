// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2008 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_trilinos_solver_h
#define dealii_trilinos_solver_h


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_TRILINOS

#  include <deal.II/base/template_constraints.h>

#  include <deal.II/lac/exceptions.h>
#  include <deal.II/lac/la_parallel_vector.h>
#  include <deal.II/lac/solver_control.h>
#  include <deal.II/lac/vector.h>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS

// for AztecOO solvers
#  include <Amesos.h>
#  include <AztecOO.h>
#  include <Epetra_LinearProblem.h>
#  include <Epetra_Operator.h>

// for Belos solvers
#  ifdef DEAL_II_TRILINOS_WITH_BELOS
#    include <BelosBlockCGSolMgr.hpp>
#    include <BelosBlockGmresSolMgr.hpp>
#    include <BelosEpetraAdapter.hpp>
#    include <BelosIteration.hpp>
#    include <BelosMultiVec.hpp>
#    include <BelosOperator.hpp>
#    include <BelosSolverManager.hpp>
#  endif

DEAL_II_ENABLE_EXTRA_DIAGNOSTICS


#  include <memory>


DEAL_II_NAMESPACE_OPEN

namespace TrilinosWrappers
{
  // forward declarations
#  ifndef DOXYGEN
  class SparseMatrix;
  class PreconditionBase;
#  endif


  /**
   * Base class for solver classes using the Trilinos solvers. Since solvers
   * in Trilinos are selected based on flags passed to a generic solver
   * object, basically all the actual solver calls happen in this class, and
   * derived classes simply set the right flags to select one solver or
   * another, or to set certain parameters for individual solvers. For a
   * general discussion on the Trilinos solver package AztecOO, we refer to
   * the <a
   * href="https://trilinos.org/docs/dev/packages/aztecoo/doc/html/index.html">AztecOO
   * user guide</a>.
   *
   * This solver class can also be used as a standalone class, where the
   * respective Krylov method is set via the flag <tt>solver_name</tt>. This
   * can be done at runtime (e.g., when parsing the solver from a
   * ParameterList) and is similar to the deal.II class SolverSelector.
   *
   * @ingroup TrilinosWrappers
   */
  class SolverBase
  {
  public:
    /**
     * Enumeration object that is set in the constructor of the derived
     * classes and tells Trilinos which solver to use. This option can also be
     * set in the user program, so one might use this base class instead of
     * one of the specialized derived classes when the solver should be set at
     * runtime. Currently enabled options are:
     */
    enum SolverName
    {
      /**
       * Use the conjugate gradient (CG) algorithm.
       */
      cg,
      /**
       * Use the conjugate gradient squared (CGS) algorithm.
       */
      cgs,
      /**
       * Use the generalized minimum residual (GMRES) algorithm.
       */
      gmres,
      /**
       * Use the biconjugate gradient stabilized (BICGStab) algorithm.
       */
      bicgstab,
      /**
       * Use the transpose-free quasi-minimal residual (TFQMR) method.
       */
      tfqmr
    } solver_name;

    /**
     * Standardized data struct to pipe additional data to the solver.
     */

    struct AdditionalData
    {
      /**
       * Set the additional data field to the desired output format and puts
       * the restart parameter in case the derived class is GMRES.
       *
       * TODO: Find a better way for setting the GMRES restart parameter since
       * it is quite inelegant to set a specific option of one solver in the
       * base class for all solvers.
       */
      explicit AdditionalData(const bool         output_solver_details = false,
                              const unsigned int gmres_restart_parameter = 30);

      /**
       * Enables/disables the output of solver details (residual in each
       * iterations etc.).
       */
      const bool output_solver_details;

      /**
       * Restart parameter for GMRES solver.
       */
      const unsigned int gmres_restart_parameter;
    };

    /**
     * Constructor. Takes the solver control object and creates the solver.
     */
    SolverBase(SolverControl        &cn,
               const AdditionalData &data = AdditionalData());

    /**
     * Second constructor. This constructor takes an enum object that
     * specifies the solver name and sets the appropriate Krylov method.
     */
    SolverBase(const enum SolverName solver_name,
               SolverControl        &cn,
               const AdditionalData &data = AdditionalData());

    /**
     * Destructor.
     */
    virtual ~SolverBase() = default;

    /**
     * Solve the linear system <tt>Ax=b</tt>. Depending on the information
     * provided by derived classes and the object passed as a preconditioner,
     * one of the linear solvers and preconditioners of Trilinos is chosen.
     */
    void
    solve(const SparseMatrix     &A,
          MPI::Vector            &x,
          const MPI::Vector      &b,
          const PreconditionBase &preconditioner);

    /**
     * Solve the linear system <tt>Ax=b</tt> where <tt>A</tt> is an operator.
     * This function can be used for matrix free computation. Depending on the
     * information provided by derived classes and the object passed as a
     * preconditioner, one of the linear solvers and preconditioners of
     * Trilinos is chosen.
     */
    void
    solve(const Epetra_Operator  &A,
          MPI::Vector            &x,
          const MPI::Vector      &b,
          const PreconditionBase &preconditioner);

    /**
     * Solve the linear system <tt>Ax=b</tt> where both <tt>A</tt> and its
     * @p preconditioner are an operator.
     * This function can be used when both <tt>A</tt> and the @p preconditioner
     * are LinearOperators derived from a TrilinosPayload.
     * Depending on the information provided by derived classes and the object
     * passed as a preconditioner, one of the linear solvers and preconditioners
     * of Trilinos is chosen.
     */
    void
    solve(const Epetra_Operator &A,
          MPI::Vector           &x,
          const MPI::Vector     &b,
          const Epetra_Operator &preconditioner);

    /**
     * Solve the linear system <tt>Ax=b</tt> where <tt>A</tt> is an operator,
     * and the vectors @p x and @p b are native Trilinos vector types.
     * This function can be used when <tt>A</tt> is a LinearOperators derived
     * from a TrilinosPayload.
     * Depending on the information provided by derived classes and the object
     * passed as a preconditioner, one of the linear solvers and preconditioners
     * of Trilinos is chosen.
     */
    void
    solve(const Epetra_Operator    &A,
          Epetra_MultiVector       &x,
          const Epetra_MultiVector &b,
          const PreconditionBase   &preconditioner);

    /**
     * Solve the linear system <tt>Ax=b</tt> where both <tt>A</tt> and its
     * @p preconditioner are an operator, and the vectors @p x and @p b are
     * native Trilinos vector types.
     * This function can be used when both <tt>A</tt> and the @p preconditioner
     * are LinearOperators derived from a TrilinosPayload.
     * Depending on the information provided by derived classes and the object
     * passed as a preconditioner, one of the linear solvers and preconditioners
     * of Trilinos is chosen.
     */
    void
    solve(const Epetra_Operator    &A,
          Epetra_MultiVector       &x,
          const Epetra_MultiVector &b,
          const Epetra_Operator    &preconditioner);



    /**
     * Solve the linear system <tt>Ax=b</tt>. Depending on the information
     * provided by derived classes and the object passed as a preconditioner,
     * one of the linear solvers and preconditioners of Trilinos is chosen.
     * This class works with matrices according to the TrilinosWrappers
     * format, but can take deal.II vectors as argument. Since deal.II are
     * serial vectors (not distributed), this function does only what you
     * expect in case the matrix is locally owned. Otherwise, an exception
     * will be thrown.
     */
    void
    solve(const SparseMatrix           &A,
          dealii::Vector<double>       &x,
          const dealii::Vector<double> &b,
          const PreconditionBase       &preconditioner);

    /**
     * Solve the linear system <tt>Ax=b</tt> where <tt>A</tt> is an operator.
     * This function can be used for matrix free computations. Depending on
     * the information provided by derived classes and the object passed as a
     * preconditioner, one of the linear solvers and preconditioners of
     * Trilinos is chosen. This class works with matrices according to the
     * TrilinosWrappers format, but can take deal.II vectors as argument.
     * Since deal.II are serial vectors (not distributed), this function does
     * only what you expect in case the matrix is locally owned. Otherwise, an
     * exception will be thrown.
     */
    void
    solve(Epetra_Operator              &A,
          dealii::Vector<double>       &x,
          const dealii::Vector<double> &b,
          const PreconditionBase       &preconditioner);

    /**
     * Solve the linear system <tt>Ax=b</tt> for deal.II's parallel
     * distributed vectors. Depending on the information provided by derived
     * classes and the object passed as a preconditioner, one of the linear
     * solvers and preconditioners of Trilinos is chosen.
     */
    void
    solve(const SparseMatrix                                       &A,
          dealii::LinearAlgebra::distributed::Vector<double>       &x,
          const dealii::LinearAlgebra::distributed::Vector<double> &b,
          const PreconditionBase &preconditioner);

    /**
     * Solve the linear system <tt>Ax=b</tt> where <tt>A</tt> is an operator.
     * This function can be used for matrix free computation. Depending on the
     * information provided by derived classes and the object passed as a
     * preconditioner, one of the linear solvers and preconditioners of
     * Trilinos is chosen.
     */
    void
    solve(Epetra_Operator                                          &A,
          dealii::LinearAlgebra::distributed::Vector<double>       &x,
          const dealii::LinearAlgebra::distributed::Vector<double> &b,
          const PreconditionBase &preconditioner);


    /**
     * Access to object that controls convergence.
     */
    SolverControl &
    control() const;

    /**
     * Exception
     */
    DeclException1(ExcTrilinosError,
                   int,
                   << "An error with error number " << arg1
                   << " occurred while calling a Trilinos function");

  protected:
    /**
     * Reference to the object that controls convergence of the iterative
     * solver. In fact, for these Trilinos wrappers, Trilinos does so itself,
     * but we copy the data from this object before starting the solution
     * process, and copy the data back into it afterwards.
     */
    SolverControl &solver_control;

  private:
    /**
     * The solve function is used to set properly the Epetra_LinearProblem,
     * once it is done this function solves the linear problem.
     */
    template <typename Preconditioner>
    void
    do_solve(const Preconditioner &preconditioner);

    /**
     * A function that sets the preconditioner that the solver will apply
     */
    template <typename Preconditioner>
    void
    set_preconditioner(AztecOO &solver, const Preconditioner &preconditioner);

    /**
     * A structure that collects the Trilinos sparse matrix, the right hand
     * side vector and the solution vector, which is passed down to the
     * Trilinos solver.
     */
    std::unique_ptr<Epetra_LinearProblem> linear_problem;

    /**
     * A structure that contains a Trilinos object that can query the linear
     * solver and determine whether the convergence criterion have been met.
     */
    std::unique_ptr<AztecOO_StatusTest> status_test;

    /**
     * A structure that contains the Trilinos solver and preconditioner
     * objects.
     */
    AztecOO solver;

    /**
     * Store a copy of the flags for this particular solver.
     */
    const AdditionalData additional_data;
  };


  // provide a declaration for two explicit specializations
  template <>
  void
  SolverBase::set_preconditioner(AztecOO                &solver,
                                 const PreconditionBase &preconditioner);

  template <>
  void
  SolverBase::set_preconditioner(AztecOO               &solver,
                                 const Epetra_Operator &preconditioner);


  /**
   * An implementation of the solver interface using the Trilinos CG solver.
   *
   * @ingroup TrilinosWrappers
   */
  class SolverCG : public SolverBase
  {
  public:
    /**
     * Constructor. In contrast to deal.II's own solvers, there is no need to
     * give a vector memory object.
     *
     * The last argument takes a structure with additional, solver dependent
     * flags for tuning.
     */
    SolverCG(SolverControl &cn, const AdditionalData &data = AdditionalData());
  };



  /**
   * An implementation of the solver interface using the Trilinos CGS solver.
   *
   * @ingroup TrilinosWrappers
   */
  class SolverCGS : public SolverBase
  {
  public:
    /**
     * Constructor. In contrast to deal.II's own solvers, there is no need to
     * give a vector memory object.
     *
     * The last argument takes a structure with additional, solver dependent
     * flags for tuning.
     */
    SolverCGS(SolverControl &cn, const AdditionalData &data = AdditionalData());
  };



  /**
   * An implementation of the solver interface using the Trilinos GMRES
   * solver.
   */
  class SolverGMRES : public SolverBase
  {
  public:
    /**
     * Constructor. In contrast to deal.II's own solvers, there is no need to
     * give a vector memory object.
     *
     * The last argument takes a structure with additional, solver dependent
     * flags for tuning.
     */
    SolverGMRES(SolverControl        &cn,
                const AdditionalData &data = AdditionalData());
  };



  /**
   * An implementation of the solver interface using the Trilinos BiCGStab
   * solver.
   *
   * @ingroup TrilinosWrappers
   */
  class SolverBicgstab : public SolverBase
  {
  public:
    /**
     * Constructor. In contrast to deal.II's own solvers, there is no need to
     * give a vector memory object.
     *
     * The last argument takes a structure with additional, solver dependent
     * flags for tuning.
     */
    SolverBicgstab(SolverControl        &cn,
                   const AdditionalData &data = AdditionalData());
  };



  /**
   * An implementation of the solver interface using the Trilinos TFQMR
   * solver.
   *
   * @ingroup TrilinosWrappers
   */
  class SolverTFQMR : public SolverBase
  {
  public:
    /**
     * Constructor. In contrast to deal.II's own solvers, there is no need to
     * give a vector memory object.
     *
     * The last argument takes a structure with additional, solver dependent
     * flags for tuning.
     */
    SolverTFQMR(SolverControl        &cn,
                const AdditionalData &data = AdditionalData());
  };



  /**
   * An implementation of Trilinos direct solvers (using the Amesos package).
   * The data field AdditionalData::solver_type can be used to specify the
   * type of solver. It allows the use of built-in solvers Amesos_Klu as well
   * as third-party solvers Amesos_Superludist or Amesos_Mumps.
   *
   * For instructions on how to install Trilinos for use with direct solvers
   * other than KLU, see the link to the Trilinos installation instructions
   * linked to from the deal.II ReadMe file.
   *
   * @ingroup TrilinosWrappers
   */
  class SolverDirect : public EnableObserverPointer
  {
  public:
    /**
     * Standardized data struct to pipe additional data to the solver.
     */

    struct AdditionalData
    {
      /**
       * Set the additional data field to the desired output format.
       */
      explicit AdditionalData(const bool         output_solver_details = false,
                              const std::string &solver_type = "Amesos_Klu");

      /**
       * Enables/disables the output of solver details (residual in each
       * iterations etc.).
       */
      bool output_solver_details;

      /**
       * Set the solver type (for third party solver support of Trilinos
       * Amesos package). Possibilities are:
       * <ul>
       * <li>  "Amesos_Lapack" </li>
       * <li>  "Amesos_Scalapack" </li>
       * <li>  "Amesos_Klu" </li>
       * <li>  "Amesos_Umfpack" </li>
       * <li>  "Amesos_Pardiso" </li>
       * <li>  "Amesos_Taucs" </li>
       * <li>  "Amesos_Superlu" </li>
       * <li>  "Amesos_Superludist" </li>
       * <li>  "Amesos_Dscpack" </li>
       * <li>  "Amesos_Mumps" </li>
       * </ul>
       * Note that the availability of these solvers in deal.II depends on
       * which solvers were set when configuring Trilinos.
       */
      std::string solver_type;
    };

    /**
     * Constructor. Creates the solver without solver control object.
     */
    explicit SolverDirect(const AdditionalData &data = AdditionalData());

    /**
     * Constructor. Takes the solver control object and creates the solver.
     */
    SolverDirect(SolverControl        &cn,
                 const AdditionalData &data = AdditionalData());

    /**
     * Destructor.
     */
    virtual ~SolverDirect() = default;

    /**
     * Initializes the direct solver for the matrix <tt>A</tt> and creates a
     * factorization for it with the package chosen from the additional
     * data structure. Note that there is no need for a preconditioner
     * here and solve() is not called.
     */
    void
    initialize(const SparseMatrix &A);

    /**
     * Initializes the direct solver for the matrix <tt>A</tt> and creates a
     * factorization for it with the package chosen from the additional
     * data structure. Note that there is no need for a preconditioner
     * here and solve() is not called. Furthermore, @p data replaces the
     * data stored in this instance.
     */
    void
    initialize(const SparseMatrix &A, const AdditionalData &data);

    /**
     * Solve the linear system <tt>Ax=b</tt> based on the
     * package set in the constructor on initialize(). Note the matrix is not
     * refactored during this call.
     */
    void
    solve(MPI::Vector &x, const MPI::Vector &b);

    /**
     * Solve the linear system <tt>Ax=b</tt> based on the package set in
     * initialize() for deal.II's own parallel vectors. Note the matrix is not
     * refactored during this call.
     */
    void
    solve(dealii::LinearAlgebra::distributed::Vector<double>       &x,
          const dealii::LinearAlgebra::distributed::Vector<double> &b);

    /**
     * Solve the linear system <tt>Ax=b</tt> based on the
     * package set in the constructor or initialize(). Note the matrix is not
     * refactored during this call.
     */
    void
    vmult(MPI::Vector &x, const MPI::Vector &b) const;

    /**
     * Solve the linear system <tt>Ax=b</tt> based on the package set in
     * initialize() for deal.II's own parallel vectors. Note the matrix is not
     * refactored during this call.
     */
    void
    vmult(dealii::LinearAlgebra::distributed::Vector<double>       &x,
          const dealii::LinearAlgebra::distributed::Vector<double> &b) const;

    /**
     * Solve the linear system <tt>Ax=b</tt>. Creates a factorization of the
     * matrix with the package chosen from the additional data structure and
     * performs the solve. Note that there is no need for a preconditioner
     * here.
     */
    void
    solve(const SparseMatrix &A, MPI::Vector &x, const MPI::Vector &b);

    /**
     * Solve the linear system <tt>Ax=b</tt>. This class works with Trilinos
     * matrices, but takes deal.II serial vectors as argument. Since these
     * vectors are not distributed, this function does only what you expect in
     * case the matrix is serial (i.e., locally owned). Otherwise, an
     * exception will be thrown.
     */
    void
    solve(const SparseMatrix           &A,
          dealii::Vector<double>       &x,
          const dealii::Vector<double> &b);

    /**
     * Solve the linear system <tt>Ax=b</tt> for deal.II's own parallel
     * vectors. Creates a factorization of the matrix with the package chosen
     * from the additional data structure and performs the solve. Note that
     * there is no need for a preconditioner here.
     */
    void
    solve(const SparseMatrix                                       &A,
          dealii::LinearAlgebra::distributed::Vector<double>       &x,
          const dealii::LinearAlgebra::distributed::Vector<double> &b);

    /**
     * Solve the linear system <tt>Ax=b</tt> where A is an operator,
     * and the vectors x and b are native Trilinos vector types.
     * This function can be used when A is a LinearOperators derived
     * from a TrilinosPayload.
     */
    void
    solve(const Epetra_Operator    &A,
          Epetra_MultiVector       &x,
          const Epetra_MultiVector &b);

    /**
     * Solve the linear system <tt>AX=B</tt> where A is an operator,
     * and the X and B are FullMatrix types. The matrices are
     * converted to multi vector and the resulting problem is solved
     * in a batched way.
     */
    void
    solve(const SparseMatrix       &sparse_matrix,
          FullMatrix<double>       &solution,
          const FullMatrix<double> &rhs);

    /**
     * Access to object that controls convergence.
     */
    SolverControl &
    control() const;

    /**
     * Exception
     */
    DeclException1(ExcTrilinosError,
                   int,
                   << "An error with error number " << arg1
                   << " occurred while calling a Trilinos function");

  private:
    /**
     * Actually performs the operations for solving the linear system,
     * including the factorization and forward and backward substitution.
     */
    void
    do_solve();

    /**
     * Local dummy solver control object.
     */
    SolverControl solver_control_own;

    /**
     * Reference to the object that controls convergence of the iterative
     * solver. In fact, for these Trilinos wrappers, Trilinos does so itself,
     * but we copy the data from this object before starting the solution
     * process, and copy the data back into it afterwards.
     */
    SolverControl &solver_control;

    /**
     * A structure that collects the Trilinos sparse matrix, the right hand
     * side vector and the solution vector, which is passed down to the
     * Trilinos solver.
     */
    std::unique_ptr<Epetra_LinearProblem> linear_problem;

    /**
     * A structure that contains the Trilinos solver and preconditioner
     * objects.
     */
    std::unique_ptr<Amesos_BaseSolver> solver;

    /**
     * Store a copy of the flags for this particular solver.
     */
    AdditionalData additional_data;
  };



#  ifdef DEAL_II_TRILINOS_WITH_BELOS
  /**
   * Wrapper around the iterate solver package from the Belos
   * package
   * (https://docs.trilinos.org/latest-release/packages/belos/doc/html/index.html),
   * targeting deal.II data structures.
   */
  template <typename VectorType>
  DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
  class SolverBelos
  {
  public:
    /**
     * Enumeration object that Trilinos which solver to use.
     * Currently enabled options are:
     */
    enum SolverName
    {
      /**
       * Use the conjugate gradient (CG) algorithm.
       */
      cg,
      /**
       * Use the generalized minimum residual (GMRES) algorithm.
       */
      gmres,
    } solver_name;

    /**
     * Struct that helps to configure SolverBelos. More advanced
     * parameters are passed to the constructor SolverBelos
     * directly via a Teuchos::ParameterList.
     */
    struct AdditionalData
    {
      /**
       * Constructor.
       */
      AdditionalData(const SolverName solver_name           = SolverName::cg,
                     const bool       right_preconditioning = false)
        : solver_name(solver_name)
        , right_preconditioning(right_preconditioning)
      {}

      /**
       * Solver type;
       */
      SolverName solver_name;

      /**
       * Flag for right preconditioning.
       */
      bool right_preconditioning;
    };

    /**
     * Constructor.
     */
    SolverBelos(SolverControl                              &solver_control,
                const AdditionalData                       &additional_data,
                const Teuchos::RCP<Teuchos::ParameterList> &belos_parameters);

    /**
     * Solve the linear system <tt>Ax=b</tt> with a given preconditioner.
     */
    template <typename OperatorType, typename PreconditionerType>
    void
    solve(const OperatorType       &a,
          VectorType               &x,
          const VectorType         &b,
          const PreconditionerType &p);

  private:
    SolverControl                              &solver_control;
    const AdditionalData                        additional_data;
    const Teuchos::RCP<Teuchos::ParameterList> &belos_parameters;
  };
#  endif

} // namespace TrilinosWrappers



#  ifndef DOXYGEN

#    ifdef DEAL_II_TRILINOS_WITH_BELOS
namespace TrilinosWrappers
{
  namespace internal
  {
    /**
     * Implementation of the abstract interface
     * Belos::MultiVec for deal.II vectors. For details,
     * see
     * https://docs.trilinos.org/latest-release/packages/belos/doc/html/classBelos_1_1MultiVec.html.
     */
    template <typename VectorType>
    DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
    class MultiVecWrapper
      : public Belos::MultiVec<typename VectorType::value_type>
    {
    public:
      /**
       * Underlying value type.
       */
      using value_type = typename VectorType::value_type;

      /**
       * Indicate that specialization exists.
       */
      static bool
      this_type_is_missing_a_specialization()
      {
        return false;
      }

      /**
       * Constructor that takes a mutable deal.II vector.
       */
      MultiVecWrapper(VectorType &vector)
      {
        this->vectors.resize(1);
        this->vectors[0].reset(
          &vector,
          [](auto *) { /*Nothing to do, since vector is owned outside.*/ });
      }

      /**
       * Constructor that takes a const deal.II vector.
       */
      MultiVecWrapper(const VectorType &vector)
      {
        this->vectors.resize(1);
        this->vectors[0].reset(
          &const_cast<VectorType &>(vector),
          [](auto *) { /*Nothing to do, since vector is owned outside.*/ });
      }

      /**
       * Destructor.
       */
      virtual ~MultiVecWrapper() = default;

      /**
       * Create a new MultiVec with numvecs columns.
       */
      virtual Belos::MultiVec<value_type> *
      Clone(const int numvecs) const
      {
        auto new_multi_vec = new MultiVecWrapper<VectorType>;

        new_multi_vec->vectors.resize(numvecs);

        for (auto &vec : new_multi_vec->vectors)
          {
            vec = std::make_shared<VectorType>();

            AssertThrow(this->vectors.size() > 0, ExcInternalError());
            vec->reinit(*this->vectors[0]);
          }

        return new_multi_vec;
      }

      /**
       * Create a new MultiVec and copy contents of *this into it (deep copy).
       */
      virtual Belos::MultiVec<value_type> *
      CloneCopy() const
      {
        AssertThrow(false, ExcNotImplemented());
      }

      /**
       * Creates a new Belos::MultiVec and copies the selected contents of
       * *this into the new multivector (deep copy). The copied vectors
       * from *this are indicated by the index.size() indices in index.
       */
      virtual Belos::MultiVec<value_type> *
      CloneCopy(const std::vector<int> &index) const
      {
        auto new_multi_vec = new MultiVecWrapper<VectorType>;

        new_multi_vec->vectors.resize(index.size());

        for (unsigned int i = 0; i < index.size(); ++i)
          {
            AssertThrow(static_cast<unsigned int>(index[i]) <
                          this->vectors.size(),
                        ExcInternalError());

            new_multi_vec->vectors[i] = std::make_shared<VectorType>();

            AssertIndexRange(index[i], this->vectors.size());
            *new_multi_vec->vectors[i] = *this->vectors[index[i]];
          }

        return new_multi_vec;
      }

      /**
       * Creates a new Belos::MultiVec that shares the selected contents
       * of *this. The index of the numvecs vectors copied from *this
       * are indicated by the indices given in index.
       */
      virtual Belos::MultiVec<value_type> *
      CloneViewNonConst(const std::vector<int> &index)
      {
        auto new_multi_vec = new MultiVecWrapper<VectorType>;

        new_multi_vec->vectors.resize(index.size());

        for (unsigned int i = 0; i < index.size(); ++i)
          {
            AssertThrow(static_cast<unsigned int>(index[i]) <
                          this->vectors.size(),
                        ExcInternalError());

            new_multi_vec->vectors[i].reset(
              this->vectors[index[i]].get(),
              [](
                auto
                  *) { /*Nothing to do, since we are creating only a view.*/ });
          }

        return new_multi_vec;
      }

      /**
       * Creates a new Belos::MultiVec that shares the selected contents
       * of *this. The index of the numvecs vectors copied from *this
       * are indicated by the indices given in index.
       */
      virtual const Belos::MultiVec<value_type> *
      CloneView(const std::vector<int> &index) const
      {
        auto new_multi_vec = new MultiVecWrapper<VectorType>;

        new_multi_vec->vectors.resize(index.size());

        for (unsigned int i = 0; i < index.size(); ++i)
          {
            AssertThrow(static_cast<unsigned int>(index[i]) <
                          this->vectors.size(),
                        ExcInternalError());

            new_multi_vec->vectors[i].reset(
              this->vectors[index[i]].get(),
              [](
                auto
                  *) { /*Nothing to do, since we are creating only a view.*/ });
          }

        return new_multi_vec;
      }

      /**
       * The number of rows in the multivector.
       */
      virtual std::ptrdiff_t
      GetGlobalLength() const
      {
        AssertThrow(this->vectors.size() > 0, ExcInternalError());

        for (unsigned int i = 1; i < this->vectors.size(); ++i)
          AssertDimension(this->vectors[0]->size(), this->vectors[i]->size());

        return this->vectors[0]->size();
      }

      /**
       * The number of vectors (i.e., columns) in the multivector.
       */
      virtual int
      GetNumberVecs() const
      {
        return vectors.size();
      }

      /**
       * Update *this with alpha * A * B + beta * (*this).
       */
      virtual void
      MvTimesMatAddMv(const value_type                                   alpha,
                      const Belos::MultiVec<value_type>                 &A_,
                      const Teuchos::SerialDenseMatrix<int, value_type> &B,
                      const value_type                                   beta)
      {
        const auto &A = try_to_get_underlying_vector(A_);

        const unsigned int n_rows = B.numRows();
        const unsigned int n_cols = B.numCols();

        AssertThrow(n_rows == static_cast<unsigned int>(A.GetNumberVecs()),
                    ExcInternalError());
        AssertThrow(n_cols == static_cast<unsigned int>(this->GetNumberVecs()),
                    ExcInternalError());

        for (unsigned int i = 0; i < n_cols; ++i)
          (*this->vectors[i]) *= beta;

        for (unsigned int i = 0; i < n_cols; ++i)
          for (unsigned int j = 0; j < n_rows; ++j)
            this->vectors[i]->add(alpha * B(j, i), *A.vectors[j]);
      }

      /**
       * Replace *this with alpha * A + beta * B.
       */
      virtual void
      MvAddMv(const value_type                   alpha,
              const Belos::MultiVec<value_type> &A_,
              const value_type                   beta,
              const Belos::MultiVec<value_type> &B_)
      {
        const auto &A = try_to_get_underlying_vector(A_);
        const auto &B = try_to_get_underlying_vector(B_);

        AssertThrow(this->vectors.size() == A.vectors.size(),
                    ExcInternalError());
        AssertThrow(this->vectors.size() == B.vectors.size(),
                    ExcInternalError());

        for (unsigned int i = 0; i < this->vectors.size(); ++i)
          {
            this->vectors[i]->equ(alpha, *A.vectors[i]);
            this->vectors[i]->add(beta, *B.vectors[i]);
          }
      }

      /**
       * Scale each element of the vectors in *this with alpha.
       */
      virtual void
      MvScale(const value_type alpha)
      {
        for (unsigned int i = 0; i < this->vectors.size(); ++i)
          (*this->vectors[i]) *= alpha;
      }

      /**
       * Scale each element of the i-th vector in *this with alpha[i].
       */
      virtual void
      MvScale(const std::vector<value_type> & /*alpha*/)
      {
        AssertThrow(false, ExcNotImplemented());
      }

      /**
       * Compute a dense matrix B through the matrix-matrix multiply
       * alpha * A^T * (*this).
       */
      virtual void
      MvTransMv(const value_type                             alpha,
                const Belos::MultiVec<value_type>           &A_,
                Teuchos::SerialDenseMatrix<int, value_type> &B) const
      {
        const auto &A = try_to_get_underlying_vector(A_);

        const unsigned int n_rows = B.numRows();
        const unsigned int n_cols = B.numCols();

        AssertThrow(n_rows == static_cast<unsigned int>(A.GetNumberVecs()),
                    ExcInternalError());
        AssertThrow(n_cols == static_cast<unsigned int>(this->GetNumberVecs()),
                    ExcInternalError());

        for (unsigned int i = 0; i < n_rows; ++i)
          for (unsigned int j = 0; j < n_cols; ++j)
            B(i, j) = alpha * ((*A.vectors[i]) * (*this->vectors[j]));
      }

      /**
       * Compute the dot product of each column of *this with the
       * corresponding column of A.
       */
      virtual void
      MvDot(const Belos::MultiVec<value_type> &A_,
            std::vector<value_type>           &b) const
      {
        const auto &A = try_to_get_underlying_vector(A_);

        AssertThrow(this->vectors.size() == A.vectors.size(),
                    ExcInternalError());
        AssertThrow(this->vectors.size() == b.size(), ExcInternalError());

        for (unsigned int i = 0; i < this->vectors.size(); ++i)
          b[i] = (*this->vectors[i]) * (*A.vectors[i]);
      }

      /**
       * Compute the norm of each vector in *this.
       */
      virtual void
      MvNorm(
        std::vector<typename Teuchos::ScalarTraits<value_type>::magnitudeType>
                       &normvec,
        Belos::NormType type = Belos::TwoNorm) const
      {
        AssertThrow(type == Belos::TwoNorm, ExcNotImplemented());
        AssertThrow(this->vectors.size() == normvec.size(), ExcInternalError());

        for (unsigned int i = 0; i < this->vectors.size(); ++i)
          normvec[i] = this->vectors[i]->l2_norm();
      }

      /**
       * Copy the vectors in A to a set of vectors in *this.
       */
      virtual void
      SetBlock(const Belos::MultiVec<value_type> & /*A*/,
               const std::vector<int> & /*index*/)
      {
        AssertThrow(false, ExcNotImplemented());
      }

      /**
       * Fill all the vectors in *this with random numbers.
       */
      virtual void
      MvRandom()
      {
        AssertThrow(false, ExcNotImplemented());
      }

      /**
       * Replace each element of the vectors in *this with alpha.
       */
      virtual void
      MvInit(const value_type /*alpha*/)
      {
        AssertThrow(false, ExcNotImplemented());
      }

      /**
       * Print *this multivector to the os output stream.
       */
      virtual void
      MvPrint(std::ostream & /*os*/) const
      {
        AssertThrow(false, ExcNotImplemented());
      }

      /**
       * Get underlying vector.
       */
      VectorType &
      genericVector()
      {
        AssertThrow(GetNumberVecs() == 1, ExcNotImplemented());

        return *vectors[0];
      }

      /**
       * Get underlying vector (const version).
       */
      const VectorType &
      genericVector() const
      {
        AssertThrow(GetNumberVecs() == 1, ExcNotImplemented());

        return *vectors[0];
      }

      /**
       * Cast Belos::MultiVec to MultiVecWrapper.
       */
      static MultiVecWrapper<VectorType> &
      try_to_get_underlying_vector(Belos::MultiVec<value_type> &vec_in)
      {
        auto vec = dynamic_cast<MultiVecWrapper<VectorType> *>(&vec_in);

        AssertThrow(vec, ExcInternalError());

        return *vec;
      }

      /**
       * Cast Belos::MultiVec to MultiVecWrapper (const version).
       */
      const static MultiVecWrapper<VectorType> &
      try_to_get_underlying_vector(const Belos::MultiVec<value_type> &vec_in)
      {
        auto vec = dynamic_cast<const MultiVecWrapper<VectorType> *>(&vec_in);

        AssertThrow(vec, ExcInternalError());

        return *vec;
      }


#      ifdef HAVE_BELOS_TSQR
      virtual void
      factorExplicit(Belos::MultiVec<value_type> &,
                     Teuchos::SerialDenseMatrix<int, value_type> &,
                     const bool = false)
      {
        DEAL_II_NOT_IMPLEMENTED();
      }

      virtual int
      revealRank(
        Teuchos::SerialDenseMatrix<int, value_type> &,
        const typename Teuchos::ScalarTraits<value_type>::magnitudeType &)
      {
        DEAL_II_NOT_IMPLEMENTED();
      }

#      endif

    private:
      /**
       * Underlying vectors.
       */
      std::vector<std::shared_ptr<VectorType>> vectors;

      /**
       * Default constructor. Only used internally to create new MultiVecWrapper
       * instances.
       */
      MultiVecWrapper() = default;
    };

    /**
     * Implementation of the abstract interface
     * Belos::Operator for deal.II vectors and deal.II
     * operators/preconditioners. For details, see
     * https://docs.trilinos.org/latest-release/packages/belos/doc/html/classBelos_1_1Operator.html.
     */
    template <typename OperatorType, typename VectorType>
    DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
    class OperatorWrapper
      : public Belos::Operator<typename VectorType::value_type>
    {
    public:
      /**
       * Underlying value type.
       */
      using value_type = typename VectorType::value_type;

      /**
       * Indicate that specialization exists.
       */
      static bool
      this_type_is_missing_a_specialization()
      {
        return false;
      }

      /**
       * Constructor.
       */
      OperatorWrapper(const OperatorType &op)
        : op(op)
      {}

      /**
       * Destructor.
       */
      virtual ~OperatorWrapper() = default;

      /**
       * Apply the operator to x, putting the result in y.
       */
      virtual void
      Apply(const Belos::MultiVec<value_type> &x,
            Belos::MultiVec<value_type>       &y,
            Belos::ETrans                      trans = Belos::NOTRANS) const
      {
        // TODO: check for Tvmult
        AssertThrow(trans == Belos::NOTRANS, ExcNotImplemented());

        op.vmult(MultiVecWrapper<VectorType>::try_to_get_underlying_vector(y)
                   .genericVector(),
                 MultiVecWrapper<VectorType>::try_to_get_underlying_vector(x)
                   .genericVector());
      }

      /**
       * Whether this operator implements applying the transpose.
       */
      virtual bool
      HasApplyTranspose() const
      {
        // TODO: check for Tvmult
        return false;
      }

    private:
      /**
       * Underlying operator.
       */
      const OperatorType &op;
    };

  } // namespace internal



  template <typename VectorType>
  DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
  SolverBelos<VectorType>::SolverBelos(
    SolverControl                              &solver_control,
    const AdditionalData                       &additional_data,
    const Teuchos::RCP<Teuchos::ParameterList> &belos_parameters)
    : solver_control(solver_control)
    , additional_data(additional_data)
    , belos_parameters(belos_parameters)
  {}



  template <typename VectorType>
  DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
  template <typename OperatorType, typename PreconditionerType>
  void SolverBelos<VectorType>::solve(const OperatorType       &A_dealii,
                                      VectorType               &x_dealii,
                                      const VectorType         &b_dealii,
                                      const PreconditionerType &P_dealii)
  {
    using value_type = typename VectorType::value_type;

    using MV = Belos::MultiVec<value_type>;
    using OP = Belos::Operator<value_type>;

    Teuchos::RCP<OP> A = Teuchos::rcp(
      new internal::OperatorWrapper<OperatorType, VectorType>(A_dealii));
    Teuchos::RCP<OP> P = Teuchos::rcp(
      new internal::OperatorWrapper<PreconditionerType, VectorType>(P_dealii));
    Teuchos::RCP<MV> X =
      Teuchos::rcp(new internal::MultiVecWrapper<VectorType>(x_dealii));
    Teuchos::RCP<MV> B =
      Teuchos::rcp(new internal::MultiVecWrapper<VectorType>(b_dealii));

    Teuchos::RCP<Belos::LinearProblem<value_type, MV, OP>> problem =
      Teuchos::rcp(new Belos::LinearProblem<value_type, MV, OP>(A, X, B));

    if (additional_data.right_preconditioning == false)
      problem->setLeftPrec(P);
    else
      problem->setRightPrec(P);

    const bool problem_set = problem->setProblem();
    AssertThrow(problem_set, ExcInternalError());

    // compute initial residal
    VectorType r;
    r.reinit(x_dealii, true);
    A_dealii.vmult(r, x_dealii);
    r.sadd(-1., 1., b_dealii);
    const auto norm_0 = r.l2_norm();

    if (solver_control.check(0, norm_0) != SolverControl::iterate)
      return;

    double relative_tolerance_to_be_achieved =
      solver_control.tolerance() / norm_0;
    const unsigned int max_steps = solver_control.max_steps();

    if (const auto *reduction_control =
          dynamic_cast<ReductionControl *>(&solver_control))
      relative_tolerance_to_be_achieved =
        std::max(relative_tolerance_to_be_achieved,
                 reduction_control->reduction());

    Teuchos::RCP<Teuchos::ParameterList> belos_parameters_copy(
      Teuchos::rcp(new Teuchos::ParameterList(*belos_parameters)));

    belos_parameters_copy->set("Convergence Tolerance",
                               relative_tolerance_to_be_achieved);
    belos_parameters_copy->set("Maximum Iterations",
                               static_cast<int>(max_steps));

    Teuchos::RCP<Belos::SolverManager<value_type, MV, OP>> solver;

    if (additional_data.solver_name == SolverName::cg)
      solver = Teuchos::rcp(
        new Belos::BlockCGSolMgr<value_type, MV, OP>(problem,
                                                     belos_parameters_copy));
    else if (additional_data.solver_name == SolverName::gmres)
      solver = Teuchos::rcp(
        new Belos::BlockGmresSolMgr<value_type, MV, OP>(problem,
                                                        belos_parameters_copy));
    else
      AssertThrow(false, ExcNotImplemented());

    const auto flag = solver->solve();

    solver_control.check(solver->getNumIters(), solver->achievedTol() * norm_0);

    AssertThrow(flag == Belos::ReturnType::Converged ||
                  ((dynamic_cast<IterationNumberControl *>(&solver_control) !=
                    nullptr) &&
                   (solver_control.last_step() == max_steps)),
                SolverControl::NoConvergence(solver_control.last_step(),
                                             solver_control.last_value()));
  }

} // namespace TrilinosWrappers
#    endif

#  endif

DEAL_II_NAMESPACE_CLOSE

#else

// Make sure the scripts that create the C++20 module input files have
// something to latch on if the preprocessor #ifdef above would
// otherwise lead to an empty content of the file.
DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_TRILINOS

#endif
