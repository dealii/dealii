// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2014 by the deal.II authors
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

#ifndef __deal2__trilinos_solver_h
#define __deal2__trilinos_solver_h


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_TRILINOS

#  include <deal.II/base/std_cxx11/shared_ptr.h>
#  include <deal.II/lac/exceptions.h>
#  include <deal.II/lac/solver_control.h>
#  include <deal.II/lac/vector.h>
#  include <deal.II/lac/parallel_vector.h>

#  include <Epetra_LinearProblem.h>
#  include <AztecOO.h>
#  include <Epetra_Operator.h>
#  include <Amesos.h>

DEAL_II_NAMESPACE_OPEN

namespace TrilinosWrappers
{
  // forward declarations
  class SparseMatrix;
  class VectorBase;
  class PreconditionBase;


  /**
   * Base class for solver classes using the Trilinos solvers. Since
   * solvers in Trilinos are selected based on flags passed to a generic
   * solver object, basically all the actual solver calls happen in this
   * class, and derived classes simply set the right flags to select one
   * solver or another, or to set certain parameters for individual
   * solvers. For a general discussion on the Trilinos solver package
   * AztecOO, we refer to the <a href =
   * "http://trilinos.sandia.gov/packages/aztecoo/AztecOOUserGuide.pdf">AztecOO
   * user guide</a>.
   *
   * This solver class can also be used as a standalone class, where the
   * respective Krylov method is set via the flag
   * <tt>solver_name</tt>. This can be done at runtime (e.g., when
   * parsing the solver from a ParameterList) and is similar to the
   * deal.II class SolverSelector.
   *
   * @ingroup TrilinosWrappers
   * @author Martin Kronbichler, 2008, 2009
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
    enum SolverName {cg, cgs, gmres, bicgstab, tfqmr} solver_name;

    /**
     * Standardized data struct to
     * pipe additional data to the
     * solver.
     */

    struct AdditionalData
    {
      /**
       * Sets the additional data field to the desired output format and puts
       * the restart parameter in case the derived class is GMRES.
       *
       * TODO: Find a better way for setting the GMRES restart parameter since
       * it is quite inelegant to set a specific option of one solver in the
       * base class for all solvers.
       */
      AdditionalData (const bool         output_solver_details   = false,
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
    SolverBase (SolverControl  &cn);

    /**
     * Second constructor. This constructor takes an enum object that
     * specifies the solver name and sets the appropriate Krylov method.
     */
    SolverBase (const enum SolverName  solver_name,
                SolverControl         &cn);

    /**
     * Destructor.
     */
    virtual ~SolverBase ();

    /**
     * Solve the linear system <tt>Ax=b</tt>. Depending on the information
     * provided by derived classes and the object passed as a preconditioner,
     * one of the linear solvers and preconditioners of Trilinos is chosen.
     */
    void
    solve (const SparseMatrix     &A,
           VectorBase             &x,
           const VectorBase       &b,
           const PreconditionBase &preconditioner);

    /**
     * Solve the linear system <tt>Ax=b</tt> where <tt>A</tt> is an
     * operator. This function can be used for matrix free
     * computation. Depending on the information provided by derived classes
     * and the object passed as a preconditioner, one of the linear solvers
     * and preconditioners of Trilinos is chosen.
     */
    void
    solve (Epetra_Operator        &A,
           VectorBase             &x,
           const VectorBase       &b,
           const PreconditionBase &preconditioner);

    /**
     * Solve the linear system <tt>Ax=b</tt>. Depending on the information
     * provided by derived classes and the object passed as a preconditioner,
     * one of the linear solvers and preconditioners of Trilinos is
     * chosen. This class works with matrices according to the
     * TrilinosWrappers format, but can take deal.II vectors as
     * argument. Since deal.II are serial vectors (not distributed), this
     * function does only what you expect in case the matrix is locally
     * owned. Otherwise, an exception will be thrown.
     */
    void
    solve (const SparseMatrix           &A,
           dealii::Vector<double>       &x,
           const dealii::Vector<double> &b,
           const PreconditionBase       &preconditioner);

    /**
     * Solve the linear system <tt>Ax=b</tt> where <tt>A</tt> is an
     * operator. This function can be used for matrix free
     * computations. Depending on the information provided by derived classes
     * and the object passed as a preconditioner, one of the linear solvers
     * and preconditioners of Trilinos is chosen. This class works with
     * matrices according to the TrilinosWrappers format, but can take deal.II
     * vectors as argument. Since deal.II are serial vectors (not
     * distributed), this function does only what you expect in case the
     * matrix is locally owned. Otherwise, an exception will be thrown.
     */
    void
    solve (Epetra_Operator              &A,
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
    solve (const SparseMatrix                                  &A,
           dealii::parallel::distributed::Vector<double>       &x,
           const dealii::parallel::distributed::Vector<double> &b,
           const PreconditionBase                              &preconditioner);

    /**
     * Solve the linear system <tt>Ax=b</tt> where <tt>A</tt> is an
     * operator. This function can be used for matrix free
     * computation. Depending on the information provided by derived classes
     * and the object passed as a preconditioner, one of the linear solvers
     * and preconditioners of Trilinos is chosen.
     */
    void
    solve (Epetra_Operator                                     &A,
           dealii::parallel::distributed::Vector<double>       &x,
           const dealii::parallel::distributed::Vector<double> &b,
           const PreconditionBase                              &preconditioner);


    /**
     * Access to object that controls convergence.
     */
    SolverControl &control() const;

    /**
     * Exception
     */
    DeclException1 (ExcTrilinosError,
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
    void do_solve(const PreconditionBase &preconditioner);

    /**
     * A structure that collects the Trilinos sparse matrix, the right hand
     * side vector and the solution vector, which is passed down to the
     * Trilinos solver.
     */
    std_cxx11::shared_ptr<Epetra_LinearProblem> linear_problem;

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



  /**
   * An implementation of the solver interface using the Trilinos CG
   * solver.
   *
   * @ingroup TrilinosWrappers
   * @author Martin Kronbichler, 2008
   */
  class SolverCG : public SolverBase
  {
  public:
    /**
     * Standardized data struct to pipe additional data to the solver.
     */

    struct AdditionalData
    {
      /**
       * Sets the additional data field to the desired output format.
       */
      AdditionalData (const bool output_solver_details = false);

      /**
       * Enables/disables the output of solver details (residual in each
       * iterations etc.).
       */
      bool output_solver_details;
    };

    /**
     * Constructor. In contrast to deal.II's own solvers, there is no need to
     * give a vector memory object.
     *
     * The last argument takes a structure with additional, solver dependent
     * flags for tuning.
     */
    SolverCG (SolverControl        &cn,
              const AdditionalData &data = AdditionalData());

  protected:
    /**
     * Store a copy of the flags for this particular solver.
     */
    const AdditionalData additional_data;
  };



  /**
   * An implementation of the solver interface using the Trilinos CGS
   * solver.
   *
   * @ingroup TrilinosWrappers
   * @author Martin Kronbichler, 2008
   */
  class SolverCGS : public SolverBase
  {
  public:
    /**
     * Standardized data struct to pipe additional data to the solver.
     */
    struct AdditionalData
    {
      /**
       * Sets the additional data field to the desired output format.
       */
      AdditionalData (const bool output_solver_details = false);

      /**
       * Enables/disables the output of solver details (residual in each
       * iterations etc.).
       */
      bool output_solver_details;
    };

    /**
     * Constructor. In contrast to deal.II's own solvers, there is no need to
     * give a vector memory object.
     *
     * The last argument takes a structure with additional, solver dependent
     * flags for tuning.
     */
    SolverCGS (SolverControl        &cn,
               const AdditionalData &data = AdditionalData());

  protected:
    /**
     * Store a copy of the flags for this particular solver.
     */
    const AdditionalData additional_data;
  };



  /**
   * An implementation of the solver interface using the Trilinos GMRES
   * solver.
   *
   * @author Martin Kronbichler, 2008
   */
  class SolverGMRES : public SolverBase
  {
  public:
    /**
     * Standardized data struct to pipe additional data to the solver.
     */
    struct AdditionalData
    {
      /**
       * Constructor. By default, set the number of temporary vectors to 30,
       * i.e. do a restart every 30 iterations.
       */
      AdditionalData (const bool         output_solver_details = false,
                      const unsigned int restart_parameter = 30);

      /**
       * Enables/disables the output of solver details (residual in each
       * iterations etc.).
       */
      bool output_solver_details;

      /**
       * Maximum number of tmp vectors.
       */
      unsigned int restart_parameter;
    };

    /**
     * Constructor. In contrast to deal.II's own solvers, there is no need to
     * give a vector memory object.
     *
     * The last argument takes a structure with additional, solver dependent
     * flags for tuning.
     */
    SolverGMRES (SolverControl        &cn,
                 const AdditionalData &data = AdditionalData());

  protected:
    /**
     * Store a copy of the flags for this particular solver.
     */
    const AdditionalData additional_data;
  };



  /**
   * An implementation of the solver interface using the Trilinos BiCGStab
   * solver.
   *
   * @ingroup TrilinosWrappers
   * @author Martin Kronbichler, 2008
   */
  class SolverBicgstab : public SolverBase
  {
  public:
    /**
     * Standardized data struct to pipe additional data to the solver.
     */
    struct AdditionalData
    {
      /**
       * Sets the additional data field to the desired output format.
       */
      AdditionalData (const bool output_solver_details = false);

      /**
       * Enables/disables the output of solver details (residual in each
       * iterations etc.).
       */
      bool output_solver_details;
    };

    /**
     * Constructor. In contrast to deal.II's own solvers, there is no need to
     * give a vector memory object.
     *
     * The last argument takes a structure with additional, solver dependent
     * flags for tuning.
     */
    SolverBicgstab (SolverControl        &cn,
                    const AdditionalData &data = AdditionalData());

  protected:
    /**
     * Store a copy of the flags for this particular solver.
     */
    const AdditionalData additional_data;
  };



  /**
   * An implementation of the solver interface using the Trilinos TFQMR
   * solver.
   *
   * @ingroup TrilinosWrappers
   * @author Martin Kronbichler, 2008
   */
  class SolverTFQMR : public SolverBase
  {
  public:
    /**
     * Standardized data struct to pipe additional data to the solver.
     */
    struct AdditionalData
    {
      /**
       * Sets the additional data field to the desired output format.
       */
      AdditionalData (const bool output_solver_details = false);

      /**
       * Enables/disables the output of solver details (residual in each
       * iterations etc.).
       */
      bool output_solver_details;
    };

    /**
     * Constructor. In contrast to deal.II's own solvers, there is no need to
     * give a vector memory object.
     *
     * The last argument takes a structure with additional, solver dependent
     * flags for tuning.
     */
    SolverTFQMR (SolverControl        &cn,
                 const AdditionalData &data = AdditionalData());

  protected:
    /**
     * Store a copy of the flags for this particular solver.
     */
    const AdditionalData additional_data;
  };



  /**
   * An implementation of Trilinos direct solvers (using the Amesos
   * package). The data field AdditionalData::solver_type can be used to
   * specify the type of solver. It allows the use of built-in solvers
   * Amesos_Klu as well as third-party solvers Amesos_Superludist or
   * Amesos_Mumps.
   *
   * For instructions on how to install Trilinos for use with
   * direct solvers other than KLU, see the link to the Trilinos
   * installation instructions linked to from the deal.II ReadMe file.
   *
   * @ingroup TrilinosWrappers
   * @author Martin Kronbichler, 2009, Uwe K&ouml;cher, 2014
   */
  class SolverDirect
  {
  public:

    /**
     * Standardized data struct to pipe additional data to the solver.
     */

    struct AdditionalData
    {
      /**
       * Sets the additional data field to the desired output format.
       */
      AdditionalData (const bool output_solver_details = false,
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
       *   <li>  "Amesos_Lapack" </li>
       *   <li>  "Amesos_Scalapack" </li>
       *   <li>  "Amesos_Klu" </li>
       *   <li>  "Amesos_Umfpack" </li>
       *   <li>  "Amesos_Pardiso" </li>
       *   <li>  "Amesos_Taucs" </li>
       *   <li>  "Amesos_Superlu" </li>
       *   <li>  "Amesos_Superludist" </li>
       *   <li>  "Amesos_Dscpack" </li>
       *   <li>  "Amesos_Mumps" </li>
       * </ul>
       * Note that the availability of these solvers in deal.II depends on
       * which solvers were set when configuring Trilinos.
       */
      std::string solver_type;
    };

    /**
     * Constructor. Takes the solver control object and creates the solver.
     */
    SolverDirect (SolverControl  &cn,
                  const AdditionalData &data = AdditionalData());

    /**
     * Destructor.
     */
    virtual ~SolverDirect ();

    /**
     * Solve the linear system <tt>Ax=b</tt>. Creates a factorization of the
     * matrix with the package chosen from the additional data structure and
     * performs the solve. Note that there is no need for a preconditioner
     * here.
     */
    void
    solve (const SparseMatrix     &A,
           VectorBase             &x,
           const VectorBase       &b);

    /**
     * Solve the linear system <tt>Ax=b</tt>. This class works with Trilinos
     * matrices, but takes deal.II serial vectors as argument. Since these
     * vectors are not distributed, this function does only what you expect in
     * case the matrix is serial (i.e., locally owned). Otherwise, an
     * exception will be thrown.
     */
    void
    solve (const SparseMatrix           &A,
           dealii::Vector<double>       &x,
           const dealii::Vector<double> &b);

    /**
     * Solve the linear system <tt>Ax=b</tt> for deal.II's own parallel
     * vectors. Creates a factorization of the matrix with the package chosen
     * from the additional data structure and performs the solve. Note that
     * there is no need for a preconditioner here.
     */
    void
    solve (const SparseMatrix                                  &A,
           dealii::parallel::distributed::Vector<double>       &x,
           const dealii::parallel::distributed::Vector<double> &b);

    /**
     * Access to object that controls
     * convergence.
     */
    SolverControl &control() const;

    /**
     * Exception
     */
    DeclException1 (ExcTrilinosError,
                    int,
                    << "An error with error number " << arg1
                    << " occurred while calling a Trilinos function");

  private:
    /**
     * Actually performs the operations for solving the linear system,
     * including the factorization and forward and backward substitution.
     */
    void do_solve();

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
    std_cxx11::shared_ptr<Epetra_LinearProblem> linear_problem;

    /**
     * A structure that contains the Trilinos solver and preconditioner
     * objects.
     */
    std_cxx11::shared_ptr<Amesos_BaseSolver> solver;

    /**
     * Store a copy of the flags for this particular solver.
     */
    const AdditionalData additional_data;
  };

}

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_TRILINOS

/*----------------------------   trilinos_solver.h     ---------------------------*/

#endif
/*----------------------------   trilinos_solver.h     ---------------------------*/
