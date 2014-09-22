// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2013 by the deal.II authors
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

#ifndef __deal2__petsc_solver_h
#define __deal2__petsc_solver_h


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_PETSC

#  include <deal.II/lac/exceptions.h>
#  include <deal.II/lac/solver_control.h>
#  include <deal.II/base/std_cxx11/shared_ptr.h>

#  include <petscksp.h>

DEAL_II_NAMESPACE_OPEN


namespace PETScWrappers
{
  // forward declarations
  class MatrixBase;
  class VectorBase;
  class PreconditionerBase;


  /**
   * Base class for solver classes using the PETSc solvers. Since solvers in
   * PETSc are selected based on flags passed to a generic solver object,
   * basically all the actual solver calls happen in this class, and derived
   * classes simply set the right flags to select one solver or another, or to
   * set certain parameters for individual solvers.
   *
   * Optionally, the user can create a solver derived from the
   * SolverBase class and can set the default arguments necessary to
   * solve the linear system of equations with SolverControl. These
   * default options can be overridden by specifying command line
   * arguments of the form @p -ksp_*. For example,
   * @p -ksp_monitor_true_residual prints out true residual norm
   * (unpreconditioned) at each iteration and @p -ksp_view provides
   * information about the linear solver and the preconditioner used in
   * the current context. The type of the solver can also be changed
   * during runtime by specifying @p -ksp_type {richardson, cg, gmres,
   * fgmres, ..} to dynamically test the optimal solver along with a
   * suitable preconditioner set using @p -pc_type {jacobi, bjacobi,
   * ilu, lu, ..}. There are several other command line options
   * available to modify the behavior of the PETSc linear solver and can
   * be obtained from the <a
   * href="http://www.mcs.anl.gov/petsc">documentation and manual
   * pages</a>.
   *
   * @note Repeated calls to solve() on a solver object with a Preconditioner
   * must be used with care. The preconditioner is initialized in the first call
   * to solve() and subsequent calls reuse the solver and preconditioner
   * object. This is done for performance reasons. The solver and preconditioner
   * can be reset by calling reset().
   *
   * One of the gotchas of PETSc is that -- in particular in MPI mode -- it
   * often does not produce very helpful error messages. In order to save
   * other users some time in searching a hard to track down error, here is
   * one situation and the error message one gets there:
   * when you don't specify an MPI communicator to your solver's constructor. In
   * this case, you will get an error of the following form from each of your
   * parallel processes:
   * @verbatim
   *   [1]PETSC ERROR: PCSetVector() line 1173 in src/ksp/pc/interface/precon.c
   *   [1]PETSC ERROR:   Arguments must have same communicators!
   *   [1]PETSC ERROR:   Different communicators in the two objects: Argument # 1 and 2!
   *   [1]PETSC ERROR: KSPSetUp() line 195 in src/ksp/ksp/interface/itfunc.c
   * @endverbatim
   *
   * This error, on which one can spend a very long time figuring out
   * what exactly goes wrong, results from not specifying an MPI
   * communicator. Note that the communicator @em must match that of the
   * matrix and all vectors in the linear system which we want to
   * solve. Aggravating the situation is the fact that the default
   * argument to the solver classes, @p PETSC_COMM_SELF, is the
   * appropriate argument for the sequential case (which is why it is
   * the default argument), so this error only shows up in parallel
   * mode.
   *
   * @ingroup PETScWrappers
   * @author Wolfgang Bangerth, 2004
   */
  class SolverBase
  {
  public:
    /**
     * Constructor. Takes the solver
     * control object and the MPI
     * communicator over which parallel
     * computations are to happen.
     *
     * Note that the communicator used here
     * must match the communicator used in
     * the system matrix, solution, and
     * right hand side object of the solve
     * to be done with this
     * solver. Otherwise, PETSc will
     * generate hard to track down errors,
     * see the documentation of the
     * SolverBase class.
     */
    SolverBase (SolverControl  &cn,
                const MPI_Comm &mpi_communicator);

    /**
     * Destructor.
     */
    virtual ~SolverBase ();

    /**
     * Solve the linear system
     * <tt>Ax=b</tt>. Depending on the
     * information provided by derived
     * classes and the object passed as a
     * preconditioner, one of the linear
     * solvers and preconditioners of PETSc
     * is chosen.  Repeated calls to
     * solve() do not reconstruct the
     * preconditioner for performance
     * reasons. See class Documentation.
     */
    void
    solve (const MatrixBase         &A,
           VectorBase               &x,
           const VectorBase         &b,
           const PreconditionerBase &preconditioner);


    /**
     * Resets the contained preconditioner
     * and solver object. See class
     * description for more details.
     */
    virtual void reset();


    /**
      * Sets a prefix name for the solver
      * object. Useful when customizing the
      * PETSc KSP object with command-line
      * options.
      */
    void set_prefix(const std::string &prefix);


    /**
     * Access to object that controls
     * convergence.
     */
    SolverControl &control() const;

    /**
     * Exception
     */
    DeclException1 (ExcPETScError,
                    int,
                    << "An error with error number " << arg1
                    << " occurred while calling a PETSc function");

  protected:

    /**
     * Reference to the object that
     * controls convergence of the
     * iterative solver. In fact, for these
     * PETSc wrappers, PETSc does so
     * itself, but we copy the data from
     * this object before starting the
     * solution process, and copy the data
     * back into it afterwards.
     */
    SolverControl &solver_control;

    /**
     * Copy of the MPI communicator object
     * to be used for the solver.
     */
    const MPI_Comm mpi_communicator;

    /**
     * Function that takes a Krylov
     * Subspace Solver context object, and
     * sets the type of solver that is
     * requested by the derived class.
     */
    virtual void set_solver_type (KSP &ksp) const = 0;

    /**
     * Solver prefix name to qualify options
     * specific to the PETSc KSP object in the
     * current context.
     * Note: A hyphen (-) must NOT be given
     * at the beginning of the prefix name.
     * The first character of all runtime
     * options is AUTOMATICALLY the hyphen.
     */
    std::string prefix_name;

  private:
    /**
     * A function that is used in PETSc as
     * a callback to check on
     * convergence. It takes the
     * information provided from PETSc and
     * checks it against deal.II's own
     * SolverControl objects to see if
     * convergence has been reached.
     */
    static
    PetscErrorCode convergence_test (KSP                 ksp,
                                     const PetscInt      iteration,
                                     const PetscReal     residual_norm,
                                     KSPConvergedReason *reason,
                                     void               *solver_control);

    /**
     * A structure that contains the PETSc
     * solver and preconditioner
     * objects. This object is preserved
     * between subsequent calls to the
     * solver if the same preconditioner is
     * used as in the previous solver
     * step. This may save some computation
     * time, if setting up a preconditioner
     * is expensive, such as in the case of
     * an ILU for example.
     *
     * The actual declaration of this class
     * is complicated by the fact that
     * PETSc changed its solver interface
     * completely and incompatibly between
     * versions 2.1.6 and 2.2.0 :-(
     *
     * Objects of this type are explicitly
     * created, but are destroyed when the
     * surrounding solver object goes out
     * of scope, or when we assign a new
     * value to the pointer to this
     * object. The respective *Destroy
     * functions are therefore written into
     * the destructor of this object, even
     * though the object does not have a
     * constructor.
     */
    struct SolverData
    {
      /**
       * Destructor
       */
      ~SolverData ();

      /**
       * Objects for Krylov subspace
       * solvers and preconditioners.
       */
      KSP  ksp;
      PC   pc;
    };

    /**
     * Pointer to an object that stores the
     * solver context. This is recreated in
     * the main solver routine if
     * necessary.
     */
    std_cxx11::shared_ptr<SolverData> solver_data;
  };



  /**
   * An implementation of the solver interface using the PETSc Richardson
   * solver.
   *
   * @ingroup PETScWrappers
   * @author Wolfgang Bangerth, 2004
   */
  class SolverRichardson : public SolverBase
  {
  public:
    /**
     * Standardized data struct to
     * pipe additional data to the
     * solver.
     */
    struct AdditionalData
    {
      /**
       * Constructor. By default,
       * set the damping parameter
       * to one.
       */
      AdditionalData (const double omega = 1);

      /**
       * Relaxation parameter.
       */
      double omega;
    };

    /**
     * Constructor. In contrast to
     * deal.II's own solvers, there is no
     * need to give a vector memory
     * object. However, PETSc solvers want
     * to have an MPI communicator context
     * over which computations are
     * parallelized. By default,
     * @p PETSC_COMM_SELF is used here,
     * but you can change this. Note that
     * for single processor (non-MPI)
     * versions, this parameter does not
     * have any effect.
     *
     * The last argument takes a structure
     * with additional, solver dependent
     * flags for tuning.
     *
     * Note that the communicator used here
     * must match the communicator used in
     * the system matrix, solution, and
     * right hand side object of the solve
     * to be done with this
     * solver. Otherwise, PETSc will
     * generate hard to track down errors,
     * see the documentation of the
     * SolverBase class.
     */
    SolverRichardson (SolverControl        &cn,
                      const MPI_Comm       &mpi_communicator = PETSC_COMM_SELF,
                      const AdditionalData &data = AdditionalData());

  protected:
    /**
     * Store a copy of the flags for this
     * particular solver.
     */
    const AdditionalData additional_data;

    /**
     * Function that takes a Krylov
     * Subspace Solver context object, and
     * sets the type of solver that is
     *appropriate for this class.
     */
    virtual void set_solver_type (KSP &ksp) const;
  };



  /**
   * An implementation of the solver interface using the PETSc Chebychev
   * solver.
   *
   * @ingroup PETScWrappers
   * @author Wolfgang Bangerth, 2004
   */
  class SolverChebychev : public SolverBase
  {
  public:
    /**
     * Standardized data struct to
     * pipe additional data to the
     * solver.
     */
    struct AdditionalData
    {};

    /**
     * Constructor. In contrast to
     * deal.II's own solvers, there is no
     * need to give a vector memory
     * object. However, PETSc solvers want
     * to have an MPI communicator context
     * over which computations are
     * parallelized. By default,
     * @p PETSC_COMM_SELF is used here,
     * but you can change this. Note that
     * for single processor (non-MPI)
     * versions, this parameter does not
     * have any effect.
     *
     * The last argument takes a structure
     * with additional, solver dependent
     * flags for tuning.
     *
     * Note that the communicator used here
     * must match the communicator used in
     * the system matrix, solution, and
     * right hand side object of the solve
     * to be done with this
     * solver. Otherwise, PETSc will
     * generate hard to track down errors,
     * see the documentation of the
     * SolverBase class.
     */
    SolverChebychev (SolverControl        &cn,
                     const MPI_Comm       &mpi_communicator = PETSC_COMM_SELF,
                     const AdditionalData &data = AdditionalData());

  protected:
    /**
     * Store a copy of the flags for this
     * particular solver.
     */
    const AdditionalData additional_data;

    /**
     * Function that takes a Krylov
     * Subspace Solver context object, and
     * sets the type of solver that is
     *appropriate for this class.
     */
    virtual void set_solver_type (KSP &ksp) const;
  };



  /**
   * An implementation of the solver interface using the PETSc CG
   * solver.
   *
   * @ingroup PETScWrappers
   * @author Wolfgang Bangerth, 2004
   */
  class SolverCG : public SolverBase
  {
  public:
    /**
     * Standardized data struct to
     * pipe additional data to the
     * solver.
     */
    struct AdditionalData
    {};

    /**
     * Constructor. In contrast to
     * deal.II's own solvers, there is no
     * need to give a vector memory
     * object. However, PETSc solvers want
     * to have an MPI communicator context
     * over which computations are
     * parallelized. By default,
     * @p PETSC_COMM_SELF is used here,
     * but you can change this. Note that
     * for single processor (non-MPI)
     * versions, this parameter does not
     * have any effect.
     *
     * The last argument takes a structure
     * with additional, solver dependent
     * flags for tuning.
     *
     * Note that the communicator used here
     * must match the communicator used in
     * the system matrix, solution, and
     * right hand side object of the solve
     * to be done with this
     * solver. Otherwise, PETSc will
     * generate hard to track down errors,
     * see the documentation of the
     * SolverBase class.
     */
    SolverCG (SolverControl        &cn,
              const MPI_Comm       &mpi_communicator = PETSC_COMM_SELF,
              const AdditionalData &data = AdditionalData());

  protected:
    /**
     * Store a copy of the flags for this
     * particular solver.
     */
    const AdditionalData additional_data;

    /**
     * Function that takes a Krylov
     * Subspace Solver context object, and
     * sets the type of solver that is
     *appropriate for this class.
     */
    virtual void set_solver_type (KSP &ksp) const;
  };



  /**
   * An implementation of the solver interface using the PETSc BiCG
   * solver.
   *
   * @ingroup PETScWrappers
   * @author Wolfgang Bangerth, 2004
   */
  class SolverBiCG : public SolverBase
  {
  public:
    /**
     * Standardized data struct to
     * pipe additional data to the
     * solver.
     */
    struct AdditionalData
    {};

    /**
     * Constructor. In contrast to
     * deal.II's own solvers, there is no
     * need to give a vector memory
     * object. However, PETSc solvers want
     * to have an MPI communicator context
     * over which computations are
     * parallelized. By default,
     * @p PETSC_COMM_SELF is used here,
     * but you can change this. Note that
     * for single processor (non-MPI)
     * versions, this parameter does not
     * have any effect.
     *
     * The last argument takes a structure
     * with additional, solver dependent
     * flags for tuning.
     *
     * Note that the communicator used here
     * must match the communicator used in
     * the system matrix, solution, and
     * right hand side object of the solve
     * to be done with this
     * solver. Otherwise, PETSc will
     * generate hard to track down errors,
     * see the documentation of the
     * SolverBase class.
     */
    SolverBiCG (SolverControl        &cn,
                const MPI_Comm       &mpi_communicator = PETSC_COMM_SELF,
                const AdditionalData &data = AdditionalData());

  protected:
    /**
     * Store a copy of the flags for this
     * particular solver.
     */
    const AdditionalData additional_data;

    /**
     * Function that takes a Krylov
     * Subspace Solver context object, and
     * sets the type of solver that is
     *appropriate for this class.
     */
    virtual void set_solver_type (KSP &ksp) const;
  };



  /**
   * An implementation of the solver interface using the PETSc GMRES
   * solver.
   *
   * @author Wolfgang Bangerth, 2004
   */
  class SolverGMRES : public SolverBase
  {
  public:
    /**
     * Standardized data struct to
     * pipe additional data to the
     * solver.
     */
    struct AdditionalData
    {
      /**
       * Constructor. By default, set the
       * number of temporary vectors to
       * 30, i.e. do a restart every 30
       * iterations.
       */
      AdditionalData (const unsigned int restart_parameter = 30,
                      const bool right_preconditioning = false);

      /**
       * Maximum number of
       * tmp vectors.
       */
      unsigned int restart_parameter;

      /**
       * Flag for right
       * preconditioning.
       */
      bool right_preconditioning;
    };

    /**
     * Constructor. In contrast to
     * deal.II's own solvers, there is no
     * need to give a vector memory
     * object. However, PETSc solvers want
     * to have an MPI communicator context
     * over which computations are
     * parallelized. By default,
     * @p PETSC_COMM_SELF is used here,
     * but you can change this. Note that
     * for single processor (non-MPI)
     * versions, this parameter does not
     * have any effect.
     *
     * The last argument takes a structure
     * with additional, solver dependent
     * flags for tuning.
     *
     * Note that the communicator used here
     * must match the communicator used in
     * the system matrix, solution, and
     * right hand side object of the solve
     * to be done with this
     * solver. Otherwise, PETSc will
     * generate hard to track down errors,
     * see the documentation of the
     * SolverBase class.
     */
    SolverGMRES (SolverControl        &cn,
                 const MPI_Comm       &mpi_communicator = PETSC_COMM_SELF,
                 const AdditionalData &data = AdditionalData());

  protected:
    /**
     * Store a copy of the flags for this
     * particular solver.
     */
    const AdditionalData additional_data;

    /**
     * Function that takes a Krylov
     * Subspace Solver context object, and
     * sets the type of solver that is
     *appropriate for this class.
     */
    virtual void set_solver_type (KSP &ksp) const;
  };



  /**
   * An implementation of the solver interface using the PETSc BiCGStab
   * solver.
   *
   * @ingroup PETScWrappers
   * @author Wolfgang Bangerth, 2004
   */
  class SolverBicgstab : public SolverBase
  {
  public:
    /**
     * Standardized data struct to
     * pipe additional data to the
     * solver.
     */
    struct AdditionalData
    {};

    /**
     * Constructor. In contrast to
     * deal.II's own solvers, there is no
     * need to give a vector memory
     * object. However, PETSc solvers want
     * to have an MPI communicator context
     * over which computations are
     * parallelized. By default,
     * @p PETSC_COMM_SELF is used here,
     * but you can change this. Note that
     * for single processor (non-MPI)
     * versions, this parameter does not
     * have any effect.
     *
     * The last argument takes a structure
     * with additional, solver dependent
     * flags for tuning.
     *
     * Note that the communicator used here
     * must match the communicator used in
     * the system matrix, solution, and
     * right hand side object of the solve
     * to be done with this
     * solver. Otherwise, PETSc will
     * generate hard to track down errors,
     * see the documentation of the
     * SolverBase class.
     */
    SolverBicgstab (SolverControl        &cn,
                    const MPI_Comm       &mpi_communicator = PETSC_COMM_SELF,
                    const AdditionalData &data = AdditionalData());

  protected:
    /**
     * Store a copy of the flags for this
     * particular solver.
     */
    const AdditionalData additional_data;

    /**
     * Function that takes a Krylov
     * Subspace Solver context object, and
     * sets the type of solver that is
     *appropriate for this class.
     */
    virtual void set_solver_type (KSP &ksp) const;
  };

  /**
   * An implementation of the solver interface using the PETSc CG Squared
   * solver.
   *
   * @ingroup PETScWrappers
   * @author Wolfgang Bangerth, 2004
   */
  class SolverCGS : public SolverBase
  {
  public:
    /**
     * Standardized data struct to
     * pipe additional data to the
     * solver.
     */
    struct AdditionalData
    {};

    /**
     * Constructor. In contrast to
     * deal.II's own solvers, there is no
     * need to give a vector memory
     * object. However, PETSc solvers want
     * to have an MPI communicator context
     * over which computations are
     * parallelized. By default,
     * @p PETSC_COMM_SELF is used here,
     * but you can change this. Note that
     * for single processor (non-MPI)
     * versions, this parameter does not
     * have any effect.
     *
     * The last argument takes a structure
     * with additional, solver dependent
     * flags for tuning.
     *
     * Note that the communicator used here
     * must match the communicator used in
     * the system matrix, solution, and
     * right hand side object of the solve
     * to be done with this
     * solver. Otherwise, PETSc will
     * generate hard to track down errors,
     * see the documentation of the
     * SolverBase class.
     */
    SolverCGS (SolverControl        &cn,
               const MPI_Comm       &mpi_communicator = PETSC_COMM_SELF,
               const AdditionalData &data = AdditionalData());

  protected:
    /**
     * Store a copy of the flags for this
     * particular solver.
     */
    const AdditionalData additional_data;

    /**
     * Function that takes a Krylov
     * Subspace Solver context object, and
     * sets the type of solver that is
     *appropriate for this class.
     */
    virtual void set_solver_type (KSP &ksp) const;
  };



  /**
   * An implementation of the solver interface using the PETSc TFQMR
   * solver.
   *
   * @ingroup PETScWrappers
   * @author Wolfgang Bangerth, 2004
   */
  class SolverTFQMR : public SolverBase
  {
  public:
    /**
     * Standardized data struct to
     * pipe additional data to the
     * solver.
     */
    struct AdditionalData
    {};

    /**
     * Constructor. In contrast to
     * deal.II's own solvers, there is no
     * need to give a vector memory
     * object. However, PETSc solvers want
     * to have an MPI communicator context
     * over which computations are
     * parallelized. By default,
     * @p PETSC_COMM_SELF is used here,
     * but you can change this. Note that
     * for single processor (non-MPI)
     * versions, this parameter does not
     * have any effect.
     *
     * The last argument takes a structure
     * with additional, solver dependent
     * flags for tuning.
     *
     * Note that the communicator used here
     * must match the communicator used in
     * the system matrix, solution, and
     * right hand side object of the solve
     * to be done with this
     * solver. Otherwise, PETSc will
     * generate hard to track down errors,
     * see the documentation of the
     * SolverBase class.
     */
    SolverTFQMR (SolverControl        &cn,
                 const MPI_Comm       &mpi_communicator = PETSC_COMM_SELF,
                 const AdditionalData &data = AdditionalData());

  protected:
    /**
     * Store a copy of the flags for this
     * particular solver.
     */
    const AdditionalData additional_data;

    /**
     * Function that takes a Krylov
     * Subspace Solver context object, and
     * sets the type of solver that is
     *appropriate for this class.
     */
    virtual void set_solver_type (KSP &ksp) const;
  };




  /**
   * An implementation of the solver interface using the PETSc TFQMR-2 solver
   * (called TCQMR in PETSc). Note that this solver had a serious bug in
   * versions up to and including PETSc 2.1.6, in that it did not check
   * convergence and always returned an error code. Thus, this class will abort
   * with an error indicating failure to converge with PETSc 2.1.6 and
   * prior. This should be fixed in later versions of PETSc, though.
   *
   * @ingroup PETScWrappers
   * @author Wolfgang Bangerth, 2004
   */
  class SolverTCQMR : public SolverBase
  {
  public:
    /**
     * Standardized data struct to
     * pipe additional data to the
     * solver.
     */
    struct AdditionalData
    {};

    /**
     * Constructor. In contrast to
     * deal.II's own solvers, there is no
     * need to give a vector memory
     * object. However, PETSc solvers want
     * to have an MPI communicator context
     * over which computations are
     * parallelized. By default,
     * @p PETSC_COMM_SELF is used here,
     * but you can change this. Note that
     * for single processor (non-MPI)
     * versions, this parameter does not
     * have any effect.
     *
     * The last argument takes a structure
     * with additional, solver dependent
     * flags for tuning.
     *
     * Note that the communicator used here
     * must match the communicator used in
     * the system matrix, solution, and
     * right hand side object of the solve
     * to be done with this
     * solver. Otherwise, PETSc will
     * generate hard to track down errors,
     * see the documentation of the
     * SolverBase class.
     */
    SolverTCQMR (SolverControl        &cn,
                 const MPI_Comm       &mpi_communicator = PETSC_COMM_SELF,
                 const AdditionalData &data = AdditionalData());

  protected:
    /**
     * Store a copy of the flags for this
     * particular solver.
     */
    const AdditionalData additional_data;

    /**
     * Function that takes a Krylov
     * Subspace Solver context object, and
     * sets the type of solver that is
     *appropriate for this class.
     */
    virtual void set_solver_type (KSP &ksp) const;
  };



  /**
   * An implementation of the solver interface using the PETSc CR
   * solver.
   *
   * @ingroup PETScWrappers
   * @author Wolfgang Bangerth, 2004
   */
  class SolverCR : public SolverBase
  {
  public:
    /**
     * Standardized data struct to
     * pipe additional data to the
     * solver.
     */
    struct AdditionalData
    {};

    /**
     * Constructor. In contrast to
     * deal.II's own solvers, there is no
     * need to give a vector memory
     * object. However, PETSc solvers want
     * to have an MPI communicator context
     * over which computations are
     * parallelized. By default,
     * @p PETSC_COMM_SELF is used here,
     * but you can change this. Note that
     * for single processor (non-MPI)
     * versions, this parameter does not
     * have any effect.
     *
     * The last argument takes a structure
     * with additional, solver dependent
     * flags for tuning.
     *
     * Note that the communicator used here
     * must match the communicator used in
     * the system matrix, solution, and
     * right hand side object of the solve
     * to be done with this
     * solver. Otherwise, PETSc will
     * generate hard to track down errors,
     * see the documentation of the
     * SolverBase class.
     */
    SolverCR (SolverControl        &cn,
              const MPI_Comm       &mpi_communicator = PETSC_COMM_SELF,
              const AdditionalData &data = AdditionalData());

  protected:
    /**
     * Store a copy of the flags for this
     * particular solver.
     */
    const AdditionalData additional_data;

    /**
     * Function that takes a Krylov
     * Subspace Solver context object, and
     * sets the type of solver that is
     *appropriate for this class.
     */
    virtual void set_solver_type (KSP &ksp) const;
  };



  /**
   * An implementation of the solver interface using the PETSc Least Squares
   * solver.
   *
   * @ingroup PETScWrappers
   * @author Wolfgang Bangerth, 2004
   */
  class SolverLSQR : public SolverBase
  {
  public:
    /**
     * Standardized data struct to
     * pipe additional data to the
     * solver.
     */
    struct AdditionalData
    {};

    /**
     * Constructor. In contrast to
     * deal.II's own solvers, there is no
     * need to give a vector memory
     * object. However, PETSc solvers want
     * to have an MPI communicator context
     * over which computations are
     * parallelized. By default,
     * @p PETSC_COMM_SELF is used here,
     * but you can change this. Note that
     * for single processor (non-MPI)
     * versions, this parameter does not
     * have any effect.
     *
     * The last argument takes a structure
     * with additional, solver dependent
     * flags for tuning.
     *
     * Note that the communicator used here
     * must match the communicator used in
     * the system matrix, solution, and
     * right hand side object of the solve
     * to be done with this
     * solver. Otherwise, PETSc will
     * generate hard to track down errors,
     * see the documentation of the
     * SolverBase class.
     */
    SolverLSQR (SolverControl        &cn,
                const MPI_Comm       &mpi_communicator = PETSC_COMM_SELF,
                const AdditionalData &data = AdditionalData());

  protected:
    /**
     * Store a copy of the flags for this
     * particular solver.
     */
    const AdditionalData additional_data;

    /**
     * Function that takes a Krylov
     * Subspace Solver context object, and
     * sets the type of solver that is
     *appropriate for this class.
     */
    virtual void set_solver_type (KSP &ksp) const;
  };


  /**
   * An implementation of the solver interface using the PETSc PREONLY
   * solver. Actually this is NOT a real solution algorithm. solve() only
   * applies the preconditioner once and returns immediately. Its only purpose
   * is to provide a solver object, when the preconditioner should be used as a
   * real solver. It is very useful in conjunction with the complete LU
   * decomposition preconditioner <tt> PreconditionLU </tt>, which in
   * conjunction with this solver class becomes a direct solver.
   *
   * @ingroup PETScWrappers
   * @author Wolfgang Bangerth, 2004, Oliver Kayser-Herold, 2004
   */
  class SolverPreOnly : public SolverBase
  {
  public:
    /**
     * Standardized data struct to
     * pipe additional data to the
     * solver.
     */
    struct AdditionalData
    {};

    /**
     * Constructor. In contrast to
     * deal.II's own solvers, there is no
     * need to give a vector memory
     * object. However, PETSc solvers want
     * to have an MPI communicator context
     * over which computations are
     * parallelized. By default,
     * @p PETSC_COMM_SELF is used here,
     * but you can change this. Note that
     * for single processor (non-MPI)
     * versions, this parameter does not
     * have any effect.
     *
     * The last argument takes a structure
     * with additional, solver dependent
     * flags for tuning.
     *
     * Note that the communicator used here
     * must match the communicator used in
     * the system matrix, solution, and
     * right hand side object of the solve
     * to be done with this
     * solver. Otherwise, PETSc will
     * generate hard to track down errors,
     * see the documentation of the
     * SolverBase class.
     */
    SolverPreOnly (SolverControl        &cn,
                   const MPI_Comm       &mpi_communicator = PETSC_COMM_SELF,
                   const AdditionalData &data = AdditionalData());

  protected:
    /**
     * Store a copy of the flags for this
     * particular solver.
     */
    const AdditionalData additional_data;

    /**
     * Function that takes a Krylov
     * Subspace Solver context object, and
     * sets the type of solver that is
     * appropriate for this class.
     */
    virtual void set_solver_type (KSP &ksp) const;
  };

  /**
   * An implementation of the solver interface using the sparse direct MUMPS
   * solver through PETSc. This class has the usual interface of all other
   * solver classes but it is of course different in that it doesn't implement
   * an iterative solver. As a consequence, things like the SolverControl object
   * have no particular meaning here.
   *
   * MUMPS allows to make use of symmetry in this matrix. In this class this is
   * made possible by the set_symmetric_mode() function. If your matrix is
   * symmetric, you can use this class as follows:
   * @code
   *    SolverControl cn;
   *    PETScWrappers::SparseDirectMUMPS solver(cn, mpi_communicator);
   *    solver.set_symmetric_mode(true);
   *    solver.solve(system_matrix, solution, system_rhs);
   * @endcode
   *
   * @note The class internally calls KSPSetFromOptions thus you are
   * able to use all the PETSc parameters for MATSOLVERMUMPS package.
   * See http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Mat/MATSOLVERMUMPS.html
   *
   * @ingroup PETScWrappers
   * @author Daniel Brauss, Alexander Grayver, 2012
   */
  class SparseDirectMUMPS : public SolverBase
  {
  public:
    /**
     * Standardized data structure
     * to pipe additional data to
     * the solver.
     */
    struct AdditionalData
    {};
    /**
     * Constructor
     */
    SparseDirectMUMPS (SolverControl        &cn,
                       const MPI_Comm       &mpi_communicator = PETSC_COMM_SELF,
                       const AdditionalData &data = AdditionalData());

    /**
     * The method to solve the
     * linear system.
     */
    void solve (const MatrixBase &A,
                VectorBase       &x,
                const VectorBase &b);

    /**
    * The method allows to take advantage
    * if the system matrix is symmetric by
    * using LDL^T decomposition instead of
    * more expensive LU. The argument
    * indicates whether the matrix is
    * symmetric or not.
    */
    void set_symmetric_mode (const bool flag);

  protected:
    /**
     * Store a copy of flags for this
     * particular solver.
     */
    const AdditionalData additional_data;

    virtual void set_solver_type (KSP &ksp) const;

  private:
    /**
     * A function that is used in PETSc
     * as a callback to check convergence.
     * It takes the information provided
     * from PETSc and checks it against
     * deal.II's own SolverControl objects
     * to see if convergence has been reached.
     */
    static
    PetscErrorCode convergence_test (KSP                ksp,
                                     const PetscInt     iteration,
                                     const PetscReal    residual_norm,
                                     KSPConvergedReason *reason,
                                     void               *solver_control);

    /**
     * A structure that contains the
     * PETSc solver and preconditioner
     * objects.  Since the solve member
     * function in the base is not used
     * here, the private SolverData struct
     * located in the base could not be used
     * either.
     */
    struct SolverDataMUMPS
    {
      /**
       * Destructor
       */
      ~SolverDataMUMPS ();

      KSP ksp;
      PC  pc;
    };

    std_cxx11::shared_ptr<SolverDataMUMPS> solver_data;

    /**
     * Flag specifies whether matrix
     * being factorized is symmetric
     * or not. It influences the type
     * of the used preconditioner
     * (PCLU or PCCHOLESKY)
     */
    bool symmetric_mode;
  };
}

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PETSC

/*----------------------------   petsc_solver.h     ---------------------------*/

#endif
/*----------------------------   petsc_solver.h     ---------------------------*/
