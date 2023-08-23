// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#ifndef dealii_petsc_solver_h
#define dealii_petsc_solver_h


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_PETSC

#  include <deal.II/base/smartpointer.h>

#  include <deal.II/lac/exceptions.h>
#  include <deal.II/lac/solver_control.h>

#  include <petscksp.h>

#  ifdef DEAL_II_WITH_SLEPC
#    include <deal.II/lac/slepc_spectral_transformation.h>
#  endif

DEAL_II_NAMESPACE_OPEN

// Forward declarations
#  ifndef DOXYGEN
#    ifdef DEAL_II_WITH_SLEPC
namespace SLEPcWrappers
{
  // forward declarations
  class TransformationBase;
} // namespace SLEPcWrappers
#    endif
#  endif

namespace PETScWrappers
{
  // forward declarations
#  ifndef DOXYGEN
  class MatrixBase;
  class VectorBase;
  class PreconditionBase;
#  endif


  /**
   * Base class for solver classes using the PETSc solvers. Since solvers in
   * PETSc are selected based on flags passed to a generic solver object,
   * basically all the actual solver calls happen in this class, and derived
   * classes simply set the right flags to select one solver or another, or to
   * set certain parameters for individual solvers.
   *
   * Optionally, the user can create a solver derived from the SolverBase
   * class and can set the default arguments necessary to solve the linear
   * system of equations with SolverControl. These default options can be
   * overridden by specifying command line arguments of the form @p -ksp_*.
   * For example, @p -ksp_monitor_true_residual prints out true residual norm
   * (unpreconditioned) at each iteration and @p -ksp_view provides
   * information about the linear solver and the preconditioner used in the
   * current context. The type of the solver can also be changed during
   * runtime by specifying @p -ksp_type {richardson, cg, gmres, fgmres, ..} to
   * dynamically test the optimal solver along with a suitable preconditioner
   * set using @p -pc_type {jacobi, bjacobi, ilu, lu, ..}. There are several
   * other command line options available to modify the behavior of the PETSc
   * linear solver and can be obtained from the <a
   * href="http://www.mcs.anl.gov/petsc">documentation and manual pages</a>.
   *
   * @note Command line options relative to convergence tolerances are ignored
   * when the solver is instantiated with a SolverControl.
   *
   * @note Repeated calls to solve() on a solver object with a Preconditioner
   * must be used with care. The preconditioner is initialized in the first
   * call to solve() and subsequent calls reuse the solver and preconditioner
   * object. This is done for performance reasons. The solver and
   * preconditioner can be reset by calling reset() or by using the command
   * line option "-ksp_reuse_preconditioner false".
   *
   * @ingroup PETScWrappers
   */
  class SolverBase
  {
  public:
    /**
     * Default constructor without SolverControl. The resulting solver will
     * use PETSc's default convergence testing routines.
     */
    SolverBase();

    /**
     * Constructor with a SolverControl instance.
     */
    SolverBase(SolverControl &cn);

    /**
     * Destructor.
     */
    virtual ~SolverBase();

    /**
     * Solve the linear system <tt>Ax=b</tt>. Depending on the information
     * provided by derived classes and the object passed as a preconditioner,
     * one of the linear solvers and preconditioners of PETSc is chosen.
     * Repeated calls to solve() do not reconstruct the preconditioner for
     * performance reasons. See class Documentation.
     */
    void
    solve(const MatrixBase       &A,
          VectorBase             &x,
          const VectorBase       &b,
          const PreconditionBase &preconditioner);

    /**
     * Resets the contained preconditioner and solver object. See class
     * description for more details.
     */
    virtual void
    reset();

    /**
     * Sets a prefix name for the solver object. Useful when customizing the
     * PETSc KSP object with command-line options.
     */
    void
    set_prefix(const std::string &prefix);

    /**
     * Access to object that controls convergence.
     */
    SolverControl &
    control() const;

    /**
     * initialize the solver with the preconditioner. This function is
     * intended for use with SLEPc spectral transformation class.
     */
    void
    initialize(const PreconditionBase &preconditioner);

    /**
     * Return the PETSc KSP object.
     */
    KSP
    petsc_ksp();

    /**
     * Conversion operator to gain access to the underlying PETSc type. If you
     * do this, you cut this class off some information it may need, so this
     * conversion operator should only be used if you know what you do.
     */
    operator KSP() const;

  protected:
    /**
     * The PETSc KSP object.
     */
    KSP ksp;

    /**
     * Reference to the object that controls convergence of the iterative
     * solver. In fact, for these PETSc wrappers, PETSc does so itself, but we
     * copy the data from this object before starting the solution process,
     * and copy the data back into it afterwards.
     */
    SmartPointer<SolverControl, SolverBase> solver_control;

    /**
     * Utility to create the KSP object and attach convergence test.
     */
    void
    initialize_ksp_with_comm(const MPI_Comm comm);

    /**
     * %Function that takes a Krylov Subspace Solver context object, and sets
     * the type of solver that is requested by the derived class.
     */
    virtual void
    set_solver_type(KSP &ksp) const;

    /**
     * Utility to use deal.II convergence testing.
     *
     * This call changes the convergence criterion when the instance of the
     * class has a SolverControl object associated.
     */
    void
    perhaps_set_convergence_test() const;

    /**
     * Solver prefix name to qualify options specific to the PETSc KSP object
     * in the current context. Note: A hyphen (-) must NOT be given at the
     * beginning of the prefix name. The first character of all runtime
     * options is AUTOMATICALLY the hyphen.
     */
    std::string prefix_name;

  private:
    /**
     * A function that is used in PETSc as a callback to check on convergence.
     * It takes the information provided from PETSc and checks it against
     * deal.II's own SolverControl objects to see if convergence has been
     * reached.
     */
    static PetscErrorCode
    convergence_test(KSP                 ksp,
                     const PetscInt      iteration,
                     const PetscReal     residual_norm,
                     KSPConvergedReason *reason,
                     void               *solver_control);


#  ifdef DEAL_II_WITH_SLEPC
    // Make the transformation class a friend, since it needs to set the KSP
    // solver.
    friend class SLEPcWrappers::TransformationBase;
#  endif
  };



  /**
   * An implementation of the solver interface using the PETSc Richardson
   * solver.
   *
   * @ingroup PETScWrappers
   */
  class SolverRichardson : public SolverBase
  {
  public:
    /**
     * Standardized data struct to pipe additional data to the solver.
     */
    struct AdditionalData
    {
      /**
       * Constructor. By default, set the damping parameter to one.
       */
      explicit AdditionalData(const double omega = 1);

      /**
       * Relaxation parameter.
       */
      double omega;
    };

    /**
     * Constructor. In contrast to deal.II's own solvers, there is no need to
     * give a vector memory object.
     *
     * The last argument takes a structure with additional, solver dependent
     * flags for tuning.
     */
    SolverRichardson(SolverControl        &cn,
                     const AdditionalData &data = AdditionalData());

    /**
     * Constructor. This constructor is deprecated and ignores the MPI
     * communicator argument. Use the other constructor instead.
     *
     * @deprecated
     */
    DEAL_II_DEPRECATED
    SolverRichardson(SolverControl        &cn,
                     const MPI_Comm        mpi_communicator,
                     const AdditionalData &data = AdditionalData());

  protected:
    /**
     * Store a copy of the flags for this particular solver.
     */
    const AdditionalData additional_data;

    /**
     * %Function that takes a Krylov Subspace Solver context object, and sets
     * the type of solver that is appropriate for this class.
     */
    virtual void
    set_solver_type(KSP &ksp) const override;
  };



  /**
   * An implementation of the solver interface using the PETSc Chebyshev (or,
   * prior version 3.3, Chebychev) solver.
   *
   * @ingroup PETScWrappers
   */
  class SolverChebychev : public SolverBase
  {
  public:
    /**
     * Standardized data struct to pipe additional data to the solver.
     */
    struct AdditionalData
    {};

    /**
     * Constructor. In contrast to deal.II's own solvers, there is no need to
     * give a vector memory object.
     *
     * The last argument takes a structure with additional, solver dependent
     * flags for tuning.
     */
    SolverChebychev(SolverControl        &cn,
                    const AdditionalData &data = AdditionalData());

    /**
     * Constructor. This constructor is deprecated and ignores the MPI
     * communicator argument. Use the other constructor instead.
     *
     * @deprecated
     */
    DEAL_II_DEPRECATED
    SolverChebychev(SolverControl        &cn,
                    const MPI_Comm        mpi_communicator,
                    const AdditionalData &data = AdditionalData());

  protected:
    /**
     * Store a copy of the flags for this particular solver.
     */
    const AdditionalData additional_data;

    /**
     * %Function that takes a Krylov Subspace Solver context object, and sets
     * the type of solver that is appropriate for this class.
     */
    virtual void
    set_solver_type(KSP &ksp) const override;
  };



  /**
   * An implementation of the solver interface using the PETSc CG solver.
   *
   * @ingroup PETScWrappers
   */
  class SolverCG : public SolverBase
  {
  public:
    /**
     * Standardized data struct to pipe additional data to the solver.
     */
    struct AdditionalData
    {};

    /**
     * Constructor. In contrast to deal.II's own solvers, there is no need to
     * give a vector memory object.
     *
     * The last argument takes a structure with additional, solver dependent
     * flags for tuning.
     */
    SolverCG(SolverControl &cn, const AdditionalData &data = AdditionalData());

    /**
     * Constructor. This constructor is deprecated and ignores the MPI
     * communicator argument. Use the other constructor instead.
     *
     * @deprecated
     */
    DEAL_II_DEPRECATED
    SolverCG(SolverControl        &cn,
             const MPI_Comm        mpi_communicator,
             const AdditionalData &data = AdditionalData());

  protected:
    /**
     * Store a copy of the flags for this particular solver.
     */
    const AdditionalData additional_data;

    /**
     * %Function that takes a Krylov Subspace Solver context object, and sets
     * the type of solver that is appropriate for this class.
     */
    virtual void
    set_solver_type(KSP &ksp) const override;
  };



  /**
   * An implementation of the solver interface using the PETSc BiCG solver.
   *
   * @ingroup PETScWrappers
   */
  class SolverBiCG : public SolverBase
  {
  public:
    /**
     * Standardized data struct to pipe additional data to the solver.
     */
    struct AdditionalData
    {};

    /**
     * Constructor. In contrast to deal.II's own solvers, there is no need to
     * give a vector memory object.
     *
     * The last argument takes a structure with additional, solver dependent
     * flags for tuning.
     */
    SolverBiCG(SolverControl        &cn,
               const AdditionalData &data = AdditionalData());

    /**
     * Constructor. This constructor is deprecated and ignores the MPI
     * communicator argument. Use the other constructor instead.
     *
     * @deprecated
     */
    DEAL_II_DEPRECATED
    SolverBiCG(SolverControl        &cn,
               const MPI_Comm        mpi_communicator,
               const AdditionalData &data = AdditionalData());

  protected:
    /**
     * Store a copy of the flags for this particular solver.
     */
    const AdditionalData additional_data;

    /**
     * %Function that takes a Krylov Subspace Solver context object, and sets
     * the type of solver that is appropriate for this class.
     */
    virtual void
    set_solver_type(KSP &ksp) const override;
  };



  /**
   * An implementation of the solver interface using the PETSc GMRES solver.
   *
   * @ingroup PETScWrappers
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
      AdditionalData(const unsigned int restart_parameter     = 30,
                     const bool         right_preconditioning = false);

      /**
       * Maximum number of tmp vectors.
       */
      unsigned int restart_parameter;

      /**
       * Flag for right preconditioning.
       */
      bool right_preconditioning;
    };

    /**
     * Constructor. In contrast to deal.II's own solvers, there is no need to
     * give a vector memory object.
     *
     * The last argument takes a structure with additional, solver dependent
     * flags for tuning.
     */
    SolverGMRES(SolverControl        &cn,
                const AdditionalData &data = AdditionalData());

    /**
     * Constructor. This constructor is deprecated and ignores the MPI
     * communicator argument. Use the other constructor instead.
     *
     * @deprecated
     */
    DEAL_II_DEPRECATED
    SolverGMRES(SolverControl        &cn,
                const MPI_Comm        mpi_communicator,
                const AdditionalData &data = AdditionalData());

  protected:
    /**
     * Store a copy of the flags for this particular solver.
     */
    const AdditionalData additional_data;

    /**
     * %Function that takes a Krylov Subspace Solver context object, and sets
     * the type of solver that is appropriate for this class.
     */
    virtual void
    set_solver_type(KSP &ksp) const override;
  };



  /**
   * An implementation of the solver interface using the PETSc BiCGStab
   * solver.
   *
   * @ingroup PETScWrappers
   */
  class SolverBicgstab : public SolverBase
  {
  public:
    /**
     * Standardized data struct to pipe additional data to the solver.
     */
    struct AdditionalData
    {};

    /**
     * Constructor. In contrast to deal.II's own solvers, there is no need to
     * give a vector memory object.
     *
     * The last argument takes a structure with additional, solver dependent
     * flags for tuning.
     */
    SolverBicgstab(SolverControl        &cn,
                   const AdditionalData &data = AdditionalData());

    /**
     * Constructor. This constructor is deprecated and ignores the MPI
     * communicator argument. Use the other constructor instead.
     *
     * @deprecated
     */
    DEAL_II_DEPRECATED
    SolverBicgstab(SolverControl        &cn,
                   const MPI_Comm        mpi_communicator,
                   const AdditionalData &data = AdditionalData());

  protected:
    /**
     * Store a copy of the flags for this particular solver.
     */
    const AdditionalData additional_data;

    /**
     * %Function that takes a Krylov Subspace Solver context object, and sets
     * the type of solver that is appropriate for this class.
     */
    virtual void
    set_solver_type(KSP &ksp) const override;
  };



  /**
   * An implementation of the solver interface using the PETSc CG Squared
   * solver.
   *
   * @ingroup PETScWrappers
   */
  class SolverCGS : public SolverBase
  {
  public:
    /**
     * Standardized data struct to pipe additional data to the solver.
     */
    struct AdditionalData
    {};

    /**
     * Constructor. In contrast to deal.II's own solvers, there is no need to
     * give a vector memory object.
     *
     * The last argument takes a structure with additional, solver dependent
     * flags for tuning.
     */
    SolverCGS(SolverControl &cn, const AdditionalData &data = AdditionalData());

    /**
     * Constructor. This constructor is deprecated and ignores the MPI
     * communicator argument. Use the other constructor instead.
     *
     * @deprecated
     */
    DEAL_II_DEPRECATED
    SolverCGS(SolverControl        &cn,
              const MPI_Comm        mpi_communicator,
              const AdditionalData &data = AdditionalData());

  protected:
    /**
     * Store a copy of the flags for this particular solver.
     */
    const AdditionalData additional_data;

    /**
     * %Function that takes a Krylov Subspace Solver context object, and sets
     * the type of solver that is appropriate for this class.
     */
    virtual void
    set_solver_type(KSP &ksp) const override;
  };



  /**
   * An implementation of the solver interface using the PETSc TFQMR solver.
   *
   * @ingroup PETScWrappers
   */
  class SolverTFQMR : public SolverBase
  {
  public:
    /**
     * Standardized data struct to pipe additional data to the solver.
     */
    struct AdditionalData
    {};

    /**
     * Constructor. In contrast to deal.II's own solvers, there is no need to
     * give a vector memory object.
     *
     * The last argument takes a structure with additional, solver dependent
     * flags for tuning.
     */
    SolverTFQMR(SolverControl        &cn,
                const AdditionalData &data = AdditionalData());

    /**
     * Constructor. This constructor is deprecated and ignores the MPI
     * communicator argument. Use the other constructor instead.
     *
     * @deprecated
     */
    DEAL_II_DEPRECATED
    SolverTFQMR(SolverControl        &cn,
                const MPI_Comm        mpi_communicator,
                const AdditionalData &data = AdditionalData());

  protected:
    /**
     * Store a copy of the flags for this particular solver.
     */
    const AdditionalData additional_data;

    /**
     * %Function that takes a Krylov Subspace Solver context object, and sets
     * the type of solver that is appropriate for this class.
     */
    virtual void
    set_solver_type(KSP &ksp) const override;
  };



  /**
   * An implementation of the solver interface using the PETSc TFQMR-2 solver
   * (called TCQMR in PETSc). Note that this solver had a serious bug in
   * versions up to and including PETSc 2.1.6, in that it did not check
   * convergence and always returned an error code. Thus, this class will
   * abort with an error indicating failure to converge with PETSc 2.1.6 and
   * prior. This should be fixed in later versions of PETSc, though.
   *
   * @ingroup PETScWrappers
   */
  class SolverTCQMR : public SolverBase
  {
  public:
    /**
     * Standardized data struct to pipe additional data to the solver.
     */
    struct AdditionalData
    {};

    /**
     * Constructor. In contrast to deal.II's own solvers, there is no need to
     * give a vector memory object.
     *
     * The last argument takes a structure with additional, solver dependent
     * flags for tuning.
     */
    SolverTCQMR(SolverControl        &cn,
                const AdditionalData &data = AdditionalData());

    /**
     * Constructor. This constructor is deprecated and ignores the MPI
     * communicator argument. Use the other constructor instead.
     *
     * @deprecated
     */
    DEAL_II_DEPRECATED
    SolverTCQMR(SolverControl        &cn,
                const MPI_Comm        mpi_communicator,
                const AdditionalData &data = AdditionalData());

  protected:
    /**
     * Store a copy of the flags for this particular solver.
     */
    const AdditionalData additional_data;

    /**
     * %Function that takes a Krylov Subspace Solver context object, and sets
     * the type of solver that is appropriate for this class.
     */
    virtual void
    set_solver_type(KSP &ksp) const override;
  };



  /**
   * An implementation of the solver interface using the PETSc CR solver.
   *
   * @ingroup PETScWrappers
   */
  class SolverCR : public SolverBase
  {
  public:
    /**
     * Standardized data struct to pipe additional data to the solver.
     */
    struct AdditionalData
    {};

    /**
     * Constructor. In contrast to deal.II's own solvers, there is no need to
     * give a vector memory object.
     *
     * The last argument takes a structure with additional, solver dependent
     * flags for tuning.
     */
    SolverCR(SolverControl &cn, const AdditionalData &data = AdditionalData());

    /**
     * Constructor. This constructor is deprecated and ignores the MPI
     * communicator argument. Use the other constructor instead.
     *
     * @deprecated
     */
    DEAL_II_DEPRECATED
    SolverCR(SolverControl        &cn,
             const MPI_Comm        mpi_communicator,
             const AdditionalData &data = AdditionalData());

  protected:
    /**
     * Store a copy of the flags for this particular solver.
     */
    const AdditionalData additional_data;

    /**
     * %Function that takes a Krylov Subspace Solver context object, and sets
     * the type of solver that is appropriate for this class.
     */
    virtual void
    set_solver_type(KSP &ksp) const override;
  };



  /**
   * An implementation of the solver interface using the PETSc Least Squares
   * solver.
   *
   * @ingroup PETScWrappers
   */
  class SolverLSQR : public SolverBase
  {
  public:
    /**
     * Standardized data struct to pipe additional data to the solver.
     */
    struct AdditionalData
    {};

    /**
     * Constructor. In contrast to deal.II's own solvers, there is no need to
     * give a vector memory object.
     *
     * The last argument takes a structure with additional, solver dependent
     * flags for tuning.
     */
    SolverLSQR(SolverControl        &cn,
               const AdditionalData &data = AdditionalData());

    /**
     * Constructor. This constructor is deprecated and ignores the MPI
     * communicator argument. Use the other constructor instead.
     *
     * @deprecated
     */
    DEAL_II_DEPRECATED
    SolverLSQR(SolverControl        &cn,
               const MPI_Comm        mpi_communicator,
               const AdditionalData &data = AdditionalData());

  protected:
    /**
     * Store a copy of the flags for this particular solver.
     */
    const AdditionalData additional_data;

    /**
     * %Function that takes a Krylov Subspace Solver context object, and sets
     * the type of solver that is appropriate for this class.
     */
    virtual void
    set_solver_type(KSP &ksp) const override;
  };


  /**
   * An implementation of the solver interface using the PETSc PREONLY solver.
   * Actually this is NOT a real solution algorithm. solve() only applies the
   * preconditioner once and returns immediately. Its only purpose is to
   * provide a solver object, when the preconditioner should be used as a real
   * solver. It is very useful in conjunction with the complete LU
   * decomposition preconditioner <tt> PreconditionLU </tt>, which in
   * conjunction with this solver class becomes a direct solver.
   *
   * @ingroup PETScWrappers
   */
  class SolverPreOnly : public SolverBase
  {
  public:
    /**
     * Standardized data struct to pipe additional data to the solver.
     */
    struct AdditionalData
    {};

    /**
     * Constructor. In contrast to deal.II's own solvers, there is no need to
     * give a vector memory object.
     *
     * The last argument takes a structure with additional, solver dependent
     * flags for tuning.
     */
    SolverPreOnly(SolverControl        &cn,
                  const AdditionalData &data = AdditionalData());

    /**
     * Constructor. This constructor is deprecated and ignores the MPI
     * communicator argument. Use the other constructor instead.
     *
     * @deprecated
     */
    DEAL_II_DEPRECATED
    SolverPreOnly(SolverControl        &cn,
                  const MPI_Comm        mpi_communicator,
                  const AdditionalData &data = AdditionalData());

  protected:
    /**
     * Store a copy of the flags for this particular solver.
     */
    const AdditionalData additional_data;

    /**
     * %Function that takes a Krylov Subspace Solver context object, and sets
     * the type of solver that is appropriate for this class.
     */
    virtual void
    set_solver_type(KSP &ksp) const override;
  };

  /**
   * An implementation of the solver interface using the sparse direct MUMPS
   * solver through PETSc. This class has the usual interface of all other
   * solver classes but it is of course different in that it doesn't implement
   * an iterative solver. As a consequence, things like the SolverControl
   * object have no particular meaning here.
   *
   * MUMPS allows to make use of symmetry in this matrix. In this class this
   * is made possible by the set_symmetric_mode() function. If your matrix is
   * symmetric, you can use this class as follows:
   * @code
   *   SolverControl cn;
   *   PETScWrappers::SparseDirectMUMPS solver(cn, mpi_communicator);
   *   solver.set_symmetric_mode(true);
   *   solver.solve(system_matrix, solution, system_rhs);
   * @endcode
   *
   * @note The class internally calls KSPSetFromOptions thus you are able to
   * use all the PETSc parameters for MATSOLVERMUMPS package. See
   * http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Mat/MATSOLVERMUMPS.html
   *
   * @ingroup PETScWrappers
   */
  class SparseDirectMUMPS : public SolverBase
  {
  public:
    /**
     * Standardized data structure to pipe additional data to the solver.
     */
    struct AdditionalData
    {};

    /**
     * Constructor.
     */
    SparseDirectMUMPS(SolverControl        &cn,
                      const AdditionalData &data = AdditionalData());

    /**
     * Constructor. This constructor is deprecated and ignores the MPI
     * communicator argument. Use the other constructor instead.
     *
     * @deprecated
     */
    DEAL_II_DEPRECATED
    SparseDirectMUMPS(SolverControl        &cn,
                      const MPI_Comm        mpi_communicator,
                      const AdditionalData &data = AdditionalData());

    /**
     * The method to solve the linear system.
     */
    void
    solve(const MatrixBase &A, VectorBase &x, const VectorBase &b);

    /**
     * The method allows to take advantage if the system matrix is symmetric
     * by using LDL^T decomposition instead of more expensive LU. The argument
     * indicates whether the matrix is symmetric or not.
     */
    void
    set_symmetric_mode(const bool flag);

  protected:
    /**
     * Store a copy of flags for this particular solver.
     */
    const AdditionalData additional_data;

    virtual void
    set_solver_type(KSP &ksp) const override;

    /**
     * Flag specifies whether matrix being factorized is symmetric or not. It
     * influences the type of the used preconditioner (PCLU or PCCHOLESKY)
     */
    bool symmetric_mode;
  };
} // namespace PETScWrappers

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PETSC

#endif
