// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2019 by the deal.II authors
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


#include <deal.II/base/logstream.h>

#include <deal.II/lac/petsc_solver.h>

#ifdef DEAL_II_WITH_PETSC

#  include <deal.II/lac/exceptions.h>
#  include <deal.II/lac/petsc_compatibility.h>
#  include <deal.II/lac/petsc_matrix_base.h>
#  include <deal.II/lac/petsc_precondition.h>
#  include <deal.II/lac/petsc_vector_base.h>

#  include <petscversion.h>

#  include <cmath>
#  include <memory>

DEAL_II_NAMESPACE_OPEN

namespace PETScWrappers
{
  SolverBase::SolverData::~SolverData()
  {
    destroy_krylov_solver(ksp);
  }



  SolverBase::SolverBase(SolverControl &cn, const MPI_Comm &mpi_communicator)
    : solver_control(cn)
    , mpi_communicator(mpi_communicator)
  {}



  void
  SolverBase::solve(const MatrixBase &      A,
                    VectorBase &            x,
                    const VectorBase &      b,
                    const PreconditionBase &preconditioner)
  {
    /*
      TODO: PETSc duplicates communicators, so this does not work (you put
    MPI_COMM_SELF in, but get something other out when you ask PETSc for the
    communicator. This mainly fails due to the MatrixFree classes, that can not
    ask PETSc for a communicator. //Timo Heister
    Assert(A.get_mpi_communicator()==mpi_communicator, ExcMessage("PETSc Solver
    and Matrix need to use the same MPI_Comm."));
    Assert(x.get_mpi_communicator()==mpi_communicator, ExcMessage("PETSc Solver
    and Vector need to use the same MPI_Comm."));
    Assert(b.get_mpi_communicator()==mpi_communicator, ExcMessage("PETSc Solver
    and Vector need to use the same MPI_Comm."));
    */

    // first create a solver object if this
    // is necessary
    if (solver_data.get() == nullptr)
      {
        solver_data = std::make_unique<SolverData>();

        PetscErrorCode ierr = KSPCreate(mpi_communicator, &solver_data->ksp);
        AssertThrow(ierr == 0, ExcPETScError(ierr));

        // let derived classes set the solver
        // type, and the preconditioning
        // object set the type of
        // preconditioner
        set_solver_type(solver_data->ksp);

        ierr = KSPSetPC(solver_data->ksp, preconditioner.get_pc());
        AssertThrow(ierr == 0, ExcPETScError(ierr));

        // make sure the preconditioner has an associated matrix set
        const Mat B = preconditioner;
        AssertThrow(B != nullptr,
                    ExcMessage("PETSc preconditioner should have an "
                               "associated matrix set to be used in solver."));

        // setting the preconditioner overwrites the used matrices.
        // hence, we need to set the matrices after the preconditioner.
#  if DEAL_II_PETSC_VERSION_LT(3, 5, 0)
        // the last argument is irrelevant here,
        // since we use the solver only once anyway
        ierr = KSPSetOperators(solver_data->ksp,
                               A,
                               preconditioner,
                               SAME_PRECONDITIONER);
#  else
        ierr = KSPSetOperators(solver_data->ksp, A, preconditioner);
#  endif
        AssertThrow(ierr == 0, ExcPETScError(ierr));

        // then a convergence monitor
        // function. that function simply
        // checks with the solver_control
        // object we have in this object for
        // convergence
        ierr = KSPSetConvergenceTest(solver_data->ksp,
                                     &convergence_test,
                                     reinterpret_cast<void *>(&solver_control),
                                     nullptr);
        AssertThrow(ierr == 0, ExcPETScError(ierr));
      }

    // set the command line option prefix name
    PetscErrorCode ierr =
      KSPSetOptionsPrefix(solver_data->ksp, prefix_name.c_str());
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    // set the command line options provided
    // by the user to override the defaults
    ierr = KSPSetFromOptions(solver_data->ksp);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    // then do the real work: set up solver
    // internal data and solve the
    // system.
    ierr = KSPSetUp(solver_data->ksp);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    ierr = KSPSolve(solver_data->ksp, b, x);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    // do not destroy solver object
    //    solver_data.reset ();

    // in case of failure: throw
    // exception
    if (solver_control.last_check() != SolverControl::success)
      AssertThrow(false,
                  SolverControl::NoConvergence(solver_control.last_step(),
                                               solver_control.last_value()));
    // otherwise exit as normal
  }


  void
  SolverBase::set_prefix(const std::string &prefix)
  {
    prefix_name = prefix;
  }


  void
  SolverBase::reset()
  {
    solver_data.reset();
  }


  SolverControl &
  SolverBase::control() const
  {
    return solver_control;
  }


  int
  SolverBase::convergence_test(KSP /*ksp*/,
                               const PetscInt      iteration,
                               const PetscReal     residual_norm,
                               KSPConvergedReason *reason,
                               void *              solver_control_x)
  {
    SolverControl &solver_control =
      *reinterpret_cast<SolverControl *>(solver_control_x);

    const SolverControl::State state =
      solver_control.check(iteration, residual_norm);

    switch (state)
      {
        case ::dealii::SolverControl::iterate:
          *reason = KSP_CONVERGED_ITERATING;
          break;

        case ::dealii::SolverControl::success:
          *reason = static_cast<KSPConvergedReason>(1);
          break;

        case ::dealii::SolverControl::failure:
          if (solver_control.last_step() > solver_control.max_steps())
            *reason = KSP_DIVERGED_ITS;
          else
            *reason = KSP_DIVERGED_DTOL;
          break;

        default:
          Assert(false, ExcNotImplemented());
      }

    // return without failure
    return 0;
  }

  void
  SolverBase::initialize(const PreconditionBase &preconditioner)
  {
    PetscErrorCode ierr;

    solver_data = std::make_unique<SolverData>();

    ierr = KSPCreate(mpi_communicator, &solver_data->ksp);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    // let derived classes set the solver
    // type, and the preconditioning
    // object set the type of
    // preconditioner
    set_solver_type(solver_data->ksp);

    ierr = KSPSetPC(solver_data->ksp, preconditioner.get_pc());
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    // then a convergence monitor
    // function. that function simply
    // checks with the solver_control
    // object we have in this object for
    // convergence
    ierr = KSPSetConvergenceTest(solver_data->ksp,
                                 &convergence_test,
                                 reinterpret_cast<void *>(&solver_control),
                                 nullptr);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    // set the command line options provided
    // by the user to override the defaults
    ierr = KSPSetFromOptions(solver_data->ksp);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }



  /* ---------------------- SolverRichardson ------------------------ */

  SolverRichardson::AdditionalData::AdditionalData(const double omega)
    : omega(omega)
  {}



  SolverRichardson::SolverRichardson(SolverControl &       cn,
                                     const MPI_Comm &      mpi_communicator,
                                     const AdditionalData &data)
    : SolverBase(cn, mpi_communicator)
    , additional_data(data)
  {}


  void
  SolverRichardson::set_solver_type(KSP &ksp) const
  {
    PetscErrorCode ierr = KSPSetType(ksp, KSPRICHARDSON);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    // set the damping factor from the data
    ierr = KSPRichardsonSetScale(ksp, additional_data.omega);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    // in the deal.II solvers, we always
    // honor the initial guess in the
    // solution vector. do so here as well:
    ierr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    // Hand over the absolute
    // tolerance and the maximum
    // iteration number to the PETSc
    // convergence criterion. The
    // custom deal.II SolverControl
    // object is ignored by the PETSc
    // Richardson method (when no
    // PETSc monitoring is present),
    // since in this case PETSc
    // uses a faster version of
    // the Richardson iteration,
    // where no residual is
    // available.
    ierr = KSPSetTolerances(ksp,
                            PETSC_DEFAULT,
                            this->solver_control.tolerance(),
                            PETSC_DEFAULT,
                            this->solver_control.max_steps() + 1);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }


  /* ---------------------- SolverChebychev ------------------------ */

  SolverChebychev::SolverChebychev(SolverControl &       cn,
                                   const MPI_Comm &      mpi_communicator,
                                   const AdditionalData &data)
    : SolverBase(cn, mpi_communicator)
    , additional_data(data)
  {}


  void
  SolverChebychev::set_solver_type(KSP &ksp) const
  {
    PetscErrorCode ierr = KSPSetType(ksp, KSPCHEBYSHEV);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    // in the deal.II solvers, we always
    // honor the initial guess in the
    // solution vector. do so here as well:
    ierr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }


  /* ---------------------- SolverCG ------------------------ */

  SolverCG::SolverCG(SolverControl &       cn,
                     const MPI_Comm &      mpi_communicator,
                     const AdditionalData &data)
    : SolverBase(cn, mpi_communicator)
    , additional_data(data)
  {}


  void
  SolverCG::set_solver_type(KSP &ksp) const
  {
    PetscErrorCode ierr = KSPSetType(ksp, KSPCG);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    // in the deal.II solvers, we always
    // honor the initial guess in the
    // solution vector. do so here as well:
    ierr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }


  /* ---------------------- SolverBiCG ------------------------ */

  SolverBiCG::SolverBiCG(SolverControl &       cn,
                         const MPI_Comm &      mpi_communicator,
                         const AdditionalData &data)
    : SolverBase(cn, mpi_communicator)
    , additional_data(data)
  {}


  void
  SolverBiCG::set_solver_type(KSP &ksp) const
  {
    PetscErrorCode ierr = KSPSetType(ksp, KSPBICG);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    // in the deal.II solvers, we always
    // honor the initial guess in the
    // solution vector. do so here as well:
    ierr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }


  /* ---------------------- SolverGMRES ------------------------ */

  SolverGMRES::AdditionalData::AdditionalData(
    const unsigned int restart_parameter,
    const bool         right_preconditioning)
    : restart_parameter(restart_parameter)
    , right_preconditioning(right_preconditioning)
  {}



  SolverGMRES::SolverGMRES(SolverControl &       cn,
                           const MPI_Comm &      mpi_communicator,
                           const AdditionalData &data)
    : SolverBase(cn, mpi_communicator)
    , additional_data(data)
  {}


  void
  SolverGMRES::set_solver_type(KSP &ksp) const
  {
    PetscErrorCode ierr = KSPSetType(ksp, KSPGMRES);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    // set the restart parameter from the
    // data. we would like to use the simple
    // code that is commented out, but this
    // leads to nasty warning and error
    // messages due to some stupidity on
    // PETSc's side: KSPGMRESSetRestart is
    // implemented as a macro in which return
    // statements are hidden. This may work
    // if people strictly follow the PETSc
    // coding style of always having
    // functions return an integer error
    // code, but the present function isn't
    // like this.
    /*
        ierr = KSPGMRESSetRestart (ksp, additional_data.restart_parameter);
        AssertThrow (ierr == 0, ExcPETScError(ierr));
    */
    // so rather expand their macros by hand,
    // and do some equally nasty stuff that at
    // least doesn't yield warnings...
    int (*fun_ptr)(KSP, int);
    ierr = PetscObjectQueryFunction(reinterpret_cast<PetscObject>(ksp),
                                    "KSPGMRESSetRestart_C",
                                    reinterpret_cast<void (**)()>(&fun_ptr));
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    ierr = (*fun_ptr)(ksp, additional_data.restart_parameter);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    // Set preconditioning side to
    // right
    if (additional_data.right_preconditioning)
      {
        ierr = KSPSetPCSide(ksp, PC_RIGHT);
        AssertThrow(ierr == 0, ExcPETScError(ierr));
      }

    // in the deal.II solvers, we always
    // honor the initial guess in the
    // solution vector. do so here as well:
    ierr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }


  /* ---------------------- SolverBicgstab ------------------------ */

  SolverBicgstab::SolverBicgstab(SolverControl &       cn,
                                 const MPI_Comm &      mpi_communicator,
                                 const AdditionalData &data)
    : SolverBase(cn, mpi_communicator)
    , additional_data(data)
  {}


  void
  SolverBicgstab::set_solver_type(KSP &ksp) const
  {
    PetscErrorCode ierr = KSPSetType(ksp, KSPBCGS);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    // in the deal.II solvers, we always
    // honor the initial guess in the
    // solution vector. do so here as well:
    ierr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }


  /* ---------------------- SolverCGS ------------------------ */

  SolverCGS::SolverCGS(SolverControl &       cn,
                       const MPI_Comm &      mpi_communicator,
                       const AdditionalData &data)
    : SolverBase(cn, mpi_communicator)
    , additional_data(data)
  {}


  void
  SolverCGS::set_solver_type(KSP &ksp) const
  {
    PetscErrorCode ierr = KSPSetType(ksp, KSPCGS);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    // in the deal.II solvers, we always
    // honor the initial guess in the
    // solution vector. do so here as well:
    ierr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }


  /* ---------------------- SolverTFQMR ------------------------ */

  SolverTFQMR::SolverTFQMR(SolverControl &       cn,
                           const MPI_Comm &      mpi_communicator,
                           const AdditionalData &data)
    : SolverBase(cn, mpi_communicator)
    , additional_data(data)
  {}


  void
  SolverTFQMR::set_solver_type(KSP &ksp) const
  {
    PetscErrorCode ierr = KSPSetType(ksp, KSPTFQMR);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    // in the deal.II solvers, we always
    // honor the initial guess in the
    // solution vector. do so here as well:
    ierr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }


  /* ---------------------- SolverTCQMR ------------------------ */

  SolverTCQMR::SolverTCQMR(SolverControl &       cn,
                           const MPI_Comm &      mpi_communicator,
                           const AdditionalData &data)
    : SolverBase(cn, mpi_communicator)
    , additional_data(data)
  {}


  void
  SolverTCQMR::set_solver_type(KSP &ksp) const
  {
    PetscErrorCode ierr = KSPSetType(ksp, KSPTCQMR);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    // in the deal.II solvers, we always
    // honor the initial guess in the
    // solution vector. do so here as well:
    ierr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }


  /* ---------------------- SolverCR ------------------------ */

  SolverCR::SolverCR(SolverControl &       cn,
                     const MPI_Comm &      mpi_communicator,
                     const AdditionalData &data)
    : SolverBase(cn, mpi_communicator)
    , additional_data(data)
  {}


  void
  SolverCR::set_solver_type(KSP &ksp) const
  {
    PetscErrorCode ierr = KSPSetType(ksp, KSPCR);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    // in the deal.II solvers, we always
    // honor the initial guess in the
    // solution vector. do so here as well:
    ierr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }


  /* ---------------------- SolverLSQR ------------------------ */

  SolverLSQR::SolverLSQR(SolverControl &       cn,
                         const MPI_Comm &      mpi_communicator,
                         const AdditionalData &data)
    : SolverBase(cn, mpi_communicator)
    , additional_data(data)
  {}


  void
  SolverLSQR::set_solver_type(KSP &ksp) const
  {
    PetscErrorCode ierr = KSPSetType(ksp, KSPLSQR);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    // in the deal.II solvers, we always
    // honor the initial guess in the
    // solution vector. do so here as well:
    ierr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }


  /* ---------------------- SolverPreOnly ------------------------ */

  SolverPreOnly::SolverPreOnly(SolverControl &       cn,
                               const MPI_Comm &      mpi_communicator,
                               const AdditionalData &data)
    : SolverBase(cn, mpi_communicator)
    , additional_data(data)
  {}


  void
  SolverPreOnly::set_solver_type(KSP &ksp) const
  {
    PetscErrorCode ierr = KSPSetType(ksp, KSPPREONLY);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    // The KSPPREONLY solver of
    // PETSc never calls the convergence
    // monitor, which leads to failure
    // even when everything was ok.
    // Therefore the SolverControl status
    // is set to some nice values, which
    // guarantee a nice result at the end
    // of the solution process.
    solver_control.check(1, 0.0);

    // Using the PREONLY solver with
    // a nonzero initial guess leads
    // PETSc to produce some error messages.
    ierr = KSPSetInitialGuessNonzero(ksp, PETSC_FALSE);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }


  /* ---------------------- SparseDirectMUMPS------------------------ */

  SparseDirectMUMPS::SolverDataMUMPS::~SolverDataMUMPS()
  {
    destroy_krylov_solver(ksp);
    // the 'pc' object is owned by the 'ksp' object, and consequently
    // does not have to be destroyed explicitly here
  }


  SparseDirectMUMPS::SparseDirectMUMPS(SolverControl &       cn,
                                       const MPI_Comm &      mpi_communicator,
                                       const AdditionalData &data)
    : SolverBase(cn, mpi_communicator)
    , additional_data(data)
    , symmetric_mode(false)
  {}


  void
  SparseDirectMUMPS::set_solver_type(KSP &ksp) const
  {
    /*
     * KSPPREONLY implements a stub method that applies only the
     * preconditioner.  Its use is due to SparseDirectMUMPS being a direct
     * (rather than iterative) solver
     */
    PetscErrorCode ierr = KSPSetType(ksp, KSPPREONLY);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    /*
     * The KSPPREONLY solver of PETSc never calls the convergence monitor,
     * which leads to failure even when everything was ok. Therefore, the
     * SolverControl status is set to some nice values, which guarantee a
     * nice result at the end of the solution process.
     */
    solver_control.check(1, 0.0);

    /*
     * Using a PREONLY solver with a nonzero initial guess leads PETSc to
     * produce some error messages.
     */
    ierr = KSPSetInitialGuessNonzero(ksp, PETSC_FALSE);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }

  void
  SparseDirectMUMPS::solve(const MatrixBase &A,
                           VectorBase &      x,
                           const VectorBase &b)
  {
#  ifdef DEAL_II_PETSC_WITH_MUMPS
    /*
     * factorization matrix to be obtained from MUMPS
     */
    Mat F;

    /*
     * setting MUMPS integer control parameters ICNTL to be passed to
     * MUMPS.  Setting entry 7 of MUMPS ICNTL array (of size 40) to a value
     * of 2. This sets use of Approximate Minimum Fill (AMF)
     */
    PetscInt ival = 2, icntl = 7;
    /*
     * number of iterations to solution (should be 1) for a direct solver
     */
    PetscInt its;
    /*
     * norm of residual
     */
    PetscReal rnorm;

    /*
     * creating a solver object if this is necessary
     */
    if (solver_data == nullptr)
      {
        solver_data = std::make_unique<SolverDataMUMPS>();

        /*
         * creates the default KSP context and puts it in the location
         * solver_data->ksp
         */
        PetscErrorCode ierr = KSPCreate(mpi_communicator, &solver_data->ksp);
        AssertThrow(ierr == 0, ExcPETScError(ierr));

        /*
         * set the matrices involved. the last argument is irrelevant here,
         * since we use the solver only once anyway
         */
#    if DEAL_II_PETSC_VERSION_LT(3, 5, 0)
        ierr =
          KSPSetOperators(solver_data->ksp, A, A, DIFFERENT_NONZERO_PATTERN);
#    else
        ierr = KSPSetOperators(solver_data->ksp, A, A);
#    endif
        AssertThrow(ierr == 0, ExcPETScError(ierr));

        /*
         * setting the solver type
         */
        set_solver_type(solver_data->ksp);

        /*
         * getting the associated preconditioner context
         */
        ierr = KSPGetPC(solver_data->ksp, &solver_data->pc);
        AssertThrow(ierr == 0, ExcPETScError(ierr));

        /*
         * build PETSc PC for particular PCLU or PCCHOLESKY preconditioner
         * depending on whether the symmetric mode has been set
         */
        if (symmetric_mode)
          ierr = PCSetType(solver_data->pc, PCCHOLESKY);
        else
          ierr = PCSetType(solver_data->pc, PCLU);
        AssertThrow(ierr == 0, ExcPETScError(ierr));

        /*
         * convergence monitor function that checks with the solver_control
         * object for convergence
         */
        ierr = KSPSetConvergenceTest(solver_data->ksp,
                                     &convergence_test,
                                     reinterpret_cast<void *>(&solver_control),
                                     PETSC_NULL);
        AssertThrow(ierr == 0, ExcPETScError(ierr));

        /*
         * set the software that is to be used to perform the lu
         * factorization here we start to see differences with the base
         * class solve function
         */
#    if DEAL_II_PETSC_VERSION_LT(3, 9, 0)
        ierr = PCFactorSetMatSolverPackage(solver_data->pc, MATSOLVERMUMPS);
#    else
        ierr = PCFactorSetMatSolverType(solver_data->pc, MATSOLVERMUMPS);
#    endif
        AssertThrow(ierr == 0, ExcPETScError(ierr));

        /*
         * set up the package to call for the factorization
         */
#    if DEAL_II_PETSC_VERSION_LT(3, 9, 0)
        ierr = PCFactorSetUpMatSolverPackage(solver_data->pc);
#    else
        ierr = PCFactorSetUpMatSolverType(solver_data->pc);
#    endif
        AssertThrow(ierr == 0, ExcPETScError(ierr));

        /*
         * get the factored matrix F from the preconditioner context.  This
         * routine is valid only for LU, ILU, Cholesky, and incomplete
         * Cholesky
         */
        ierr = PCFactorGetMatrix(solver_data->pc, &F);
        AssertThrow(ierr == 0, ExcPETScError(ierr));

        /*
         * Passing the control parameters to MUMPS
         */
        ierr = MatMumpsSetIcntl(F, icntl, ival);
        AssertThrow(ierr == 0, ExcPETScError(ierr));

        /*
         * set the command line option prefix name
         */
        ierr = KSPSetOptionsPrefix(solver_data->ksp, prefix_name.c_str());
        AssertThrow(ierr == 0, ExcPETScError(ierr));

        /*
         * set the command line options provided by the user to override
         * the defaults
         */
        ierr = KSPSetFromOptions(solver_data->ksp);
        AssertThrow(ierr == 0, ExcPETScError(ierr));
      }

    /*
     * solve the linear system
     */
    PetscErrorCode ierr = KSPSolve(solver_data->ksp, b, x);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    /*
     * in case of failure throw exception
     */
    if (solver_control.last_check() != SolverControl::success)
      {
        AssertThrow(false,
                    SolverControl::NoConvergence(solver_control.last_step(),
                                                 solver_control.last_value()));
      }
    else
      {
        /*
         * obtain convergence information. obtain the number of iterations
         * and residual norm
         */
        ierr = KSPGetIterationNumber(solver_data->ksp, &its);
        AssertThrow(ierr == 0, ExcPETScError(ierr));
        ierr = KSPGetResidualNorm(solver_data->ksp, &rnorm);
        AssertThrow(ierr == 0, ExcPETScError(ierr));
      }

#  else // DEAL_II_PETSC_WITH_MUMPS
    Assert(
      false,
      ExcMessage(
        "Your PETSc installation does not include a copy of "
        "the MUMPS package necessary for this solver. You will need to configure "
        "PETSc so that it includes MUMPS, recompile it, and then re-configure "
        "and recompile deal.II as well."));

    // Cast to void to silence compiler warnings
    (void)A;
    (void)x;
    (void)b;
#  endif
  }

  PetscErrorCode
  SparseDirectMUMPS::convergence_test(KSP /*ksp*/,
                                      const PetscInt      iteration,
                                      const PetscReal     residual_norm,
                                      KSPConvergedReason *reason,
                                      void *              solver_control_x)
  {
    SolverControl &solver_control =
      *reinterpret_cast<SolverControl *>(solver_control_x);

    const SolverControl::State state =
      solver_control.check(iteration, residual_norm);

    switch (state)
      {
        case ::dealii::SolverControl::iterate:
          *reason = KSP_CONVERGED_ITERATING;
          break;

        case ::dealii::SolverControl::success:
          *reason = static_cast<KSPConvergedReason>(1);
          break;

        case ::dealii::SolverControl::failure:
          if (solver_control.last_step() > solver_control.max_steps())
            *reason = KSP_DIVERGED_ITS;
          else
            *reason = KSP_DIVERGED_DTOL;
          break;

        default:
          Assert(false, ExcNotImplemented());
      }

    return 0;
  }

  void
  SparseDirectMUMPS::set_symmetric_mode(const bool flag)
  {
    symmetric_mode = flag;
  }

} // namespace PETScWrappers

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PETSC
