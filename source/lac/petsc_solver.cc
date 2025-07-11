// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2004 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/logstream.h>

#include <deal.II/lac/petsc_solver.h>

#ifdef DEAL_II_WITH_PETSC

#  include <deal.II/lac/exceptions.h>
#  include <deal.II/lac/petsc_compatibility.h>
#  include <deal.II/lac/petsc_matrix_base.h>
#  include <deal.II/lac/petsc_precondition.h>
#  include <deal.II/lac/petsc_vector_base.h>

// Shorthand notation for PETSc error codes.
#  define AssertPETSc(code)                          \
    do                                               \
      {                                              \
        PetscErrorCode ierr = (code);                \
        AssertThrow(ierr == 0, ExcPETScError(ierr)); \
      }                                              \
    while (false)

DEAL_II_NAMESPACE_OPEN

namespace PETScWrappers
{
  SolverBase::SolverBase()
    : ksp(nullptr)
    , solver_control(nullptr)
  {}


  SolverBase::SolverBase(SolverControl &cn)
    : ksp(nullptr)
    , solver_control(&cn)
  {}



  void
  SolverBase::set_solver_type(KSP &) const
  {}



  SolverBase::~SolverBase()
  {
    AssertPETSc(KSPDestroy(&ksp));
  }



  KSP
  SolverBase::petsc_ksp()
  {
    return ksp;
  }



  SolverBase::operator KSP() const
  {
    return ksp;
  }



  void
  SolverBase::solve(const MatrixBase       &A,
                    VectorBase             &x,
                    const VectorBase       &b,
                    const PreconditionBase &preconditioner)
  {
    // first create a solver object if this
    // is necessary
    if (ksp == nullptr)
      {
        initialize_ksp_with_comm(A.get_mpi_communicator());

        // let derived classes set the solver
        // type, and the preconditioning
        // object set the type of
        // preconditioner
        set_solver_type(ksp);

        AssertPETSc(KSPSetPC(ksp, preconditioner.get_pc()));

        /*
         * by default we set up the preconditioner only once.
         * this can be overridden by command line.
         */
        AssertPETSc(KSPSetReusePreconditioner(ksp, PETSC_TRUE));
      }

    // setting the preconditioner overwrites the used matrices.
    // hence, we need to set the matrices after the preconditioner.
    Mat B;
    AssertPETSc(KSPGetOperators(ksp, nullptr, &B));
    AssertPETSc(KSPSetOperators(ksp, A, B));

    // set the command line option prefix name
    AssertPETSc(KSPSetOptionsPrefix(ksp, prefix_name.c_str()));

    // set the command line options provided
    // by the user to override the defaults
    AssertPETSc(KSPSetFromOptions(ksp));

    // then do the real work: set up solver
    // internal data and solve the
    // system.
    AssertPETSc(KSPSetUp(ksp));

    AssertPETSc(KSPSolve(ksp, b, x));

    // in case of failure: throw
    // exception
    if (solver_control &&
        solver_control->last_check() != SolverControl::success)
      AssertThrow(false,
                  SolverControl::NoConvergence(solver_control->last_step(),
                                               solver_control->last_value()));
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
    AssertPETSc(KSPDestroy(&ksp));
  }


  SolverControl &
  SolverBase::control() const
  {
    AssertThrow(
      solver_control,
      ExcMessage(
        "You need to create the solver with a SolverControl object if you want to call the function that returns it."));
    return *solver_control;
  }


  PetscErrorCode
  SolverBase::convergence_test(KSP /*ksp*/,
                               const PetscInt      iteration,
                               const PetscReal     residual_norm,
                               KSPConvergedReason *reason,
                               void               *solver_control_x)
  {
    PetscFunctionBeginUser;
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
          *reason = KSP_CONVERGED_RTOL;
          break;

        case ::dealii::SolverControl::failure:
          if (solver_control.last_step() > solver_control.max_steps())
            *reason = KSP_DIVERGED_ITS;
          else
            *reason = KSP_DIVERGED_DTOL;
          break;

        default:
          DEAL_II_NOT_IMPLEMENTED();
      }

    PetscFunctionReturn(PETSC_SUCCESS);
  }



  void
  SolverBase::initialize_ksp_with_comm(const MPI_Comm comm)
  {
    // Create the PETSc KSP object
    AssertPETSc(KSPCreate(comm, &ksp));

    // then a convergence monitor
    // function that simply
    // checks with the solver_control
    // object we have in this object for
    // convergence
    perhaps_set_convergence_test();
  }



  void
  SolverBase::perhaps_set_convergence_test() const
  {
    if (ksp && solver_control)
      AssertPETSc(
        KSPSetConvergenceTest(ksp, &convergence_test, solver_control, nullptr));
  }


  void
  SolverBase::initialize(const PreconditionBase &preconditioner)
  {
    initialize_ksp_with_comm(preconditioner.get_mpi_communicator());

    // let derived classes set the solver
    // type, and the preconditioning
    // object set the type of
    // preconditioner
    set_solver_type(ksp);

    // set the command line options provided
    // by the user to override the defaults
    AssertPETSc(KSPSetFromOptions(ksp));
  }



  /* ---------------------- SolverRichardson ------------------------ */

  SolverRichardson::AdditionalData::AdditionalData(const double omega)
    : omega(omega)
  {}



  SolverRichardson::SolverRichardson(SolverControl        &cn,
                                     const AdditionalData &data)
    : SolverBase(cn)
    , additional_data(data)
  {}



  SolverRichardson::SolverRichardson(SolverControl &cn,
                                     const MPI_Comm,
                                     const AdditionalData &data)
    : SolverRichardson(cn, data)
  {}


  void
  SolverRichardson::set_solver_type(KSP &ksp) const
  {
    AssertPETSc(KSPSetType(ksp, KSPRICHARDSON));

    // set the damping factor from the data
    AssertPETSc(KSPRichardsonSetScale(ksp, additional_data.omega));

    // in the deal.II solvers, we always
    // honor the initial guess in the
    // solution vector. do so here as well:
    AssertPETSc(KSPSetInitialGuessNonzero(ksp, PETSC_TRUE));

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
    AssertPETSc(KSPSetTolerances(ksp,
                                 PETSC_DEFAULT,
                                 this->solver_control->tolerance(),
                                 PETSC_DEFAULT,
                                 this->solver_control->max_steps() + 1));
  }


  /* ---------------------- SolverChebychev ------------------------ */

  SolverChebychev::SolverChebychev(SolverControl        &cn,
                                   const AdditionalData &data)
    : SolverBase(cn)
    , additional_data(data)
  {}


  SolverChebychev::SolverChebychev(SolverControl &cn,
                                   const MPI_Comm,
                                   const AdditionalData &data)
    : SolverChebychev(cn, data)
  {}


  void
  SolverChebychev::set_solver_type(KSP &ksp) const
  {
    AssertPETSc(KSPSetType(ksp, KSPCHEBYSHEV));

    // in the deal.II solvers, we always
    // honor the initial guess in the
    // solution vector. do so here as well:
    AssertPETSc(KSPSetInitialGuessNonzero(ksp, PETSC_TRUE));
  }


  /* ---------------------- SolverCG ------------------------ */

  SolverCG::SolverCG(SolverControl &cn, const AdditionalData &data)
    : SolverBase(cn)
    , additional_data(data)
  {}


  SolverCG::SolverCG(SolverControl &cn,
                     const MPI_Comm,
                     const AdditionalData &data)
    : SolverCG(cn, data)
  {}


  void
  SolverCG::set_solver_type(KSP &ksp) const
  {
    AssertPETSc(KSPSetType(ksp, KSPCG));

    // in the deal.II solvers, we always
    // honor the initial guess in the
    // solution vector. do so here as well:
    AssertPETSc(KSPSetInitialGuessNonzero(ksp, PETSC_TRUE));
  }


  /* ---------------------- SolverBiCG ------------------------ */

  SolverBiCG::SolverBiCG(SolverControl &cn, const AdditionalData &data)
    : SolverBase(cn)
    , additional_data(data)
  {}


  SolverBiCG::SolverBiCG(SolverControl &cn,
                         const MPI_Comm,
                         const AdditionalData &data)
    : SolverBiCG(cn, data)
  {}


  void
  SolverBiCG::set_solver_type(KSP &ksp) const
  {
    AssertPETSc(KSPSetType(ksp, KSPBICG));

    // in the deal.II solvers, we always
    // honor the initial guess in the
    // solution vector. do so here as well:
    AssertPETSc(KSPSetInitialGuessNonzero(ksp, PETSC_TRUE));
  }


  /* ---------------------- SolverGMRES ------------------------ */

  SolverGMRES::AdditionalData::AdditionalData(
    const unsigned int restart_parameter,
    const bool         right_preconditioning)
    : restart_parameter(restart_parameter)
    , right_preconditioning(right_preconditioning)
  {}



  SolverGMRES::SolverGMRES(SolverControl &cn, const AdditionalData &data)
    : SolverBase(cn)
    , additional_data(data)
  {}


  SolverGMRES::SolverGMRES(SolverControl &cn,
                           const MPI_Comm,
                           const AdditionalData &data)
    : SolverGMRES(cn, data)
  {}


  void
  SolverGMRES::set_solver_type(KSP &ksp) const
  {
    AssertPETSc(KSPSetType(ksp, KSPGMRES));

    AssertPETSc(KSPGMRESSetRestart(ksp, additional_data.restart_parameter));

    // Set preconditioning side to right
    if (additional_data.right_preconditioning)
      {
        AssertPETSc(KSPSetPCSide(ksp, PC_RIGHT));
      }

    // in the deal.II solvers, we always
    // honor the initial guess in the
    // solution vector. do so here as well:
    AssertPETSc(KSPSetInitialGuessNonzero(ksp, PETSC_TRUE));
  }


  /* ---------------------- SolverBicgstab ------------------------ */

  SolverBicgstab::SolverBicgstab(SolverControl &cn, const AdditionalData &data)
    : SolverBase(cn)
    , additional_data(data)
  {}


  SolverBicgstab::SolverBicgstab(SolverControl &cn,
                                 const MPI_Comm,
                                 const AdditionalData &data)
    : SolverBicgstab(cn, data)
  {}


  void
  SolverBicgstab::set_solver_type(KSP &ksp) const
  {
    AssertPETSc(KSPSetType(ksp, KSPBCGS));

    // in the deal.II solvers, we always
    // honor the initial guess in the
    // solution vector. do so here as well:
    AssertPETSc(KSPSetInitialGuessNonzero(ksp, PETSC_TRUE));
  }


  /* ---------------------- SolverCGS ------------------------ */

  SolverCGS::SolverCGS(SolverControl &cn, const AdditionalData &data)
    : SolverBase(cn)
    , additional_data(data)
  {}


  SolverCGS::SolverCGS(SolverControl &cn,
                       const MPI_Comm,
                       const AdditionalData &data)
    : SolverCGS(cn, data)
  {}


  void
  SolverCGS::set_solver_type(KSP &ksp) const
  {
    AssertPETSc(KSPSetType(ksp, KSPCGS));

    // in the deal.II solvers, we always
    // honor the initial guess in the
    // solution vector. do so here as well:
    AssertPETSc(KSPSetInitialGuessNonzero(ksp, PETSC_TRUE));
  }


  /* ---------------------- SolverTFQMR ------------------------ */

  SolverTFQMR::SolverTFQMR(SolverControl &cn, const AdditionalData &data)
    : SolverBase(cn)
    , additional_data(data)
  {}


  SolverTFQMR::SolverTFQMR(SolverControl &cn,
                           const MPI_Comm,
                           const AdditionalData &data)
    : SolverTFQMR(cn, data)
  {}


  void
  SolverTFQMR::set_solver_type(KSP &ksp) const
  {
    AssertPETSc(KSPSetType(ksp, KSPTFQMR));

    // in the deal.II solvers, we always
    // honor the initial guess in the
    // solution vector. do so here as well:
    AssertPETSc(KSPSetInitialGuessNonzero(ksp, PETSC_TRUE));
  }


  /* ---------------------- SolverTCQMR ------------------------ */

  SolverTCQMR::SolverTCQMR(SolverControl &cn, const AdditionalData &data)
    : SolverBase(cn)
    , additional_data(data)
  {}


  SolverTCQMR::SolverTCQMR(SolverControl &cn,
                           const MPI_Comm,
                           const AdditionalData &data)
    : SolverTCQMR(cn, data)
  {}


  void
  SolverTCQMR::set_solver_type(KSP &ksp) const
  {
    AssertPETSc(KSPSetType(ksp, KSPTCQMR));

    // in the deal.II solvers, we always
    // honor the initial guess in the
    // solution vector. do so here as well:
    AssertPETSc(KSPSetInitialGuessNonzero(ksp, PETSC_TRUE));
  }


  /* ---------------------- SolverCR ------------------------ */

  SolverCR::SolverCR(SolverControl &cn, const AdditionalData &data)
    : SolverBase(cn)
    , additional_data(data)
  {}


  SolverCR::SolverCR(SolverControl &cn,
                     const MPI_Comm,
                     const AdditionalData &data)
    : SolverCR(cn, data)
  {}


  void
  SolverCR::set_solver_type(KSP &ksp) const
  {
    AssertPETSc(KSPSetType(ksp, KSPCR));

    // in the deal.II solvers, we always
    // honor the initial guess in the
    // solution vector. do so here as well:
    AssertPETSc(KSPSetInitialGuessNonzero(ksp, PETSC_TRUE));
  }


  /* ---------------------- SolverLSQR ------------------------ */

  SolverLSQR::SolverLSQR(SolverControl &cn, const AdditionalData &data)
    : SolverBase(cn)
    , additional_data(data)
  {}



  SolverLSQR::SolverLSQR(SolverControl &cn,
                         const MPI_Comm,
                         const AdditionalData &data)
    : SolverLSQR(cn, data)
  {}



  void
  SolverLSQR::set_solver_type(KSP &ksp) const
  {
    AssertPETSc(KSPSetType(ksp, KSPLSQR));

    // in the deal.II solvers, we always
    // honor the initial guess in the
    // solution vector. do so here as well:
    AssertPETSc(KSPSetInitialGuessNonzero(ksp, PETSC_TRUE));

    // The KSPLSQR implementation overwrites the user-defined
    // convergence test at creation (i.e. KSPSetType) time.
    // This is probably a bad design decision in PETSc.
    // Anyway, here we make sure we use our own convergence
    // test.
    perhaps_set_convergence_test();
  }


  /* ---------------------- SolverPreOnly ------------------------ */

  SolverPreOnly::SolverPreOnly(SolverControl &cn, const AdditionalData &data)
    : SolverBase(cn)
    , additional_data(data)
  {}


  SolverPreOnly::SolverPreOnly(SolverControl &cn,
                               const MPI_Comm,
                               const AdditionalData &data)
    : SolverPreOnly(cn, data)
  {}


  void
  SolverPreOnly::set_solver_type(KSP &ksp) const
  {
    AssertPETSc(KSPSetType(ksp, KSPPREONLY));

    // The KSPPREONLY solver of
    // PETSc never calls the convergence
    // monitor, which leads to failure
    // even when everything was ok.
    // Therefore the SolverControl status
    // is set to some nice values, which
    // guarantee a nice result at the end
    // of the solution process.
    solver_control->check(1, 0.0);

    // Using the PREONLY solver with
    // a nonzero initial guess leads
    // PETSc to produce some error messages.
    AssertPETSc(KSPSetInitialGuessNonzero(ksp, PETSC_FALSE));
  }


  /* ---------------------- SparseDirectMUMPS------------------------ */

  SparseDirectMUMPS::SparseDirectMUMPS(SolverControl        &cn,
                                       const AdditionalData &data)
    : SolverBase(cn)
    , additional_data(data)
    , symmetric_mode(false)
  {}



  SparseDirectMUMPS::SparseDirectMUMPS(SolverControl &cn,
                                       const MPI_Comm,
                                       const AdditionalData &data)
    : SparseDirectMUMPS(cn, data)
  {}



  void
  SparseDirectMUMPS::set_solver_type(KSP &ksp) const
  {
    /*
     * KSPPREONLY implements a stub method that applies only the
     * preconditioner.  Its use is due to SparseDirectMUMPS being a direct
     * (rather than iterative) solver
     */
    AssertPETSc(KSPSetType(ksp, KSPPREONLY));

    /*
     * The KSPPREONLY solver of PETSc never calls the convergence monitor,
     * which leads to failure even when everything was ok. Therefore, the
     * SolverControl status is set to some nice values, which guarantee a
     * nice result at the end of the solution process.
     */
    solver_control->check(1, 0.0);

    /*
     * Using a PREONLY solver with a nonzero initial guess leads PETSc to
     * produce some error messages.
     */
    AssertPETSc(KSPSetInitialGuessNonzero(ksp, PETSC_FALSE));
  }

  void
  SparseDirectMUMPS::solve(const MatrixBase &A,
                           VectorBase       &x,
                           const VectorBase &b)
  {
#  ifdef DEAL_II_PETSC_WITH_MUMPS
    /*
     * creating a solver object if this is necessary
     */
    if (ksp == nullptr)
      {
        initialize_ksp_with_comm(A.get_mpi_communicator());

        /*
         * setting the solver type
         */
        set_solver_type(ksp);

        /*
         * set the matrices involved. the last argument is irrelevant here,
         * since we use the solver only once anyway
         */
        AssertPETSc(KSPSetOperators(ksp, A, A));

        /*
         * getting the associated preconditioner context
         */
        PC pc;
        AssertPETSc(KSPGetPC(ksp, &pc));

        /*
         * build PETSc PC for particular PCLU or PCCHOLESKY preconditioner
         * depending on whether the symmetric mode has been set
         */
        if (symmetric_mode)
          AssertPETSc(PCSetType(pc, PCCHOLESKY));
        else
          AssertPETSc(PCSetType(pc, PCLU));

          /*
           * set the software that is to be used to perform the lu
           * factorization here we start to see differences with the base
           * class solve function
           */
#    if DEAL_II_PETSC_VERSION_LT(3, 9, 0)
        AssertPETSc(PCFactorSetMatSolverPackage(pc, MATSOLVERMUMPS));
#    else
        AssertPETSc(PCFactorSetMatSolverType(pc, MATSOLVERMUMPS));
#    endif

        /*
         * set up the package to call for the factorization
         */
#    if DEAL_II_PETSC_VERSION_LT(3, 9, 0)
        AssertPETSc(PCFactorSetUpMatSolverPackage(pc));
#    else
        AssertPETSc(PCFactorSetUpMatSolverType(pc));
#    endif

        /*
         * get the factored matrix F from the preconditioner context.
         */
        Mat F;
        AssertPETSc(PCFactorGetMatrix(pc, &F));

        /*
         * pass control parameters to MUMPS.
         * Setting entry 7 of MUMPS ICNTL array to a value
         * of 2. This sets use of Approximate Minimum Fill (AMF)
         */
        AssertPETSc(MatMumpsSetIcntl(F, 7, 2));

        /*
         * by default we set up the preconditioner only once.
         * this can be overridden by command line.
         */
        AssertPETSc(KSPSetReusePreconditioner(ksp, PETSC_TRUE));
      }

    /*
     * set the matrices involved. the last argument is irrelevant here,
     * since we use the solver only once anyway
     */
    AssertPETSc(KSPSetOperators(ksp, A, A));

    /*
     * set the command line option prefix name
     */
    AssertPETSc(KSPSetOptionsPrefix(ksp, prefix_name.c_str()));

    /*
     * set the command line options provided by the user to override
     * the defaults
     */
    AssertPETSc(KSPSetFromOptions(ksp));

    /*
     * solve the linear system
     */
    AssertPETSc(KSPSolve(ksp, b, x));

    /*
     * in case of failure throw exception
     */
    if (solver_control &&
        solver_control->last_check() != SolverControl::success)
      {
        AssertThrow(false,
                    SolverControl::NoConvergence(solver_control->last_step(),
                                                 solver_control->last_value()));
      }

#  else // DEAL_II_PETSC_WITH_MUMPS
    Assert(
      false,
      ExcMessage(
        "Your PETSc installation does not include a copy of "
        "the MUMPS package necessary for this solver. You will need to configure "
        "PETSc so that it includes MUMPS, recompile it, and then re-configure "
        "and recompile deal.II as well."));
    (void)A;
    (void)x;
    (void)b;
#  endif
  }



  void
  SparseDirectMUMPS::set_symmetric_mode(const bool matrix_is_symmetric)
  {
    symmetric_mode = matrix_is_symmetric;
  }

} // namespace PETScWrappers

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PETSC
