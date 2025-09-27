// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2009 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/lac/slepc_solver.h>

#ifdef DEAL_II_WITH_SLEPC

#  include <deal.II/lac/petsc_matrix_base.h>
#  include <deal.II/lac/petsc_vector_base.h>
#  include <deal.II/lac/slepc_spectral_transformation.h>

#  include <petscversion.h>

#  include <slepcversion.h>

#  include <cmath>
#  include <vector>

DEAL_II_NAMESPACE_OPEN

namespace SLEPcWrappers
{
  SolverBase::SolverBase(SolverControl &cn, const MPI_Comm mpi_communicator)
    : solver_control(cn)
    , mpi_communicator(mpi_communicator)
    , reason(EPS_CONVERGED_ITERATING)
  {
    // create eigensolver context
    PetscErrorCode ierr = EPSCreate(mpi_communicator, &eps);
    AssertThrow(ierr == 0, ExcSLEPcError(ierr));

    // hand over the absolute tolerance and the maximum number of
    // iteration steps to the SLEPc convergence criterion.
    ierr = EPSSetTolerances(eps,
                            this->solver_control.tolerance(),
                            this->solver_control.max_steps());
    AssertThrow(ierr == 0, ExcSLEPcError(ierr));

    // default values:
    set_which_eigenpairs(EPS_LARGEST_MAGNITUDE);
    set_problem_type(EPS_GNHEP);

    // TODO:
    // By default, EPS initializes the starting vector or the initial subspace
    // randomly.
  }



  SolverBase::~SolverBase()
  {
    if (eps != nullptr)
      {
        // Destroy the solver object.
        const PetscErrorCode ierr = EPSDestroy(&eps);
        AssertNothrow(ierr == 0, ExcSLEPcError(ierr));
      }
  }



  void
  SolverBase::set_matrices(const PETScWrappers::MatrixBase &A)
  {
    // standard eigenspectrum problem
    const PetscErrorCode ierr = EPSSetOperators(eps, A, nullptr);
    AssertThrow(ierr == 0, ExcSLEPcError(ierr));
  }



  void
  SolverBase::set_matrices(const PETScWrappers::MatrixBase &A,
                           const PETScWrappers::MatrixBase &B)
  {
    // generalized eigenspectrum problem
    const PetscErrorCode ierr = EPSSetOperators(eps, A, B);
    AssertThrow(ierr == 0, ExcSLEPcError(ierr));
  }



  void
  SolverBase::set_transformation(
    SLEPcWrappers::TransformationBase &transformation)
  {
    // set transformation type if any
    // STSetShift is called inside
    PetscErrorCode ierr = EPSSetST(eps, transformation.st);
    AssertThrow(ierr == 0, SolverBase::ExcSLEPcError(ierr));

#  if DEAL_II_SLEPC_VERSION_GTE(3, 8, 0)
    // see
    // https://lists.mcs.anl.gov/mailman/htdig/petsc-users/2017-October/033649.html
    // From 3.8.0 onward, SLEPc insists that when looking for smallest
    // eigenvalues with shift-and-invert, users should (a) set target,
    // (b) use EPS_TARGET_MAGNITUDE. The former, however, needs to be
    // applied to the 'eps' object and not the spectral transformation.
    if (SLEPcWrappers::TransformationShiftInvert *sinv =
          dynamic_cast<SLEPcWrappers::TransformationShiftInvert *>(
            &transformation))
      {
        ierr = EPSSetTarget(eps, sinv->additional_data.shift_parameter);
        AssertThrow(ierr == 0, SolverBase::ExcSLEPcError(ierr));
      }
#  endif
  }



  void
  SolverBase::set_target_eigenvalue(const PetscScalar &this_target)
  {
    // set target eigenvalues to solve for
    // in all transformation except STSHIFT there is a direct connection between
    // the target and the shift, read more on p41 of SLEPc manual.
    const PetscErrorCode ierr = EPSSetTarget(eps, this_target);
    AssertThrow(ierr == 0, ExcSLEPcError(ierr));
  }



  void
  SolverBase::set_which_eigenpairs(const EPSWhich eps_which)
  {
    // set which portion of the eigenspectrum to solve for
    const PetscErrorCode ierr = EPSSetWhichEigenpairs(eps, eps_which);
    AssertThrow(ierr == 0, ExcSLEPcError(ierr));
  }



  void
  SolverBase::set_problem_type(const EPSProblemType eps_problem)
  {
    const PetscErrorCode ierr = EPSSetProblemType(eps, eps_problem);
    AssertThrow(ierr == 0, ExcSLEPcError(ierr));
  }



  void
  SolverBase::solve(const unsigned int n_eigenpairs, unsigned int *n_converged)
  {
    // set number of eigenvectors to compute
    PetscErrorCode ierr =
      EPSSetDimensions(eps, n_eigenpairs, PETSC_DECIDE, PETSC_DECIDE);
    AssertThrow(ierr == 0, ExcSLEPcError(ierr));

    // set the solve options to the eigenvalue problem solver context
    ierr = EPSSetFromOptions(eps);
    AssertThrow(ierr == 0, ExcSLEPcError(ierr));

    // TODO breaks step-36
    // force Krylov solver to use true residual instead of an estimate.
    // EPSSetTrueResidual(solver_data->eps, PETSC_TRUE);
    // AssertThrow (ierr == 0, ExcSLEPcError(ierr));

    // Set convergence test to be absolute
    ierr = EPSSetConvergenceTest(eps, EPS_CONV_ABS);
    AssertThrow(ierr == 0, ExcSLEPcError(ierr));

    // TODO Set the convergence test function
    // ierr = EPSSetConvergenceTestFunction (solver_data->eps,
    // &convergence_test,
    //              reinterpret_cast<void *>(&solver_control));
    // AssertThrow (ierr == 0, ExcSLEPcError(ierr));

    // solve the eigensystem
    ierr = EPSSolve(eps);
    AssertThrow(ierr == 0, ExcSLEPcError(ierr));

    // Get number of converged eigenstates. We need to go around with a
    // temporary variable once because the function wants to have a
    // PetscInt as second argument whereas the `n_converged` argument
    // to this function is just an unsigned int.
    {
      PetscInt petsc_n_converged = *n_converged;
      ierr                       = EPSGetConverged(eps, &petsc_n_converged);
      AssertThrow(ierr == 0, ExcSLEPcError(ierr));
      *n_converged = petsc_n_converged;
    }

    PetscInt  n_iterations  = 0;
    PetscReal residual_norm = 0;

    // @todo Investigate elaborating on some of this to act on the
    // complete eigenspectrum
    {
      // get the number of solver iterations
      ierr = EPSGetIterationNumber(eps, &n_iterations);
      AssertThrow(ierr == 0, ExcSLEPcError(ierr));

      // get the maximum of residual norm among converged eigenvectors.
      for (unsigned int i = 0; i < *n_converged; ++i)
        {
          double residual_norm_i = 0.0;
          // EPSComputeError (or, in older versions of SLEPc,
          // EPSComputeResidualNorm) uses an L2-norm and is not consistent
          // with the stopping criterion used during the solution process (see
          // the SLEPC manual, section 2.5). However, the norm that gives error
          // bounds (Saad, 1992, ch3) is (for Hermitian problems)
          // | \lambda - \widehat\lambda | <= ||r||_2
          //
          // Similarly, EPSComputeRelativeError may not be consistent with the
          // stopping criterion used in the solution process.
          //
          // EPSGetErrorEstimate is (according to the SLEPc manual) consistent
          // with the residual norm used during the solution process. However,
          // it is not guaranteed to be derived from the residual even when
          // EPSSetTrueResidual is set: see the discussion in the thread
          //
          // https://lists.mcs.anl.gov/pipermail/petsc-users/2014-November/023509.html
          //
          // for more information.
#  if DEAL_II_SLEPC_VERSION_GTE(3, 6, 0)
          ierr = EPSComputeError(eps, i, EPS_ERROR_ABSOLUTE, &residual_norm_i);
#  else
          ierr = EPSComputeResidualNorm(eps, i, &residual_norm_i);
#  endif

          AssertThrow(ierr == 0, ExcSLEPcError(ierr));
          residual_norm = std::max(residual_norm, residual_norm_i);
        }

      // check the solver state
      const SolverControl::State state =
        solver_control.check(n_iterations, residual_norm);

      // get the solver state according to SLEPc
      get_solver_state(state);

      // as SLEPc uses different stopping criteria, we have to omit this step.
      // This can be checked only in conjunction with EPSGetErrorEstimate.
      // and in case of failure: throw exception
      // if (solver_control.last_check () != SolverControl::success)
      //   AssertThrow(false, SolverControl::NoConvergence
      //   (solver_control.last_step(),
      //                                                    solver_control.last_value()));
    }
  }



  void
  SolverBase::get_eigenpair(const unsigned int         index,
                            PetscScalar               &eigenvalues,
                            PETScWrappers::VectorBase &eigenvectors)
  {
    // get converged eigenpair
    const PetscErrorCode ierr =
      EPSGetEigenpair(eps, index, &eigenvalues, nullptr, eigenvectors, nullptr);
    AssertThrow(ierr == 0, ExcSLEPcError(ierr));
  }



  void
  SolverBase::get_eigenpair(const unsigned int         index,
                            double                    &real_eigenvalues,
                            double                    &imag_eigenvalues,
                            PETScWrappers::VectorBase &real_eigenvectors,
                            PETScWrappers::VectorBase &imag_eigenvectors)
  {
#  ifndef PETSC_USE_COMPLEX
    // get converged eigenpair
    const PetscErrorCode ierr = EPSGetEigenpair(eps,
                                                index,
                                                &real_eigenvalues,
                                                &imag_eigenvalues,
                                                real_eigenvectors,
                                                imag_eigenvectors);
    AssertThrow(ierr == 0, ExcSLEPcError(ierr));
#  else
    Assert(
      (false),
      ExcMessage(
        "Your PETSc/SLEPc installation was configured with scalar-type complex "
        "but this function is not defined for complex types."));

    // Cast to void to silence compiler warnings
    (void)index;
    (void)real_eigenvalues;
    (void)imag_eigenvalues;
    (void)real_eigenvectors;
    (void)imag_eigenvectors;
#  endif
  }



  void
  SolverBase::get_solver_state(const SolverControl::State state)
  {
    switch (state)
      {
        case ::dealii::SolverControl::iterate:
          reason = EPS_CONVERGED_ITERATING;
          break;

        case ::dealii::SolverControl::success:
          reason = static_cast<EPSConvergedReason>(1);
          break;

        case ::dealii::SolverControl::failure:
          if (solver_control.last_step() > solver_control.max_steps())
            reason = EPS_DIVERGED_ITS;
          else
            reason = EPS_DIVERGED_BREAKDOWN;
          break;

        default:
          DEAL_II_NOT_IMPLEMENTED();
      }
  }



  /* ---------------------- SolverControls ----------------------- */
  SolverControl &
  SolverBase::control() const
  {
    return solver_control;
  }



  int
  SolverBase::convergence_test(
    EPS /*eps             */,
    PetscScalar /*real_eigenvalue */,
    PetscScalar /*imag_eigenvalue */,
    PetscReal /*residual norm associated to the eigenpair   */,
    PetscReal * /*(output) computed error estimate */,
    void * /*solver_control_x*/)
  {
    // If the error estimate returned by the convergence test function is less
    // than the tolerance, then the eigenvalue is accepted as converged.
    // This function is undefined (future reference only).

    // return without failure.
    return 0;
  }



  /* ---------------------- SolverKrylovSchur ------------------------ */
  SolverKrylovSchur::SolverKrylovSchur(SolverControl        &cn,
                                       const MPI_Comm        mpi_communicator,
                                       const AdditionalData &data)
    : SolverBase(cn, mpi_communicator)
    , additional_data(data)
  {
    const PetscErrorCode ierr =
      EPSSetType(eps, const_cast<char *>(EPSKRYLOVSCHUR));
    AssertThrow(ierr == 0, ExcSLEPcError(ierr));
  }



  /* ---------------------- SolverArnoldi ------------------------ */
  SolverArnoldi::AdditionalData::AdditionalData(
    const bool delayed_reorthogonalization)
    : delayed_reorthogonalization(delayed_reorthogonalization)
  {}



  SolverArnoldi::SolverArnoldi(SolverControl        &cn,
                               const MPI_Comm        mpi_communicator,
                               const AdditionalData &data)
    : SolverBase(cn, mpi_communicator)
    , additional_data(data)
  {
    PetscErrorCode ierr = EPSSetType(eps, const_cast<char *>(EPSARNOLDI));
    AssertThrow(ierr == 0, ExcSLEPcError(ierr));

    // if requested, set delayed reorthogonalization in the Arnoldi
    // iteration.
    if (additional_data.delayed_reorthogonalization)
      {
        ierr = EPSArnoldiSetDelayed(eps, PETSC_TRUE);
        AssertThrow(ierr == 0, ExcSLEPcError(ierr));
      }
  }



  /* ---------------------- Lanczos ------------------------ */
  SolverLanczos::AdditionalData::AdditionalData(const EPSLanczosReorthogType r)
    : reorthog(r)
  {}



  SolverLanczos::SolverLanczos(SolverControl        &cn,
                               const MPI_Comm        mpi_communicator,
                               const AdditionalData &data)
    : SolverBase(cn, mpi_communicator)
    , additional_data(data)
  {
    PetscErrorCode ierr = EPSSetType(eps, const_cast<char *>(EPSLANCZOS));
    AssertThrow(ierr == 0, ExcSLEPcError(ierr));

    ierr = EPSLanczosSetReorthog(eps, additional_data.reorthog);
    AssertThrow(ierr == 0, ExcSLEPcError(ierr));
  }



  /* ----------------------- Power ------------------------- */
  SolverPower::SolverPower(SolverControl        &cn,
                           const MPI_Comm        mpi_communicator,
                           const AdditionalData &data)
    : SolverBase(cn, mpi_communicator)
    , additional_data(data)
  {
    PetscErrorCode ierr = EPSSetType(eps, const_cast<char *>(EPSPOWER));
    AssertThrow(ierr == 0, ExcSLEPcError(ierr));
  }



  /* ---------------- Generalized Davidson ----------------- */
  SolverGeneralizedDavidson::AdditionalData::AdditionalData(
    bool double_expansion)
    : double_expansion(double_expansion)
  {}



  SolverGeneralizedDavidson::SolverGeneralizedDavidson(
    SolverControl        &cn,
    const MPI_Comm        mpi_communicator,
    const AdditionalData &data)
    : SolverBase(cn, mpi_communicator)
    , additional_data(data)
  {
    PetscErrorCode ierr = EPSSetType(eps, const_cast<char *>(EPSGD));
    AssertThrow(ierr == 0, ExcSLEPcError(ierr));

    if (additional_data.double_expansion)
      {
        ierr = EPSGDSetDoubleExpansion(eps, PETSC_TRUE);
        AssertThrow(ierr == 0, ExcSLEPcError(ierr));
      }
  }



  /* ------------------ Jacobi Davidson -------------------- */
  SolverJacobiDavidson::SolverJacobiDavidson(SolverControl &cn,
                                             const MPI_Comm mpi_communicator,
                                             const AdditionalData &data)
    : SolverBase(cn, mpi_communicator)
    , additional_data(data)
  {
    const PetscErrorCode ierr = EPSSetType(eps, const_cast<char *>(EPSJD));
    AssertThrow(ierr == 0, ExcSLEPcError(ierr));
  }



  /* ---------------------- LAPACK ------------------------- */
  SolverLAPACK::SolverLAPACK(SolverControl        &cn,
                             const MPI_Comm        mpi_communicator,
                             const AdditionalData &data)
    : SolverBase(cn, mpi_communicator)
    , additional_data(data)
  {
    // 'Tis overwhelmingly likely that PETSc/SLEPc *always* has
    // BLAS/LAPACK, but let's be defensive.
#  if PETSC_HAVE_BLASLAPACK
    const PetscErrorCode ierr = EPSSetType(eps, const_cast<char *>(EPSLAPACK));
    AssertThrow(ierr == 0, ExcSLEPcError(ierr));
#  else
    Assert(
      (false),
      ExcMessage(
        "Your PETSc/SLEPc installation was not configured with BLAS/LAPACK "
        "but this is needed to use the LAPACK solver."));
#  endif
  }
} // namespace SLEPcWrappers

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_SLEPC
