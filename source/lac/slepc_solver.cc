// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2014 by the deal.II authors
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

#include <deal.II/lac/slepc_solver.h>

#ifdef DEAL_II_WITH_SLEPC

#  include <deal.II/lac/petsc_matrix_base.h>
#  include <deal.II/lac/petsc_vector_base.h>
#  include <deal.II/lac/petsc_vector.h>
#  include <deal.II/lac/slepc_spectral_transformation.h>

#  include <cmath>
#  include <vector>

#  include <petscversion.h>
#  include <slepcversion.h>

DEAL_II_NAMESPACE_OPEN

namespace SLEPcWrappers
{

  SolverBase::SolverData::~SolverData ()
  {
    // Destroy the solver object.
#if DEAL_II_PETSC_VERSION_LT(3,2,0)
    int ierr = EPSDestroy (eps);
#else
    int ierr = EPSDestroy (&eps);
#endif
    AssertThrow (ierr == 0, ExcSLEPcError(ierr));
  }

  SolverBase::SolverBase (SolverControl  &cn,
                          const MPI_Comm &mpi_communicator)
    :
    solver_control (cn),
    mpi_communicator (mpi_communicator),
    set_which (EPS_LARGEST_MAGNITUDE),
    set_problem (EPS_GNHEP),
    opA (NULL),
    opB (NULL),
    transformation (NULL)
  {}

  SolverBase::~SolverBase ()
  {}

  void
  SolverBase::set_matrices (const PETScWrappers::MatrixBase &A)
  {
    // standard eigenspectrum problem
    opA = &A;
    opB = NULL;
  }

  void
  SolverBase::set_matrices (const PETScWrappers::MatrixBase &A,
                            const PETScWrappers::MatrixBase &B)
  {
    // generalized eigenspectrum problem
    opA = &A;
    opB = &B;
  }

  void
  SolverBase::set_initial_vector (const PETScWrappers::VectorBase &this_initial_vector)
  {
    initial_space.resize(1);
    initial_space[0] = this_initial_vector;
  }

  void
  SolverBase::set_transformation (SLEPcWrappers::TransformationBase &this_transformation)
  {
    transformation = &this_transformation;
  }

  void
  SolverBase::set_target_eigenvalue (const PetscScalar &this_target)
  {
    target_eigenvalue.reset (new PetscScalar(this_target));
  }

  void
  SolverBase::set_which_eigenpairs (const EPSWhich eps_which)
  {
    set_which = eps_which;
  }

  void
  SolverBase::set_problem_type (const EPSProblemType eps_problem)
  {
    set_problem = eps_problem;
  }

  void
  SolverBase::solve (const unsigned int  n_eigenpairs,
                     unsigned int       *n_converged)
  {
    int ierr;

    // create a solver object if this is necessary
    if (solver_data.get() == 0)
      {
        // reset solver dtaa
        solver_data.reset (new SolverData());

        // create eigensolver context and set operators
        ierr = EPSCreate (mpi_communicator, &solver_data->eps);
        AssertThrow (ierr == 0, ExcSLEPcError(ierr));

        // set eigenspectrum problem type (general/standard)
        AssertThrow (opA, ExcSLEPcWrappersUsageError());
        if (opB)
          ierr = EPSSetOperators (solver_data->eps, *opA, *opB);
        else
          ierr = EPSSetOperators (solver_data->eps, *opA, PETSC_NULL);
        AssertThrow (ierr == 0, ExcSLEPcError(ierr));

        // set runtime options
        set_solver_type (solver_data->eps);
      }

    // set the problem type
    ierr = EPSSetProblemType (solver_data->eps, set_problem);
    AssertThrow (ierr == 0, ExcSLEPcError(ierr));

    // set the initial vector(s) if any
    if (initial_space.size() != 0)
      {

#if DEAL_II_PETSC_VERSION_LT(3,1,0)
        ierr = EPSSetInitialVector (solver_data->eps, &initial_space[0]);
#else
        ierr = EPSSetInitialSpace (solver_data->eps, initial_space.size(), &initial_space[0]);
#endif

        AssertThrow (ierr == 0, ExcSLEPcError(ierr));
      }

    // if a spectral transformation is to be used, set the
    // transformation and target the wanted eigenvalues
    if (transformation)
      {
        // set transformation type if any
        // STSetShift is called inside
        transformation->set_context (solver_data->eps);
      }

    // set target eigenvalues to solve for
    // in all transformation except STSHIFT there is a direct connection between
    // the target and the shift, read more on p41 of SLEPc manual.
    if (target_eigenvalue.get() != 0)
      {
        int ierr = EPSSetTarget (solver_data->eps, *(target_eigenvalue.get()) );
        AssertThrow (ierr == 0, ExcSLEPcError(ierr));
      }

    // set which portion of the eigenspectrum to solve for
    ierr = EPSSetWhichEigenpairs (solver_data->eps, set_which);
    AssertThrow (ierr == 0, ExcSLEPcError(ierr));

    // set number of eigenvectors to compute
    ierr = EPSSetDimensions (solver_data->eps, n_eigenpairs,
                             PETSC_DECIDE, PETSC_DECIDE);
    AssertThrow (ierr == 0, ExcSLEPcError(ierr));

    // set the solve options to the eigenvalue problem solver context
    ierr = EPSSetFromOptions (solver_data->eps);
    AssertThrow (ierr == 0, ExcSLEPcError(ierr));

    // TODO breaks step-36
    // force solvers to use true residual
    //EPSSetTrueResidual(solver_data->eps, PETSC_TRUE);
    //AssertThrow (ierr == 0, ExcSLEPcError(ierr));

    // Set convergence test to be absolute
    ierr = EPSSetConvergenceTest (solver_data->eps, EPS_CONV_ABS);
    AssertThrow (ierr == 0, ExcSLEPcError(ierr));

    // Set the convergence test function
    // ierr = EPSSetConvergenceTestFunction (solver_data->eps, &convergence_test,
    //              reinterpret_cast<void *>(&solver_control));
    // AssertThrow (ierr == 0, ExcSLEPcError(ierr));

    // solve the eigensystem
    ierr = EPSSolve (solver_data->eps);
    AssertThrow (ierr == 0, ExcSLEPcError(ierr));

    // get number of converged eigenstates
    ierr = EPSGetConverged (solver_data->eps,
                            reinterpret_cast<PetscInt *>(n_converged));
    AssertThrow (ierr == 0, ExcSLEPcError(ierr));

    PetscInt n_iterations   = 0;
    PetscReal residual_norm = 0;

    // @todo Investigate elaborating on some of this to act on the
    // complete eigenspectrum
    {
      // get the number of solver iterations
      ierr = EPSGetIterationNumber (solver_data->eps, &n_iterations);
      AssertThrow (ierr == 0, ExcSLEPcError(ierr));

      // get the maximum of residual norm among converged eigenvectors.
      for (unsigned int i = 0; i < *n_converged; i++)
        {
          double residual_norm_i = 0.0;
          // EPSComputeResidualNorm is L2-norm and is not consistent with the stopping criteria
          // used during the solution process.
          // Yet, this is the norm which gives error bounds (Saad, 1992, ch3):
          //   | \lambda - \widehat\lambda | <= ||r||_2
          ierr = EPSComputeResidualNorm (solver_data->eps, i, &residual_norm_i);

          // EPSComputeRelativeError may not be consistent with the stopping criteria
          // used during the solution process. Given EPS_CONV_ABS set above,
          // this can be either the l2 norm or the mass-matrix induced norm
          // when EPS_GHEP is set.
          // ierr = EPSComputeRelativeError (solver_data->eps, i, &residual_norm_i);

          // EPSGetErrorEstimate is consistent with the residual norm
          // used during the solution process. However, it is not guaranteed to
          // be derived from the residual even when EPSSetTrueResidual is set.
          // ierr = EPSGetErrorEstimate (solver_data->eps, i, &residual_norm_i);

          AssertThrow (ierr == 0, ExcSLEPcError(ierr));
          residual_norm = std::max (residual_norm, residual_norm_i);
        }

      // check the solver state
      const SolverControl::State state
        = solver_control.check (n_iterations, residual_norm);

      // get the solver state according to SLEPc
      get_solver_state (state);

      // as SLEPc uses different stopping criteria, we have to omit this step.
      // This can be checked only in conjunction with EPSGetErrorEstimate.
      // and in case of failure: throw exception
      // if (solver_control.last_check () != SolverControl::success)
      //   AssertThrow(false, SolverControl::NoConvergence (solver_control.last_step(),
      //                                                    solver_control.last_value()));
    }
  }

  void
  SolverBase::get_eigenpair (const unsigned int            index,
                             PetscScalar               &eigenvalues,
                             PETScWrappers::VectorBase &eigenvectors)
  {
    AssertThrow (solver_data.get() != 0, ExcSLEPcWrappersUsageError());

    // get converged eigenpair
    int ierr = EPSGetEigenpair (solver_data->eps, index,
                                &eigenvalues, PETSC_NULL,
                                eigenvectors, PETSC_NULL);
    AssertThrow (ierr == 0, ExcSLEPcError(ierr));
  }


  void
  SolverBase::get_eigenpair (const unsigned int         index,
                             double                    &real_eigenvalues,
                             double                    &imag_eigenvalues,
                             PETScWrappers::VectorBase &real_eigenvectors,
                             PETScWrappers::VectorBase &imag_eigenvectors)
  {
#ifndef PETSC_USE_COMPLEX
    AssertThrow (solver_data.get() != 0, ExcSLEPcWrappersUsageError());

    // get converged eigenpair
    int ierr = EPSGetEigenpair (solver_data->eps, index,
                                &real_eigenvalues, &imag_eigenvalues,
                                real_eigenvectors, imag_eigenvectors);
    AssertThrow (ierr == 0, ExcSLEPcError(ierr));
#else
    Assert ((false),
            ExcMessage ("Your PETSc/SLEPc installation was configured with scalar-type complex "
                        "but this function is not defined for complex types."));
#endif
  }


  void
  SolverBase::reset ()
  {
    AssertThrow (solver_data.get() != 0, ExcSLEPcWrappersUsageError());

    // destroy solver object.
    solver_data.reset ();
  }

  EPS *
  SolverBase::get_eps ()
  {
    if (solver_data.get () == 0)
      return NULL;

    return &solver_data->eps;
  }

  void
  SolverBase::get_solver_state (const SolverControl::State state)
  {
    switch (state)
      {
      case ::dealii::SolverControl::iterate:
        solver_data->reason = EPS_CONVERGED_ITERATING;
        break;

      case ::dealii::SolverControl::success:
        solver_data->reason = static_cast<EPSConvergedReason>(1);
        break;

      case ::dealii::SolverControl::failure:
        if (solver_control.last_step() > solver_control.max_steps())
          solver_data->reason = EPS_DIVERGED_ITS;
        else
          solver_data->reason = EPS_DIVERGED_BREAKDOWN;
        break;

      default:
        Assert (false, ExcNotImplemented());
      }
  }

  /* ---------------------- SolverControls ----------------------- */
  SolverControl &
  SolverBase::control () const
  {
    return solver_control;
  }

  int
  SolverBase::convergence_test (EPS          /*eps             */,
                                PetscScalar  /*real_eigenvalue */,
                                PetscScalar  /*imag_eigenvalue */,
                                PetscReal    /*residual_norm   */,
                                PetscReal   */*estimated_error */,
                                void        */*solver_control_x*/)
  {
    // This function is undefined (future reference only).

    // return without failure.
    return 0;
  }

  /* ---------------------- SolverKrylovSchur ------------------------ */
  SolverKrylovSchur::SolverKrylovSchur (SolverControl        &cn,
                                        const MPI_Comm       &mpi_communicator,
                                        const AdditionalData &data)
    :
    SolverBase (cn, mpi_communicator),
    additional_data (data)
  {}

  void
  SolverKrylovSchur::set_solver_type (EPS &eps) const
  {
    int ierr;
    ierr = EPSSetType (eps, const_cast<char *>(EPSKRYLOVSCHUR));
    AssertThrow (ierr == 0, ExcSLEPcError(ierr));

    // hand over the absolute tolerance and the maximum number of
    // iteration steps to the SLEPc convergence criterion.
    ierr = EPSSetTolerances(eps, this->solver_control.tolerance(),
                            this->solver_control.max_steps());
    AssertThrow (ierr == 0, ExcSLEPcError(ierr));
  }

  /* ---------------------- SolverArnoldi ------------------------ */
  SolverArnoldi::AdditionalData::
  AdditionalData (const bool delayed_reorthogonalization)
    :
    delayed_reorthogonalization (delayed_reorthogonalization)
  {}

  SolverArnoldi::SolverArnoldi (SolverControl        &cn,
                                const MPI_Comm       &mpi_communicator,
                                const AdditionalData &data)
    :
    SolverBase (cn, mpi_communicator),
    additional_data (data)
  {}

  void
  SolverArnoldi::set_solver_type (EPS &eps) const
  {
    int ierr;
    ierr = EPSSetType (eps, const_cast<char *>(EPSARNOLDI));
    AssertThrow (ierr == 0, ExcSLEPcError(ierr));

    // hand over the absolute tolerance and the maximum number of
    // iteration steps to the SLEPc convergence criterion.
    ierr = EPSSetTolerances(eps, this->solver_control.tolerance(),
                            this->solver_control.max_steps());
    AssertThrow (ierr == 0, ExcSLEPcError(ierr));

    // if requested, set delayed reorthogonalization in the Arnoldi
    // iteration.
    if (additional_data.delayed_reorthogonalization)
      {
        ierr = EPSArnoldiSetDelayed (eps, PETSC_TRUE);
        AssertThrow (ierr == 0, ExcSLEPcError(ierr));
      }
  }

  /* ---------------------- Lanczos ------------------------ */
  SolverLanczos::AdditionalData::
  AdditionalData(const EPSLanczosReorthogType r)
    : reorthog(r)
  {}

  SolverLanczos::SolverLanczos (SolverControl        &cn,
                                const MPI_Comm       &mpi_communicator,
                                const AdditionalData &data)
    :
    SolverBase (cn, mpi_communicator),
    additional_data (data)
  {}

  void
  SolverLanczos::set_solver_type (EPS &eps) const
  {
    int ierr;
    ierr = EPSSetType (eps, const_cast<char *>(EPSLANCZOS));
    AssertThrow (ierr == 0, ExcSLEPcError(ierr));

    // hand over the absolute tolerance and the maximum number of
    // iteration steps to the SLEPc convergence criterion.
    ierr = EPSSetTolerances (eps, this->solver_control.tolerance(),
                             this->solver_control.max_steps());
    AssertThrow (ierr == 0, ExcSLEPcError(ierr));

    ierr = EPSLanczosSetReorthog(eps,additional_data.reorthog);
    AssertThrow (ierr == 0, ExcSLEPcError(ierr));
  }

  /* ----------------------- Power ------------------------- */
  SolverPower::SolverPower (SolverControl        &cn,
                            const MPI_Comm       &mpi_communicator,
                            const AdditionalData &data)
    :
    SolverBase (cn, mpi_communicator),
    additional_data (data)
  {}

  void
  SolverPower::set_solver_type (EPS &eps) const
  {
    int ierr;
    ierr = EPSSetType (eps, const_cast<char *>(EPSPOWER));
    AssertThrow (ierr == 0, ExcSLEPcError(ierr));

    // hand over the absolute tolerance and the maximum number of
    // iteration steps to the SLEPc convergence criterion.
    ierr = EPSSetTolerances (eps, this->solver_control.tolerance(),
                             this->solver_control.max_steps());
    AssertThrow (ierr == 0, ExcSLEPcError(ierr));
  }

  /* ---------------- Generalized Davidson ----------------- */
  SolverGeneralizedDavidson::AdditionalData::
  AdditionalData(bool double_expansion)
    :  double_expansion(double_expansion)
  {}

  SolverGeneralizedDavidson::SolverGeneralizedDavidson (SolverControl        &cn,
                                                        const MPI_Comm       &mpi_communicator,
                                                        const AdditionalData &data)
    :
    SolverBase (cn, mpi_communicator),
    additional_data (data)
  {}

  void
  SolverGeneralizedDavidson::set_solver_type (EPS &eps) const
  {
#if DEAL_II_PETSC_VERSION_GTE(3,1,0)
    int ierr;
    ierr = EPSSetType (eps, const_cast<char *>(EPSGD));
    AssertThrow (ierr == 0, ExcSLEPcError(ierr));

    // hand over the absolute tolerance and the maximum number of
    // iteration steps to the SLEPc convergence criterion.
    ierr = EPSSetTolerances (eps, this->solver_control.tolerance(),
                             this->solver_control.max_steps());
    AssertThrow (ierr == 0, ExcSLEPcError(ierr));

    if (additional_data.double_expansion)
      {
        ierr = EPSGDSetDoubleExpansion (eps, PETSC_TRUE);
        AssertThrow (ierr == 0, ExcSLEPcError(ierr));
      }
#else
    // Suppress compiler warnings about unused parameters.
    (void) eps;

    // PETSc/SLEPc version must be > 3.1.0.
    Assert ((false),
            ExcMessage ("Your SLEPc installation does not include a copy of the "
                        "Generalized Davidson solver. A SLEPc version > 3.1.0 is required."));
#endif
  }

  /* ------------------ Jacobi Davidson -------------------- */
  SolverJacobiDavidson::SolverJacobiDavidson (SolverControl        &cn,
                                              const MPI_Comm       &mpi_communicator,
                                              const AdditionalData &data)
    :
    SolverBase (cn, mpi_communicator),
    additional_data (data)
  {}

  void
  SolverJacobiDavidson::set_solver_type (EPS &eps) const
  {
#if DEAL_II_PETSC_VERSION_GTE(3,1,0)
    int ierr;
    ierr = EPSSetType (eps, const_cast<char *>(EPSJD));
    AssertThrow (ierr == 0, ExcSLEPcError(ierr));

    // hand over the absolute tolerance and the maximum number of
    // iteration steps to the SLEPc convergence criterion.
    ierr = EPSSetTolerances (eps, this->solver_control.tolerance(),
                             this->solver_control.max_steps());
    AssertThrow (ierr == 0, ExcSLEPcError(ierr));
#else
    // Suppress compiler warnings about unused parameters.
    (void) eps;

    // PETSc/SLEPc version must be > 3.1.0.
    Assert ((false),
            ExcMessage ("Your SLEPc installation does not include a copy of the "
                        "Jacobi-Davidson solver. A SLEPc version > 3.1.0 is required."));
#endif
  }


  /* ---------------------- LAPACK ------------------------- */
  SolverLAPACK::SolverLAPACK (SolverControl        &cn,
                              const MPI_Comm       &mpi_communicator,
                              const AdditionalData &data)
    :
    SolverBase (cn, mpi_communicator),
    additional_data (data)
  {}

  void
  SolverLAPACK::set_solver_type (EPS &eps) const
  {
    // 'Tis overwhelmingly likely that PETSc/SLEPc *always* has
    // BLAS/LAPACK, but let's be defensive.
#if PETSC_HAVE_BLASLAPACK
    int ierr;
    ierr = EPSSetType (eps, const_cast<char *>(EPSLAPACK));
    AssertThrow (ierr == 0, ExcSLEPcError(ierr));

    // hand over the absolute tolerance and the maximum number of
    // iteration steps to the SLEPc convergence criterion.
    ierr = EPSSetTolerances (eps, this->solver_control.tolerance(),
                             this->solver_control.max_steps());
    AssertThrow (ierr == 0, ExcSLEPcError(ierr));
#else
    // Suppress compiler warnings about unused parameters.
    (void) eps;

    Assert ((false),
            ExcMessage ("Your PETSc/SLEPc installation was not configured with BLAS/LAPACK "
                        "but this is needed to use the LAPACK solver."));
#endif
  }

}

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_SLEPC

