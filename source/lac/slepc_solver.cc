// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2013 by the deal.II authors
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
    target_eigenvalue (0.),
    set_which (EPS_LARGEST_MAGNITUDE),
    set_problem (EPS_GNHEP),
    opA (NULL),
    opB (NULL),
    initial_vector (NULL),
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
    initial_vector = (&this_initial_vector);
  }

  void
  SolverBase::set_transformation (SLEPcWrappers::TransformationBase &this_transformation)
  {
    transformation = &this_transformation;
  }

  void
  SolverBase::set_target_eigenvalue (const PetscScalar &this_target)
  {
    target_eigenvalue = this_target;
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
    if (initial_vector && initial_vector->size() != 0)
      {

#if DEAL_II_PETSC_VERSION_LT(3,1,0)
        ierr = EPSSetInitialVector (solver_data->eps, *initial_vector);
#else
        Vec this_vector = *initial_vector;
        ierr = EPSSetInitialSpace (solver_data->eps, 1, &this_vector);
#endif

        AssertThrow (ierr == 0, ExcSLEPcError(ierr));
      }

    // if a spectral transformation is to be used, set the
    // transformation and target the wanted eigenvalues
    if (transformation)
      {
        // set transformation type if any
        transformation->set_context (solver_data->eps);

        // set target eigenvalues to solve for
        ierr = EPSSetTarget (solver_data->eps, target_eigenvalue);
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
    PetscReal residual_norm = 1.e300;

    // @todo Investigate elaborating on some of this to act on the
    // complete eigenspectrum
    {
      // get the number of solver iterations
      ierr = EPSGetIterationNumber (solver_data->eps, &n_iterations);
      AssertThrow (ierr == 0, ExcSLEPcError(ierr));

      // get the residual norm of the most extreme eigenvalue if and
      // only if at least one eigenvector has converged.
      if ((*n_converged)>0)
        {
          ierr = EPSComputeResidualNorm (solver_data->eps, 0, &residual_norm);
          AssertThrow (ierr == 0, ExcSLEPcError(ierr));
        }

      // check the solver state
      const SolverControl::State state
        = solver_control.check (n_iterations, residual_norm);

      // get the solver state according to SLEPc
      get_solver_state (state);

      // and in case of failure: throw exception
      if (solver_control.last_check () != SolverControl::success)
        AssertThrow(false, SolverControl::NoConvergence (solver_control.last_step(),
                                                         solver_control.last_value()));
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

