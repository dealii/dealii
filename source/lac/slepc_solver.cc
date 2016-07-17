// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2016 by the deal.II authors
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

  SolverBase::SolverBase (SolverControl  &cn,
                          const MPI_Comm &mpi_communicator)
    :
    solver_control (cn),
    mpi_communicator (mpi_communicator)
  {
    // create eigensolver context
    int ierr = EPSCreate (mpi_communicator, &eps);
    AssertThrow (ierr == 0, ExcSLEPcError(ierr));

    // hand over the absolute tolerance and the maximum number of
    // iteration steps to the SLEPc convergence criterion.
    ierr = EPSSetTolerances(eps, this->solver_control.tolerance(),
                            this->solver_control.max_steps());
    AssertThrow (ierr == 0, ExcSLEPcError(ierr));

    // default values:
    set_which_eigenpairs(EPS_LARGEST_MAGNITUDE);
    set_problem_type(EPS_GNHEP);

    // TODO:
    // By default, EPS initializes the starting vector or the initial subspace randomly.
  }

  SolverBase::~SolverBase ()
  {
    if (eps != NULL)
      {
        // Destroy the solver object.
#if DEAL_II_PETSC_VERSION_LT(3,2,0)
        int ierr = EPSDestroy (eps);
#else
        int ierr = EPSDestroy (&eps);
#endif
        AssertThrow (ierr == 0, ExcSLEPcError(ierr));
      }
  }

  void
  SolverBase::set_matrices (const PETScWrappers::MatrixBase &A)
  {
    // standard eigenspectrum problem
    int ierr = EPSSetOperators (eps, A, PETSC_NULL);
    AssertThrow (ierr == 0, ExcSLEPcError(ierr));
  }

  void
  SolverBase::set_matrices (const PETScWrappers::MatrixBase &A,
                            const PETScWrappers::MatrixBase &B)
  {
    // generalized eigenspectrum problem
    int ierr = EPSSetOperators (eps, A, B);
    AssertThrow (ierr == 0, ExcSLEPcError(ierr));
  }

  void
  SolverBase::set_transformation (SLEPcWrappers::TransformationBase &transformation)
  {
    // set transformation type if any
    // STSetShift is called inside
    int ierr = EPSSetST(eps,transformation.st);
    AssertThrow (ierr == 0, SolverBase::ExcSLEPcError(ierr));
  }

  void
  SolverBase::set_initial_vector (const PETScWrappers::VectorBase &this_initial_vector)
  {
    Assert(this_initial_vector.l2_norm()>0.0,
           ExcMessage("Initial vector should be nonzero."));

    int ierr;
    Vec vec = this_initial_vector;
#if DEAL_II_PETSC_VERSION_LT(3,1,0)
    ierr = EPSSetInitialVector (eps, &vec);
#else
    ierr = EPSSetInitialSpace (eps, 1, &vec);
#endif
    AssertThrow (ierr == 0, ExcSLEPcError(ierr));
  }

  void
  SolverBase::set_target_eigenvalue (const PetscScalar &this_target)
  {
    // set target eigenvalues to solve for
    // in all transformation except STSHIFT there is a direct connection between
    // the target and the shift, read more on p41 of SLEPc manual.
    int ierr = EPSSetTarget (eps, this_target );
    AssertThrow (ierr == 0, ExcSLEPcError(ierr));
  }

  void
  SolverBase::set_which_eigenpairs (const EPSWhich eps_which)
  {
    // set which portion of the eigenspectrum to solve for
    int ierr = EPSSetWhichEigenpairs (eps, eps_which);
    AssertThrow (ierr == 0, ExcSLEPcError(ierr));
  }

  void
  SolverBase::set_problem_type (const EPSProblemType eps_problem)
  {
    int ierr = EPSSetProblemType (eps, eps_problem);
    AssertThrow (ierr == 0, ExcSLEPcError(ierr));
  }

  void
  SolverBase::solve (const unsigned int  n_eigenpairs,
                     unsigned int       *n_converged)
  {
    int ierr;

    // set number of eigenvectors to compute
    ierr = EPSSetDimensions (eps, n_eigenpairs,
                             PETSC_DECIDE, PETSC_DECIDE);
    AssertThrow (ierr == 0, ExcSLEPcError(ierr));

    // set the solve options to the eigenvalue problem solver context
    ierr = EPSSetFromOptions (eps);
    AssertThrow (ierr == 0, ExcSLEPcError(ierr));

    // TODO breaks step-36
    // force Krylov solver to use true residual instead of an estimate.
    //EPSSetTrueResidual(solver_data->eps, PETSC_TRUE);
    //AssertThrow (ierr == 0, ExcSLEPcError(ierr));

    // Set convergence test to be absolute
    ierr = EPSSetConvergenceTest (eps, EPS_CONV_ABS);
    AssertThrow (ierr == 0, ExcSLEPcError(ierr));

    // TODO Set the convergence test function
    // ierr = EPSSetConvergenceTestFunction (solver_data->eps, &convergence_test,
    //              reinterpret_cast<void *>(&solver_control));
    // AssertThrow (ierr == 0, ExcSLEPcError(ierr));

    // solve the eigensystem
    ierr = EPSSolve (eps);
    AssertThrow (ierr == 0, ExcSLEPcError(ierr));

    // get number of converged eigenstates
    ierr = EPSGetConverged (eps,
                            reinterpret_cast<PetscInt *>(n_converged));
    AssertThrow (ierr == 0, ExcSLEPcError(ierr));

    PetscInt n_iterations   = 0;
    PetscReal residual_norm = 0;

    // @todo Investigate elaborating on some of this to act on the
    // complete eigenspectrum
    {
      // get the number of solver iterations
      ierr = EPSGetIterationNumber (eps, &n_iterations);
      AssertThrow (ierr == 0, ExcSLEPcError(ierr));

      // get the maximum of residual norm among converged eigenvectors.
      for (unsigned int i = 0; i < *n_converged; i++)
        {
          double residual_norm_i = 0.0;
          // EPSComputeResidualNorm is L2-norm and is not consistent with the stopping criteria
          // used during the solution process.
          // Yet, this is the norm which gives error bounds (Saad, 1992, ch3):
          //   | \lambda - \widehat\lambda | <= ||r||_2
          ierr = EPSComputeResidualNorm (eps, i, &residual_norm_i);

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
    // get converged eigenpair
    int ierr = EPSGetEigenpair (eps, index,
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
    // get converged eigenpair
    int ierr = EPSGetEigenpair (eps, index,
                                &real_eigenvalues, &imag_eigenvalues,
                                real_eigenvectors, imag_eigenvectors);
    AssertThrow (ierr == 0, ExcSLEPcError(ierr));
#else
    Assert ((false),
            ExcMessage ("Your PETSc/SLEPc installation was configured with scalar-type complex "
                        "but this function is not defined for complex types."));

    // Cast to void to silence compiler warnings
    (void) index;
    (void) real_eigenvalues;
    (void) imag_eigenvalues;
    (void) real_eigenvectors;
    (void) imag_eigenvectors;
#endif
  }

  void
  SolverBase::get_solver_state (const SolverControl::State state)
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
                                PetscReal    /*residual norm associated to the eigenpair   */,
                                PetscReal   */*(output) computed error estimate */,
                                void        */*solver_control_x*/)
  {
    // If the error estimate returned by the convergence test function is less
    // than the tolerance, then the eigenvalue is accepted as converged.
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
  {
    int ierr = EPSSetType (eps, const_cast<char *>(EPSKRYLOVSCHUR));
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
  {
    int ierr = EPSSetType (eps, const_cast<char *>(EPSARNOLDI));
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
  {
    int ierr = EPSSetType (eps, const_cast<char *>(EPSLANCZOS));
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
  {
    int ierr = EPSSetType (eps, const_cast<char *>(EPSPOWER));
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
  {
#if DEAL_II_PETSC_VERSION_GTE(3,1,0)
    int ierr = EPSSetType (eps, const_cast<char *>(EPSGD));
    AssertThrow (ierr == 0, ExcSLEPcError(ierr));

    if (additional_data.double_expansion)
      {
        ierr = EPSGDSetDoubleExpansion (eps, PETSC_TRUE);
        AssertThrow (ierr == 0, ExcSLEPcError(ierr));
      }
#else
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
  {
#if DEAL_II_PETSC_VERSION_GTE(3,1,0)
    int ierr;
    ierr = EPSSetType (eps, const_cast<char *>(EPSJD));
    AssertThrow (ierr == 0, ExcSLEPcError(ierr));
#else
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
  {
    // 'Tis overwhelmingly likely that PETSc/SLEPc *always* has
    // BLAS/LAPACK, but let's be defensive.
#if PETSC_HAVE_BLASLAPACK
    int ierr;
    ierr = EPSSetType (eps, const_cast<char *>(EPSLAPACK));
    AssertThrow (ierr == 0, ExcSLEPcError(ierr));
#else
    Assert ((false),
            ExcMessage ("Your PETSc/SLEPc installation was not configured with BLAS/LAPACK "
                        "but this is needed to use the LAPACK solver."));
#endif
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_SLEPC

