
#include <lac/petsc_matrix_base.h>
#include <lac/petsc_vector_base.h>
#include <lac/petsc_vector.h>
#include <lac/slepc_solver.h>
#include <lac/slepc_spectral_transformation.h>

#include <cmath>
#include <vector>

#ifdef DEAL_II_USE_SLEPC

#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR < 2)
#include <petscsles.h>
#endif
#include <petscversion.h>

DEAL_II_NAMESPACE_OPEN

namespace SLEPcWrappers
{

  SolverBase::SolverData::~SolverData ()
  {
                                   // Destroy the solver object.
    int ierr = EPSDestroy (eps);
    AssertThrow (ierr == 0, ExcSLEPcError(ierr)); 
  }
  
  SolverBase::SolverBase (SolverControl  &cn,
                          const MPI_Comm &mpi_communicator)
    :
    solver_control (cn),
    mpi_communicator (mpi_communicator),
    set_which (EPS_LARGEST_MAGNITUDE),
    opA (NULL), opB (NULL),
    ini_vec (NULL),
    transform (NULL)
  {
  }
  
  SolverBase::~SolverBase ()
  {
    if( solver_data != 0 )
      solver_data.reset ();
  }

  void
  SolverBase::set_matrices (const PETScWrappers::MatrixBase &A)
  {
    opA = &A;
    opB = NULL;
  }

  void
  SolverBase::set_matrices (const PETScWrappers::MatrixBase &A,
			    const PETScWrappers::MatrixBase &B)
  {
    opA = &A;
    opB = &B;
  }

  void
  SolverBase::set_initial_vector (const PETScWrappers::VectorBase &initial_vec)
  {
    ini_vec = &initial_vec;
  }

  void
  SolverBase::set_transformation (SLEPcWrappers::TransformationBase &trans)
  {
    transform = &trans;
  }

  void
  SolverBase::set_which_eigenpairs (const EPSWhich eps_which)
  {
    set_which = eps_which;
  }

  void
  SolverBase::solve (const unsigned int n_eigenvectors, unsigned int *n_converged)
  {
    int ierr;
    
    AssertThrow (solver_data == 0, ExcSLEPcWrappersUsageError());    
    solver_data.reset (new SolverData());
    
                                    // create eigensolver context and
                                    // set operators.
    ierr = EPSCreate (mpi_communicator, &solver_data->eps);
    AssertThrow (ierr == 0, ExcSLEPcError(ierr));

    AssertThrow (opA, ExcSLEPcWrappersUsageError());
    if (opB)
      ierr = EPSSetOperators (solver_data->eps, *opA, *opB);
    else
      ierr = EPSSetOperators (solver_data->eps, *opA, PETSC_NULL);
    AssertThrow (ierr == 0, ExcSLEPcError(ierr));
    
    if (ini_vec && ini_vec->size() != 0) 
      {
	ierr = EPSSetInitialVector(solver_data->eps, *ini_vec);
	AssertThrow (ierr == 0, ExcSLEPcError(ierr));
      }
    
    if (transform)
      transform->set_context(solver_data->eps);
    
                                    // set runtime options.
    set_solver_type (solver_data->eps);

    ierr = EPSSetWhichEigenpairs (solver_data->eps, set_which);
    AssertThrow (ierr == 0, ExcSLEPcError(ierr));

                                    // set number of eigenvectors to
                                    // compute
    ierr = EPSSetDimensions (solver_data->eps, n_eigenvectors, 
			     PETSC_DECIDE, PETSC_DECIDE);
    AssertThrow (ierr == 0, ExcSLEPcError(ierr));

    ierr = EPSSetFromOptions (solver_data->eps);
    AssertThrow (ierr == 0, ExcSLEPcError(ierr));

                                    // solve the eigensystem
    ierr = EPSSolve (solver_data->eps);
    AssertThrow (ierr == 0, ExcSLEPcError(ierr));

                                    // get number of converged
                                    // eigenstates
    ierr = EPSGetConverged (solver_data->eps, 
#ifdef PETSC_USE_64BIT_INDICES
			    reinterpret_cast<PetscInt *>(n_converged));
#else
			    reinterpret_cast<int *>(n_converged));
#endif
    AssertThrow (ierr == 0, ExcSLEPcError(ierr));
  }

  void
  SolverBase::get_eigenpair (const unsigned int index, 
			     double &kr,
			     PETScWrappers::VectorBase &vr) 
  {
    AssertThrow (solver_data != 0, ExcSLEPcWrappersUsageError());    

                                    // get converged eigenpair
    int ierr = EPSGetEigenpair(solver_data->eps, index, 
			       &kr, PETSC_NULL, vr, PETSC_NULL);
    AssertThrow (ierr == 0, ExcSLEPcError(ierr));
  }


  void
  SolverBase::reset ()
  {
    AssertThrow (solver_data != 0, ExcSLEPcWrappersUsageError());

                                    // destroy solver object.
    solver_data.reset ();
  }

  EPS *
  SolverBase::get_EPS ()
  {
    if( solver_data == 0 )
      return NULL;
    return &solver_data->eps;
  }
    
  /* ---------------------- SolverControls ----------------------- */
  SolverControl &
  SolverBase::control () const
  {
    return solver_control;
  }

  int
  SolverBase::convergence_test (EPS                 /*eps*/,
                                const int           iteration,
                                const PetscScalar   residual_norm,
                                EPSConvergedReason *reason,
                                void               *solver_control_x)
  {
    SolverControl &solver_control 
      = *reinterpret_cast<SolverControl*>(solver_control_x);
    
    const SolverControl::State state
      = solver_control.check (iteration, residual_norm);
    
    switch (state)
      {
      case ::dealii::SolverControl::iterate:
	*reason = EPS_CONVERGED_ITERATING;
	break;
	
      case ::dealii::SolverControl::success:
	*reason = static_cast<EPSConvergedReason>(1);
	break;
	
      case ::dealii::SolverControl::failure:
	if (solver_control.last_step() > solver_control.max_steps())
	  *reason = EPS_DIVERGED_ITS;
	break;
	
      default:
	Assert (false, ExcNotImplemented());
      }
    
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

                                    // hand over the absolute
                                    // tolerance and the maximum
                                    // number of iteration steps to
                                    // the SLEPc convergence
                                    // criterion.
    ierr = EPSSetTolerances(eps, this->solver_control.tolerance(),
			    this->solver_control.max_steps());
    AssertThrow (ierr == 0, ExcSLEPcError(ierr));
  }
  
  /* ---------------------- SolverArnoldi ------------------------ */
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

                                    // hand over the absolute
                                    // tolerance and the maximum
                                    // number of iteration steps to
                                    // the SLEPc convergence
                                    // criterion.
    ierr = EPSSetTolerances(eps, this->solver_control.tolerance(),
			    this->solver_control.max_steps());
    AssertThrow (ierr == 0, ExcSLEPcError(ierr));
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

                                    // hand over the absolute
                                    // tolerance and the maximum
                                    // number of iteration steps to
                                    // the SLEPc convergence
                                    // criterion.
    ierr = EPSSetTolerances(eps, this->solver_control.tolerance(),
			    this->solver_control.max_steps());
    AssertThrow (ierr == 0, ExcSLEPcError(ierr));
  }

}

DEAL_II_NAMESPACE_CLOSE

#else
// On gcc2.95 on Alpha OSF1, the native assembler does not like empty
// files, so provide some dummy code
namespace { void dummy () {} }
#endif // DEAL_II_USE_SLEPC

