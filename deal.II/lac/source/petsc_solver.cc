//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <lac/petsc_solver.h>
#include <lac/petsc_matrix_base.h>
#include <lac/petsc_vector_base.h>
#include <lac/petsc_precondition.h>

#include <cmath>

#ifdef DEAL_II_USE_PETSC

#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR < 2)
#  include <petscsles.h>
#endif
#include <petscversion.h>


namespace PETScWrappers
{

  SolverBase::SolverData::~SolverData ()
  {
                                     // destroy the solver object
    int ierr = KSPDestroy (ksp);
    AssertThrow (ierr == 0, ExcPETScError(ierr)); 
  
                                     // and destroy the solver object if we
                                     // are in old PETSc mode
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR < 2)
    ierr = SLESDestroy (sles);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
#endif
  }

  
  
  SolverBase::SolverBase (SolverControl  &cn,
                          const MPI_Comm &mpi_communicator)
                  :
                  solver_control (cn),
                  mpi_communicator (mpi_communicator)
  {}

  

  SolverBase::~SolverBase ()
  {}

  

  void
  SolverBase::solve (const MatrixBase         &A,
                     VectorBase               &x,
                     const VectorBase         &b,
                     const PreconditionerBase &preconditioner)
  {
    int ierr;

                                     // PETSc unfortunately changed its
                                     // interfaces completely and incompatibly
                                     // between 2.1.6 and 2.2.0. This forces
                                     // us to work around this by offering two
                                     // completely different
                                     // implementations. sigh...
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR < 2)
                                     // first create a solver object
    SLES sles;
    ierr = SLESCreate (mpi_communicator, &sles);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
    
                                     // set the matrices involved. the
                                     // last argument is irrelevant here,
                                     // since we use the solver only once
                                     // anyway
    ierr = SLESSetOperators (sles, A, preconditioner,
                             SAME_PRECONDITIONER);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
    
                                     // let derived classes set the solver
                                     // type, and the preconditioning object
                                     // set the type of preconditioner
    KSP ksp;
    ierr = SLESGetKSP (sles, &ksp);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
    
    set_solver_type (ksp);
    
    PC pc;
    ierr = SLESGetPC (sles, &pc);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
    
    preconditioner.set_preconditioner_type (pc);
  
                                     // in the deal.II solvers, we always
                                     // honor the initial guess in the
                                     // solution vector. do so here as well:
    KSPSetInitialGuessNonzero (ksp, PETSC_TRUE);
  
                                     // then a convergence monitor
                                     // function. that function simply checks
                                     // with the solver_control object we have
                                     // in this object for convergence
    KSPSetConvergenceTest (ksp, &convergence_test,
                           reinterpret_cast<void *>(&solver_control));
      
                                     // then do the real work: set up solver
                                     // internal data and solve the
                                     // system. this could be joined in one
                                     // operation, but it is recommended this
                                     // way to be able to separate statistic
                                     // output if requested
    int iterations = 0;
    ierr = SLESSetUp (sles, b, x);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
      
    ierr = SLESSolve (sles, b, x, &iterations);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

#else
                                     // and here is the PETSc post 2.1.6
                                     // version: SLES objects have been
                                     // completely abandoned, so use KSP and
                                     // PC only

                                     // first create a solver object if this
                                     // is necessary
    if (solver_data == 0)
      {
        solver_data.reset (new SolverData());
        
        ierr = KSPCreate (mpi_communicator, &solver_data->ksp);
        AssertThrow (ierr == 0, ExcPETScError(ierr));

                                         // set the matrices involved. the
                                         // last argument is irrelevant here,
                                         // since we use the solver only once
                                         // anyway
        ierr = KSPSetOperators (solver_data->ksp, A, preconditioner,
                                SAME_PRECONDITIONER);
        AssertThrow (ierr == 0, ExcPETScError(ierr));

                                         // let derived classes set the solver
                                         // type, and the preconditioning
                                         // object set the type of
                                         // preconditioner
        set_solver_type (solver_data->ksp);
    
        ierr = KSPGetPC (solver_data->ksp, &solver_data->pc);
        AssertThrow (ierr == 0, ExcPETScError(ierr));

        preconditioner.set_preconditioner_type (solver_data->pc);

                                         // in the deal.II solvers, we always
                                         // honor the initial guess in the
                                         // solution vector. do so here as
                                         // well:
        KSPSetInitialGuessNonzero (solver_data->ksp, PETSC_TRUE);

                                         // then a convergence monitor
                                         // function. that function simply
                                         // checks with the solver_control
                                         // object we have in this object for
                                         // convergence
        KSPSetConvergenceTest (solver_data->ksp, &convergence_test,
                               reinterpret_cast<void *>(&solver_control));
      }
    
                                     // then do the real work: set up solver
                                     // internal data and solve the
                                     // system.
    ierr = KSPSetRhs(solver_data->ksp,b);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
    ierr = KSPSetSolution(solver_data->ksp, x);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
    ierr = KSPSetUp (solver_data->ksp);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
    
    ierr = KSPSolve (solver_data->ksp);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

#endif
                                     // destroy solver object
    solver_data.reset ();

                                     // in case of failure: throw
                                     // exception
    if (solver_control.last_check() != SolverControl::success)
      throw SolverControl::NoConvergence (solver_control.last_step(),
                                          solver_control.last_value());
                                     // otherwise exit as normal
  }



  SolverControl &
  SolverBase::control() const
  {
    return solver_control;
  }
  


  int
  SolverBase::convergence_test (KSP                 /*ksp*/,
                                const int           iteration,
                                const PetscScalar   residual_norm,
                                KSPConvergedReason *reason,
                                void               *solver_control_x)
  {
    SolverControl &solver_control = *reinterpret_cast<SolverControl*>(solver_control_x);
    
    const SolverControl::State state
        = solver_control.check (iteration, residual_norm);
    
    switch (state)
      {
        case ::SolverControl::iterate:
              *reason = KSP_CONVERGED_ITERATING;
              break;
              
        case ::SolverControl::success:
              *reason = static_cast<KSPConvergedReason>(1);
              break;
              
        case ::SolverControl::failure:
              if (solver_control.last_step() > solver_control.max_steps())
                *reason = KSP_DIVERGED_ITS;
              else
                *reason = KSP_DIVERGED_DTOL;
              break;
              
        default:
              Assert (false, ExcNotImplemented());
      }

                                     // return without failure
    return 0;
  }
  
    

/* ---------------------- SolverRichardson ------------------------ */


  SolverRichardson::AdditionalData::
  AdditionalData (const double omega)
                  :
                  omega (omega)
  {}

  
  
  SolverRichardson::SolverRichardson (SolverControl        &cn,
                                      const MPI_Comm       &mpi_communicator,
                                      const AdditionalData &data)
                  :
                  SolverBase (cn, mpi_communicator),
                  additional_data (data)
  {}


  void
  SolverRichardson::set_solver_type (KSP &ksp) const
  {
                                     // set the type of solver. work around a
                                     // problem in PETSc 2.1.6, where it asks
                                     // for a char*, even though KSPXXXX is of
                                     // type const char*
    int ierr;
    ierr = KSPSetType (ksp, const_cast<char *>(KSPRICHARDSON));
    AssertThrow (ierr == 0, ExcPETScError(ierr));

                                     // set the damping factor from the data
    ierr = KSPRichardsonSetScale (ksp, additional_data.omega);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }
  

/* ---------------------- SolverChebychev ------------------------ */

  SolverChebychev::SolverChebychev (SolverControl        &cn,
                                    const MPI_Comm       &mpi_communicator,
                                    const AdditionalData &data)
                  :
                  SolverBase (cn, mpi_communicator),
                  additional_data (data)
  {}


  void
  SolverChebychev::set_solver_type (KSP &ksp) const
  {
                                     // set the type of solver. work around a
                                     // problem in PETSc 2.1.6, where it asks
                                     // for a char*, even though KSPXXXX is of
                                     // type const char*
    int ierr;
    ierr = KSPSetType (ksp, const_cast<char *>(KSPCHEBYCHEV));
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }

  
/* ---------------------- SolverCG ------------------------ */

  SolverCG::SolverCG (SolverControl        &cn,
                      const MPI_Comm       &mpi_communicator,
                      const AdditionalData &data)
                  :
                  SolverBase (cn, mpi_communicator),
                  additional_data (data)
  {}


  void
  SolverCG::set_solver_type (KSP &ksp) const
  {
                                     // set the type of solver. work around a
                                     // problem in PETSc 2.1.6, where it asks
                                     // for a char*, even though KSPXXXX is of
                                     // type const char*
    int ierr;
    ierr = KSPSetType (ksp, const_cast<char *>(KSPCG));
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }
  

/* ---------------------- SolverBiCG ------------------------ */

  SolverBiCG::SolverBiCG (SolverControl        &cn,
                          const MPI_Comm       &mpi_communicator,
                          const AdditionalData &data)
                  :
                  SolverBase (cn, mpi_communicator),
                  additional_data (data)
  {}


  void
  SolverBiCG::set_solver_type (KSP &ksp) const
  {
                                     // set the type of solver. work around a
                                     // problem in PETSc 2.1.6, where it asks
                                     // for a char*, even though KSPXXXX is of
                                     // type const char*
    int ierr;
    ierr = KSPSetType (ksp, const_cast<char *>(KSPBICG));
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }

  
/* ---------------------- SolverGMRES ------------------------ */

  SolverGMRES::AdditionalData::
  AdditionalData (const unsigned int restart_parameter)
                  :
                  restart_parameter (restart_parameter)
  {}

  
  
  SolverGMRES::SolverGMRES (SolverControl        &cn,
                            const MPI_Comm       &mpi_communicator,
                            const AdditionalData &data)
                  :
                  SolverBase (cn, mpi_communicator),
                  additional_data (data)
  {}


  void
  SolverGMRES::set_solver_type (KSP &ksp) const
  {
                                     // set the type of solver. work around a
                                     // problem in PETSc 2.1.6, where it asks
                                     // for a char*, even though KSPXXXX is of
                                     // type const char*
    int ierr;
    ierr = KSPSetType (ksp, const_cast<char *>(KSPGMRES));
    AssertThrow (ierr == 0, ExcPETScError(ierr));

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
                                     // and do some equally nast stuff that at
                                     // least doesn't yield warnings...
    int (*fun_ptr)(KSP,int); 
    ierr = PetscObjectQueryFunction((PetscObject)(ksp),
                                    "KSPGMRESSetRestart_C",
                                    (PetscVoidFunction)&fun_ptr);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    ierr = (*fun_ptr)(ksp,additional_data.restart_parameter);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }
  

/* ---------------------- SolverBicgstab ------------------------ */

  SolverBicgstab::SolverBicgstab (SolverControl        &cn,
                                  const MPI_Comm       &mpi_communicator,
                                  const AdditionalData &data)
                  :
                  SolverBase (cn, mpi_communicator),
                  additional_data (data)
  {}


  void
  SolverBicgstab::set_solver_type (KSP &ksp) const
  {
                                     // set the type of solver. work around a
                                     // problem in PETSc 2.1.6, where it asks
                                     // for a char*, even though KSPXXXX is of
                                     // type const char*
    int ierr;
    ierr = KSPSetType (ksp, const_cast<char *>(KSPBCGS));
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }

  
/* ---------------------- SolverCGS ------------------------ */

  SolverCGS::SolverCGS (SolverControl        &cn,
                        const MPI_Comm       &mpi_communicator,
                        const AdditionalData &data)
                  :
                  SolverBase (cn, mpi_communicator),
                  additional_data (data)
  {}


  void
  SolverCGS::set_solver_type (KSP &ksp) const
  {
                                     // set the type of solver. work around a
                                     // problem in PETSc 2.1.6, where it asks
                                     // for a char*, even though KSPXXXX is of
                                     // type const char*
    int ierr;
    ierr = KSPSetType (ksp, const_cast<char *>(KSPCGS));
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }
  

/* ---------------------- SolverTFQMR ------------------------ */

  SolverTFQMR::SolverTFQMR (SolverControl        &cn,
                            const MPI_Comm       &mpi_communicator,
                            const AdditionalData &data)
                  :
                  SolverBase (cn, mpi_communicator),
                  additional_data (data)
  {}


  void
  SolverTFQMR::set_solver_type (KSP &ksp) const
  {
                                     // set the type of solver. work around a
                                     // problem in PETSc 2.1.6, where it asks
                                     // for a char*, even though KSPXXXX is of
                                     // type const char*
    int ierr;
    ierr = KSPSetType (ksp, const_cast<char *>(KSPTFQMR));
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }

  
/* ---------------------- SolverTCQMR ------------------------ */

  SolverTCQMR::SolverTCQMR (SolverControl        &cn,
                            const MPI_Comm       &mpi_communicator,
                            const AdditionalData &data)
                  :
                  SolverBase (cn, mpi_communicator),
                  additional_data (data)
  {}


  void
  SolverTCQMR::set_solver_type (KSP &ksp) const
  {
                                     // set the type of solver. work around a
                                     // problem in PETSc 2.1.6, where it asks
                                     // for a char*, even though KSPXXXX is of
                                     // type const char*
    int ierr;
    ierr = KSPSetType (ksp, const_cast<char *>(KSPTCQMR));
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }

  
/* ---------------------- SolverCR ------------------------ */

  SolverCR::SolverCR (SolverControl        &cn,
                      const MPI_Comm       &mpi_communicator,
                      const AdditionalData &data)
                  :
                  SolverBase (cn, mpi_communicator),
                  additional_data (data)
  {}


  void
  SolverCR::set_solver_type (KSP &ksp) const
  {
                                     // set the type of solver. work around a
                                     // problem in PETSc 2.1.6, where it asks
                                     // for a char*, even though KSPXXXX is of
                                     // type const char*
    int ierr;
    ierr = KSPSetType (ksp, const_cast<char *>(KSPCR));
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }

  
/* ---------------------- SolverLSQR ------------------------ */

  SolverLSQR::SolverLSQR (SolverControl        &cn,
                          const MPI_Comm       &mpi_communicator,
                          const AdditionalData &data)
                  :
                  SolverBase (cn, mpi_communicator),
                  additional_data (data)
  {}


  void
  SolverLSQR::set_solver_type (KSP &ksp) const
  {
                                     // set the type of solver. work around a
                                     // problem in PETSc 2.1.6, where it asks
                                     // for a char*, even though KSPXXXX is of
                                     // type const char*
    int ierr;
    ierr = KSPSetType (ksp, const_cast<char *>(KSPLSQR));
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }

  
}

#else
// On gcc2.95 on Alpha OSF1, the native assembler does not like empty
// files, so provide some dummy code
namespace { void dummy () {} }
#endif // DEAL_II_USE_PETSC
