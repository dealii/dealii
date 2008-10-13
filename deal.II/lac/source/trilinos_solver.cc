//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <lac/trilinos_solver.h>
#include <lac/trilinos_sparse_matrix.h>
#include <lac/trilinos_vector_base.h>
#include <lac/trilinos_precondition.h>

#include <cmath>

#ifdef DEAL_II_USE_TRILINOS

#include <Epetra_Operator.h>

DEAL_II_NAMESPACE_OPEN

namespace TrilinosWrappers
{

  SolverBase::AdditionalData::AdditionalData (const unsigned int gmres_restart_parameter)
                  :
                  gmres_restart_parameter (gmres_restart_parameter)
  {}

  
  
  SolverBase::SolverBase (SolverControl  &cn)
                  :
		  solver_name    (gmres),
                  solver_control (cn)
  {}

  

  SolverBase::SolverBase (const enum SolverBase::SolverName  solver_name,
			  SolverControl                     &cn)
                  :
		  solver_name    (solver_name),
                  solver_control (cn)
  {}

  

  SolverBase::~SolverBase ()
  {}

  

  SolverControl &
  SolverBase::control() const
  {
    return solver_control;
  }



  void
  SolverBase::solve (const SparseMatrix     &A,
                     VectorBase             &x,
                     const VectorBase       &b,
                     const PreconditionBase &preconditioner)
  {
    int ierr;
    
    linear_problem.reset();

					// We need an
					// Epetra_LinearProblem object
					// to let the AztecOO solver
					// know about the matrix and
					// vectors.
    linear_problem = std::auto_ptr<Epetra_LinearProblem> (
			      new Epetra_LinearProblem(&*(A.matrix), &*x.vector,
						       &*b.vector));

					// Next we can allocate the
					// AztecOO solver...
    solver.SetProblem(*linear_problem);

					// ... and we can specify the
					// solver to be used.
    switch (solver_name)
      {
      case cg:
	solver.SetAztecOption(AZ_solver, AZ_cg);
	break;
      case cgs:
	solver.SetAztecOption(AZ_solver, AZ_cgs);
	break;
      case gmres:
	solver.SetAztecOption(AZ_solver, AZ_gmres);
	solver.SetAztecOption(AZ_kspace, additional_data.gmres_restart_parameter);
	break;
      case bicgstab:
	solver.SetAztecOption(AZ_solver, AZ_bicgstab);
	break;
      case tfqmr:
	solver.SetAztecOption(AZ_solver, AZ_tfqmr);
	break;
      default:
	Assert (false, ExcNotImplemented());
      }

					// Introduce the
					// preconditioner, ...
    ierr = solver.SetPrecOperator (const_cast<Epetra_Operator*>
				     (&*preconditioner.preconditioner));
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

					// ... set some options, ...
    solver.SetAztecOption (AZ_output, AZ_none);
    solver.SetAztecOption (AZ_conv, AZ_noscaled);

					// ... and then solve!
    ierr = solver.Iterate (solver_control.max_steps(), 
			   solver_control.tolerance());
    AssertThrow (ierr >= 0, ExcTrilinosError(ierr));

					// Finally, let the deal.II
					// SolverControl object know
					// what has happened. If the
					// solve succeeded, the status
					// of the solver control will
					// turn into
					// SolverControl::success.
    solver_control.check (solver.NumIters(), solver.TrueResidual());

    if (solver_control.last_check() != SolverControl::success)
      throw SolverControl::NoConvergence (solver_control.last_step(),
                                          solver_control.last_value());
  }
  
    


  
/* ---------------------- SolverCG ------------------------ */

  SolverCG::SolverCG (SolverControl        &cn,
                      const AdditionalData &data)
                  :
                  SolverBase (cn),
                  additional_data (data)
  {
    solver_name = cg;
  }
  

/* ---------------------- SolverGMRES ------------------------ */

  SolverGMRES::AdditionalData::
  AdditionalData (const unsigned int restart_parameter)
                  :
                  restart_parameter (restart_parameter)
  {}

  
  
  SolverGMRES::SolverGMRES (SolverControl        &cn,
                            const AdditionalData &data)
                  :
                  SolverBase (cn),
                  additional_data (data.restart_parameter)
  {
    solver_name = gmres;
  }
  

/* ---------------------- SolverBicgstab ------------------------ */

  SolverBicgstab::SolverBicgstab (SolverControl        &cn,
                                  const AdditionalData &data)
                  :
                  SolverBase (cn),
                  additional_data (data)
  {
    solver_name = bicgstab;
  }

  
/* ---------------------- SolverCGS ------------------------ */

  SolverCGS::SolverCGS (SolverControl        &cn,
                        const AdditionalData &data)
                  :
                  SolverBase (cn),
                  additional_data (data)
  {
    solver_name = cgs;
  }
  

/* ---------------------- SolverTFQMR ------------------------ */

  SolverTFQMR::SolverTFQMR (SolverControl        &cn,
                            const AdditionalData &data)
                  :
                  SolverBase (cn),
                  additional_data (data)
  {
    solver_name = tfqmr;
  }

}

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_USE_PETSC
