//----------------------------  petsc_solver.h  ---------------------------
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
//----------------------------  petsc_solver.h  ---------------------------
#ifndef __deal2__petsc_solver_h
#define __deal2__petsc_solver_h

#include <base/config.h>
#include <base/exceptions.h>
#include <lac/solver_control.h>

#ifdef DEAL_II_USE_PETSC

#include <petscksp.h>


/*! @addtogroup PETSc
 *@{
 */

namespace PETScWrappers
{
                                   // forward declarations
  class MatrixBase;
  class VectorBase;
  class PreconditionerBase;
  
  
/**
 * Base class for solver classes using the PETSc solvers. Since solvers in
 * PETSc are selected based on flags passed to a generic solver object,
 * basically all the actual solver calls happen in this class, and derived
 * classes simply set the right flags to select one solver or another, or to
 * set certain parameters for individual solvers.
 *
 * @author Wolfgang Bangerth, 2004
 */
  class SolverBase
  {
    public:
                                       /**
                                        * Constructor. Takes the solver
                                        * control object and the MPI
                                        * communicator over which parallel
                                        * computations are to happen.
                                        *
                                        * Note that the communicator used here
                                        * must match the communicator used in
                                        * the system matrix, solution, and
                                        * right hand side object of the solve
                                        * to be done with this
                                        * solver. Otherwise, PETSc will
                                        * generate hard to track down errors.
                                        */
      SolverBase (SolverControl &cn,
                  MPI_Comm      &mpi_communicator);
      
                                       /**
                                        * Destructor.
                                        */
      virtual ~SolverBase ();

                                       /**
                                        * Solve the linear system
                                        * <tt>Ax=b</tt>. Depending on the
                                        * information provided by derived
                                        * classes and the object passed as a
                                        * preconditioner, one of the linear
                                        * solvers and preconditioners of PETSc
                                        * is chosen.
                                        */
      void
      solve (const MatrixBase         &A,
             VectorBase               &x,
             const VectorBase         &b,
             const PreconditionerBase &preconditioner) const;


                                       /**
                                        * Access to object that controls
                                        * convergence.
                                        */
      SolverControl & control() const;

                                       /**
                                        * Exception
                                        */
      DeclException1 (ExcPETScError,
                      int,
                      << "An error with error number " << arg1
                      << " occured while calling a PETSc function");

    protected:

                                       /**
                                        * Reference to the object that
                                        * controls convergence of the
                                        * iterative solver. In fact, for these
                                        * PETSc wrappers, PETSc does so
                                        * itself, but we copy the data from
                                        * this object before starting the
                                        * solution process, and copy the data
                                        * back into it afterwards.
                                        */
      SolverControl &solver_control;

                                       /**
                                        * Copy of the MPI communicator object
                                        * to be used for the solver.
                                        */
      MPI_Comm mpi_communicator;

                                       /**
                                        * Function that takes a Krylov
                                        * Subspace Solver context object, and
                                        * sets the type of solver that is
                                        * requested by the derived class.
                                        */
      virtual void set_solver_type (KSP &ksp) const = 0;

    private:
                                       /**
                                        * A function that is used in PETSc as
                                        * a callback to check on
                                        * convergence. It takes the
                                        * information provided from PETSc and
                                        * checks it against deal.II's own
                                        * SolverControl objects to see if
                                        * convergence has been reached.
                                        */
      static
      int
      convergence_test (KSP                 ksp,
                        const int           iteration,
                        const double        residual_norm,
                        KSPConvergedReason *reason,
                        void               *solver_control);
  };



/**
 * An implementation of the solver interface using the PETSc Richardson
 * solver.
 *
 * @author Wolfgang Bangerth, 2004
 */
  class SolverRichardson : public SolverBase
  {
    public:
                                       /**
                                        * Standardized data struct to
                                        * pipe additional data to the
                                        * solver.
                                        */
      struct AdditionalData 
      {
                                           /**
                                            * Constructor. By default,
                                            * set the damping parameter
                                            * to one.
                                            */
          AdditionalData (const double omega = 1);
	
                                           /**
                                            * Relaxation parameter.
                                            */
          double omega;
      };
	
                                       /**
                                        * Constructor. In contrast to
                                        * deal.II's own solvers, there is no
                                        * need to give a vector memory
                                        * object. However, PETSc solvers want
                                        * to have an MPI communicator context
                                        * over which computations are
                                        * parallelized. By default,
                                        * @p PETSC_COMM_SELF is used here,
                                        * but you can change this. Note that
                                        * for single processor (non-MPI)
                                        * versions, this parameter does not
                                        * have any effect.
                                        *
                                        * The last argument takes a structure
                                        * with additional, solver dependent
                                        * flags for tuning. We specify the
                                        * (default) value to the constructor
                                        * call in this default argument
                                        * because otherwise gcc 2.95 generates
                                        * a compiler fault.
                                        *
                                        * Note that the communicator used here
                                        * must match the communicator used in
                                        * the system matrix, solution, and
                                        * right hand side object of the solve
                                        * to be done with this
                                        * solver. Otherwise, PETSc will
                                        * generate hard to track down errors.
                                        */
      SolverRichardson (SolverControl        &cn,
                        MPI_Comm             &mpi_communicator = PETSC_COMM_SELF,
                        const AdditionalData &data = AdditionalData(1));

    protected:
                                       /**
                                        * Store a copy of the flags for this
                                        * particular solver.
                                        */
      const AdditionalData additional_data;

                                       /**
                                        * Function that takes a Krylov
                                        * Subspace Solver context object, and
                                        * sets the type of solver that is
                                        *appropriate for this class.
                                        */
      virtual void set_solver_type (KSP &ksp) const;
  };



/**
 * An implementation of the solver interface using the PETSc Chebychev
 * solver.
 *
 * @author Wolfgang Bangerth, 2004
 */
  class SolverChebychev : public SolverBase
  {
    public:
                                       /**
                                        * Standardized data struct to
                                        * pipe additional data to the
                                        * solver.
                                        */
      struct AdditionalData
      {};
	
                                       /**
                                        * Constructor. In contrast to
                                        * deal.II's own solvers, there is no
                                        * need to give a vector memory
                                        * object. However, PETSc solvers want
                                        * to have an MPI communicator context
                                        * over which computations are
                                        * parallelized. By default,
                                        * @p PETSC_COMM_SELF is used here,
                                        * but you can change this. Note that
                                        * for single processor (non-MPI)
                                        * versions, this parameter does not
                                        * have any effect.
                                        *
                                        * The last argument takes a structure
                                        * with additional, solver dependent
                                        * flags for tuning.
                                        *
                                        * Note that the communicator used here
                                        * must match the communicator used in
                                        * the system matrix, solution, and
                                        * right hand side object of the solve
                                        * to be done with this
                                        * solver. Otherwise, PETSc will
                                        * generate hard to track down errors.
                                        */
      SolverChebychev (SolverControl        &cn,
                       MPI_Comm             &mpi_communicator = PETSC_COMM_SELF,
                       const AdditionalData &data = AdditionalData());

    protected:
                                       /**
                                        * Store a copy of the flags for this
                                        * particular solver.
                                        */
      const AdditionalData additional_data;

                                       /**
                                        * Function that takes a Krylov
                                        * Subspace Solver context object, and
                                        * sets the type of solver that is
                                        *appropriate for this class.
                                        */
      virtual void set_solver_type (KSP &ksp) const;
  };



/**
 * An implementation of the solver interface using the PETSc CG
 * solver.
 *
 * @author Wolfgang Bangerth, 2004
 */
  class SolverCG : public SolverBase
  {
    public:
                                       /**
                                        * Standardized data struct to
                                        * pipe additional data to the
                                        * solver.
                                        */
      struct AdditionalData
      {};
	
                                       /**
                                        * Constructor. In contrast to
                                        * deal.II's own solvers, there is no
                                        * need to give a vector memory
                                        * object. However, PETSc solvers want
                                        * to have an MPI communicator context
                                        * over which computations are
                                        * parallelized. By default,
                                        * @p PETSC_COMM_SELF is used here,
                                        * but you can change this. Note that
                                        * for single processor (non-MPI)
                                        * versions, this parameter does not
                                        * have any effect.
                                        *
                                        * The last argument takes a structure
                                        * with additional, solver dependent
                                        * flags for tuning.
                                        *
                                        * Note that the communicator used here
                                        * must match the communicator used in
                                        * the system matrix, solution, and
                                        * right hand side object of the solve
                                        * to be done with this
                                        * solver. Otherwise, PETSc will
                                        * generate hard to track down errors.
                                        */
      SolverCG (SolverControl        &cn,
                MPI_Comm             &mpi_communicator = PETSC_COMM_SELF,
                const AdditionalData &data = AdditionalData());

    protected:
                                       /**
                                        * Store a copy of the flags for this
                                        * particular solver.
                                        */
      const AdditionalData additional_data;

                                       /**
                                        * Function that takes a Krylov
                                        * Subspace Solver context object, and
                                        * sets the type of solver that is
                                        *appropriate for this class.
                                        */
      virtual void set_solver_type (KSP &ksp) const;
  };



/**
 * An implementation of the solver interface using the PETSc BiCG
 * solver.
 *
 * @author Wolfgang Bangerth, 2004
 */
  class SolverBiCG : public SolverBase
  {
    public:
                                       /**
                                        * Standardized data struct to
                                        * pipe additional data to the
                                        * solver.
                                        */
      struct AdditionalData
      {};
	
                                       /**
                                        * Constructor. In contrast to
                                        * deal.II's own solvers, there is no
                                        * need to give a vector memory
                                        * object. However, PETSc solvers want
                                        * to have an MPI communicator context
                                        * over which computations are
                                        * parallelized. By default,
                                        * @p PETSC_COMM_SELF is used here,
                                        * but you can change this. Note that
                                        * for single processor (non-MPI)
                                        * versions, this parameter does not
                                        * have any effect.
                                        *
                                        * The last argument takes a structure
                                        * with additional, solver dependent
                                        * flags for tuning.
                                        *
                                        * Note that the communicator used here
                                        * must match the communicator used in
                                        * the system matrix, solution, and
                                        * right hand side object of the solve
                                        * to be done with this
                                        * solver. Otherwise, PETSc will
                                        * generate hard to track down errors.
                                        */
      SolverBiCG (SolverControl        &cn,
                  MPI_Comm             &mpi_communicator = PETSC_COMM_SELF,
                  const AdditionalData &data = AdditionalData());

    protected:
                                       /**
                                        * Store a copy of the flags for this
                                        * particular solver.
                                        */
      const AdditionalData additional_data;

                                       /**
                                        * Function that takes a Krylov
                                        * Subspace Solver context object, and
                                        * sets the type of solver that is
                                        *appropriate for this class.
                                        */
      virtual void set_solver_type (KSP &ksp) const;
  };



/**
 * An implementation of the solver interface using the PETSc GMRES
 * solver.
 *
 * @author Wolfgang Bangerth, 2004
 */
  class SolverGMRES : public SolverBase
  {
    public:
                                       /**
                                        * Standardized data struct to
                                        * pipe additional data to the
                                        * solver.
                                        */
      struct AdditionalData
      {
                                           /**
                                            * Constructor. By default, set the
                                            * number of temporary vectors to
                                            * 30, i.e. do a restart every 30
                                            * iterations.
                                            */
          AdditionalData (const unsigned int restart_parameter = 30);
	
                                           /**
                                            * Maximum number of
                                            * tmp vectors.
                                            */
          unsigned int restart_parameter;
      };
	
                                       /**
                                        * Constructor. In contrast to
                                        * deal.II's own solvers, there is no
                                        * need to give a vector memory
                                        * object. However, PETSc solvers want
                                        * to have an MPI communicator context
                                        * over which computations are
                                        * parallelized. By default,
                                        * @p PETSC_COMM_SELF is used here,
                                        * but you can change this. Note that
                                        * for single processor (non-MPI)
                                        * versions, this parameter does not
                                        * have any effect.
                                        *
                                        * The last argument takes a structure
                                        * with additional, solver dependent
                                        * flags for tuning. We specify the
                                        * (default) value to the constructor
                                        * call in this default argument
                                        * because otherwise gcc 2.95 generates
                                        * a compiler fault.
                                        *
                                        * Note that the communicator used here
                                        * must match the communicator used in
                                        * the system matrix, solution, and
                                        * right hand side object of the solve
                                        * to be done with this
                                        * solver. Otherwise, PETSc will
                                        * generate hard to track down errors.
                                        */
      SolverGMRES (SolverControl        &cn,
                   MPI_Comm             &mpi_communicator = PETSC_COMM_SELF,
                   const AdditionalData &data = AdditionalData(30));

    protected:
                                       /**
                                        * Store a copy of the flags for this
                                        * particular solver.
                                        */
      const AdditionalData additional_data;

                                       /**
                                        * Function that takes a Krylov
                                        * Subspace Solver context object, and
                                        * sets the type of solver that is
                                        *appropriate for this class.
                                        */
      virtual void set_solver_type (KSP &ksp) const;
  };



/**
 * An implementation of the solver interface using the PETSc BiCGStab
 * solver.
 *
 * @author Wolfgang Bangerth, 2004
 */
  class SolverBicgstab : public SolverBase
  {
    public:
                                       /**
                                        * Standardized data struct to
                                        * pipe additional data to the
                                        * solver.
                                        */
      struct AdditionalData
      {};
	
                                       /**
                                        * Constructor. In contrast to
                                        * deal.II's own solvers, there is no
                                        * need to give a vector memory
                                        * object. However, PETSc solvers want
                                        * to have an MPI communicator context
                                        * over which computations are
                                        * parallelized. By default,
                                        * @p PETSC_COMM_SELF is used here,
                                        * but you can change this. Note that
                                        * for single processor (non-MPI)
                                        * versions, this parameter does not
                                        * have any effect.
                                        *
                                        * The last argument takes a structure
                                        * with additional, solver dependent
                                        * flags for tuning.
                                        *
                                        * Note that the communicator used here
                                        * must match the communicator used in
                                        * the system matrix, solution, and
                                        * right hand side object of the solve
                                        * to be done with this
                                        * solver. Otherwise, PETSc will
                                        * generate hard to track down errors.
                                        */
      SolverBicgstab (SolverControl        &cn,
                      MPI_Comm             &mpi_communicator = PETSC_COMM_SELF,
                      const AdditionalData &data = AdditionalData());

    protected:
                                       /**
                                        * Store a copy of the flags for this
                                        * particular solver.
                                        */
      const AdditionalData additional_data;

                                       /**
                                        * Function that takes a Krylov
                                        * Subspace Solver context object, and
                                        * sets the type of solver that is
                                        *appropriate for this class.
                                        */
      virtual void set_solver_type (KSP &ksp) const;
  };



/**
 * An implementation of the solver interface using the PETSc CG Squared
 * solver.
 *
 * @author Wolfgang Bangerth, 2004
 */
  class SolverCGS : public SolverBase
  {
    public:
                                       /**
                                        * Standardized data struct to
                                        * pipe additional data to the
                                        * solver.
                                        */
      struct AdditionalData
      {};
	
                                       /**
                                        * Constructor. In contrast to
                                        * deal.II's own solvers, there is no
                                        * need to give a vector memory
                                        * object. However, PETSc solvers want
                                        * to have an MPI communicator context
                                        * over which computations are
                                        * parallelized. By default,
                                        * @p PETSC_COMM_SELF is used here,
                                        * but you can change this. Note that
                                        * for single processor (non-MPI)
                                        * versions, this parameter does not
                                        * have any effect.
                                        *
                                        * The last argument takes a structure
                                        * with additional, solver dependent
                                        * flags for tuning.
                                        *
                                        * Note that the communicator used here
                                        * must match the communicator used in
                                        * the system matrix, solution, and
                                        * right hand side object of the solve
                                        * to be done with this
                                        * solver. Otherwise, PETSc will
                                        * generate hard to track down errors.
                                        */
      SolverCGS (SolverControl        &cn,
                 MPI_Comm             &mpi_communicator = PETSC_COMM_SELF,
                 const AdditionalData &data = AdditionalData());

    protected:
                                       /**
                                        * Store a copy of the flags for this
                                        * particular solver.
                                        */
      const AdditionalData additional_data;

                                       /**
                                        * Function that takes a Krylov
                                        * Subspace Solver context object, and
                                        * sets the type of solver that is
                                        *appropriate for this class.
                                        */
      virtual void set_solver_type (KSP &ksp) const;
  };



/**
 * An implementation of the solver interface using the PETSc TFQMR
 * solver.
 *
 * @author Wolfgang Bangerth, 2004
 */
  class SolverTFQMR : public SolverBase
  {
    public:
                                       /**
                                        * Standardized data struct to
                                        * pipe additional data to the
                                        * solver.
                                        */
      struct AdditionalData
      {};
	
                                       /**
                                        * Constructor. In contrast to
                                        * deal.II's own solvers, there is no
                                        * need to give a vector memory
                                        * object. However, PETSc solvers want
                                        * to have an MPI communicator context
                                        * over which computations are
                                        * parallelized. By default,
                                        * @p PETSC_COMM_SELF is used here,
                                        * but you can change this. Note that
                                        * for single processor (non-MPI)
                                        * versions, this parameter does not
                                        * have any effect.
                                        *
                                        * The last argument takes a structure
                                        * with additional, solver dependent
                                        * flags for tuning.
                                        *
                                        * Note that the communicator used here
                                        * must match the communicator used in
                                        * the system matrix, solution, and
                                        * right hand side object of the solve
                                        * to be done with this
                                        * solver. Otherwise, PETSc will
                                        * generate hard to track down errors.
                                        */
      SolverTFQMR (SolverControl        &cn,
                   MPI_Comm             &mpi_communicator = PETSC_COMM_SELF,
                   const AdditionalData &data = AdditionalData());

    protected:
                                       /**
                                        * Store a copy of the flags for this
                                        * particular solver.
                                        */
      const AdditionalData additional_data;

                                       /**
                                        * Function that takes a Krylov
                                        * Subspace Solver context object, and
                                        * sets the type of solver that is
                                        *appropriate for this class.
                                        */
      virtual void set_solver_type (KSP &ksp) const;
  };




/**
 * An implementation of the solver interface using the PETSc TFQMR-2 solver
 * (called TCQMR in PETSc). Note that this solver had a serious bug in
 * versions up to and including PETSc 2.1.6, in that it did not check
 * convergence and always returned an error code. Thus, this class will abort
 * with an error indicating failure to converge with PETSc 2.1.6 and
 * prior. This should be fixed in later versions of PETSc, though.
 *
 * @author Wolfgang Bangerth, 2004
 */
  class SolverTCQMR : public SolverBase
  {
    public:
                                       /**
                                        * Standardized data struct to
                                        * pipe additional data to the
                                        * solver.
                                        */
      struct AdditionalData
      {};
	
                                       /**
                                        * Constructor. In contrast to
                                        * deal.II's own solvers, there is no
                                        * need to give a vector memory
                                        * object. However, PETSc solvers want
                                        * to have an MPI communicator context
                                        * over which computations are
                                        * parallelized. By default,
                                        * @p PETSC_COMM_SELF is used here,
                                        * but you can change this. Note that
                                        * for single processor (non-MPI)
                                        * versions, this parameter does not
                                        * have any effect.
                                        *
                                        * The last argument takes a structure
                                        * with additional, solver dependent
                                        * flags for tuning.
                                        *
                                        * Note that the communicator used here
                                        * must match the communicator used in
                                        * the system matrix, solution, and
                                        * right hand side object of the solve
                                        * to be done with this
                                        * solver. Otherwise, PETSc will
                                        * generate hard to track down errors.
                                        */
      SolverTCQMR (SolverControl        &cn,
                   MPI_Comm             &mpi_communicator = PETSC_COMM_SELF,
                   const AdditionalData &data = AdditionalData());

    protected:
                                       /**
                                        * Store a copy of the flags for this
                                        * particular solver.
                                        */
      const AdditionalData additional_data;

                                       /**
                                        * Function that takes a Krylov
                                        * Subspace Solver context object, and
                                        * sets the type of solver that is
                                        *appropriate for this class.
                                        */
      virtual void set_solver_type (KSP &ksp) const;
  };



/**
 * An implementation of the solver interface using the PETSc CR
 * solver.
 *
 * @author Wolfgang Bangerth, 2004
 */
  class SolverCR : public SolverBase
  {
    public:
                                       /**
                                        * Standardized data struct to
                                        * pipe additional data to the
                                        * solver.
                                        */
      struct AdditionalData
      {};
	
                                       /**
                                        * Constructor. In contrast to
                                        * deal.II's own solvers, there is no
                                        * need to give a vector memory
                                        * object. However, PETSc solvers want
                                        * to have an MPI communicator context
                                        * over which computations are
                                        * parallelized. By default,
                                        * @p PETSC_COMM_SELF is used here,
                                        * but you can change this. Note that
                                        * for single processor (non-MPI)
                                        * versions, this parameter does not
                                        * have any effect.
                                        *
                                        * The last argument takes a structure
                                        * with additional, solver dependent
                                        * flags for tuning.
                                        *
                                        * Note that the communicator used here
                                        * must match the communicator used in
                                        * the system matrix, solution, and
                                        * right hand side object of the solve
                                        * to be done with this
                                        * solver. Otherwise, PETSc will
                                        * generate hard to track down errors.
                                        */
      SolverCR (SolverControl        &cn,
                MPI_Comm             &mpi_communicator = PETSC_COMM_SELF,
                const AdditionalData &data = AdditionalData());

    protected:
                                       /**
                                        * Store a copy of the flags for this
                                        * particular solver.
                                        */
      const AdditionalData additional_data;

                                       /**
                                        * Function that takes a Krylov
                                        * Subspace Solver context object, and
                                        * sets the type of solver that is
                                        *appropriate for this class.
                                        */
      virtual void set_solver_type (KSP &ksp) const;
  };



/**
 * An implementation of the solver interface using the PETSc Least Squares
 * solver.
 *
 * @author Wolfgang Bangerth, 2004
 */
  class SolverLSQR : public SolverBase
  {
    public:
                                       /**
                                        * Standardized data struct to
                                        * pipe additional data to the
                                        * solver.
                                        */
      struct AdditionalData
      {};
	
                                       /**
                                        * Constructor. In contrast to
                                        * deal.II's own solvers, there is no
                                        * need to give a vector memory
                                        * object. However, PETSc solvers want
                                        * to have an MPI communicator context
                                        * over which computations are
                                        * parallelized. By default,
                                        * @p PETSC_COMM_SELF is used here,
                                        * but you can change this. Note that
                                        * for single processor (non-MPI)
                                        * versions, this parameter does not
                                        * have any effect.
                                        *
                                        * The last argument takes a structure
                                        * with additional, solver dependent
                                        * flags for tuning.
                                        *
                                        * Note that the communicator used here
                                        * must match the communicator used in
                                        * the system matrix, solution, and
                                        * right hand side object of the solve
                                        * to be done with this
                                        * solver. Otherwise, PETSc will
                                        * generate hard to track down errors.
                                        */
      SolverLSQR (SolverControl        &cn,
                  MPI_Comm             &mpi_communicator = PETSC_COMM_SELF,
                  const AdditionalData &data = AdditionalData());

    protected:
                                       /**
                                        * Store a copy of the flags for this
                                        * particular solver.
                                        */
      const AdditionalData additional_data;

                                       /**
                                        * Function that takes a Krylov
                                        * Subspace Solver context object, and
                                        * sets the type of solver that is
                                        *appropriate for this class.
                                        */
      virtual void set_solver_type (KSP &ksp) const;
  };
  
}

/*@}*/

#endif // DEAL_II_USE_PETSC

/*----------------------------   petsc_solver.h     ---------------------------*/

#endif
/*----------------------------   petsc_solver.h     ---------------------------*/
