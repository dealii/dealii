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
#include <lac/exceptions.h>
#include <lac/solver_control.h>

#ifdef DEAL_II_USE_PETSC

#include <boost/shared_ptr.hpp>
#include <petscksp.h>


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
 * One of the gotchas of PETSc is that -- in particular in MPI mode -- it
 * often does not produce very helpful error messages. In order to save
 * other users some time in searching a hard to track down error, here is
 * one situation and the error message one gets there:
 * when you don't specify an MPI communicator to your solver's constructor. In
 * this case, you will get an error of the following form from each of your
 * parallel processes:
 * @verbatim
 *   [1]PETSC ERROR: PCSetVector() line 1173 in src/ksp/pc/interface/precon.c
 *   [1]PETSC ERROR:   Arguments must have same communicators!
 *   [1]PETSC ERROR:   Different communicators in the two objects: Argument # 1 and 2!
 *   [1]PETSC ERROR: KSPSetUp() line 195 in src/ksp/ksp/interface/itfunc.c
 * @endverbatim
 *
 * This error, on which one can spend a very long time figuring out what
 * exactly goes wrong, results from not specifying an MPI communicator. Note
 * that the communicator @em must match that of the matrix and all vectors in
 * the linear system which we want to solve. Aggravating the situation is the
 * fact that the default argument to the solver classes, @p PETSC_COMM_SELF,
 * is the appropriate argument for the sequential case (which is why it is the
 * default argument), so this error only shows up in parallel mode.
 * 
 * @ingroup PETScWrappers
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
                                        * generate hard to track down errors,
                                        * see the documentation of the
                                        * SolverBase class.
                                        */
      SolverBase (SolverControl  &cn,
                  const MPI_Comm &mpi_communicator);
      
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
             const PreconditionerBase &preconditioner);


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
      const MPI_Comm mpi_communicator;

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
                        const PetscScalar   residual_norm,
                        KSPConvergedReason *reason,
                        void               *solver_control);

                                       /**
                                        * A structure that contains the PETSc
                                        * solver and preconditioner
                                        * objects. This object is preserved
                                        * between subsequent calls to the
                                        * solver if the same preconditioner is
                                        * used as in the previous solver
                                        * step. This may save some computation
                                        * time, if setting up a preconditioner
                                        * is expensive, such as in the case of
                                        * an ILU for example.
                                        *
                                        * The actual declaration of this class
                                        * is complicated by the fact that
                                        * PETSc changed its solver interface
                                        * completely and incompatibly between
                                        * versions 2.1.6 and 2.2.0 :-(
                                        *
                                        * Objects of this type are explicitly
                                        * created, but are destroyed when the
                                        * surrounding solver object goes out
                                        * of scope, or when we assign a new
                                        * value to the pointer to this
                                        * object. The respective *Destroy
                                        * functions are therefore written into
                                        * the destructor of this object, even
                                        * though the object does not have a
                                        * constructor.
                                        */
      struct SolverData
      {
                                           /**
                                            * Destructor
                                            */
          ~SolverData ();
          
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR < 2)
                                           /**
                                            * A PETSc solver object.
                                            */
          SLES sles;
#endif
                                           /**
                                            * Objects for Krylov subspace
                                            * solvers and preconditioners.
                                            */
          KSP  ksp;
          PC   pc;
      };

                                       /**
                                        * Pointer to an object that stores the
                                        * solver context. This is recreated in
                                        * the main solver routine if
                                        * necessary.
                                        */
      boost::shared_ptr<SolverData> solver_data;
  };



/**
 * An implementation of the solver interface using the PETSc Richardson
 * solver.
 *
 * @ingroup PETScWrappers
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
                                        * generate hard to track down errors,
                                        * see the documentation of the
                                        * SolverBase class.
                                        */
      SolverRichardson (SolverControl        &cn,
                        const MPI_Comm       &mpi_communicator = PETSC_COMM_SELF,
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
 * @ingroup PETScWrappers
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
                                        * generate hard to track down errors,
                                        * see the documentation of the
                                        * SolverBase class.
                                        */
      SolverChebychev (SolverControl        &cn,
                       const MPI_Comm       &mpi_communicator = PETSC_COMM_SELF,
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
 * @ingroup PETScWrappers
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
                                        * generate hard to track down errors,
                                        * see the documentation of the
                                        * SolverBase class.
                                        */
      SolverCG (SolverControl        &cn,
                const MPI_Comm       &mpi_communicator = PETSC_COMM_SELF,
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
 * @ingroup PETScWrappers
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
                                        * generate hard to track down errors,
                                        * see the documentation of the
                                        * SolverBase class.
                                        */
      SolverBiCG (SolverControl        &cn,
                  const MPI_Comm       &mpi_communicator = PETSC_COMM_SELF,
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
          AdditionalData (const unsigned int restart_parameter = 30,
			  const bool right_preconditioning = false);
	
                                           /**
                                            * Maximum number of
                                            * tmp vectors.
                                            */
          unsigned int restart_parameter;

					   /**
					    * Flag for right
					    * preconditioning.
					    */
	  bool right_preconditioning;
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
                                        * generate hard to track down errors,
                                        * see the documentation of the
                                        * SolverBase class.
                                        */
      SolverGMRES (SolverControl        &cn,
                   const MPI_Comm       &mpi_communicator = PETSC_COMM_SELF,
                   const AdditionalData &data = AdditionalData(30,false));

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
 * @ingroup PETScWrappers
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
                                        * generate hard to track down errors,
                                        * see the documentation of the
                                        * SolverBase class.
                                        */
      SolverBicgstab (SolverControl        &cn,
                      const MPI_Comm       &mpi_communicator = PETSC_COMM_SELF,
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
 * @ingroup PETScWrappers
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
                                        * generate hard to track down errors,
                                        * see the documentation of the
                                        * SolverBase class.
                                        */
      SolverCGS (SolverControl        &cn,
                 const MPI_Comm       &mpi_communicator = PETSC_COMM_SELF,
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
 * @ingroup PETScWrappers
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
                                        * generate hard to track down errors,
                                        * see the documentation of the
                                        * SolverBase class.
                                        */
      SolverTFQMR (SolverControl        &cn,
                   const MPI_Comm       &mpi_communicator = PETSC_COMM_SELF,
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
 * @ingroup PETScWrappers
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
                                        * generate hard to track down errors,
                                        * see the documentation of the
                                        * SolverBase class.
                                        */
      SolverTCQMR (SolverControl        &cn,
                   const MPI_Comm       &mpi_communicator = PETSC_COMM_SELF,
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
 * @ingroup PETScWrappers
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
                                        * generate hard to track down errors,
                                        * see the documentation of the
                                        * SolverBase class.
                                        */
      SolverCR (SolverControl        &cn,
                const MPI_Comm       &mpi_communicator = PETSC_COMM_SELF,
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
 * @ingroup PETScWrappers
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
                                        * generate hard to track down errors,
                                        * see the documentation of the
                                        * SolverBase class.
                                        */
      SolverLSQR (SolverControl        &cn,
                  const MPI_Comm       &mpi_communicator = PETSC_COMM_SELF,
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
 * An implementation of the solver interface using the PETSc PREONLY
 * solver. Actually this is NOT a real solution algorithm. Its only
 * purpose is to provide a solver object, when the preconditioner
 * should be used as real solver. It is very useful in conjunction with
 * the complete LU decomposition preconditioner <tt> PreconditionLU </tt>, 
 * which in conjunction with this solver class becomes a direct solver.
 *
 * @ingroup PETScWrappers
 * @author Wolfgang Bangerth, 2004, Oliver Kayser-Herold, 2004
 */
  class SolverPreOnly : public SolverBase
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
                                        * generate hard to track down errors,
                                        * see the documentation of the
                                        * SolverBase class.
                                        */
      SolverPreOnly (SolverControl        &cn,
		     const MPI_Comm       &mpi_communicator = PETSC_COMM_SELF,
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
                                        * appropriate for this class.
                                        */
      virtual void set_solver_type (KSP &ksp) const;
  };

}

#endif // DEAL_II_USE_PETSC

/*----------------------------   petsc_solver.h     ---------------------------*/

#endif
/*----------------------------   petsc_solver.h     ---------------------------*/
