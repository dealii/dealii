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
#ifndef __deal2__trilinos_block_solver_h
#define __deal2__trilinos_block_solver_h


#include <base/config.h>

#ifdef DEAL_II_USE_TRILINOS

#  include <lac/exceptions.h>
#  include <lac/solver_control.h>

DEAL_II_NAMESPACE_OPEN


namespace TrilinosWrappers
{
                                   // forward declarations
  class BlockSparseMatrix;
  class BlockVector;
  namespace MPI
  {
    class BlockVector;
  }
  class SparseMatrix;
  class PreconditionBase;
  class PreconditionBlockBase;
  
  
/**
 * Base class for solver classes using the Trilinos solvers on block
 * matrices based on the Trilinos abstract solver interface Thyra.
 *
 * @ingroup TrilinosWrappers
 * @author Martin Kronbichler, 2008
 */
  class SolverBlockBase
  {
    public:

				       /**
					* Enumeration object that is
					* set in the constructor of
					* the derived classes and
					* tells Trilinos which solver
					* to use. This option can also
					* be set in the user program,
					* so one might use this base
					* class instead of one of the
					* specialized derived classes
					* when the solver should be
					* set at runtime. Currently
					* enabled options are:
					*/
      enum SolverBlockName {cg, gmres} solver_name;

                                       /**
                                        * Standardized data struct to
                                        * pipe additional data to the
                                        * solver.
                                        */
      struct AdditionalData
      {
				       /**
					* Restart parameter in case
					* the derived class is
					* GMRES.
					* 
					* TODO: Find a better way for
					* doing this.
					*/
	AdditionalData (const unsigned int gmres_restart_parameter = 30,
			const bool         right_preconditioning   = false,
			const bool         output_details          = false);

				       /**
					* The number of vectors that
					* should be used in the GMRES
					* algorithm before the solver
					* is restarted.
					*/
	unsigned int gmres_restart_parameter;

				       /**
					* Flag that turns right
					* preconditioning on.
					*/
	bool right_preconditioning;

				       /**
					* A flag to determine whether
					* solver details should be
					* written to screen.
					*/
	bool output_details;
      };

                                       /**
                                        * Constructor. Takes the
                                        * solver control object and
                                        * creates the solver.
                                        */
      SolverBlockBase (SolverControl  &cn);

                                       /**
                                        * Second constructor. This
                                        * constructor takes an enum
                                        * object that specifies the
                                        * solver name and sets the
                                        * appropriate Krylov
                                        * method.
                                        */
      SolverBlockBase (const enum SolverBlockName  solver_name,
		       SolverControl              &cn);

                                       /**
                                        * Destructor.
                                        */
      virtual ~SolverBlockBase ();

                                       /**
                                        * Solve the linear system
                                        * <tt>Ax=b</tt>. Depending on
                                        * the information provided by
                                        * derived classes and the
                                        * object passed as a
                                        * preconditioner, one of the
                                        * linear solvers and
                                        * preconditioners of Trilinos
                                        * is chosen. This interface is
                                        * designed for Trilinos
                                        * running in parallel.
                                        */
      void
      solve (const BlockSparseMatrix     &A,
             MPI::BlockVector            &x,
             const MPI::BlockVector      &b,
             const PreconditionBlockBase &preconditioner);

                                       /**
                                        * Solve the linear system
                                        * <tt>Ax=b</tt>. Depending on
                                        * the information provided by
                                        * derived classes and the
                                        * object passed as a
                                        * preconditioner, one of the
                                        * linear solvers and
                                        * preconditioners of Trilinos
                                        * is chosen. This interface is
                                        * designed for Trilinos
                                        * running in serial.
                                        */
      void
      solve (const BlockSparseMatrix     &A,
             BlockVector                 &x,
             const BlockVector           &b,
             const PreconditionBlockBase &preconditioner);

                                       /**
                                        * Access to object that controls
                                        * convergence.
                                        */
      SolverControl & control() const;

                                       /**
                                        * Exception
                                        */
      DeclException1 (ExcTrilinosError,
                      int,
                      << "An error with error number " << arg1
                      << " occured while calling a Trilinos function");

                                       /**
                                        * Exception
                                        */
      DeclException2 (ExcOverlappingMaps,
                      std::string, std::string,
                      << "The Epetra map in the " << arg1 << " " << arg2
                      << " is overlapping between the processors, which"
		      << " is forbidden when doing a solve.");

                                       /**
					* Exception.
					*/
      DeclException1 (ExcNonMatchingMaps,
		      std::string,
		      << "The Epetra map in the vector "
		      << arg1
		      << " does not match the one used in the matrix. "
		      << "Check vector and matrix setup.");

    protected:

                                       /**
                                        * Reference to the object that
                                        * controls convergence of the
                                        * iterative solver. In fact,
                                        * for these Trilinos wrappers,
                                        * Trilinos does so itself, but
                                        * we copy the data from this
                                        * object before starting the
                                        * solution process, and copy
                                        * the data back into it
                                        * afterwards.
                                        */
      SolverControl &solver_control;

    protected:

                                       /**
                                        * Store a copy of the flags for this
                                        * particular solver.
                                        */
      AdditionalData additional_data;

  };



/**
 * An implementation of the solver interface using the Trilinos CG
 * solver for block matrices and vectors.
 *
 * @ingroup TrilinosWrappers
 * @author Martin Kronbichler, 2008
 */
  class SolverBlockCG : public SolverBlockBase
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
					* Restart parameter in case
					* the derived class is
					* GMRES.
					* 
					* TODO: Find a better way for
					* doing this.
					*/
	AdditionalData (const bool         right_preconditioning   = false,
			const bool         output_details          = false);

				       /**
					* Flag that turns right
					* preconditioning on.
					*/
	bool right_preconditioning;

				       /**
					* A flag to determine whether
					* solver details should be
					* written to screen.
					*/
	bool output_details;
      };

                                       /**
                                        * Constructor. In contrast to
                                        * deal.II's own solvers, there is no
                                        * need to give a vector memory
                                        * object.
					*
                                        * The last argument takes a structure
                                        * with additional, solver dependent
                                        * flags for tuning.
                                        */
      SolverBlockCG (SolverControl        &cn,
		     const AdditionalData &data = AdditionalData());

    protected:

                                       /**
                                        * Store a copy of the flags for this
                                        * particular solver.
                                        */
      const AdditionalData data;
  };



/**
 * An implementation of the solver interface using the Trilinos GMRES
 * solver for block matrices and vectors.
 *
 * @author Martin Kronbichler, 2008
 */
  class SolverBlockGMRES : public SolverBlockBase
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
					* Restart parameter in case
					* the derived class is
					* GMRES.
					* 
					* TODO: Find a better way for
					* doing this.
					*/
	AdditionalData (const unsigned int gmres_restart_parameter = 30,
			const bool         right_preconditioning   = false,
			const bool         output_details          = false);

				       /**
					* The number of vectors that
					* should be used in the GMRES
					* algorithm before the solver
					* is restarted.
					*/
	unsigned int gmres_restart_parameter;

				       /**
					* Flag that turns right
					* preconditioning on.
					*/
	bool right_preconditioning;

				       /**
					* A flag to determine whether
					* solver details should be
					* written to screen.
					*/
	bool output_details;
      };

                                       /**
                                        * Constructor. In contrast to
                                        * deal.II's own solvers, there is no
                                        * need to give a vector memory
                                        * object.
                                        *
                                        * The last argument takes a structure
                                        * with additional, solver dependent
                                        * flags for tuning.
                                        */
      SolverBlockGMRES (SolverControl        &cn,
			const AdditionalData &data = AdditionalData());

    protected:

                                       /**
                                        * Store a copy of the flags for this
                                        * particular solver.
                                        */
      const AdditionalData data;
  };


}

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_USE_TRILINOS

/*----------------------   trilinos_solver_block.h     -----------------------*/

#endif
/*----------------------   trilinos_solver_block.h     -----------------------*/
