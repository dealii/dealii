//----------------------------  petsc_precondition.h  ---------------------------
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
//----------------------------  petsc_precondition.h  ---------------------------
#ifndef __deal2__petsc_precondition_h
#define __deal2__petsc_precondition_h

#include <base/config.h>
#include <lac/exceptions.h>

#ifdef DEAL_II_USE_PETSC

#include <petscpc.h>


namespace PETScWrappers
{
                                   // forward declarations
  class MatrixBase;
  class VectorBase;
  class SolverBase;
  
  
/**
 * Base class for preconditioner classes using the PETSc functionality. The
 * classes in this hierarchy don't do a whole lot, except for providing a
 * function that sets the preconditioner and certain parameters on the
 * preconditioning context of the solver. These classes are basically here
 * only to allow a similar interface than already used for the deal.II solver
 * and preconditioner classes.
 *
 * @ingroup PETScWrappers
 * @author Wolfgang Bangerth, 2004
 */
  class PreconditionerBase
  {
    public:
                                       /**
                                        * Constructor. Take a pointer to the
                                        * matrix from which the preconditioner
                                        * shall be constructed.
                                        */
      PreconditionerBase (const MatrixBase &matrix);
      
                                       /**
                                        * Destructor.
                                        */
      virtual ~PreconditionerBase ();


    protected:
                                       /**
                                        * A pointer to the matrix that acts as
                                        * a preconditioner.
                                        */
      const Mat matrix;
      
                                       /**
                                        * Conversion operator to get a
                                        * representation of the matrix that
                                        * represents this preconditioner. We
                                        * use this inside the actual solver,
                                        * where we need to pass this matrix to
                                        * the PETSc solvers.
                                        */
      operator const Mat () const;

                                       /**
                                        * Function that takes a Krylov
                                        * Subspace Preconditioner context
                                        * object, and sets the type of
                                        * preconditioner that is requested by
                                        * the derived class.
                                        */
      virtual void set_preconditioner_type (PC &pc) const = 0;

                                       /**
                                        * Make the solver class a friend,
                                        * since it needs to call the
                                        * conversion operator.
                                        */
      friend class SolverBase;
  };



/**
 * A class that implements the interface to use the PETSc Jacobi
 * preconditioner.
 *
 * @ingroup PETScWrappers
 * @author Wolfgang Bangerth, 2004
 */
  class PreconditionJacobi : public PreconditionerBase
  {
    public:
                                       /**
                                        * Standardized data struct to
                                        * pipe additional flags to the
                                        * preconditioner.
                                        */      
      struct AdditionalData
      {};

                                       /**
                                        * Constructor. Take the matrix which
                                        * is used to form the preconditioner,
                                        * and additional flags if there are
                                        * any.
                                        */
      PreconditionJacobi (const MatrixBase     &matrix,
                          const AdditionalData &additional_data = AdditionalData());
      
    protected:
                                       /**
                                        * Store a copy of the flags for this
                                        * particular preconditioner.
                                        */
      const AdditionalData additional_data;

                                       /**
                                        * Function that takes a Krylov
                                        * Subspace Preconditioner context
                                        * object, and sets the type of
                                        * preconditioner that is appropriate
                                        * for the present class.
                                        */
      virtual void set_preconditioner_type (PC &pc) const;      
  };

      

/**
 * A class that implements the interface to use the PETSc Block Jacobi
 * preconditioner. The blocking structure of the matrix is determined by the
 * association of degrees of freedom to the individual processors in an
 * MPI-parallel job. If you use this preconditioner on a sequential job (or an
 * MPI job with only one process) then the entire matrix is the only block.
 *
 * By default, PETSc uses an ILU(0) decomposition of each diagonal block of
 * the matrix for preconditioning. This can be changed, as is explained in the
 * relevant section of the PETSc manual, but is not implemented here.
 *
 * @ingroup PETScWrappers
 * @author Wolfgang Bangerth, 2004
 */
  class PreconditionBlockJacobi : public PreconditionerBase
  {
    public:
                                       /**
                                        * Standardized data struct to
                                        * pipe additional flags to the
                                        * preconditioner.
                                        */      
      struct AdditionalData
      {};

                                       /**
                                        * Constructor. Take the matrix which
                                        * is used to form the preconditioner,
                                        * and additional flags if there are
                                        * any.
                                        */
      PreconditionBlockJacobi (const MatrixBase     &matrix,
                               const AdditionalData &additional_data = AdditionalData());
      
    protected:
                                       /**
                                        * Store a copy of the flags for this
                                        * particular preconditioner.
                                        */
      const AdditionalData additional_data;

                                       /**
                                        * Function that takes a Krylov
                                        * Subspace Preconditioner context
                                        * object, and sets the type of
                                        * preconditioner that is appropriate
                                        * for the present class.
                                        */
      virtual void set_preconditioner_type (PC &pc) const;      
  };

      

/**
 * A class that implements the interface to use the PETSc SOR
 * preconditioner.
 *
 * @ingroup PETScWrappers
 * @author Wolfgang Bangerth, 2004
 */
  class PreconditionSOR : public PreconditionerBase
  {
    public:
                                       /**
                                        * Standardized data struct to
                                        * pipe additional flags to the
                                        * preconditioner.
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
                                        * Constructor. Take the matrix which
                                        * is used to form the preconditioner,
                                        * and additional flags if there are
                                        * any.
                                        *
                                        * We specify the (default) value to
                                        * the constructor call in the default
                                        * argument because otherwise gcc 2.95
                                        * generates a compiler fault.
                                        */
      PreconditionSOR (const MatrixBase     &matrix,
                       const AdditionalData &additional_data = AdditionalData(1));
      
    protected:
                                       /**
                                        * Store a copy of the flags for this
                                        * particular preconditioner.
                                        */
      const AdditionalData additional_data;

                                       /**
                                        * Function that takes a Krylov
                                        * Subspace Preconditioner context
                                        * object, and sets the type of
                                        * preconditioner that is appropriate
                                        * for the present class.
                                        */
      virtual void set_preconditioner_type (PC &pc) const;      
  };
      


/**
 * A class that implements the interface to use the PETSc SSOR
 * preconditioner.
 *
 * @ingroup PETScWrappers
 * @author Wolfgang Bangerth, 2004
 */
  class PreconditionSSOR : public PreconditionerBase
  {
    public:
                                       /**
                                        * Standardized data struct to
                                        * pipe additional flags to the
                                        * preconditioner.
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
                                        * Constructor. Take the matrix which
                                        * is used to form the preconditioner,
                                        * and additional flags if there are
                                        * any.
                                        *
                                        * We specify the (default) value to
                                        * the constructor call in the default
                                        * argument because otherwise gcc 2.95
                                        * generates a compiler fault.
                                        */
      PreconditionSSOR (const MatrixBase     &matrix,
                        const AdditionalData &additional_data = AdditionalData(1));
      
    protected:
                                       /**
                                        * Store a copy of the flags for this
                                        * particular preconditioner.
                                        */
      const AdditionalData additional_data;

                                       /**
                                        * Function that takes a Krylov
                                        * Subspace Preconditioner context
                                        * object, and sets the type of
                                        * preconditioner that is appropriate
                                        * for the present class.
                                        */
      virtual void set_preconditioner_type (PC &pc) const;      
  };
      


/**
 * A class that implements the interface to use the PETSc Eisenstat
 * preconditioner.
 *
 * @ingroup PETScWrappers
 * @author Wolfgang Bangerth, 2004
 */
  class PreconditionEisenstat : public PreconditionerBase
  {
    public:
                                       /**
                                        * Standardized data struct to
                                        * pipe additional flags to the
                                        * preconditioner.
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
                                        * Constructor. Take the matrix which
                                        * is used to form the preconditioner,
                                        * and additional flags if there are
                                        * any.
                                        *
                                        * We specify the (default) value to
                                        * the constructor call in the default
                                        * argument because otherwise gcc 2.95
                                        * generates a compiler fault.
                                        */
      PreconditionEisenstat (const MatrixBase     &matrix,
                             const AdditionalData &additional_data = AdditionalData(1));
      
    protected:
                                       /**
                                        * Store a copy of the flags for this
                                        * particular preconditioner.
                                        */
      const AdditionalData additional_data;

                                       /**
                                        * Function that takes a Krylov
                                        * Subspace Preconditioner context
                                        * object, and sets the type of
                                        * preconditioner that is appropriate
                                        * for the present class.
                                        */
      virtual void set_preconditioner_type (PC &pc) const;      
  };
      


/**
 * A class that implements the interface to use the PETSc Incomplete Cholesky
 * preconditioner.
 *
 * @ingroup PETScWrappers
 * @author Wolfgang Bangerth, 2004
 */
  class PreconditionICC : public PreconditionerBase
  {
    public:
                                       /**
                                        * Standardized data struct to
                                        * pipe additional flags to the
                                        * preconditioner.
                                        */      
      struct AdditionalData
      {
                                           /**
                                            * Constructor. By default,
                                            * set the fill-in parameter
                                            * to zero.
                                            */
          AdditionalData (const unsigned int levels = 0);
	
                                           /**
                                            * Fill-in parameter.
                                            */
          unsigned int levels;
      };

                                       /**
                                        * Constructor. Take the matrix which
                                        * is used to form the preconditioner,
                                        * and additional flags if there are
                                        * any.
                                        *
                                        * We specify the (default) value to
                                        * the constructor call in the default
                                        * argument because otherwise gcc 2.95
                                        * generates a compiler fault.
                                        */
      PreconditionICC (const MatrixBase     &matrix,
                       const AdditionalData &additional_data = AdditionalData(0));
      
    protected:
                                       /**
                                        * Store a copy of the flags for this
                                        * particular preconditioner.
                                        */
      const AdditionalData additional_data;

                                       /**
                                        * Function that takes a Krylov
                                        * Subspace Preconditioner context
                                        * object, and sets the type of
                                        * preconditioner that is appropriate
                                        * for the present class.
                                        */
      virtual void set_preconditioner_type (PC &pc) const;      
  };


  
/**
 * A class that implements the interface to use the PETSc ILU
 * preconditioner.
 *
 * @ingroup PETScWrappers
 * @author Wolfgang Bangerth, 2004
 */
  class PreconditionILU : public PreconditionerBase
  {
    public:
                                       /**
                                        * Standardized data struct to
                                        * pipe additional flags to the
                                        * preconditioner.
                                        */      
      struct AdditionalData
      {
                                           /**
                                            * Constructor. By default,
                                            * set the fill-in parameter
                                            * to zero.
                                            */
          AdditionalData (const unsigned int levels = 0);
	
                                           /**
                                            * Fill-in parameter.
                                            */
          unsigned int levels;
      };

                                       /**
                                        * Constructor. Take the matrix which
                                        * is used to form the preconditioner,
                                        * and additional flags if there are
                                        * any.
                                        *
                                        * We specify the (default) value to
                                        * the constructor call in the default
                                        * argument because otherwise gcc 2.95
                                        * generates a compiler fault.
                                        */
      PreconditionILU (const MatrixBase     &matrix,
                       const AdditionalData &additional_data = AdditionalData(0));
      
    protected:
                                       /**
                                        * Store a copy of the flags for this
                                        * particular preconditioner.
                                        */
      const AdditionalData additional_data;

                                       /**
                                        * Function that takes a Krylov
                                        * Subspace Preconditioner context
                                        * object, and sets the type of
                                        * preconditioner that is appropriate
                                        * for the present class.
                                        */
      virtual void set_preconditioner_type (PC &pc) const;      
  };



/**
 * A class that implements the interface to use the PETSc LU
 * preconditioner. The LU decomposition is only implemented for single
 * processor machines. It should provide a convenient interface to
 * another direct solver.
 *
 * @ingroup PETScWrappers
 * @author Oliver Kayser-Herold, 2004
 */
  class PreconditionLU : public PreconditionerBase
  {
    public:
                                       /**
                                        * Standardized data struct to
                                        * pipe additional flags to the
                                        * preconditioner.
                                        */      
      struct AdditionalData
      {
                                           /**
                                            * Constructor. (Default values
					    * taken from function PCCreate_LU
					    * of the PetSC lib.)
                                            */
          AdditionalData (const double pivoting = 1.e-6,
			  const double zero_pivot = 1.e-12,
			  const double damping = 0.0);

                                           /**
                                            * Determines, when Pivoting is
					    * done during LU decomposition.
					    * 0.0 indicates no pivoting,
					    * and 1.0 complete pivoting.
					    * Confer PetSC manual for more
					    * details.
                                            */
          double pivoting;

                                           /**
                                            * Size at which smaller pivots
					    * are declared to be zero.
					    * Confer PetSC manual for more
					    * details.
					    */
	  double zero_pivot;					    

                                           /**
                                            * This quantity is added to the
					    * diagonal of the matrix during
					    * factorisation.
					    */
	  double damping;
      };

                                       /**
                                        * Constructor. Take the matrix which
                                        * is used to form the preconditioner,
                                        * and additional flags if there are
                                        * any.
                                        *
                                        * We specify the (default) value to
                                        * the constructor call in the default
                                        * argument because otherwise gcc 2.95
                                        * generates a compiler fault.
                                        */
      PreconditionLU (const MatrixBase     &matrix,
		      const AdditionalData &additional_data = 
		      AdditionalData(1.e-6, 1.e-12, 0.0));
      
    protected:
                                       /**
                                        * Store a copy of the flags for this
                                        * particular preconditioner.
                                        */
      const AdditionalData additional_data;

                                       /**
                                        * Function that takes a Krylov
                                        * Subspace Preconditioner context
                                        * object, and sets the type of
                                        * preconditioner that is appropriate
                                        * for the present class.
                                        */
      virtual void set_preconditioner_type (PC &pc) const;      
  };

}


#endif // DEAL_II_USE_PETSC

/*----------------------------   petsc_precondition.h     ---------------------------*/

#endif
/*----------------------------   petsc_precondition.h     ---------------------------*/
