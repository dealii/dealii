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
#ifndef __deal2__trilinos_precondition_h
#define __deal2__trilinos_precondition_h


#include <base/config.h>
#include <base/subscriptor.h>


#ifdef DEAL_II_USE_TRILINOS

#include <Teuchos_RefCountPtr.hpp>

// forward declarations
class Ifpack_Preconditioner;

DEAL_II_NAMESPACE_OPEN

/*! @addtogroup TrilinosWrappers
 *@{
 */

namespace TrilinosWrappers
{

  class VectorBase;
  class SparseMatrix;
  class SolverBase;

/**
 * The base class for all preconditioners based on Trilinos sparse
 * matrices.
 *
 * @ingroup TrilinosWrappers
 * @ingroup Preconditioners
 * @author Martin Kronbichler, 2008
 */
  class PreconditionBase : public Subscriptor
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
                                        * Constructor. Does not do
                                        * anything. The
                                        * <tt>initialize</tt> function
                                        * of the derived classes will
                                        * have to create the
                                        * preconditioner from a given
                                        * sparse matrix.
                                        */
      PreconditionBase ();
      
				       /**
					* Apply the preconditioner.
					*/
      void vmult (VectorBase       &dst,
		  const VectorBase &src) const;

                                       /**
					* Exception.
					*/
      DeclException1 (ExcNonMatchingMaps,
		      std::string,
		      << "The sparse matrix that the preconditioner is based "
		      << "on a map that is not compatible to the one in vector "
		      << arg1
		      << ". Check preconditioner and matrix setup.");

      friend class SolverBase;

    protected:
				       /**
					* This is a pointer to the
					* preconditioner object.
					*/
        Teuchos::RefCountPtr<Ifpack_Preconditioner> preconditioner;

  };

  
/**
 * A wrapper class for a (pointwise) Jacobi preconditioner for
 * Trilinos matrices. This preconditioner works both in serial and in
 * parallel, depending on the matrix it is based on.
 *
 * The AdditionalData data structure allows to set preconditioner
 * options. For the Jacobi preconditioner, these options are the
 * damping parameter <tt>omega</tt> and a <tt>min_diagonal</tt>
 * argument that can be used to make the preconditioner work even if
 * the matrix contains some zero elements on the diagonal. The default
 * settings are 1 for the damping parameter and zero for the diagonal
 * augmentation.
 *
 * @ingroup TrilinosWrappers 
 * @ingroup Preconditioners 
 * @author Martin Kronbichler, 2008
 */
  class PreconditionJacobi : public PreconditionBase
  {
    public:

                                       /**
                                        * Standardized data struct to
                                        * pipe additional flags to the
                                        * preconditioner. The
                                        * parameter <tt>omega</tt>
                                        * specifies the relaxation
                                        * parameter in the Jacobi
                                        * preconditioner. The
                                        * parameter
                                        * <tt>min_diagonal</tt> can be
                                        * used to make the application
                                        * of the preconditioner also
                                        * possible when some diagonal
                                        * elements are zero. In a
                                        * default application this
                                        * would mean that we divide by
                                        * zero, so by setting the
                                        * parameter
                                        * <tt>min_diagonal</tt> to a
                                        * small nonzero value the SOR
                                        * will work on a matrix that
                                        * is not too far away from the
                                        * one we want to
                                        * treat.
                                        */      
      struct AdditionalData
      {
                                       /**
					* Constructor. By default, set
					* the damping parameter to
					* one, and do not modify the
					* diagonal.
					*/
	AdditionalData (const double       omega = 1,
			const double       min_diagonal = 0);
	
                                       /**
					* Relaxation parameter and
					* minimal diagonal value.
					*/
	double omega, min_diagonal;
      };

                                       /**
					* Take the sparse matrix the
                                        * preconditioner object should
                                        * be built of, and additional
                                        * flags (damping parameter,
                                        * etc.) if there are any.
					*/
      void initialize (const SparseMatrix   &matrix,
		       const AdditionalData &additional_data = AdditionalData());
  };
  

  
  
/**
 * A wrapper class for a (pointwise) SSOR preconditioner for Trilinos
 * matrices. This preconditioner works both in serial and in parallel,
 * depending on the matrix it is based on.
 *
 * The AdditionalData data structure allows to set preconditioner
 * options. For the SSOR preconditioner, these options are the
 * damping/relaxation parameter <tt>omega</tt>, a
 * <tt>min_diagonal</tt> argument that can be used to make the
 * preconditioner work even if the matrix contains some zero elements
 * on the diagonal, and a parameter <tt>overlap</tt> that determines
 * if and how much overlap there should be between the matrix
 * partitions on the various MPI processes. The default settings are 1
 * for the relaxation parameter, 0 for the diagonal augmentation and 0
 * for the overlap.
 *
 * Note that a parallel applicatoin of the SSOR preconditioner is
 * actually a block-Jacobi preconditioner with block size equal to the
 * local matrix size. Spoken more technically, this parallel operation
 * is an <a
 * href="http://en.wikipedia.org/wiki/Additive_Schwarz_method">additive
 * Schwarz method</a> with an SSOR <em>approximate solve</em> as inner
 * solver, based on the outer parallel partitioning.
 *
 * @ingroup TrilinosWrappers
 * @ingroup Preconditioners
 * @author Wolfgang Bangerth, 2008
 */
  class PreconditionSSOR : public PreconditionBase
  {
    public:
 
                                       /**
                                        * Standardized data struct to
                                        * pipe additional flags to the
                                        * preconditioner. The
                                        * parameter <tt>omega</tt>
                                        * specifies the relaxation
                                        * parameter in the SSOR
                                        * preconditioner. The
                                        * parameter
                                        * <tt>min_diagonal</tt> can be
                                        * used to make the application
                                        * of the preconditioner also
                                        * possible when some diagonal
                                        * elements are zero. In a
                                        * default application this
                                        * would mean that we divide by
                                        * zero, so by setting the
                                        * parameter
                                        * <tt>min_diagonal</tt> to a
                                        * small nonzero value the SOR
                                        * will work on a matrix that
                                        * is not too far away from the
                                        * one we want to
                                        * treat. Finally,
                                        * <tt>overlap</tt> governs the
                                        * overlap of the partitions
                                        * when the preconditioner runs
                                        * in parallel, forming a
                                        * so-called additive Schwarz
                                        * preconditioner.
                                        */      
      struct AdditionalData
      {
                                       /**
					* Constructor. By default, set
					* the damping parameter to
					* one, we do not modify the
					* diagonal, and there is no
					* overlap (i.e. in parallel,
					* we run a BlockJacobi
					* preconditioner, where each
					* block is inverted
					* approximately by an SSOR.
					*/
	AdditionalData (const double       omega = 1,
			const double       min_diagonal = 0,
			const unsigned int overlap = 0);
	
                                       /**
					* Relaxation parameter,
					* minimal diagonal element,
					* and overlap.
					*/
	double omega, min_diagonal;
	unsigned int overlap; 
      };

                                       /**
					* Take the sparse matrix the
                                        * preconditioner object should
                                        * be built of, and additional
                                        * flags (damping parameter,
                                        * overlap in parallel
                                        * computations, etc.) if there
                                        * are any.
					*/
      void initialize (const SparseMatrix   &matrix,
		       const AdditionalData &additional_data = AdditionalData());
  };
  

  
  
/**
 * A wrapper class for a (pointwise) SOR preconditioner for Trilinos
 * matrices. This preconditioner works both in serial and in parallel,
 * depending on the matrix it is based on.
 *
 * The AdditionalData data structure allows to set preconditioner
 * options. For the SOR preconditioner, these options are the
 * damping/relaxation parameter <tt>omega</tt>, a
 * <tt>min_diagonal</tt> argument that can be used to make the
 * preconditioner work even if the matrix contains some zero elements
 * on the diagonal, and a parameter <tt>overlap</tt> that determines
 * if and how much overlap there should be between the matrix
 * partitions on the various MPI processes. The default settings are 1
 * for the relaxation parameter, 0 for the diagonal augmentation and 0
 * for the overlap. 
 *
 * Note that a parallel applicatoin of the SOR preconditioner is
 * actually a block-Jacobi preconditioner with block size equal to the
 * local matrix size. Spoken more technically, this parallel operation
 * is an <a
 * href="http://en.wikipedia.org/wiki/Additive_Schwarz_method">additive
 * Schwarz method</a> with an SOR <em>approximate solve</em> as inner
 * solver, based on the outer parallel partitioning.
 *
 * @ingroup TrilinosWrappers
 * @ingroup Preconditioners
 * @author Martin Kronbichler, 2008
 */
  class PreconditionSOR : public PreconditionBase
  {
    public:

                                       /**
                                        * Standardized data struct to
                                        * pipe additional flags to the
                                        * preconditioner. The
                                        * parameter <tt>omega</tt>
                                        * specifies the relaxation
                                        * parameter in the SOR
                                        * preconditioner. The
                                        * parameter
                                        * <tt>min_diagonal</tt> can be
                                        * used to make the application
                                        * of the preconditioner also
                                        * possible when some diagonal
                                        * elements are zero. In a
                                        * default application this
                                        * would mean that we divide by
                                        * zero, so by setting the
                                        * parameter
                                        * <tt>min_diagonal</tt> to a
                                        * small nonzero value the SOR
                                        * will work on a matrix that
                                        * is not too far away from the
                                        * one we want to
                                        * treat. Finally,
                                        * <tt>overlap</tt> governs the
                                        * overlap of the partitions
                                        * when the preconditioner runs
                                        * in parallel, forming a
                                        * so-called additive Schwarz
                                        * preconditioner.
                                        */      
      struct AdditionalData
      {
                                       /**
					* Constructor. By default, set
					* the damping parameter to
					* one, we do not modify the
					* diagonal, and there is no
					* overlap (i.e. in parallel,
					* we run a BlockJacobi
					* preconditioner, where each
					* block is inverted
					* approximately by an SOR.
					*/
	AdditionalData (const double       omega = 1,
			const double       min_diagonal = 0,
			const unsigned int overlap = 0);
	
                                       /**
					* Relaxation parameter,
					* minimal diagonal element,
					* and overlap.
					*/
	double omega, min_diagonal;
	unsigned int overlap; 
      };

                                       /**
					* Take the sparse matrix the
                                        * preconditioner object should
                                        * be built of, and additional
                                        * flags (damping parameter,
                                        * overlap in parallel
                                        * computations etc.) if there
                                        * are any.
					*/
      void initialize (const SparseMatrix   &matrix,
		       const AdditionalData &additional_data = AdditionalData());
  };



/**
 * A wrapper class for an incomplete Cholesky factorization (IC)
 * preconditioner for @em symmetric Trilinos matrices. This
 * preconditioner works both in serial and in parallel, depending on
 * the matrix it is based on. In general, an incomplete factorization
 * does not take all fill-in elements that would appear in a full
 * factorization (that is the basis for a direct solve). Trilinos
 * allows to set the amount of fill-in elements, governed by the
 * additional data argument <tt>ic_fill</tt>, so one can gradually
 * choose between a factorization on the sparse matrix structure only
 * (<tt>ic_fill=0</tt>) to a full factorization (<tt>ic_fill</tt> in
 * the range of 10 to 50, depending on the spatial dimension of the
 * PDE problem and the degree of the finite element basis functions;
 * generally, more required fill-in elements require this parameter to
 * be set to a higher integer value).
 *
 * The AdditionalData data structure allows to set preconditioner
 * options. Besides the fill-in argument, these options are some
 * options for perturbations (see the documentation of the
 * AdditionalData structure for details), and a parameter
 * <tt>overlap</tt> that determines if and how much overlap there
 * should be between the matrix partitions on the various MPI
 * processes. The default settings are 1 for the relaxation parameter,
 * 0 for the diagonal augmentation and 0 for the overlap.
 *
 * Note that a parallel applicatoin of the IC preconditioner is
 * actually a block-Jacobi preconditioner with block size equal to the
 * local matrix size. Spoken more technically, this parallel operation
 * is an <a
 * href="http://en.wikipedia.org/wiki/Additive_Schwarz_method">additive
 * Schwarz method</a> with an IC <em>approximate solve</em> as inner
 * solver, based on the (outer) parallel partitioning.
 *
 * @ingroup TrilinosWrappers
 * @ingroup Preconditioners
 * @author Martin Kronbichler, 2008
 */
  class PreconditionIC : public PreconditionBase
  {
    public:
                                       /**
                                        * Standardized data struct to
                                        * pipe additional parameters
                                        * to the preconditioner. The
                                        * Trilinos IC decomposition
                                        * allows for some fill-in, so
                                        * it actually is a threshold
                                        * incomplete Cholesky
                                        * factorization. The amount of
                                        * fill-in, and hence, the
                                        * amount of memory used by
                                        * this preconditioner, is
                                        * controlled by the parameter
                                        * <tt>ic_fill</tt>, which
                                        * specifies this as a
                                        * double. When forming the
                                        * preconditioner, for certain
                                        * problems bad conditioning
                                        * (or just bad luck) can cause
                                        * the preconditioner to be
                                        * very poorly
                                        * conditioned. Hence it can
                                        * help to add diagonal
                                        * perturbations to the
                                        * original matrix and form the
                                        * preconditioner for this
                                        * slightly better
                                        * matrix. <tt>ic_atol</tt> is
                                        * an absolute perturbation
                                        * that is added to the
                                        * diagonal before forming the
                                        * prec, and <tt>ic_rtol</tt>
                                        * is a scaling factor $rtol
                                        * \geq 1$. The last parameter
                                        * specifies the overlap of the
                                        * partitions when the
                                        * preconditioner runs in
                                        * parallel.
                                        */      
      struct AdditionalData
      {
                                       /**
					* Constructor. By default, set
					* the drop tolerance to 0, the
					* level of extra fill-ins is
					* set to be zero (just use the
					* matrix structure, do not
					* generate any additional
					* fill-in), the tolerance
					* level are 0 and 1,
					* respectively, and the
					* overlap in case of a
					* parallel execution is
					* zero. This overlap in a
					* block-application of the IC
					* in the parallel case makes
					* the preconditioner a
					* so-called additive Schwarz
					* preconditioner.
					*/
	AdditionalData (const unsigned int ic_fill = 0,
			const double       ic_atol = 0.,
			const double       ic_rtol = 1.,
			const unsigned int overlap = 0);
	
                                       /**
					* IC parameters and overlap.
					*/
	unsigned int ic_fill;
	double ic_atol, ic_rtol;
	unsigned int overlap;
      };

                                       /**
                                        * Initialize function. Takes
                                        * the matrix the
                                        * preconditioner should be
                                        * computed of, and additional
                                        * flags if there are any.
                                        */
      void initialize (const SparseMatrix   &matrix,
		       const AdditionalData &additional_data = AdditionalData());
  };



/**
 * A wrapper class for an incomplete LU factorization (ILU)
 * preconditioner for Trilinos matrices. This preconditioner works
 * both in serial and in parallel, depending on the matrix it is based
 * on. In general, an incomplete factorization does not take all
 * fill-in elements that would appear in a full factorization (that is
 * the basis for a direct solve). Trilinos allows to set the amount of
 * fill-in elements, governed by the additional data argument
 * <tt>ilu_fill</tt>, so one can gradually choose between a
 * factorization on the sparse matrix structure only
 * (<tt>ilu_fill=0</tt>) to a full factorization (<tt>ilu_fill</tt> in
 * the range of 10 to 50, depending on the spatial dimension of the
 * PDE problem and the degree of the finite element basis functions;
 * generally, more required fill-in elements require this parameter to
 * be set to a higher integer value).
 *
 * The AdditionalData data structure allows to set preconditioner
 * options. Besides the fill-in argument, these options are some
 * options for perturbations (see the documentation of the
 * AdditionalData structure for details), and a parameter
 * <tt>overlap</tt> that determines if and how much overlap there
 * should be between the matrix partitions on the various MPI
 * processes. The default settings are 1 for the relaxation parameter,
 * 0 for the diagonal augmentation and 0 for the overlap.
 *
 * Note that a parallel applicatoin of the ILU preconditioner is
 * actually a block-Jacobi preconditioner with block size equal to the
 * local matrix size. Spoken more technically, this parallel operation
 * is an <a
 * href="http://en.wikipedia.org/wiki/Additive_Schwarz_method">additive
 * Schwarz method</a> with an ILU <em>approximate solve</em> as inner
 * solver, based on the (outer) parallel partitioning.
 *
 * @ingroup TrilinosWrappers
 * @ingroup Preconditioners
 * @author Martin Kronbichler, 2008
 */
  class PreconditionILU : public PreconditionBase
  {
    public:
                                       /**
                                        * Standardized data struct to
                                        * pipe additional parameters
                                        * to the preconditioner. The
                                        * Trilinos ILU decomposition
                                        * allows for some fill-in, so
                                        * it actually is a threshold
                                        * incomplete LU
                                        * factorization. The amount of
                                        * fill-in, and hence, the
                                        * amount of memory used by
                                        * this preconditioner, is
                                        * controlled by the parameter
                                        * <tt>ilu_fill</tt>, which
                                        * specifies this as a
                                        * double. When forming the
                                        * preconditioner, for certain
                                        * problems bad conditioning
                                        * (or just bad luck) can cause
                                        * the preconditioner to be
                                        * very poorly
                                        * conditioned. Hence it can
                                        * help to add diagonal
                                        * perturbations to the
                                        * original matrix and form the
                                        * preconditioner for this
                                        * slightly better
                                        * matrix. <tt>ilu_atol</tt> is
                                        * an absolute perturbation
                                        * that is added to the
                                        * diagonal before forming the
                                        * prec, and <tt>ilu_rtol</tt>
                                        * is a scaling factor $rtol
                                        * \geq 1$. The last parameter
                                        * specifies the overlap of the
                                        * partitions when the
                                        * preconditioner runs in
                                        * parallel.
                                        */      
      struct AdditionalData
      {
                                       /**
					* Constructor. By default, the
					* level of extra fill-ins is
					* set to be zero (just use the
					* matrix structure, do not
					* generate any additional
					* fill-in), the tolerance
					* level are 0 and 1,
					* respectively, and the
					* overlap in case of a
					* parallel execution is
					* zero. This overlap in a
					* block-application of the ILU
					* in the parallel case makes
					* the preconditioner a
					* so-called additive Schwarz
					* preconditioner.
					*/
	AdditionalData (const unsigned int ilu_fill = 0,
			const double       ilu_atol = 0.,
			const double       ilu_rtol = 1.,
			const unsigned int overlap  = 0);
	
                                       /**
					* ILU parameters and overlap.
					*/
	unsigned int ilu_fill;
	double ilu_atol, ilu_rtol;
	unsigned int overlap;
      };

                                       /**
                                        * Initialize function. Takes
                                        * the matrix which is used to
                                        * form the preconditioner, and
                                        * additional flags if there
                                        * are any.
                                        */
      void initialize (const SparseMatrix   &matrix,
		       const AdditionalData &additional_data = AdditionalData());
  };
  

}

/*@}*/


DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_USE_TRILINOS

/*----------------------------   trilinos_precondition.h     ---------------------------*/

#endif
/*----------------------------   trilinos_precondition.h     ---------------------------*/
