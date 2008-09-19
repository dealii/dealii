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
  
/**  
 * @ingroup TrilinosWrappers
 * @ingroup Preconditioners
 * @author Martin Kronbichler, 2008
 */
  class PreconditionJacobi : public Subscriptor
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
                                        * Constructor. Does not do
                                        * anything. The
                                        * <tt>initialize</tt> function
                                        * will have to set create the
                                        * partitioner.
                                        */
      PreconditionJacobi ();

                                       /**
					* Take the matrix which is
                                        * used to form the
                                        * preconditioner, and
                                        * additional flags if there
                                        * are any.
					*/
      void initialize (const SparseMatrix   &matrix,
		       const AdditionalData &additional_data = AdditionalData());
      
				       /**
					* Apply the preconditioner.
					*/
      void vmult (VectorBase       &dst,
		  const VectorBase &src) const;

    private:
      Teuchos::RefCountPtr<Ifpack_Preconditioner> preconditioner;
      //std::auto_ptr<Ifpack_Preconditioner> preconditioner;
      //boost::shared_ptr<Ifpack_Preconditioner> preconditioner;
  };
  

  
  
/**  
 * @ingroup TrilinosWrappers
 * @ingroup Preconditioners
 * @author Wolfgang Bangerth, 2008
 */
  class PreconditionSSOR : public Subscriptor
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
                                        * Constructor. Does not do
                                        * anything. The
                                        * <tt>initialize</tt> function
                                        * will have to set create the
                                        * partitioner.
                                        */
      PreconditionSSOR ();

                                       /**
					* Take the matrix which is
                                        * used to form the
                                        * preconditioner, and
                                        * additional flags if there
                                        * are any.
					*/
      void initialize (const SparseMatrix   &matrix,
		       const AdditionalData &additional_data = AdditionalData());
      
				       /**
					* Apply the preconditioner.
					*/
      void vmult (VectorBase       &dst,
		  const VectorBase &src) const;

    private:
      Teuchos::RefCountPtr<Ifpack_Preconditioner> preconditioner;
  };
  

  
  
/**  
 * @ingroup TrilinosWrappers
 * @ingroup Preconditioners
 * @author Martin Kronbichler, 2008
 */
  class PreconditionSOR : public Subscriptor
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
                                        * Constructor. Does not do
                                        * anything. The
                                        * <tt>initialize</tt> function
                                        * will have to set create the
                                        * partitioner.
                                        */
      PreconditionSOR ();

                                       /**
					* Take the matrix which is
                                        * used to form the
                                        * preconditioner, and
                                        * additional flags if there
                                        * are any.
					*/
      void initialize (const SparseMatrix   &matrix,
		       const AdditionalData &additional_data = AdditionalData());
      
				       /**
					* Apply the preconditioner.
					*/
      void vmult (VectorBase       &dst,
		  const VectorBase &src) const;

    private:
      Teuchos::RefCountPtr<Ifpack_Preconditioner> preconditioner;
  };



/**  
 * @ingroup TrilinosWrappers
 * @ingroup Preconditioners
 * @author Martin Kronbichler, 2008
 */
  class PreconditionIC : public Subscriptor
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
	AdditionalData (const double       ic_fill = 0.,
			const double       ic_atol = 0.,
			const double       ic_rtol = 1.,
			const unsigned int overlap = 0);
	
                                       /**
					* IC parameters and overlap.
					*/
	double ic_fill, ic_atol, ic_rtol;
	unsigned int overlap;
      };

                                       /**
					* (Empty) constructor.
					*/
      PreconditionIC ();

                                       /**
                                        * Initialize function. Takes
                                        * the matrix which is used to
                                        * form the preconditioner, and
                                        * additional flags if there
                                        * are any.
                                        */
      void initialize (const SparseMatrix   &matrix,
		       const AdditionalData &additional_data = AdditionalData());
      
				       /**
					* Apply the preconditioner.
					*/
      void vmult (VectorBase       &dst,
		  const VectorBase &src) const;

    private:
      Teuchos::RefCountPtr<Ifpack_Preconditioner> preconditioner;
  };



/**  
 * @ingroup TrilinosWrappers
 * @ingroup Preconditioners
 * @author Martin Kronbichler, 2008
 */
  class PreconditionILU : public Subscriptor
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
	AdditionalData (const double       ilu_fill = 0.,
			const double       ilu_atol = 0.,
			const double       ilu_rtol = 1.,
			const unsigned int overlap  = 0);
	
                                       /**
					* ILU parameters and overlap.
					*/
	double ilu_drop, ilu_fill, ilu_atol, ilu_rtol;
	unsigned int overlap;
      };

                                       /**
					* (Empty) constructor.
					*/
      PreconditionILU ();

                                       /**
                                        * Initialize function. Takes
                                        * the matrix which is used to
                                        * form the preconditioner, and
                                        * additional flags if there
                                        * are any.
                                        */
      void initialize (const SparseMatrix   &matrix,
		       const AdditionalData &additional_data = AdditionalData());
      
				       /**
					* Apply the preconditioner.
					*/
      void vmult (VectorBase       &dst,
		  const VectorBase &src) const;

    private:
      Teuchos::RefCountPtr<Ifpack_Preconditioner> preconditioner;
  };
  

}

/*@}*/


DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_USE_TRILINOS

/*----------------------------   trilinos_precondition.h     ---------------------------*/

#endif
/*----------------------------   trilinos_precondition.h     ---------------------------*/
