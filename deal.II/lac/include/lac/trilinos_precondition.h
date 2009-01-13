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
#include <lac/vector.h>
#include <lac/sparse_matrix.h>

#include <base/std_cxx0x/shared_ptr.h>


#ifdef DEAL_II_USE_TRILINOS

#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
#  include <Epetra_MpiComm.h>
#else
#  include <Epetra_SerialComm.h>
#endif
#include <Epetra_Map.h>

#include <Teuchos_RCP.hpp>
#include <Epetra_Operator.h>

// forward declarations
class Ifpack_Preconditioner;
class Ifpack_Chebyshev;
namespace ML_Epetra
{
  class MultiLevelPreconditioner;
}


DEAL_II_NAMESPACE_OPEN

/*! @addtogroup TrilinosWrappers
 *@{
 */

namespace TrilinosWrappers
{

  class VectorBase;
  class SparseMatrix;
  class BlockSparseMatrix;
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
                                        * Copy constructor.
					*/
      PreconditionBase (const PreconditionBase &);

                                       /**
                                        * Destructor.
                                        */
      ~PreconditionBase ();
      
				       /**
					* Apply the preconditioner.
					*/
      void vmult (VectorBase       &dst,
		  const VectorBase &src) const;

				       /**
					* Apply the preconditioner on
				        * deal.II data structures
				        * instead of the ones provided
				        * in the Trilinos wrapper
				        * class.
					*/
      void vmult (dealii::Vector<double>       &dst,
		  const dealii::Vector<double> &src) const;

                                       /**
					* Exception.
					*/
      DeclException1 (ExcNonMatchingMaps,
		      std::string,
		      << "The sparse matrix the preconditioner is based on "
		      << "uses a map that is not compatible to the one in vector "
		      << arg1
		      << ". Check preconditioner and matrix setup.");

      friend class SolverBase;
      friend class PreconditionStokes;

    protected:
				       /**
					* This is a pointer to the
					* preconditioner object that
					* is used when applying the
					* preconditioner.
					*/
      Teuchos::RCP<const Epetra_Operator> preconditioner;

                                       /**
					* Internal communication
					* pattern in case the matrix
					* needs to be copied from
					* deal.II format.
					*/
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
      Epetra_MpiComm     communicator;
#else
      Epetra_SerialComm  communicator;
#endif

                                       /**
					* Internal Trilinos map in
					* case the matrix needs to be
					* copied from deal.II format.
					*/
      std::auto_ptr<Epetra_Map>   map;
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
					* This specifies the
					* relaxation parameter in the
					* Jacobi preconditioner.
					*/
	double omega;

                                       /**
					* This specifies the minimum
					* value the diagonal elements
					* should have. This might be
					* necessary when the Jacobi
					* preconditioner is used on
					* matrices with zero diagonal
					* elements. In that case, a
					* straight-forward application
					* would not be possible since
					* we would divide by zero.
					*/
	double min_diagonal;
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

    private:
				       /**
					* This is a pointer to the
					* Ifpack data contained in
					* this preconditioner.
					*/
      Teuchos::RCP<Ifpack_Preconditioner> ifpack;
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
					* This specifies the (over-)
					* relaxation parameter in the
					* SSOR preconditioner.
					*/
	double omega;

                                       /**
					* This specifies the minimum
					* value the diagonal elements
					* should have. This might be
					* necessary when the SSOR
					* preconditioner is used on
					* matrices with zero diagonal
					* elements. In that case, a
					* straight-forward application
					* would not be possible since
					* we divide by the diagonal
					* element.
					*/
	double min_diagonal;

                                       /**
					* This determines how large
					* the overlap of the local
					* matrix portions on each
					* processor in a parallel
					* application should be.
					*/
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

    private:
				       /**
					* This is a pointer to the
					* Ifpack data contained in
					* this preconditioner.
					*/
      Teuchos::RCP<Ifpack_Preconditioner> ifpack;
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
					* This specifies the (over-)
					* relaxation parameter in the
					* SOR preconditioner.
					*/
	double omega;

                                       /**
					* This specifies the minimum
					* value the diagonal elements
					* should have. This might be
					* necessary when the SOR
					* preconditioner is used on
					* matrices with zero diagonal
					* elements. In that case, a
					* straight-forward application
					* would not be possible since
					* we divide by the diagonal
					* element.
					*/
	double min_diagonal;

                                       /**
					* This determines how large
					* the overlap of the local
					* matrix portions on each
					* processor in a parallel
					* application should be.
					*/
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

    private:
				       /**
					* This is a pointer to the
					* Ifpack data contained in
					* this preconditioner.
					*/
      Teuchos::RCP<Ifpack_Preconditioner> ifpack;
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
 * processes.  The default settings are 0 for the additional fill-in, 0
 * for the absolute augmentation tolerance, 1 for the relative
 * augmentation tolerance, 0 for the overlap.
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
					* This specifies the amount of
					* additional fill-in elements
					* besides the sparse matrix
					* structure. When
					* <tt>ic_fill</tt> is large,
					* this means that many
					* fill-ins will be added, so
					* that the IC preconditioner
					* comes closer to a direct
					* sparse Cholesky
					* decomposition. Note,
					* however, that this will
					* drastically increase the
					* memory requirement,
					* especially when the
					* preconditioner is used in
					* 3D.
					*/
	unsigned int ic_fill;

                                       /**
					* This specifies the amount of
					* an absolute perturbation
					* that will be added to the
					* diagonal of the matrix,
					* which sometimes can help to
					* get better preconditioners.
					*/
	double ic_atol;

                                       /**
					* This specifies the factor by
					* which the diagonal of the
					* matrix will be scaled, which
					* sometimes can help to get
					* better preconditioners.
					*/
	double ic_rtol;

                                       /**
					* This determines how large
					* the overlap of the local
					* matrix portions on each
					* processor in a parallel
					* application should be.
					*/
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

    private:
				       /**
					* This is a pointer to the
					* Ifpack data contained in
					* this preconditioner.
					*/
      Teuchos::RCP<Ifpack_Preconditioner> ifpack;
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
 * processes. The default settings are 0 for the additional fill-in, 0
 * for the absolute augmentation tolerance, 1 for the relative
 * augmentation tolerance, 0 for the overlap.
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
					* This specifies the amount of
					* additional fill-in elements
					* besides the sparse matrix
					* structure. When
					* <tt>ilu_fill</tt> is large,
					* this means that many
					* fill-ins will be added, so
					* that the ILU preconditioner
					* comes closer to a (direct)
					* sparse LU
					* decomposition. Note,
					* however, that this will
					* drastically increase the
					* memory requirement,
					* especially when the
					* preconditioner is used in
					* 3D.
					*/
	unsigned int ilu_fill;

                                       /**
					* This specifies the amount of
					* an absolute perturbation
					* that will be added to the
					* diagonal of the matrix,
					* which sometimes can help to
					* get better preconditioners.
					*/
	double ilu_atol;

                                       /**
					* This specifies the factor by
					* which the diagonal of the
					* matrix will be scaled, which
					* sometimes can help to get
					* better preconditioners.
					*/
	double ilu_rtol;

                                       /**
					* This determines how large
					* the overlap of the local
					* matrix portions on each
					* processor in a parallel
					* application should be.
					*/
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

    private:
				       /**
					* This is a pointer to the
					* Ifpack data contained in
					* this preconditioner.
					*/
      Teuchos::RCP<Ifpack_Preconditioner> ifpack;
  };






/**
 * A wrapper class for a thresholded incomplete LU factorization (ILU-T)
 * preconditioner for Trilinos matrices. This preconditioner works both in
 * serial and in parallel, depending on the matrix it is based on. In
 * general, an incomplete factorization does not take all fill-in elements
 * that would appear in a full factorization (that is the basis for a direct
 * solve). For the ILU-T precondtioner, the parameter <tt>ilut_drop</tt>
 * lets the user specify which elements should be dropped (i.e., should not
 * be part of the incomplete decomposition). Trilinos calculates first the
 * complete factorization for one row, and then skips those elements that
 * are lower than the threshold. This is the main difference to the
 * non-thresholded ILU preconditioner, where the parameter
 * <tt>ilut_fill</tt> governs the incomplete factorization structure. This
 * parameter is available here as well, but provides only some extra
 * information here.
 *
 * The AdditionalData data structure allows to set preconditioner
 * options. Besides the fill-in arguments, these options are some options
 * for perturbations (see the documentation of the AdditionalData structure
 * for details), and a parameter <tt>overlap</tt> that determines if and how
 * much overlap there should be between the matrix partitions on the various
 * MPI processes. The default settings are 0 for the additional fill-in, 0
 * for the absolute augmentation tolerance, 1 for the relative augmentation
 * tolerance, 0 for the overlap.
 *
 * Note that a parallel applicatoin of the ILU-T preconditioner is
 * actually a block-Jacobi preconditioner with block size equal to the
 * local matrix size. Spoken more technically, this parallel operation
 * is an <a
 * href="http://en.wikipedia.org/wiki/Additive_Schwarz_method">additive
 * Schwarz method</a> with an ILU <em>approximate solve</em> as inner
 * solver, based on the (outer) parallel partitioning.
 *
 * @ingroup TrilinosWrappers
 * @ingroup Preconditioners
 * @author Martin Kronbichler, 2009
 */
  class PreconditionILUT : public PreconditionBase
  {
    public:
                                       /**
                                        * Standardized data struct to pipe
                                        * additional parameters to the
                                        * preconditioner. The Trilinos ILU-T
                                        * decomposition allows for some
                                        * fill-in, so it actually is a
                                        * threshold incomplete LU
                                        * factorization. The amount of
                                        * fill-in, and hence, the amount of
                                        * memory used by this
                                        * preconditioner, is controlled by
                                        * the parameters <tt>ilut_drop</tt>
                                        * and <tt>ilut_fill</tt>, which
                                        * specifies a threshold about which
                                        * values should form the incomplete
                                        * factorization and the level of
                                        * additional fill-in. When forming
                                        * the preconditioner, for certain
                                        * problems bad conditioning (or just
                                        * bad luck) can cause the
                                        * preconditioner to be very poorly
                                        * conditioned. Hence it can help to
                                        * add diagonal perturbations to the
                                        * original matrix and form the
                                        * preconditioner for this slightly
                                        * better matrix. <tt>ilut_atol</tt>
                                        * is an absolute perturbation that
                                        * is added to the diagonal before
                                        * forming the prec, and
                                        * <tt>ilu_rtol</tt> is a scaling
                                        * factor $rtol \geq 1$. The last
                                        * parameter specifies the overlap of
                                        * the partitions when the
                                        * preconditioner runs in parallel.
                                        */      
      struct AdditionalData
      {
                                       /**
					* Constructor. By default, no
					* element will be dropped, the level
					* of extra fill-ins is set to be
					* zero (just use the matrix
					* structure, do not generate any
					* additional fill-in except the one
					* that results from non-dropping
					* large elements), the tolerance
					* level are 0 and 1, respectively,
					* and the overlap in case of a
					* parallel execution is zero. This
					* overlap in a block-application of
					* the ILU in the parallel case makes
					* the preconditioner a so-called
					* additive Schwarz preconditioner.
					*/
	AdditionalData (const double       ilut_drop = 0.,
			const unsigned int ilut_fill = 0,
			const double       ilut_atol = 0.,
			const double       ilut_rtol = 1.,
			const unsigned int overlap  = 0);

                                       /**
					* This specifies the relative size
					* of elements which should be
					* dropped when forming an incomplete
					* LU decomposition with threshold.
					*/
	double ilut_drop;

                                       /**
					* This specifies the amount of
					* additional fill-in elements
					* besides the sparse matrix
					* structure. When
					* <tt>ilu_fill</tt> is large,
					* this means that many
					* fill-ins will be added, so
					* that the ILU preconditioner
					* comes closer to a (direct)
					* sparse LU
					* decomposition. Note,
					* however, that this will
					* drastically increase the
					* memory requirement,
					* especially when the
					* preconditioner is used in
					* 3D.
					*/
	unsigned int ilut_fill;

                                       /**
					* This specifies the amount of
					* an absolute perturbation
					* that will be added to the
					* diagonal of the matrix,
					* which sometimes can help to
					* get better preconditioners.
					*/
	double ilut_atol;

                                       /**
					* This specifies the factor by
					* which the diagonal of the
					* matrix will be scaled, which
					* sometimes can help to get
					* better preconditioners.
					*/
	double ilut_rtol;

                                       /**
					* This determines how large
					* the overlap of the local
					* matrix portions on each
					* processor in a parallel
					* application should be.
					*/
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

    private:
				       /**
					* This is a pointer to the
					* Ifpack data contained in
					* this preconditioner.
					*/
      Teuchos::RCP<Ifpack_Preconditioner> ifpack;
  };



/**
 * A wrapper class for a sparse direct LU decomposition on parallel
 * blocks for Trilinos matrices. When run in serial, this corresponds
 * to a direct solve on the matrix.
 *
 * The AdditionalData data structure allows to set preconditioner
 * options.
 *
 * Note that a parallel applicatoin of the block direct solve
 * preconditioner is actually a block-Jacobi preconditioner with block
 * size equal to the local matrix size. Spoken more technically, this
 * parallel operation is an <a
 * href="http://en.wikipedia.org/wiki/Additive_Schwarz_method">additive
 * Schwarz method</a> with an <em>exact solve</em> as inner solver,
 * based on the (outer) parallel partitioning.
 *
 * @ingroup TrilinosWrappers
 * @ingroup Preconditioners
 * @author Martin Kronbichler, 2008
 */
  class PreconditionBlockwiseDirect : public PreconditionBase
  {
    public:
                                       /**
                                        * Standardized data struct to
                                        * pipe additional parameters
                                        * to the preconditioner.
                                        */      
      struct AdditionalData
      {
                                       /**
					* Constructor.
					*/
	AdditionalData (const unsigned int overlap  = 0);
	

                                       /**
					* This determines how large
					* the overlap of the local
					* matrix portions on each
					* processor in a parallel
					* application should be.
					*/
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

    private:
				       /**
					* This is a pointer to the
					* Ifpack data contained in
					* this preconditioner.
					*/
      Teuchos::RCP<Ifpack_Preconditioner> ifpack;
  };
  





/**
 * A wrapper class for a Chebyshev preconditioner for Trilinos matrices.
 *
 * The AdditionalData data structure allows to set preconditioner
 * options.
 *
 * @ingroup TrilinosWrappers
 * @ingroup Preconditioners
 * @author Martin Kronbichler, 2008
 */
  class PreconditionChebyshev : public PreconditionBase
  {
    public:
                                       /**
                                        * Standardized data struct to
                                        * pipe additional parameters
                                        * to the preconditioner.
                                        */      
      struct AdditionalData
      {
                                       /**
					* Constructor.
					*/
	AdditionalData (const unsigned int degree           = 1,
			const double       max_eigenvalue   = 10.,
			const double       eigenvalue_ratio = 30.,
			const double       min_eigenvalue   = 1.,
			const double       min_diagonal     = 1e-12,
			const bool         nonzero_starting = false);
	
                                       /**
					* This determines the degree of the
					* Chebyshev polynomial. The degree
					* of the polynomial gives the number
					* of matrix-vector products to be
					* performed for one application of
					* the vmult() operation.
					*/
	unsigned int degree;

                                       /**
					* This sets the maximum eigenvalue
					* of the matrix, which needs to be
					* set properly for appropriate
					* performance of the Chebyshev
					* preconditioner.
					*/
	double max_eigenvalue;

                                       /**
					* This sets the ratio between the
					* maximum and the minimum
					* eigenvalue.
					*/
	double eigenvalue_ratio;

                                       /**
					* This sets the minimum eigenvalue,
					* which is an optional parameter
					* only used internally for checking
					* whether we use an identity matrix.
					*/
	double min_eigenvalue;

                                       /**
					* This sets a threshold below which
					* the diagonal element will not be
					* inverted in the Chebyshev
					* algorithm.
					*/
	double min_diagonal;

                                       /**
					* When this flag is set to
					* <tt>true</tt>, it enables the
					* method <tt>vmult(dst, src)</tt> to
					* use non-zero data in the vector
					* <tt>dst</tt>, appending to it the
					* Chebyshev corrections. This can be
					* useful in some situations
					* (e.g. when used for high-frequency
					* error smoothing), but not the way
					* the solver classes expect a
					* preconditioner to work (where one
					* ignores the content in
					* <tt>dst</tt> for the
					* preconditioner application). The
					* user should really know what she
					* is doing when touching this flag.
					*/
	bool nonzero_starting;
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

    private:
				       /**
					* This is a pointer to the
					* Ifpack data contained in
					* this preconditioner.
					*/
      Teuchos::RCP<Ifpack_Chebyshev> ifpack;
  };
  


/**
 * This class implements an algebraic multigrid (AMG) preconditioner
 * based on the Trilinos ML implementation, which is a black-box
 * preconditioner that works well for many PDE-based linear problems.
 * What this class does is twofold.  When the initialize() function is
 * invoked, a ML preconditioner object is created based on the matrix
 * that we want the preconditioner to be based on. A call of the
 * respective <code>vmult</code> function does call the respective
 * operation in the Trilinos package, where it is called
 * <code>ApplyInverse</code>. Use of this class is explained in the
 * @ref step_31 "step-31" tutorial program.
 *
 * There are a few pecularities in initialize(). Since the Trilinos
 * objects we want to use are heavily dependent on Epetra objects, the
 * fundamental construction routines for vectors and matrices in
 * Trilinos, we do a copy of our deal.II preconditioner matrix to a
 * Epetra matrix. This is of course not optimal, but for the time
 * being there is no direct support for our data interface.  When
 * doing this time-consuming operation, we can still profit from the
 * fact that some of the entries in the preconditioner matrix are zero
 * and hence can be neglected.
 *
 * The implementation is able to distinguish between matrices from elliptic
 * problems and convection dominated problems. We use the standard options
 * provided by Trilinos ML for elliptic problems, except that we use a
 * Chebyshev smoother instead of a symmetric Gauss-Seidel smoother.  For
 * most elliptic problems, Chebyshev provides a better damping of high
 * frequencies (in the algebraic sense) than Gauss-Seidel (SSOR), and is
 * faster (Chebyshev requires only some matrix-vector products, whereas SSOR
 * requires substitutions which are more expensive).
 *
 * @ingroup TrilinosWrappers
 * @ingroup Preconditioners
 * @author Martin Kronbichler, 2008
 */
  class PreconditionAMG : public PreconditionBase
  {
    public:

    struct AdditionalData
      {
                                       /**
					* Constructor. By default, we
					* pretend to work on elliptic
					* problems with linear finite
					* elements on a scalar
					* equation.
					*/
	AdditionalData (const bool                             elliptic = true,
			const bool                             higher_order_elements = false,
			const double                           aggregation_threshold = 1e-4,
			const std::vector<std::vector<bool> > &constant_modes = std::vector<std::vector<bool> > (1),
			const unsigned int                     smoother_sweeps = 3,
			const unsigned int                     smoother_overlap = 0,
			const bool                             output_details = false);

				       /**
					* Determines whether the AMG
					* preconditioner should be
					* optimized for elliptic
					* problems (ML option smoothed
					* aggregation SA, using a
					* Chebyshev smoother) or for
					* non-elliptic problems (ML
					* option non-symmetric
					* smoothed aggregation NSSA,
					* smoother is SSOR with
					* underrelaxation).
					*/
	bool elliptic;

				       /**
					* Determines whether the
					* matrix that the
					* preconditioner is built upon
					* is generated from linear or
					* higher-order elements.
					*/
	bool higher_order_elements;

				       /**
					* This threshold tells the AMG
					* setup how the coarsening
					* should be performed. In the
					* AMG used by ML, all points
					* that strongly couple with
					* the tentative coarse-level
					* point form one
					* aggregate. The term
					* <em>strong coupling</em> is
					* controlled by the variable
					* <tt>aggregation_threshold</tt>,
					* meaning that all elements
					* that are not smaller than
					* <tt>aggregation_threshold</tt>
					* times the diagonal element
					* do couple strongly.
					*/
	double aggregation_threshold;

				       /**
					* Specifies the constant modes
					* (near null space) of the
					* matrix. This parameter tells
					* AMG whether we work on a
					* scalar equation (where the
					* near null space only
					* consists of ones) or on a
					* vector-valued equation.
					*/
	std::vector<std::vector<bool> > constant_modes;

				       /**
					* Determines how many sweeps of the
					* smoother should be performed. When
					* the flag <tt>elliptic</tt> is set
					* to <tt>true</tt>, i.e., for
					* elliptic or almost elliptic
					* problems, the polynomial degree of
					* the Chebyshev smoother is set to
					* <tt>smoother_sweeps</tt>. In the
					* non-elliptic case,
					* <tt>smoother_sweeps</tt> sets the
					* number of SSOR relaxation sweeps
					* for post-smoothing to be
					* performed.
					*/
	unsigned int smoother_sweeps;

				       /**
					* Determines the overlap in
					* the SSOR/Chebyshev error
					* smoother when run in
					* parallel.
					*/
	unsigned int smoother_overlap;

				       /**
					* If this flag is set to
					* <tt>true</tt>, then internal
					* information from the ML
					* preconditioner is printed to
					* screen. This can be useful
					* when debugging the
					* preconditioner.
					*/
	bool output_details;
      };

				       /**
					* Let Trilinos compute a
					* multilevel hierarchy for the
					* solution of a linear system
					* with the given matrix. The
					* function uses the matrix
					* format specified in
					* TrilinosWrappers::SparseMatrix.
					*/
      void initialize (const SparseMatrix                    &matrix,
		       const AdditionalData &additional_data = AdditionalData());

				       /**
					* Let Trilinos compute a
					* multilevel hierarchy for the
					* solution of a linear system
					* with the given matrix. This
					* function takes a deal.ii
					* matrix and copies the
					* content into a Trilinos
					* matrix, so the function can
					* be considered rather
					* inefficient.
					*/
      void initialize (const ::dealii::SparseMatrix<double> &deal_ii_sparse_matrix,
		       const AdditionalData                 &additional_data = AdditionalData(),
		       const double                          drop_tolerance = 1e-13);

				       /**
					* This function can be used
				        * for a faster recalculation
				        * of the preconditioner
				        * construction when the matrix
				        * entries underlying the
				        * preconditioner have changed,
				        * but the matrix sparsity
				        * pattern has remained the
				        * same. What this function
				        * does is taking the already
				        * generated coarsening
				        * structure, computing the AMG
				        * prolongation and restriction
				        * according to a smoothed
				        * aggregation strategy and
				        * then building the whole
				        * multilevel hiearchy. This
				        * function can be considerably
				        * faster than the initialize
				        * function, since the
				        * coarsening pattern is
				        * usually the most difficult
				        * thing to do when setting up
				        * the AMG ML preconditioner.
					*/
      void reinit ();

    private:

				       /**
					* A pointer to the
					* preconditioner object.
					*/
      Teuchos::RCP<ML_Epetra::MultiLevelPreconditioner> multilevel_operator;

				       /**
					* A copy of the deal.II matrix
					* into Trilinos format.
					*/
      std_cxx0x::shared_ptr<SparseMatrix> Matrix;
  };

}

/*@}*/


DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_USE_TRILINOS

/*----------------------------   trilinos_precondition.h     ---------------------------*/

#endif
/*----------------------------   trilinos_precondition.h     ---------------------------*/
