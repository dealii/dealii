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
#ifndef __deal2__trilinos_precondition_amg_h
#define __deal2__trilinos_precondition_amg_h


#include <base/config.h>
#include <lac/exceptions.h>
#include <lac/trilinos_vector_base.h>
#include <lac/trilinos_sparse_matrix.h>
#include <lac/vector.h>
#include <lac/sparse_matrix.h>

#include <boost/shared_ptr.hpp>


#ifdef DEAL_II_USE_TRILINOS
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
#  include <Epetra_MpiComm.h>
#else
#  include <Epetra_SerialComm.h>
#endif
#include <Teuchos_RefCountPtr.hpp>

// some forward declarations
namespace ML_Epetra
{
  class MultiLevelPreconditioner;
}
class Epetra_Map;
class Epetra_CrsMatrix;



DEAL_II_NAMESPACE_OPEN

/*! @addtogroup TrilinosWrappers
 *@{
 */

namespace TrilinosWrappers
{

					// Forward declarations.
  class SolverBase;

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
 * The implementation is able to distinguish between matrices from
 * elliptic problems and convection dominated problems. We use the
 * standard options provided by Trilinos ML for elliptic problems,
 * except that we use a Chebyshev smoother instead of a symmetric
 * Gauss-Seidel smoother.  For most elliptic problems, Chebyshev
 * provides a better damping of high frequencies (in the algebraic
 * sense) than Gauss-Seidel (SSOR).
 *
 * @ingroup TrilinosWrappers
 * @ingroup Preconditioners
 * @author Martin Kronbichler, 2008
 */
  class PreconditionAMG : public Subscriptor
  {
    public:
				       /**
					* Constructor.
					*/
      PreconditionAMG ();

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
		       const bool                             elliptic = true,
		       const bool                             higher_order_elements = false,
		       const double                           aggregation_threshold = 1e-4,
		       const std::vector<std::vector<bool> > &null_space = std::vector<std::vector<bool> > (1),
		       const bool                             output_details = false);

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
      void initialize (const ::dealii::SparseMatrix<double>  &deal_ii_sparse_matrix,
		       const bool                             elliptic = true,
		       const bool                             higher_order_elements = false,
		       const double                           aggregation_threshold = 1e-4,
		       const std::vector<std::vector<bool> > &null_space = std::vector<std::vector<bool> > (),
		       const bool                            output_details = false);

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
				        * does is to take the already
				        * generated coarsening
				        * structure, compute the AMG
				        * prolongation and restriction
				        * according to a smoothed
				        * aggregation strategy and
				        * then builds the whole
				        * multilevel hiearchy. This
				        * function can be considerably
				        * faster in that case, since
				        * the coarsening pattern is
				        * usually the most difficult
				        * thing to do when setting up
				        * the AMG ML preconditioner.
					*/
      void reinit ();
      
				       /**
					* Apply the preconditioner.
					*/
      void vmult (VectorBase       &dst,
		  const VectorBase &src) const;

				       /**
					* Do the same as before, but 
				        * now use deal.II vectors instead
				        * of the ones provided in the
				        * Trilinos wrapper class.
					*/
      void vmult (dealii::Vector<double>       &dst,
		  const dealii::Vector<double> &src) const;

                                       /**
					* Exception.
					*/
      DeclException1 (ExcNonMatchingMaps,
		      std::string,
		      << "The sparse matrix that the preconditioner is based "
		      << "on a map that is not compatible to the one in vector"
		      << arg1
		      << ". Check preconditioner and matrix setup.");

    private:

				       /**
					* A pointer to the
					* preconditioner object.
					*/
      Teuchos::RefCountPtr<ML_Epetra::MultiLevelPreconditioner> 
	multilevel_operator;

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
      boost::shared_ptr<Epetra_Map>   Map;

				       /**
					* A copy of the deal.II matrix
					* into Trilinos format.
					*/
      boost::shared_ptr<SparseMatrix> Matrix;

      friend class SolverBase;
  };
}

/*@}*/


DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_USE_TRILINOS

/*----------------------------   trilinos_precondition_amg.h     ---------------------------*/

#endif
/*----------------------------   trilinos_precondition_amg.h     ---------------------------*/
