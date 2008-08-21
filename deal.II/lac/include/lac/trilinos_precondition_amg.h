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
#include <lac/trilinos_vector.h>
#include <lac/trilinos_sparse_matrix.h>

#include <boost/shared_ptr.hpp>


#ifdef DEAL_II_USE_TRILINOS

// some forward declarations
namespace ML_Epetra
{
  class MultiLevelPreconditioner;
}
class Epetra_Map;
class Epetra_CrsMatrix;



DEAL_II_NAMESPACE_OPEN

namespace TrilinosWrappers
{
  
/**  
 * This class implements an algebraic multigrid (AMG) preconditioner
 * based on the Trilinos ML implementation.  What this class does is
 * twofold.  When the initialize() function is invoked, a ML
 * preconditioner object is created based on the matrix
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
					* Destructor.
					*/
      ~PreconditionAMG ();

				       /**
					* Let Trilinos compute a
					* multilevel hierarchy for the
					* solution of a linear system
					* with the given matrix.
					*/
      void initialize (const dealii::SparseMatrix<double> &matrix,
		       const std::vector<double>  &null_space,
		       const unsigned int          null_space_dimension,
		       const bool                  higher_order_elements,
		       const bool                  elliptic,
		       const bool                  output_details,
		       const double                drop_tolerance = 1e-13);

				       /**
					* Multiply the source vector
					* with the preconditioner
					* represented by this object,
					* and return it in the
					* destination vector.
					*/
      void vmult (dealii::Vector<double>       &dst,
		  const dealii::Vector<double> &src) const;

    private:

				       /**
					* Pointer to the Trilinos AMG object.
					*/
      boost::shared_ptr<ML_Epetra::MultiLevelPreconditioner> multigrid_operator;

      Epetra_SerialComm  communicator;
      boost::shared_ptr<Epetra_Map>       Map;

				       /**
					* A copy of the deal.II matrix
					* into Trilinos format.
					*/
      boost::shared_ptr<Epetra_CrsMatrix> Matrix;
  };
}



DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_USE_TRILINOS

/*----------------------------   trilinos_precondition_amg_base.h     ---------------------------*/

#endif
/*----------------------------   trilinos_precondition_amg_base.h     ---------------------------*/
