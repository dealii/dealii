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


#include <lac/trilinos_precondition_amg.h>

#ifdef DEAL_II_USE_TRILINOS

#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>
#include <ml_include.h>
#include <ml_MultiLevelPreconditioner.h>

DEAL_II_NAMESPACE_OPEN

namespace TrilinosWrappers
{

  PreconditionAMG::PreconditionAMG () 
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
    :
      communicator (MPI_COMM_WORLD)
#endif
  {}


  
  PreconditionAMG::~PreconditionAMG ()
  {}



  void
  PreconditionAMG::
  initialize (const SparseMatrix        &matrix,
	      const std::vector<double> &null_space,
	      const unsigned int         null_space_dimension,
	      const bool                 elliptic,
	      const bool                 higher_order_elements,
	      const bool                 output_details)
  {
    const unsigned int n_rows = matrix.m();
    
    if (!Matrix)
      {
	Matrix = boost::shared_ptr<SparseMatrix>
			  (const_cast<SparseMatrix*>(&matrix));
	Map = boost::shared_ptr<Epetra_Map>
			  (const_cast<Epetra_Map*>(&matrix.row_map));
      }

				     // Build the AMG preconditioner.
    Teuchos::ParameterList parameter_list;
  
    if (elliptic)
      {
	ML_Epetra::SetDefaults("SA",parameter_list);
	parameter_list.set("smoother: type", "Chebyshev");
	parameter_list.set("smoother: sweeps", 4);
      }
    else
      {
	ML_Epetra::SetDefaults("NSSA",parameter_list);
	parameter_list.set("aggregation: type", "Uncoupled");
	parameter_list.set("aggregation: block scaling", true);
      }
  
    if (output_details)
      parameter_list.set("ML output", 10);
    else
      parameter_list.set("ML output", 0);
  
    if (higher_order_elements)
      parameter_list.set("aggregation: type", "MIS");
  
    Assert (n_rows * null_space_dimension == null_space.size(),
	    ExcDimensionMismatch(n_rows * null_space_dimension,
				 null_space.size()));
  
    if (null_space_dimension > 1)
      {
	parameter_list.set("null space: type", "pre-computed");
	parameter_list.set("null space: dimension", int(null_space_dimension));
	parameter_list.set("null space: vectors", (double *)&null_space[0]);
      }

    multigrid_operator = boost::shared_ptr<ML_Epetra::MultiLevelPreconditioner>
			 (new ML_Epetra::MultiLevelPreconditioner(
				      *matrix.matrix, parameter_list, true));

    if (output_details)
      multigrid_operator->PrintUnused(0);
  }



  void
  PreconditionAMG::
  initialize (const dealii::SparseMatrix<double> &deal_ii_sparse_matrix,
	      const std::vector<double>  &null_space,
	      const unsigned int          null_space_dimension,
	      const bool                  elliptic,
	      const bool                  higher_order_elements,
	      const bool                  output_details)
  {
    const unsigned int n_rows = deal_ii_sparse_matrix.m();
  
				     // Init Epetra Matrix, avoid 
				     // storing the nonzero elements.

    Map.reset (new Epetra_Map(n_rows, 0, communicator));

    Matrix.reset();
    Matrix = boost::shared_ptr<SparseMatrix> (new SparseMatrix());

    Matrix->reinit (*Map, deal_ii_sparse_matrix);
    Matrix->compress();

    initialize (*Matrix, null_space, null_space_dimension, elliptic,
	        higher_order_elements, output_details);
  }
  
  
  void PreconditionAMG::
  reinit ()
  {
    multigrid_operator->ReComputePreconditioner();
  }
  
  
  
  void PreconditionAMG::vmult (Vector        &dst,
			       const Vector  &src) const
  {
    const int ierr = multigrid_operator->ApplyInverse (*src.vector, *dst.vector);

    Assert (ierr == 0, ExcTrilinosError(ierr));
  }


				   // For the implementation of the
				   // <code>vmult</code> function we
				   // note that invoking a call of 
				   // the Trilinos preconditioner 
				   // requires us to use Epetra vectors
				   // as well. It is faster
				   // to provide a view, i.e., feed 
				   // Trilinos with a pointer to the
				   // data, so we avoid copying the
				   // content of the vectors during
				   // the iteration. In the declaration
				   // of the right hand side, we need
				   // to cast the source vector (that
				   // is <code>const</code> in all deal.II 
				   // calls) to non-constant value, as
				   // this is the way Trilinos wants to
				   // have them.
  void PreconditionAMG::vmult (dealii::Vector<double>       &dst,
			       const dealii::Vector<double> &src) const
  {
    Assert (Map->SameAs(Matrix->row_map),
	    ExcMessage("The sparse matrix given to the preconditioner uses "
		       "a map that is not compatible. Check ML preconditioner "
		       "setup."));
    Assert (Map->SameAs(Matrix->col_map),
	    ExcMessage("The sparse matrix given to the preconditioner uses "
		       "a map that is not compatible. Check ML preconditioner "
		       "setup."));
    
    Epetra_Vector LHS (View, *Map, dst.begin());
    Epetra_Vector RHS (View, *Map, const_cast<double*>(src.begin()));
  
    const int res = multigrid_operator->ApplyInverse (RHS, LHS);
  
    Assert (res == 0,
	    ExcMessage ("Trilinos AMG MultiLevel preconditioner returned "
			"with an error!"));
  }

}

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_USE_TRILINOS
