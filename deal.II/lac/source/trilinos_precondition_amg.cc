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
#include <lac/vector.h>
#include <lac/sparse_matrix.h>

#include <Epetra_SerialComm.h>
#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>
#include <ml_include.h>
#include <ml_MultiLevelPreconditioner.h>



#ifdef DEAL_II_USE_TRILINOS

DEAL_II_NAMESPACE_OPEN

namespace TrilinosWrappers
{

  PreconditionAMG::PreconditionAMG ()
  {}


  PreconditionAMG::~PreconditionAMG ()
  {}

  
  void PreconditionAMG::initialize (
    const dealii::SparseMatrix<double> &matrix,
    const std::vector<double>  &null_space,
    const unsigned int          null_space_dimension,
    const bool                  elliptic,
    const bool                  higher_order_elements,
    const bool                  output_details,
    const double                drop_tolerance
  )
  {
    Assert (drop_tolerance >= 0,
	    ExcMessage ("Drop tolerance must be a non-negative number."));
  
    const unsigned int n_rows = matrix.m();
    const SparsityPattern *sparsity_pattern = &(matrix.get_sparsity_pattern());
  
				     // Init Epetra Matrix, avoid 
				     // storing the nonzero elements.
    {
      Map.reset (new Epetra_Map(n_rows, 0, communicator));
    
      std::vector<int> row_lengths (n_rows);
      for (dealii::SparseMatrix<double>::const_iterator p = matrix.begin();
	   p != matrix.end(); ++p)
	if (std::abs(p->value()) > drop_tolerance)
	  ++row_lengths[p->row()];
  
      Matrix.reset (new Epetra_CrsMatrix(Copy, *Map, &row_lengths[0], true));
  
      const unsigned int max_nonzero_entries
	= *std::max_element (row_lengths.begin(), row_lengths.end());
  
      std::vector<double> values(max_nonzero_entries, 0);
      std::vector<int> row_indices(max_nonzero_entries);
  
      for (unsigned int row=0; row<n_rows; ++row)
	{
	  unsigned int index = 0;
	  for (dealii::SparseMatrix<double>::const_iterator p = matrix.begin(row);
	       p != matrix.end(row); ++p)
	    if (std::abs(p->value()) > drop_tolerance)
	      {
		row_indices[index] = p->column();
		values[index]      = p->value();
		++index;
	      }

	  Assert (index == static_cast<unsigned int>(row_lengths[row]),
		  ExcMessage("Filtering out zeros could not "
			     "be successfully finished!"));
  
	  Matrix->InsertGlobalValues(row, row_lengths[row],
				     &values[0], &row_indices[0]);
	}
      
      Matrix->FillComplete();
    }
  
				     // Build the AMG preconditioner.
    Teuchos::ParameterList parameter_list;
  
				     // The implementation is able
				     // to distinguish between
				     // matrices from elliptic problems
				     // and convection dominated 
				     // problems. We use the standard
				     // options for elliptic problems,
				     // except that we use a 
				     // Chebyshev smoother instead
				     // of a symmetric Gauss-Seidel
				     // smoother. For most elliptic 
				     // problems, Chebyshev is better
				     // than Gauss-Seidel (SSOR).
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
			 (new ML_Epetra::MultiLevelPreconditioner(*Matrix, parameter_list, true));

    if (output_details)
      multigrid_operator->PrintUnused(0);
  }

				   // For the implementation of the
				   // <code>vmult</code> function we
				   // note that invoking a call of 
				   // the Trilinos preconditioner 
				   // requires us to use Epetra vectors
				   // as well. Luckily, it is sufficient
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
