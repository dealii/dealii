//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2004, 2005, 2006, 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <lac/trilinos_vector.h>
#include <lac/trilinos_sparse_matrix.h>

#ifdef DEAL_II_USE_TRILINOS

DEAL_II_NAMESPACE_OPEN

namespace TrilinosWrappers
{
  namespace MatrixIterators
  {
    void
    SparseMatrix::const_iterator::Accessor::
    visit_present_row ()
    {
                                       // if we are asked to visit the
                                       // past-the-end line, then simply
                                       // release all our caches and go on
                                       // with life
      if (this->a_row == matrix->m())
        {
          colnum_cache.reset ();
          value_cache.reset ();

          return;
        }
      
                                       // otherwise first flush Trilinos caches
      matrix->compress ();

                                       // get a representation of the present
                                       // row
      int ncols;
      int colnums = matrix->n();
      TrilinosScalar *values = new TrilinosScalar(colnums);
      
      int ierr;
      ierr = matrix->matrix->ExtractGlobalRowCopy((int)this->a_row, colnums,  
						  ncols, &(values[0]));
      AssertThrow (ierr == 0, ExcTrilinosError(ierr));

                                       // copy it into our caches if the line
                                       // isn't empty. if it is, then we've
                                       // done something wrong, since we
                                       // shouldn't have initialized an
                                       // iterator for an empty line (what
                                       // would it point to?)
      Assert (ncols != 0, ExcInternalError());
      colnum_cache.reset (new std::vector<unsigned int> (colnums,
                                                         colnums+ncols));
      value_cache.reset (new std::vector<TrilinosScalar> (values, values+ncols));
    }
  }


				      // The constructor is actually the
				      // only point where we have to check
                                      // whether we build a serial or
                                      // a parallel Trilinos matrix.
                                      // In the end, it does not even
                                      // matter how many threads there
                                      // are, but only if we use an
                                      // MPI compiler or a standard 
				      // compiler. So, one thread on
                                      // an MPI compiler will still get
                                      // a parallel interface.
  SparseMatrix::SparseMatrix ()
                  :
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
		  row_map (0, 0, Epetra_MpiComm(MPI_COMM_WORLD)),
#else
		  row_map (0, 0, Epetra_SerialComm()),
#endif
		  col_map (row_map),
		  matrix (std::auto_ptr<Epetra_FECrsMatrix>
				 (new Epetra_FECrsMatrix(Copy, row_map, 0))),
                  last_action (Insert)
  {}

  SparseMatrix::SparseMatrix (const Epetra_Map  &InputMap,
			      const unsigned int n_max_entries_per_row)
                  :
		  row_map (InputMap),
		  col_map (row_map),
		  matrix (std::auto_ptr<Epetra_FECrsMatrix>
				 (new Epetra_FECrsMatrix(Copy, row_map, 
					int(n_max_entries_per_row), false))),
                  last_action (Insert)
  {}

  SparseMatrix::SparseMatrix (const Epetra_Map  &InputMap,
			      const std::vector<unsigned int> &n_entries_per_row)
                  :
		  row_map (InputMap),
		  col_map (row_map),
		  matrix (std::auto_ptr<Epetra_FECrsMatrix>
		    (new Epetra_FECrsMatrix(Copy, row_map, 
		      (int*)const_cast<unsigned int*>(&(n_entries_per_row[0])),
		      false))),
                  last_action (Insert)
  {}

  SparseMatrix::SparseMatrix (const Epetra_Map  &InputRowMap,
			      const Epetra_Map  &InputColMap,
			      const unsigned int n_max_entries_per_row)
                  :
		  row_map (InputRowMap),
		  col_map (InputColMap),
		  matrix (std::auto_ptr<Epetra_FECrsMatrix>
				 (new Epetra_FECrsMatrix(Copy, row_map, col_map, 
					int(n_max_entries_per_row), false))),
                  last_action (Insert)
  {}

  SparseMatrix::SparseMatrix (const Epetra_Map  &InputRowMap,
			      const Epetra_Map  &InputColMap,
			      const std::vector<unsigned int> &n_entries_per_row)
                  :
		  row_map (InputRowMap),
		  col_map (InputColMap),
		  matrix (std::auto_ptr<Epetra_FECrsMatrix>
		    (new Epetra_FECrsMatrix(Copy, row_map, col_map, 
		      (int*)const_cast<unsigned int*>(&(n_entries_per_row[0])),
		      false))),
                  last_action (Insert)
  {}



  SparseMatrix::~SparseMatrix ()
  {
  }


  
  void
  SparseMatrix::reinit (const CompressedSparsityPattern &sparsity_pattern,
		        const unsigned int               n_max_entries_per_row)
  {

    unsigned int n_rows = sparsity_pattern.n_rows();

    Assert (matrix->NumGlobalRows() == (int)sparsity_pattern.n_rows(),
	    ExcDimensionMismatch (matrix->NumGlobalRows(),
				  sparsity_pattern.n_rows()));

    std::vector<double> values(n_max_entries_per_row, 0.);
    std::vector<int>    row_indices(n_max_entries_per_row);

    for (unsigned int row=0; row<n_rows; ++row)
      {
	const int row_length = sparsity_pattern.row_length(row);
	row_indices.resize (row_length, 0);
        values.resize (row_length, 0.);
	
	for (int col=0; col< row_length; ++col)
	  row_indices[col] = sparsity_pattern.column_number (row, col);
	
	matrix->InsertGlobalValues(row, row_length,
                                   &values[0], &row_indices[0]);
      }

				  // In the end, the matrix needs to
				  // be compressed in order to be
				  // really ready. However, that is
				  // a collective operation, so it
				  // has to be called on all processes
				  // by the user, whereas this function
				  // should only be used on one processor
				  // since our sparsity pattern data
				  // types are all serial.
  }



  void
  SparseMatrix::reinit (const CompressedSparsityPattern &sparsity_pattern)
  {
    unsigned int n_rows = sparsity_pattern.n_rows();
    
    Assert (matrix->NumGlobalRows() == (int)sparsity_pattern.n_rows(),
	    ExcDimensionMismatch (matrix->NumGlobalRows(),
				  sparsity_pattern.n_rows()));
    
				 // Trilinos seems to have a bug for
				 // rectangular matrices at this point,
				 // so do not check for consistent 
				 // column numbers here.
//    Assert (matrix->NumGlobalCols() == (int)sparsity_pattern.n_cols(),
//	    ExcDimensionMismatch (matrix->NumGlobalCols(),
//				  sparsity_pattern.n_cols()));

    std::vector<int> n_entries_per_row(n_rows);
    
    for (unsigned int row=0; row<n_rows; ++row)
      n_entries_per_row[(int)row] = sparsity_pattern.row_length(row);
    
    const unsigned int n_max_entries_per_row = *std::max_element (
		    &n_entries_per_row[0], &n_entries_per_row[n_rows-1]);

    reinit (sparsity_pattern, n_max_entries_per_row);
  }



  void
  SparseMatrix::reinit (const Epetra_Map                &input_map,
		        const CompressedSparsityPattern &sparsity_pattern)
  {

    unsigned int n_rows = sparsity_pattern.n_rows();
    matrix.reset();
    row_map = input_map;
    col_map = row_map;
    
    Assert (input_map.NumGlobalElements() == (int)sparsity_pattern.n_rows(),
	    ExcDimensionMismatch (input_map.NumGlobalElements(),
				  sparsity_pattern.n_rows()));
    Assert (input_map.NumGlobalElements() == (int)sparsity_pattern.n_cols(),
	    ExcDimensionMismatch (input_map.NumGlobalElements(),
				  sparsity_pattern.n_cols()));
    
    std::vector<int> n_entries_per_row(n_rows);
    
    for (unsigned int row=0; row<n_rows; ++row)
      n_entries_per_row[(int)row] = sparsity_pattern.row_length(row);
    
    matrix = std::auto_ptr<Epetra_FECrsMatrix> (new Epetra_FECrsMatrix (
				Copy, row_map, &n_entries_per_row[0], true));
    
    const unsigned int n_max_entries_per_row = *std::max_element (
		    &n_entries_per_row[0], &n_entries_per_row[n_rows-1]);

    reinit (sparsity_pattern, n_max_entries_per_row);
  }



  void
  SparseMatrix::reinit (const Epetra_Map                &input_row_map,
			const Epetra_Map                &input_col_map,
		        const CompressedSparsityPattern &sparsity_pattern)
  {

    unsigned int n_rows = sparsity_pattern.n_rows();
    matrix.reset();
    row_map = input_row_map;
    col_map = input_col_map;
    
    Assert (input_row_map.NumGlobalElements() == (int)sparsity_pattern.n_rows(),
	    ExcDimensionMismatch (input_row_map.NumGlobalElements(),
				  sparsity_pattern.n_rows()));
    Assert (input_col_map.NumGlobalElements() == (int)sparsity_pattern.n_cols(),
	    ExcDimensionMismatch (input_col_map.NumGlobalElements(),
				  sparsity_pattern.n_cols()));
    
    std::vector<int> n_entries_per_row(n_rows);
    
    for (unsigned int row=0; row<n_rows; ++row)
      n_entries_per_row[(int)row] = sparsity_pattern.row_length(row);
    
    matrix = std::auto_ptr<Epetra_FECrsMatrix>
	      (new Epetra_FECrsMatrix(Copy, row_map, col_map,
	                              &n_entries_per_row[0], true));
    
    const unsigned int n_max_entries_per_row = *std::max_element (
		    &n_entries_per_row[0], &n_entries_per_row[n_rows-1]);

    reinit (sparsity_pattern, n_max_entries_per_row);
  }



  void
  SparseMatrix::clear ()
  {
                                     // When we clear the matrix,
				     // reset the pointer and 
				     // generate an empty matrix.
    matrix.reset();
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
    row_map = Epetra_Map (0,0,Epetra_MpiComm(MPI_COMM_WORLD)),
#else
    row_map = Epetra_Map (0,0,Epetra_SerialComm()),
#endif

    matrix = std::auto_ptr<Epetra_FECrsMatrix> 
	      (new Epetra_FECrsMatrix(Copy, row_map, 0));
  }



  void
  SparseMatrix::compress ()
  {
                                     // flush buffers
    int ierr;
    if (row_map.SameAs(col_map))
      ierr = matrix->GlobalAssemble ();
    else
      ierr = matrix->GlobalAssemble (col_map, row_map);
    
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    ierr = matrix->OptimizeStorage ();
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



  SparseMatrix &
  SparseMatrix::operator = (const double d)
  {
    Assert (d==0, ExcScalarAssignmentOnlyForZeroValue());

    compress ();

    const int ierr = matrix->PutScalar(d);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    return *this;
  }



  void
  SparseMatrix::set (const unsigned int   i,
		     const unsigned int   j,
		     const TrilinosScalar value)
  {

    Assert (numbers::is_finite(value),
	    ExcMessage("The given value is not finite but either "
		       "infinite or Not A Number (NaN)"));

    if (last_action == Add)
      {
        int ierr;
	if (row_map.SameAs(col_map))
	  ierr = matrix->GlobalAssemble (false);
	else
	  ierr = matrix->GlobalAssemble(col_map, row_map, false);
	
        AssertThrow (ierr == 0, ExcTrilinosError(ierr));
	last_action = Insert;
      }

    int trilinos_i = i;
    int trilinos_j = j;

    const int ierr = matrix->ReplaceGlobalValues (trilinos_i, 1,
						  const_cast<double*>(&value), 
						  &trilinos_j);

    AssertThrow (ierr <= 0, ExcAccessToNonPresentElement(i,j));
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



  void
  SparseMatrix::add (const unsigned int   i,
		     const unsigned int   j,
		     const TrilinosScalar value)
  {

    Assert (numbers::is_finite(value), 
	    ExcMessage("The given value is not finite but either "
		       "infinite or Not A Number (NaN)"));

    if (last_action == Insert)
      {
        int ierr;
	if (row_map.SameAs(col_map))
	  ierr = matrix->GlobalAssemble (false);
	else
	  ierr = matrix->GlobalAssemble(col_map, row_map, false);
        
	AssertThrow (ierr == 0, ExcTrilinosError(ierr));

        last_action = Add;
      }

                                     // we have to do above actions in any
                                     // case to be consistent with the MPI
                                     // communication model (see the
                                     // comments in the documentation of
                                     // TrilinosWrappers::MPI::Vector), but we
                                     // can save some work if the addend is
                                     // zero
    if (value == 0)
      return;

    int trilinos_i = i;
    int trilinos_j = j;

    const int ierr = matrix->SumIntoGlobalValues (trilinos_i, 1, 
						  const_cast<double*>(&value), 
						  &trilinos_j);

    AssertThrow (ierr <= 0, ExcAccessToNonPresentElement(i,j));
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



  TrilinosScalar
  SparseMatrix::el (const unsigned int i,
		    const unsigned int j) const
  {
                                      // Extract local indices in
                                      // the matrix.
    int trilinos_i = matrix->LRID(i), trilinos_j = matrix->LRID(j);
    TrilinosScalar value = 0.;
    
                                      // If the data is not on the
                                      // present processor, we can't
                                      // continue.
    if ((trilinos_i == -1 ) || (trilinos_j == -1))
      {
        Assert (false, ExcAccessToNonLocalElement(i, j, local_range().first,
						  local_range().second));
      }
    else
    {
				      // Check whether the matrix 
				      // already is transformed to
				      // local indices.
      if (!matrix->Filled())
	matrix->FillComplete(true);

				      // Prepare pointers for extraction
				      // of a view of the row.
      int nnz_present = matrix->NumMyEntries(trilinos_i);
      int nnz_extracted;
      int *col_indices;
      TrilinosScalar *values;

				      // Generate the view and make
				      // sure that we have not generated
				      // an error.
      int ierr = matrix->ExtractMyRowView(trilinos_i, nnz_extracted,
					  values, col_indices);
      Assert (ierr==0, ExcTrilinosError(ierr));

      Assert (nnz_present == nnz_extracted,
	      ExcDimensionMismatch(nnz_present, nnz_extracted));

				      // Search the index where we
				      // look for the value, and then
				      // finally get it.
      int* index = std::find(&col_indices[0],&col_indices[0] + nnz_present,
			     trilinos_j);

      int position;
      if (!index)
	value = 0;
      else
	{
	  position = (int)(index - &(col_indices[0]));
	  value = values[position];
	}
    }

    return value;
  }



  TrilinosScalar
  SparseMatrix::diag_element (const unsigned int i) const
  {
    Assert (m() == n(), ExcNotQuadratic());

                                     // this doesn't seem to work any
                                     // different than any other element
    return el(i,i);
  }



  unsigned int
  SparseMatrix::m () const
  {
    int n_rows = matrix->NumGlobalRows();

    return n_rows;
  }



  unsigned int
  SparseMatrix::n () const
  {
    int n_cols = matrix -> NumGlobalCols();
    return n_cols;
  }



  unsigned int
  SparseMatrix::local_size () const
  {
    int n_rows = matrix -> NumMyRows();

    return n_rows;
  }



  std::pair<unsigned int, unsigned int>
  SparseMatrix::local_range () const
  {
    int begin, end;
    begin = matrix->RowMap().MinMyGID();
    end = matrix->RowMap().MaxMyGID();
    
    return std::make_pair (begin, end);
  }



  unsigned int
  SparseMatrix::n_nonzero_elements () const
  {
    int nnz = matrix->NumGlobalNonzeros();

    return static_cast<unsigned int>(nnz);
  }



  unsigned int
  SparseMatrix::row_length (const unsigned int row) const
  {
    Assert (row < m(), ExcInternalError());

				 // get a representation of the present
				 // row
    int ncols = -1;
    int local_row = matrix->RowMap().LID(row);

				 // on the processor who owns this
				 // row, we'll have a non-negative
				 // value.
    if (local_row >= 0)
      {
	int ierr = matrix->NumMyRowEntries (local_row, ncols);
	AssertThrow (ierr == 0, ExcTrilinosError(ierr));
      }

    return ncols;
  }


  TrilinosScalar
  SparseMatrix::l1_norm () const
  {
    if (!matrix->Filled())
      matrix->FillComplete();

    TrilinosScalar result = matrix->NormOne();

    return result;
  }
  
  

  TrilinosScalar
  SparseMatrix::linfty_norm () const
  {
    if (!matrix->Filled())
      matrix->FillComplete();

    TrilinosScalar result = matrix->NormInf();

    return result;
  }



  TrilinosScalar
  SparseMatrix::frobenius_norm () const
  {
    if (!matrix->Filled())
      matrix->FillComplete();

    TrilinosScalar result = matrix->NormFrobenius();

    return result;
  }



  SparseMatrix &
  SparseMatrix::operator *= (const TrilinosScalar a)
  {
    const int ierr = matrix->Scale (a);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    return *this;
  }



  SparseMatrix &
  SparseMatrix::operator /= (const TrilinosScalar a)
  {
    const TrilinosScalar factor = 1./a;

    const int ierr = matrix->Scale (factor);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    return *this;
  }


  void
  SparseMatrix::vmult (Vector       &dst,
		       const Vector &src) const
  {
    Assert (&src != &dst, ExcSourceEqualsDestination());

    if (!matrix->Filled())
      matrix->FillComplete();

    const int ierr = matrix->Multiply (false, *(src.vector), *(dst.vector));
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



  void
  SparseMatrix::Tvmult (Vector       &dst,
			const Vector &src) const
  {
    Assert (&src != &dst, ExcSourceEqualsDestination());

    if (!matrix->Filled())
      matrix->FillComplete();

    const int ierr = matrix->Multiply (true, *(src.vector), *(dst.vector));
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



  void
  SparseMatrix::vmult_add (Vector       &dst,
			   const Vector &src) const
  {
    Assert (&src != &dst, ExcSourceEqualsDestination());

    Vector tmp = dst;
    vmult (dst, src);
    dst += tmp;
  }



  void
  SparseMatrix::Tvmult_add (Vector       &dst,
			    const Vector &src) const
  {
    Assert (&src != &dst, ExcSourceEqualsDestination());

    Vector tmp = dst;
    vmult (dst, src);
    dst += tmp;
  }



  TrilinosScalar
  SparseMatrix::matrix_norm_square (const Vector &v) const
  {
    Vector tmp(v.map);
    vmult (tmp, v);
    return tmp*v;
  }



  TrilinosScalar
  SparseMatrix::matrix_scalar_product (const Vector &u,
				       const Vector &v) const
  {
    Vector tmp(v.map);
    vmult (tmp, v);
    return u*tmp;
  }



  TrilinosScalar
  SparseMatrix::residual (Vector       &dst,
			  const Vector &x,
			  const Vector &b) const
  {
    vmult (dst, x);
    dst -= b;
    dst *= -1;

    return dst.l2_norm();
  }



				  // TODO: Currently this only flips
                                  // a flag that tells Trilinos that
                                  // any application should be done with
                                  // the transpose. However, the matrix
                                  // structure is not reset.
  void
  SparseMatrix::transpose () 
  {
    int ierr;

    if (!matrix->UseTranspose())
      {
        ierr = matrix->SetUseTranspose (true);
        AssertThrow (ierr == 0, ExcTrilinosError(ierr));
      }
    else
      {
        ierr = matrix->SetUseTranspose (false);
        AssertThrow (ierr == 0, ExcTrilinosError(ierr));
      }
  }



  bool
  SparseMatrix::is_symmetric (const double tolerance) 
  {
    //bool truth;
    if (tolerance == 0)
      Assert (false, ExcNotImplemented());

    return false;
  }  



  bool
  SparseMatrix::is_hermitian () 
  {
    //bool truth;

    Assert (false, ExcNotImplemented());
    return false;
  }  

  void
  SparseMatrix::write_ascii ()
  {
    Assert (false, ExcNotImplemented());
  }

}

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_USE_TRILINOS
