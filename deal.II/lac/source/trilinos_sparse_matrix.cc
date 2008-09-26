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

				  // copy it into our caches if the
				  // line isn't empty. if it is, then
				  // we've done something wrong, since
				  // we shouldn't have initialized an
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
				  // whether we build a serial or a
				  // parallel Trilinos matrix.
				  // Actually, it does not even matter
				  // how many threads there are, but
				  // only if we use an MPI compiler or
				  // a standard compiler. So, even one
				  // thread on a configuration with
				  // MPI will still get a parallel
				  // interface.
  SparseMatrix::SparseMatrix ()
		  :
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
                  row_map (0, 0, Epetra_MpiComm(MPI_COMM_WORLD)),
#else
                  row_map (0, 0, Epetra_SerialComm()),
#endif
		  col_map (row_map),
		  last_action (Zero),
		  compressed (true),
		  matrix (std::auto_ptr<Epetra_FECrsMatrix>
				(new Epetra_FECrsMatrix(Copy, row_map, 0)))
  {}

  SparseMatrix::SparseMatrix (const Epetra_Map  &InputMap,
			      const unsigned int n_max_entries_per_row)
		  :
                  row_map (InputMap),
		  col_map (row_map),
		  last_action (Zero),
		  compressed (true),
		  matrix (std::auto_ptr<Epetra_FECrsMatrix>
				(new Epetra_FECrsMatrix(Copy, row_map, 
					int(n_max_entries_per_row), false)))
  {}

  SparseMatrix::SparseMatrix (const Epetra_Map                &InputMap,
			      const std::vector<unsigned int> &n_entries_per_row)
		  :
                  row_map (InputMap),
		  col_map (row_map),
		  last_action (Zero),
		  compressed (true),
		  matrix (std::auto_ptr<Epetra_FECrsMatrix>
		    (new Epetra_FECrsMatrix(Copy, row_map, 
		      (int*)const_cast<unsigned int*>(&(n_entries_per_row[0])),
					    false)))
  {}

  SparseMatrix::SparseMatrix (const Epetra_Map  &InputRowMap,
			      const Epetra_Map  &InputColMap,
			      const unsigned int n_max_entries_per_row)
		  :
                  row_map (InputRowMap),
                  col_map (InputColMap),
		  last_action (Zero),
		  compressed (true),
		  matrix (std::auto_ptr<Epetra_FECrsMatrix>
				(new Epetra_FECrsMatrix(Copy, row_map, col_map, 
					int(n_max_entries_per_row), false)))
  {}

  SparseMatrix::SparseMatrix (const Epetra_Map                &InputRowMap,
			      const Epetra_Map                &InputColMap,
			      const std::vector<unsigned int> &n_entries_per_row)
		  :
                  row_map (InputRowMap),
                  col_map (InputColMap),
		  last_action (Zero),
		  compressed (true),
		  matrix (std::auto_ptr<Epetra_FECrsMatrix>
		    (new Epetra_FECrsMatrix(Copy, row_map, col_map, 
		      (int*)const_cast<unsigned int*>(&(n_entries_per_row[0])),
					    false)))
  {}



  SparseMatrix::~SparseMatrix ()
  {
  }



  SparseMatrix &
  SparseMatrix::copy_from (const SparseMatrix &m)
  {
    *matrix = *m.matrix;
    return *this;
  }
  


  void
  SparseMatrix::reinit (const SparsityPattern &sparsity_pattern)
  {
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
    Epetra_MpiComm    trilinos_communicator (MPI_COMM_WORLD);
#else
    Epetra_SerialComm trilinos_communicator;
#endif

    const Epetra_Map rows (sparsity_pattern.n_rows(),
			   0,
			   trilinos_communicator);
    const Epetra_Map columns (sparsity_pattern.n_cols(),
			      0,
			      trilinos_communicator);

    reinit (rows, columns, sparsity_pattern);
  }



  void
  SparseMatrix::reinit (const Epetra_Map       &input_map,
			const SparsityPattern  &sparsity_pattern)
  {
    reinit (input_map, input_map, sparsity_pattern);
  }



  void
  SparseMatrix::reinit (const Epetra_Map       &input_row_map,
			const Epetra_Map       &input_col_map,
			const SparsityPattern  &sparsity_pattern)
  {
    matrix.reset();

    unsigned int n_rows = sparsity_pattern.n_rows();

    if (input_row_map.Comm().MyPID() == 0)
      {
	Assert (input_row_map.NumGlobalElements() == (int)sparsity_pattern.n_rows(),
		ExcDimensionMismatch (input_row_map.NumGlobalElements(),
				      sparsity_pattern.n_rows()));
	Assert (input_col_map.NumGlobalElements() == (int)sparsity_pattern.n_cols(),
		ExcDimensionMismatch (input_col_map.NumGlobalElements(),
				      sparsity_pattern.n_cols()));
      }

    row_map = input_row_map;
    col_map = input_col_map;

    std::vector<int> n_entries_per_row(n_rows);

    for (unsigned int row=0; row<n_rows; ++row)
      n_entries_per_row[row] = sparsity_pattern.row_length(row);

				  // TODO: There seems to be problem
				  // in Epetra when a matrix is
				  // initialized with both row and
				  // column map. Maybe find something
				  // more out about this... It could
				  // be related to the bug #4123. For
				  // the moment, just ignore the
				  // column information and generate
				  // the matrix as if it were
				  // square. The call to
				  // GlobalAssemble later will set the
				  // correct values.
    matrix = std::auto_ptr<Epetra_FECrsMatrix>
	        (new Epetra_FECrsMatrix(Copy, row_map,
					&n_entries_per_row[0], false));


				  // TODO: As of now, assume that the
				  // sparsity pattern sits at the all
				  // processors (completely), and let
				  // each processor set its rows. Since
				  // this is wasteful, a better solution
				  // should be found in the future.
    Assert (matrix->NumGlobalRows() == (int)sparsity_pattern.n_rows(),
	    ExcDimensionMismatch (matrix->NumGlobalRows(),
				  sparsity_pattern.n_rows()));
    
				  // Trilinos has a bug for rectangular
				  // matrices at this point, so do not
				  // check for consistent column numbers
				  // here.
				  //
				  // this bug is filed in the Sandia
				  // bugzilla under #4123 and should be
				  // fixed for version 9.0
//        Assert (matrix->NumGlobalCols() == (int)sparsity_pattern.n_cols(),
//		ExcDimensionMismatch (matrix->NumGlobalCols(),
//				      sparsity_pattern.n_cols()));

    std::vector<double> values;
    std::vector<int>    row_indices;
    
    for (unsigned int row=0; row<n_rows; ++row)
      if (row_map.MyGID(row))
	{
	  const int row_length = sparsity_pattern.row_length(row);
	  row_indices.resize (row_length, -1);
	  values.resize (row_length, 0.);

	  for (int col=0; col < row_length; ++col)
	    row_indices[col] = sparsity_pattern.column_number (row, col);

	  matrix->InsertGlobalValues (row, row_length,
				      &values[0], &row_indices[0]);
	}
    
    last_action = Zero;

				  // In the end, the matrix needs to
				  // be compressed in order to be
				  // really ready.
    compress();    
  }



  void
  SparseMatrix::reinit (const SparseMatrix &sparse_matrix)
  {
    matrix.reset();

    row_map = sparse_matrix.row_map;
    col_map = sparse_matrix.col_map;
    matrix = std::auto_ptr<Epetra_FECrsMatrix>(new Epetra_FECrsMatrix(
				                   *sparse_matrix.matrix));
    compressed = true;
  }



  void
  SparseMatrix::reinit (const Epetra_Map                     &input_map,
			const ::dealii::SparseMatrix<double> &dealii_sparse_matrix,
			const double                          drop_tolerance)
  {
    reinit (input_map, input_map, dealii_sparse_matrix, drop_tolerance);
  }



  void
  SparseMatrix::reinit (const Epetra_Map                     &input_row_map,
			const Epetra_Map                     &input_col_map,
			const ::dealii::SparseMatrix<double> &dealii_sparse_matrix,
			const double                          drop_tolerance)
  {
    matrix.reset();

    unsigned int n_rows = dealii_sparse_matrix.m();

    Assert (input_row_map.NumGlobalElements() == (int)n_rows,
	    ExcDimensionMismatch (input_row_map.NumGlobalElements(),
				  n_rows));
    Assert (input_col_map.NumGlobalElements() == (int)dealii_sparse_matrix.n(),
	    ExcDimensionMismatch (input_col_map.NumGlobalElements(),
				  dealii_sparse_matrix.n()));

    row_map = input_row_map;
    col_map = input_col_map;

    std::vector<int> n_entries_per_row(n_rows);

    for (unsigned int row=0; row<n_rows; ++row)
      n_entries_per_row[(int)row] = 
	  dealii_sparse_matrix.get_sparsity_pattern().row_length(row);

				  // TODO: There seems to be problem
				  // in Epetra when a matrix is
				  // initialized with both row and
				  // column map. Maybe find something
				  // more out about this... It could
				  // be related to the bug #4123. For
				  // the moment, just ignore the
				  // column information and generate
				  // the matrix as if it were
				  // square. The call to
				  // GlobalAssemble later will set the
				  // correct values.
    matrix = std::auto_ptr<Epetra_FECrsMatrix>
	        (new Epetra_FECrsMatrix(Copy, row_map,
					&n_entries_per_row[0], false));

    std::vector<double> values;
    std::vector<int> row_indices;

    for (unsigned int row=0; row<n_rows; ++row)
      {
	values.resize (n_entries_per_row[row],0.);
	row_indices.resize (n_entries_per_row[row],-1);
	
	unsigned int index = 0;
	for (dealii::SparseMatrix<double>::const_iterator  
	      p  = dealii_sparse_matrix.begin(row);
	      p != dealii_sparse_matrix.end(row); ++p)
	  if (std::fabs(p->value()) > drop_tolerance)
	    {
	      row_indices[index] = p->column();
	      values[index]      = p->value();
	      ++index;
	    }
	const int n_row_entries = index;

	matrix->InsertGlobalValues(row, n_row_entries,
				   &values[0], &row_indices[0]);
      }

    compress();
  }



  void
  SparseMatrix::clear ()
  {
				  // When we clear the matrix, reset
				  // the pointer and generate an
				  // empty matrix.
    matrix.reset();
 
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
    row_map = Epetra_Map (0, 0, Epetra_MpiComm(MPI_COMM_WORLD));
#else
    row_map = Epetra_Map (0, 0, Epetra_SerialComm());
#endif

    col_map = row_map;

    matrix = std::auto_ptr<Epetra_FECrsMatrix> 
	      (new Epetra_FECrsMatrix(Copy, row_map, 0));
    compressed = true;
  }



  void
  SparseMatrix::compress ()
  {
				  // flush buffers
    int ierr;
    ierr = matrix->GlobalAssemble (col_map, row_map, true);
    
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    ierr = matrix->OptimizeStorage ();
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    last_action = Zero;

    compressed = true;
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
	ierr = matrix->GlobalAssemble(col_map, row_map, false);
	
	AssertThrow (ierr == 0, ExcTrilinosError(ierr));
      }

    if (last_action != Insert)
      last_action = Insert;

#ifdef DEBUG
    if (in_local_range (i) == false)
      compressed = false;
#endif
      
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
	ierr = matrix->GlobalAssemble(col_map, row_map, false);

	AssertThrow (ierr == 0, ExcTrilinosError(ierr));
    }

    if (last_action != Add)
      last_action = Add;

				  // we have to do above actions in any
				  // case to be consistent with the MPI
				  // communication model (see the
				  // comments in the documentation of
				  // TrilinosWrappers::Vector), but we
				  // can save some work if the addend is
				  // zero
    if (value == 0)
      return;

#ifdef DEBUG
    if (in_local_range (i) == false)
      compressed = false;
#endif
      
    int trilinos_i = i;
    int trilinos_j = j;

    const int ierr = matrix->SumIntoGlobalValues (trilinos_i, 1, 
						  const_cast<double*>(&value), 
						  &trilinos_j);

    AssertThrow (ierr <= 0, ExcAccessToNonPresentElement(i,j));
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



  void
  SparseMatrix::clear_row (const unsigned int   row,
			   const TrilinosScalar new_diag_value)
  {
    Assert (matrix->Filled()==true,
	    ExcMessage("Matrix must be compressed before invoking clear_row."));

				  // Only do this on the rows owned
				  // locally on this processor.
    int local_row = matrix->LRID(row);
    if (local_row >= 0)
      {
	TrilinosScalar *values;
	int *col_indices;
	int num_entries;
	const int ierr = matrix->ExtractMyRowView(local_row, num_entries,
						  values, col_indices);
	
	Assert (ierr == 0,
		ExcTrilinosError(ierr));
	
	int* diag_find = std::find(col_indices,col_indices+num_entries, 
				  local_row);
	int diag_index = (int)(diag_find - col_indices);

	for (int j=0; j<num_entries; ++j)
	  if (diag_index != j || new_diag_value == 0)
	    values[j] = 0.;

	if (diag_find && std::fabs(values[diag_index]) == 0 && new_diag_value!=0)
	  values[diag_index] = new_diag_value;
      }
  }



  void
  SparseMatrix::clear_rows (const std::vector<unsigned int> &rows,
			    const TrilinosScalar             new_diag_value)
  {
    for (unsigned int row=0; row<rows.size(); ++row)
      clear_row(rows[row], new_diag_value);
  }



  TrilinosScalar
  SparseMatrix::operator() (const unsigned int i,
			    const unsigned int j) const
  {
				      // Extract local indices in
				      // the matrix.
    int trilinos_i = matrix->LRID(i), trilinos_j = matrix->LRID(j);
    TrilinosScalar value = 0.;

				      // If the data is not on the
				      // present processor, we throw
				      // an exception. This is on of
				      // the two tiny differences to
				      // the el(i,j) call, which does
				      // not throw any assertions.
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
	if (matrix->Filled() == false)
	  matrix->FillComplete(col_map, row_map, true);

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
      
	int* el_find = std::find(&col_indices[0],&col_indices[0] + nnz_present,
				 trilinos_j);
      
	int el_index = (int)(el_find - col_indices);

				        // This is actually the only
				        // difference to the el(i,j)
				        // function, which means that
				        // we throw an exception in
				        // this case instead of just
				        // returning zero for an
				        // element that is not present
				        // in the sparsity pattern.
	if (!el_find)
	  {
	    Assert (false, ExcInvalidIndex (i,j));
	  }
	else
	  value = values[el_index];

      }

    return value;
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
				      // continue. Just print out
				      // zero.

				      // TODO: Is this reasonable? Or
				      // should we retain the assert
				      // call?
    if ((trilinos_i == -1 ) || (trilinos_j == -1))
      {
	return 0.;
	//Assert (false, ExcAccessToNonLocalElement(i, j, local_range().first,
	//				  local_range().second));
      }
    else
    {
				      // Check whether the matrix 
				      // already is transformed to
				      // local indices.
      if (!matrix->Filled())
	matrix->FillComplete(col_map, row_map, true);

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
      
      int* el_find = std::find(&col_indices[0],&col_indices[0] + nnz_present,
				trilinos_j);
      
      int el_index = (int)(el_find - col_indices);

      if (!el_find)
	value = 0;
      else
	{
	  value = values[el_index];
	}
    }

    return value;
  }



  TrilinosScalar
  SparseMatrix::diag_element (const unsigned int i) const
  {
    Assert (m() == n(), ExcNotQuadratic());

				  // Trilinos doesn't seem to have a
				  // more efficient way to access the
				  // diagonal than by just using the
				  // standard el(i,j) function.
    return el(i,i);
  }



  unsigned int
  SparseMatrix::m () const
  {
    int n_rows = matrix -> NumGlobalRows();

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
    begin = matrix -> RowMap().MinMyGID();
    end = matrix -> RowMap().MaxMyGID()+1;
    
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

				  // get a representation of the
				  // present row
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
    if (matrix->Filled() == false)
      matrix->FillComplete(col_map, row_map, true);

    TrilinosScalar result = matrix->NormOne();

    return result;
  }
  
  

  TrilinosScalar
  SparseMatrix::linfty_norm () const
  {
    if (matrix->Filled() == false)
      matrix->FillComplete(col_map, row_map, true);

    TrilinosScalar result = matrix->NormInf();

    return result;
  }



  TrilinosScalar
  SparseMatrix::frobenius_norm () const
  {
    if (matrix->Filled() == false)
      matrix->FillComplete(col_map, row_map, true);

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
    Assert (a !=0, ExcDivideByZero());

    const TrilinosScalar factor = 1./a;

    const int ierr = matrix->Scale (factor);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    return *this;
  }



  void
  SparseMatrix::vmult (VectorBase       &dst,
		       const VectorBase &src) const
  {
    Assert (&src != &dst, ExcSourceEqualsDestination());
    
    Assert (col_map.SameAs(src.vector->Map()) == true,
	    ExcMessage ("Column map of matrix does not fit with vector map!"));
    Assert (row_map.SameAs(dst.vector->Map()) == true,
	    ExcMessage ("Row map of matrix does not fit with vector map!"));

    if (matrix->Filled() == false)
      matrix->FillComplete(col_map, row_map, true);

    const int ierr = matrix->Multiply (false, *(src.vector), *(dst.vector));
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



  void
  SparseMatrix::Tvmult (VectorBase       &dst,
			const VectorBase &src) const
  {
    Assert (&src != &dst, ExcSourceEqualsDestination());

    Assert (row_map.SameAs(src.vector->Map()) == true,
	    ExcMessage ("Row map of matrix does not fit with vector map!"));
    Assert (col_map.SameAs(dst.vector->Map()) == true,
	    ExcMessage ("Column map of matrix does not fit with vector map!"));

    if (matrix->Filled() == false)
      matrix->FillComplete(col_map, row_map, true);

    const int ierr = matrix->Multiply (true, *(src.vector), *(dst.vector));
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



  void
  SparseMatrix::vmult_add (VectorBase       &dst,
			   const VectorBase &src) const
  {
    Assert (&src != &dst, ExcSourceEqualsDestination());

    VectorBase tmp(dst);
    vmult (dst, src);
    dst += tmp;
  }



  void
  SparseMatrix::Tvmult_add (VectorBase       &dst,
			    const VectorBase &src) const
  {
    Assert (&src != &dst, ExcSourceEqualsDestination());

    VectorBase tmp(dst);
    vmult (dst, src);
    dst += tmp;
  }



  TrilinosScalar
  SparseMatrix::matrix_norm_square (const VectorBase &v) const
  {
    VectorBase tmp (v);
    vmult (tmp, v);
    return tmp*v;
  }



  TrilinosScalar
  SparseMatrix::matrix_scalar_product (const VectorBase &u,
				       const VectorBase &v) const
  {
    VectorBase tmp (v);
    vmult (tmp, v);
    return u*tmp;
  }



  TrilinosScalar
  SparseMatrix::residual (VectorBase       &dst,
			  const VectorBase &x,
			  const VectorBase &b) const
  {
    vmult (dst, x);
    dst -= b;
    dst *= -1.;

    return dst.l2_norm();
  }



  void
  SparseMatrix::add (const TrilinosScalar  factor,
		     const SparseMatrix   &rhs)
  {
    Assert (rhs.m() == m(), ExcDimensionMismatch (rhs.m(), m()));
    Assert (rhs.n() == n(), ExcDimensionMismatch (rhs.n(), n()));
    
				  // I bet that there must be a
                                  // better way to do this but it
                                  // has not been found: currently,
				  // we simply go through each row
				  // of the argument matrix, copy
				  // it, scale it, and add it to
				  // the current matrix. that's
				  // probably not the most
				  // efficient way to do things.
    const std::pair<unsigned int, unsigned int>
      local_range = rhs.local_range();

    unsigned int max_row_length = 0;
    for (unsigned int row=local_range.first;
	 row < local_range.second; ++row)
      max_row_length
	= std::max (max_row_length,
		    static_cast<unsigned int>(rhs.matrix->NumGlobalEntries(row)));
    
    std::vector<int>    column_indices (max_row_length);
    std::vector<double> values (max_row_length);
    
    for (unsigned int row=local_range.first;
	 row < local_range.second; ++row)
      {
	int n_entries;
	rhs.matrix->ExtractGlobalRowCopy (row, max_row_length,
					  n_entries,
					  &values[0],
					  &column_indices[0]);
	for (int i=0; i<n_entries; ++i)
	  values[i] *= factor;

	matrix->SumIntoGlobalValues (row, n_entries,
				     &values[0],
				     &column_indices[0]);
      }

    compress ();
  }
  


				  // TODO: Currently this only flips a
				  // flag that tells Trilinos that any
				  // application should be done with
				  // the transpose. However, the
				  // matrix structure is not
				  // reset. Can we leave it like this?
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
  SparseMatrix::is_symmetric (const double tolerance) const
  {
    (void)tolerance;
    
    Assert (false, ExcNotImplemented());

    return false;
  }  



  bool
  SparseMatrix::is_hermitian () const
  {
    Assert (false, ExcNotImplemented());
    return false;
  }  



  void
  SparseMatrix::write_ascii ()
  {
    Assert (false, ExcNotImplemented());
  }



				  // As of now, no particularly neat
				  // ouput is generated in case of
				  // multiple processors.
  void SparseMatrix::print (std::ostream &out) const
  {
    double * values;
    int * indices;
    int num_entries;
  
    for (int i=0; i<matrix->NumMyRows(); ++i)
      {
	matrix->ExtractMyRowView (i, num_entries, values, indices);
	for (int j=0; j<num_entries; ++j)
	  out << "(" << i << "," << indices[matrix->GRID(j)] << ") " 
	      << values[j] << std::endl;
      }
  
    AssertThrow (out, ExcIO());
  }

}

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_USE_TRILINOS
