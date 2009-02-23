//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2008, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <lac/trilinos_sparse_matrix.h>

#include <lac/trilinos_sparsity_pattern.h>
#include <lac/sparsity_pattern.h>
#include <lac/compressed_sparsity_pattern.h>
#include <lac/compressed_set_sparsity_pattern.h>
#include <lac/compressed_simple_sparsity_pattern.h>

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
                  row_map (0, 0, Epetra_MpiComm(MPI_COMM_SELF)),
#else
                  row_map (0, 0, Epetra_SerialComm()),
#endif
		  col_map (row_map),
		  last_action (Zero),
		  compressed (true),
		  matrix (std::auto_ptr<Epetra_FECrsMatrix>
			  (new Epetra_FECrsMatrix(View, row_map, col_map, 0)))
  {
    matrix->FillComplete();
  }



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
				(new Epetra_FECrsMatrix(Copy, row_map, 
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
		    (new Epetra_FECrsMatrix(Copy, row_map,
		      (int*)const_cast<unsigned int*>(&(n_entries_per_row[0])),
					    false)))
  {}



  SparseMatrix::SparseMatrix (const unsigned int m,
			      const unsigned int n,
			      const unsigned int n_max_entries_per_row)
		  :
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
                  row_map (m, 0, Epetra_MpiComm(MPI_COMM_SELF)),
                  col_map (n, 0, Epetra_MpiComm(MPI_COMM_SELF)),
#else
                  row_map (m, 0, Epetra_SerialComm()),
                  col_map (n, 0, Epetra_SerialComm()),
#endif
		  last_action (Zero),
		  compressed (true),
				   // on one processor only, we know how the
				   // columns of the matrix will be
				   // distributed (everything on one
				   // processor), so we can hand in this
				   // information to the constructor. we
				   // can't do so in parallel, where the
				   // information from columns is only
				   // available when entries have been added
		  matrix (std::auto_ptr<Epetra_FECrsMatrix>
			  (new Epetra_FECrsMatrix(Copy, row_map, col_map,
						  int(n_max_entries_per_row), 
						  false)))
  {}



  SparseMatrix::SparseMatrix (const unsigned int               m,
			      const unsigned int               n,
			      const std::vector<unsigned int> &n_entries_per_row)
		  :
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
                  row_map (m, 0, Epetra_MpiComm(MPI_COMM_SELF)),
                  col_map (n, 0, Epetra_MpiComm(MPI_COMM_SELF)),
#else
                  row_map (m, 0, Epetra_SerialComm()),
                  col_map (n, 0, Epetra_SerialComm()),
#endif
		  last_action (Zero),
		  compressed (true),
		  matrix (std::auto_ptr<Epetra_FECrsMatrix>
			  (new Epetra_FECrsMatrix(Copy, row_map, col_map, 
			   (int*)const_cast<unsigned int*>(&(n_entries_per_row[0])), 
						  false)))
  {}



  SparseMatrix::SparseMatrix (const SparsityPattern &InputSP)
		  :
                  Subscriptor(),
                  row_map (InputSP.row_map),
		  col_map (InputSP.col_map),
		  last_action (Zero),
		  compressed (true),
		  matrix (std::auto_ptr<Epetra_FECrsMatrix>
			  (new Epetra_FECrsMatrix(Copy, *InputSP.graph, false)))
  {
    Assert(InputSP.graph->Filled() == true,
	   ExcMessage("The Trilinos sparsity pattern has not been compressed."));
    compress();
  }



  SparseMatrix::SparseMatrix (const SparseMatrix &InputMatrix)
		  :
                  Subscriptor(),
                  row_map (InputMatrix.row_map),
		  col_map (InputMatrix.col_map),
		  last_action (Zero),
		  compressed (true),
		  matrix (std::auto_ptr<Epetra_FECrsMatrix>
			  (new Epetra_FECrsMatrix(*InputMatrix.matrix)))
  {}



  SparseMatrix::~SparseMatrix ()
  {}



  SparseMatrix &
  SparseMatrix::copy_from (const SparseMatrix &m)
  {
    row_map = m.row_map;
    col_map = m.col_map;
    *matrix = *m.matrix;
    compress();
    return *this;
  }
  


  template <typename SparsityType>
  void
  SparseMatrix::reinit (const SparsityType &sparsity_pattern)
  {
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
    Epetra_MpiComm    trilinos_communicator (MPI_COMM_SELF);
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



  template <typename SparsityType>
  void
  SparseMatrix::reinit (const Epetra_Map    &input_map,
			const SparsityType  &sparsity_pattern)
  {
    reinit (input_map, input_map, sparsity_pattern);
  }



  template <typename SparsityType>
  void
  SparseMatrix::reinit (const Epetra_Map    &input_row_map,
			const Epetra_Map    &input_col_map,
			const SparsityType  &sparsity_pattern)
  {
    matrix.reset();

    const unsigned int n_rows = sparsity_pattern.n_rows();

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

				  // The deal.II notation of a Sparsity
				  // pattern corresponds to the Epetra
				  // concept of a Graph. Hence, we generate
				  // a graph by copying the sparsity pattern
				  // into it, and then build up the matrix
				  // from the graph. This is considerable
				  // faster than directly filling elements
				  // into the matrix. Moreover, it consumes
				  // less memory, since the internal
				  // reordering does not need to be done on
				  // all the double values.

				  // TODO: There seems to be problem in
				  // Epetra when a Graph/matrix is
				  // initialized with both row and column
				  // map. Maybe find something more out
				  // about this... It could be related to
				  // the bug #4123. For the moment, just
				  // ignore the column information and
				  // generate the graph to the matrix as if
				  // it were square. The call to
				  // GlobalAssemble later will set the
				  // correct values.

				   // for more than one processor, need to
				   // specify only row map first and let the
				   // matrix entries decide about the column
				   // map (which says which columns are
				   // present in the matrix, not to be
				   // confused with the col_map that tells
				   // how the domain dofs of the matrix will
				   // be distributed). for only one
				   // processor, we can directly assign the
				   // columns as well.
    std::auto_ptr<Epetra_CrsGraph> graph;
    if (row_map.Comm().NumProc() > 1)
      graph = std::auto_ptr<Epetra_CrsGraph> 
	(new Epetra_CrsGraph (Copy, row_map, 
			      &n_entries_per_row[input_row_map.MinMyGID()], true));
    else
      graph = std::auto_ptr<Epetra_CrsGraph> 
	(new Epetra_CrsGraph (Copy, row_map, col_map,
			      &n_entries_per_row[input_row_map.MinMyGID()], true));

				  // TODO: As of now, assume that the
				  // sparsity pattern sits at the all
				  // processors (completely), and let
				  // each processor set its rows. Since
				  // this is wasteful, a better solution
				  // should be found in the future.
    Assert (graph->NumGlobalRows() == (int)sparsity_pattern.n_rows(),
    	    ExcDimensionMismatch (graph->NumGlobalRows(),
    				  sparsity_pattern.n_rows()));

				  // Trilinos has a bug for rectangular
				  // matrices at this point, so do not
				  // check for consistent column numbers
				  // here.
				  //
				  // this bug is filed in the Sandia
				  // bugzilla under #4123 and should be
				  // fixed for version 9.0
//        Assert (graph.NumGlobalCols() == (int)sparsity_pattern.n_cols(),
//		ExcDimensionMismatch (graph.NumGlobalCols(),
//				      sparsity_pattern.n_cols()));

    std::vector<int>   row_indices;
    
    for (unsigned int row=0; row<n_rows; ++row)
      if (row_map.MyGID(row))
	{
	  const int row_length = sparsity_pattern.row_length(row);
	  row_indices.resize (row_length, -1);

	  for (int col=0; col < row_length; ++col)
	    row_indices[col] = sparsity_pattern.column_number (row, col);

	  graph->InsertGlobalIndices (row, row_length, &row_indices[0]);
	}

				  // Now, fill the graph (sort indices, make
				  // memory contiguous, etc).
    graph->FillComplete(col_map, row_map);

				  // And now finally generate the matrix.
    matrix =  std::auto_ptr<Epetra_FECrsMatrix>
                       (new Epetra_FECrsMatrix(Copy, *graph, false));

    last_action = Zero;

				  // In the end, the matrix needs to
				  // be compressed in order to be
				  // really ready.
    compress();
  }



				   // The CompressedSetSparsityPattern
				   // class stores the columns
				   // differently, so we need to
				   // explicitly rewrite this function
				   // in that case.
  template<>
  void
  SparseMatrix::reinit (const Epetra_Map                    &input_row_map,
			const Epetra_Map                    &input_col_map,
			const CompressedSetSparsityPattern  &sparsity_pattern)
  {
    matrix.reset();

    const unsigned int n_rows = sparsity_pattern.n_rows();

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


				  // The deal.II notation of a Sparsity
				  // pattern corresponds to the Epetra
				  // concept of a Graph. Hence, we generate
				  // a graph by copying the sparsity pattern
				  // into it, and then build up the matrix
				  // from the graph. This is considerable
				  // faster than directly filling elements
				  // into the matrix.

				  // TODO: There seems to be problem in
				  // Epetra when a Graph/matrix is
				  // initialized with both row and column
				  // map. Maybe find something more out
				  // about this... It could be related to
				  // the bug #4123. For the moment, just
				  // ignore the column information and
				  // generate the graph to the matrix as if
				  // it were square. The call to
				  // GlobalAssemble later will set the
				  // correct values.
    std::auto_ptr<Epetra_CrsGraph> graph;
    if (row_map.Comm().NumProc() > 1)
      graph = std::auto_ptr<Epetra_CrsGraph> 
	(new Epetra_CrsGraph (Copy, row_map, 
			      &n_entries_per_row[input_row_map.MinMyGID()], true));
    else
      graph = std::auto_ptr<Epetra_CrsGraph> 
	(new Epetra_CrsGraph (Copy, row_map, col_map,
			      &n_entries_per_row[input_row_map.MinMyGID()], true));

				  // TODO: As of now, assume that the
				  // sparsity pattern sits at the all
				  // processors (completely), and let
				  // each processor set its rows. Since
				  // this is wasteful, a better solution
				  // should be found in the future.
    Assert (graph->NumGlobalRows() == (int)sparsity_pattern.n_rows(),
    	    ExcDimensionMismatch (graph->NumGlobalRows(),
    				  sparsity_pattern.n_rows()));

				  // Trilinos has a bug for rectangular
				  // matrices at this point, so do not
				  // check for consistent column numbers
				  // here.
				  //
				  // this bug is filed in the Sandia
				  // bugzilla under #4123 and should be
				  // fixed for version 9.0
//        Assert (graph.NumGlobalCols() == (int)sparsity_pattern.n_cols(),
//		ExcDimensionMismatch (graph.NumGlobalCols(),
//				      sparsity_pattern.n_cols()));

    std::vector<int>   row_indices;
    
    for (unsigned int row=0; row<n_rows; ++row)
      if (row_map.MyGID(row))
	{
	  const int row_length = sparsity_pattern.row_length(row);
	  row_indices.resize (row_length, -1);

	  CompressedSetSparsityPattern::row_iterator col_num = 
	    sparsity_pattern.row_begin (row);

	  for (unsigned int col = 0; 
	       col_num != sparsity_pattern.row_end (row); 
	       ++col_num, ++col)
	    row_indices[col] = *col_num;

	  graph->InsertGlobalIndices (row, row_length, &row_indices[0]);
	}

				  // Now, fill the graph (sort indices, make
				  // memory contiguous, etc).
    graph->FillComplete(col_map, row_map);

				  // And now finally generate the matrix.
    matrix =  std::auto_ptr<Epetra_FECrsMatrix>
                       (new Epetra_FECrsMatrix(Copy, *graph, false));
    
    last_action = Zero;

				  // In the end, the matrix needs to
				  // be compressed in order to be
				  // really ready.
    compress();
  }



  void
  SparseMatrix::reinit (const SparsityPattern &sparsity_pattern)
  {
    matrix.reset();

    row_map = sparsity_pattern.row_map;
    col_map = sparsity_pattern.col_map;

    Assert (sparsity_pattern.graph->Filled() == true,
	    ExcMessage("The Trilinos sparsity pattern has not been compressed"));

    matrix = std::auto_ptr<Epetra_FECrsMatrix>
      (new Epetra_FECrsMatrix(Copy, *sparsity_pattern.graph, false));

    compress();
  }



  void
  SparseMatrix::reinit (const SparseMatrix &sparse_matrix)
  {
    matrix.reset();

    row_map = sparse_matrix.row_map;
    col_map = sparse_matrix.col_map;
    matrix = std::auto_ptr<Epetra_FECrsMatrix>
      (new Epetra_FECrsMatrix(Copy, sparse_matrix.matrix->Graph(), false));

    compress();
  }



  void
  SparseMatrix::reinit (const ::dealii::SparseMatrix<double> &dealii_sparse_matrix,
			const double                          drop_tolerance)
  {
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
    Epetra_MpiComm    trilinos_communicator (MPI_COMM_SELF);
#else
    Epetra_SerialComm trilinos_communicator;
#endif

    const Epetra_Map rows (dealii_sparse_matrix.m(),
			   0,
			   trilinos_communicator);
    const Epetra_Map columns (dealii_sparse_matrix.n(),
			      0,
			      trilinos_communicator);
    reinit (rows, columns, dealii_sparse_matrix, drop_tolerance);
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

    std::vector<TrilinosScalar> values;
    std::vector<unsigned int>   row_indices;

    for (unsigned int row=0; row<n_rows; ++row)
      {
	values.resize (n_entries_per_row[row],0.);
	row_indices.resize (n_entries_per_row[row],
			    numbers::invalid_unsigned_int);
	
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

	set (row, index, &row_indices[0], &values[0], false);
      }

    compress();
  }



  void 
  SparseMatrix::reinit (const Epetra_CrsMatrix &input_matrix)
  {
    matrix.reset();

    Assert (input_matrix.Filled()==true,
	    ExcMessage("Input CrsMatrix has not called FillComplete()!"));

    row_map = input_matrix.RangeMap();
    col_map = input_matrix.DomainMap();

    const Epetra_CrsGraph *graph = &input_matrix.Graph();

    matrix = std::auto_ptr<Epetra_FECrsMatrix> 
                      (new Epetra_FECrsMatrix(Copy, *graph, false));

    matrix->FillComplete (col_map, row_map, true);

    int length;
    int *row_indices;
    TrilinosScalar *values;
    for (unsigned int row=0; row<m(); ++row)
      if (row_map.MyGID(row))
	{
	  const int local_row = row_map.LID(row);
	  input_matrix.ExtractMyRowView(local_row, length, values, row_indices);
	  matrix->ReplaceMyValues(local_row, length, values, row_indices);
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
    row_map = Epetra_Map (0, 0, Epetra_MpiComm(MPI_COMM_SELF));
#else
    row_map = Epetra_Map (0, 0, Epetra_SerialComm());
#endif

    col_map = row_map;

    matrix = std::auto_ptr<Epetra_FECrsMatrix> 
	      (new Epetra_FECrsMatrix(View, row_map, 0));

    matrix->FillComplete();

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

	if (diag_find && std::fabs(values[diag_index]) == 0.0 && 
	    new_diag_value != 0.0)
	  values[diag_index] = new_diag_value;
      }
  }



  void
  SparseMatrix::clear_rows (const std::vector<unsigned int> &rows,
			    const TrilinosScalar             new_diag_value)
  {
    compress();
    for (unsigned int row=0; row<rows.size(); ++row)
      clear_row(rows[row], new_diag_value);

				        // This function needs to be called
				        // on all processors. We change some
				        // data, so we need to flush the
				        // buffers to make sure that the
				        // right data is used.
    compress();
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
				      // an exception. This is one of
				      // the two tiny differences to
				      // the el(i,j) call, which does
				      // not throw any assertions.
    if (trilinos_i == -1)
      {
	Assert (false, ExcAccessToNonLocalElement(i, j, local_range().first,
						  local_range().second));
      }
    else
      {
				      // Check whether the matrix has
				      // already been transformed to local
				      // indices.
	if (matrix->Filled() == false)
	  matrix->GlobalAssemble(col_map, row_map, true);

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
      
	int* el_find = std::find(col_indices, col_indices + nnz_present,
				 trilinos_j);

	int local_col_index = (int)(el_find - col_indices);

				        // This is actually the only
				        // difference to the el(i,j)
				        // function, which means that
				        // we throw an exception in
				        // this case instead of just
				        // returning zero for an
				        // element that is not present
				        // in the sparsity pattern.
	if (local_col_index == nnz_present)
	  {
	    Assert (false, ExcInvalidIndex (i,j));
	  }
	else
	  value = values[local_col_index];
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
	matrix->GlobalAssemble(col_map, row_map, true);

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
      int* el_find = std::find(col_indices, col_indices + nnz_present,
			       trilinos_j);

      int local_col_index = (int)(el_find - col_indices);


				        // This is actually the only
				        // difference to the () function
				        // querying (i,j), where we throw an
				        // exception instead of just
				        // returning zero for an element
				        // that is not present in the
				        // sparsity pattern.
      if (local_col_index == nnz_present)
	value = 0;
      else
	value = values[local_col_index];
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
    unsigned int n_cols = matrix -> NumGlobalCols();
    return n_cols;
  }



  unsigned int
  SparseMatrix::local_size () const
  {
    unsigned int n_rows = matrix -> NumMyRows();

    return n_rows;
  }



  std::pair<unsigned int, unsigned int>
  SparseMatrix::local_range () const
  {
    unsigned int begin, end;
    begin = matrix -> RowMap().MinMyGID();
    end = matrix -> RowMap().MaxMyGID()+1;
    
    return std::make_pair (begin, end);
  }



  unsigned int
  SparseMatrix::n_nonzero_elements () const
  {
    unsigned int nnz = matrix->NumGlobalNonzeros();

    return nnz;
  }



  unsigned int
  SparseMatrix::row_length (const unsigned int row) const
  {
    Assert (row < m(), ExcInternalError());

				  // get a representation of the
				  // present row
    int ncols = -1;
    int local_row = matrix->LRID(row);

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

    if (matrix->Filled() == false)
      matrix->FillComplete(col_map, row_map, true);

    Assert (src.vector->Map().SameAs(matrix->DomainMap()) == true,
	    ExcMessage ("Column map of matrix does not fit with vector map!"));
    Assert (dst.vector->Map().SameAs(matrix->RangeMap()) == true,
	    ExcMessage ("Row map of matrix does not fit with vector map!"));

    const int ierr = matrix->Multiply (false, *(src.vector), *(dst.vector));
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



  void
  SparseMatrix::Tvmult (VectorBase       &dst,
			const VectorBase &src) const
  {
    Assert (&src != &dst, ExcSourceEqualsDestination());

    if (matrix->Filled() == false)
      matrix->FillComplete(col_map, row_map, true);

    Assert (src.vector->Map().SameAs(matrix->RangeMap()) == true,
	    ExcMessage ("Column map of matrix does not fit with vector map!"));
    Assert (dst.vector->Map().SameAs(matrix->DomainMap()) == true,
	    ExcMessage ("Row map of matrix does not fit with vector map!"));

    const int ierr = matrix->Multiply (true, *(src.vector), *(dst.vector));
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



  void
  SparseMatrix::vmult_add (VectorBase       &dst,
			   const VectorBase &src) const
  {
    Assert (&src != &dst, ExcSourceEqualsDestination());

				   // Choose to reinit the vector with fast
				   // argument set, which does not overwrite
				   // the content -- this is what we need
				   // since we're going to overwrite that
				   // anyway in the vmult operation.
    temp_vector.reinit(dst, true);

    vmult (temp_vector, src);
    dst += temp_vector;
  }



  void
  SparseMatrix::Tvmult_add (VectorBase       &dst,
			    const VectorBase &src) const
  {
    Assert (&src != &dst, ExcSourceEqualsDestination());

    temp_vector.reinit(dst, true);

    vmult (temp_vector, src);
    dst += temp_vector;
  }



  TrilinosScalar
  SparseMatrix::matrix_norm_square (const VectorBase &v) const
  {
    Assert (row_map.SameAs(col_map),
	    ExcDimensionMismatch(row_map.NumGlobalElements(),
				 col_map.NumGlobalElements()));

    temp_vector.reinit(v);

    vmult (temp_vector, v);
    return temp_vector*v;
  }



  TrilinosScalar
  SparseMatrix::matrix_scalar_product (const VectorBase &u,
				       const VectorBase &v) const
  {
    Assert (row_map.SameAs(col_map),
	    ExcDimensionMismatch(row_map.NumGlobalElements(),
				 col_map.NumGlobalElements()));

    temp_vector.reinit(v);

    vmult (temp_vector, v);
    return u*temp_vector;
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
    
    const std::pair<unsigned int, unsigned int>
      local_range = rhs.local_range();

    unsigned int max_row_length = 0;
    for (unsigned int row=local_range.first;
	 row < local_range.second; ++row)
      max_row_length
	= std::max (max_row_length,
		    static_cast<unsigned int>(rhs.matrix->NumGlobalEntries(row)));
    
    std::vector<int>            column_indices (max_row_length);
    std::vector<TrilinosScalar> values (max_row_length);
    int ierr;

				   // If both matrices have been transformed
				   // to local index space (in Trilinos
				   // speak: they are filled), we're having
				   // matrices based on the same indices
				   // with the same number of nonzeros
				   // (actually, we'd need sparsity pattern,
				   // but that is too expensive to check),
				   // we can extract views of the column
				   // data on both matrices and simply
				   // manipulate the values that are
				   // addressed by the pointers.
    if (matrix->Filled() == true && 
	rhs.matrix->Filled() == true &&
	this->local_range() == local_range && 
	matrix->NumMyNonzeros() == rhs.matrix->NumMyNonzeros())
      for (unsigned int row=local_range.first;
	   row < local_range.second; ++row)
	{
	  Assert (matrix->NumGlobalEntries(row) == 
		  rhs.matrix->NumGlobalEntries(row),
		  ExcDimensionMismatch(matrix->NumGlobalEntries(row),
				       rhs.matrix->NumGlobalEntries(row)));

	  const int row_local = matrix->RowMap().LID(row);
	  int n_entries, rhs_n_entries;
	  TrilinosScalar *value_ptr, *rhs_value_ptr;

				   // In debug mode, we want to check
				   // whether the indices really are the
				   // same in the calling matrix and the
				   // input matrix. The reason for doing
				   // this only in debug mode is that both
				   // extracting indices and comparing
				   // indices is relatively slow compared to
				   // just working with the values.
#ifdef DEBUG
	  int *index_ptr, *rhs_index_ptr;
	  ierr = rhs.matrix->ExtractMyRowView (row_local, rhs_n_entries, 
					       rhs_value_ptr, rhs_index_ptr);
	  Assert (ierr == 0, ExcTrilinosError(ierr));

	  ierr = matrix->ExtractMyRowView (row_local, n_entries, value_ptr,
					   index_ptr);
	  Assert (ierr == 0, ExcTrilinosError(ierr));
#else
	  rhs.matrix->ExtractMyRowView (row_local, rhs_n_entries,rhs_value_ptr);
	  matrix->ExtractMyRowView (row_local, n_entries, value_ptr);
#endif

	  AssertThrow (n_entries == rhs_n_entries,
		       ExcDimensionMismatch (n_entries, rhs_n_entries));

	  for (int i=0; i<n_entries; ++i)
	    {
	      *value_ptr++ += *rhs_value_ptr++ * factor;
#ifdef DEBUG
	      Assert (*index_ptr++ == *rhs_index_ptr++,
		      ExcInternalError());
#endif
	    }
	}
				   // If we have different sparsity patterns
				   // (expressed by a different number of
				   // nonzero elements), we have to be more
				   // careful and extract a copy of the row
				   // data, multiply it by the factor and
				   // then add it to the matrix using the
				   // respective add() function.
    else if (matrix->Filled() == true && rhs.matrix->Filled() == true &&
	     this->local_range() == local_range)
      for (unsigned int row=local_range.first;
	   row < local_range.second; ++row)
	{
	  const int row_local = matrix->RowMap().LID(row);
	  int n_entries;

	  ierr = rhs.matrix->ExtractMyRowCopy (row_local, max_row_length,
					       n_entries,
					       &values[0],
					       &column_indices[0]);
	  Assert (ierr == 0, ExcTrilinosError(ierr));

	  for (int i=0; i<n_entries; ++i)
	    values[i] *= factor;

	  TrilinosScalar *value_ptr = &values[0];

	  ierr = matrix->SumIntoMyValues (row_local, n_entries, value_ptr,
					  &column_indices[0]);
	  Assert (ierr == 0, ExcTrilinosError(ierr));
	}
    else
      {
	for (unsigned int row=local_range.first;
	     row < local_range.second; ++row)
	  {
	    int n_entries;
	    ierr = rhs.matrix->Epetra_CrsMatrix::ExtractGlobalRowCopy 
	      ((int)row, max_row_length, n_entries, &values[0], &column_indices[0]);
	    Assert (ierr == 0, ExcTrilinosError(ierr));

	    for (int i=0; i<n_entries; ++i)
	      values[i] *= factor;

	    ierr = matrix->Epetra_CrsMatrix::SumIntoGlobalValues 
	      ((int)row, n_entries, &values[0], &column_indices[0]);
	    Assert (ierr == 0, ExcTrilinosError(ierr));
	  }
	compress ();
      }
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
  void
  SparseMatrix::print (std::ostream &out) const
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




  // explicit instantiations
  //
  template void
  SparseMatrix::reinit (const dealii::SparsityPattern &);
  template void
  SparseMatrix::reinit (const CompressedSparsityPattern &);
  template void
  SparseMatrix::reinit (const CompressedSetSparsityPattern &);
  template void
  SparseMatrix::reinit (const CompressedSimpleSparsityPattern &);


  template void
  SparseMatrix::reinit (const Epetra_Map &,
			const dealii::SparsityPattern &);
  template void
  SparseMatrix::reinit (const Epetra_Map &,
			const CompressedSparsityPattern &);
  template void
  SparseMatrix::reinit (const Epetra_Map &,
			const CompressedSetSparsityPattern &);
  template void
  SparseMatrix::reinit (const Epetra_Map &,
			const CompressedSimpleSparsityPattern &);


  template void
  SparseMatrix::reinit (const Epetra_Map &,
			const Epetra_Map &,
			const dealii::SparsityPattern &);
  template void
  SparseMatrix::reinit (const Epetra_Map &,
			const Epetra_Map &,
			const CompressedSparsityPattern &);
  template void
  SparseMatrix::reinit (const Epetra_Map &,
			const Epetra_Map &,
			const CompressedSimpleSparsityPattern &);

}

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_USE_TRILINOS
