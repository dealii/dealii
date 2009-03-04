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

#include <lac/trilinos_sparsity_pattern.h>

#include <lac/sparsity_pattern.h>
#include <lac/compressed_sparsity_pattern.h>
#include <lac/compressed_set_sparsity_pattern.h>
#include <lac/compressed_simple_sparsity_pattern.h>

#ifdef DEAL_II_USE_TRILINOS

DEAL_II_NAMESPACE_OPEN

namespace TrilinosWrappers
{
  namespace SparsityPatternIterators
  {
    void
    const_iterator::Accessor::
    visit_present_row ()
    {
				  // if we are asked to visit the
				  // past-the-end line, then simply
				  // release all our caches and go on
				  // with life
      if (this->a_row == sparsity_pattern->n_rows())
	{
	  colnum_cache.reset ();

	  return;
	}
      
				  // otherwise first flush Trilinos caches
      sparsity_pattern->compress ();

				  // get a representation of the present
				  // row
      int ncols;
      int colnums = sparsity_pattern->n_cols();
      
      int ierr;
      ierr = sparsity_pattern->graph->ExtractGlobalRowCopy((int)this->a_row, 
							   colnums,
							   ncols, 
							   (int*)&(*colnum_cache)[0]);
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
  SparsityPattern::SparsityPattern ()
		  :
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
                  row_map (0, 0, Epetra_MpiComm(MPI_COMM_SELF)),
#else
                  row_map (0, 0, Epetra_SerialComm()),
#endif
		  col_map (row_map),
		  compressed (true),
		  graph (std::auto_ptr<Epetra_FECrsGraph>
			 (new Epetra_FECrsGraph(View, row_map, col_map, 0)))
  {
    graph->FillComplete();
  }

  SparsityPattern::SparsityPattern (const Epetra_Map  &InputMap,
				    const unsigned int n_entries_per_row)
		  :
                  row_map (InputMap),
		  col_map (row_map),
		  compressed (false),
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
		  graph (row_map.Comm().NumProc() > 1 ?
			 (std::auto_ptr<Epetra_FECrsGraph>
				   (new Epetra_FECrsGraph(Copy, row_map, 
							  (int)n_entries_per_row, 
							  false)))
			 :
			 (std::auto_ptr<Epetra_FECrsGraph>
			           (new Epetra_FECrsGraph(Copy, row_map, col_map, 
							  (int)n_entries_per_row, 
							  false)))
			 )
  {}
 
  SparsityPattern::SparsityPattern (const Epetra_Map                &InputMap,
				    const std::vector<unsigned int> &n_entries_per_row)
		  :
                  row_map (InputMap),
		  col_map (row_map),
		  compressed (false),
		  graph (row_map.Comm().NumProc() > 1 ?
			 (std::auto_ptr<Epetra_FECrsGraph>
			  (new Epetra_FECrsGraph(Copy, row_map, 
						 (int*)const_cast<unsigned int*>
						 (&(n_entries_per_row[row_map.MinMyGID()])),
						 false)))
			 :
			 (std::auto_ptr<Epetra_FECrsGraph>
			  (new Epetra_FECrsGraph(Copy, row_map, col_map, 
						 (int*)const_cast<unsigned int*>
						 (&(n_entries_per_row[row_map.MinMyGID()])),
						 false)))
			 )
  {}

  SparsityPattern::SparsityPattern (const Epetra_Map  &InputRowMap,
				    const Epetra_Map  &InputColMap,
				    const unsigned int n_entries_per_row)
		  :
                  row_map (InputRowMap),
                  col_map (InputColMap),
		  compressed (false),
		  graph (row_map.Comm().NumProc() > 1 ?
			 (std::auto_ptr<Epetra_FECrsGraph>
				   (new Epetra_FECrsGraph(Copy, row_map, 
							  (int)n_entries_per_row, 
							  false)))
			 :
			 (std::auto_ptr<Epetra_FECrsGraph>
			           (new Epetra_FECrsGraph(Copy, row_map, col_map, 
							  (int)n_entries_per_row, 
							  false)))
			 )
  {}

  SparsityPattern::SparsityPattern (const Epetra_Map                &InputRowMap,
				    const Epetra_Map                &InputColMap,
				    const std::vector<unsigned int> &n_entries_per_row)
		  :
                  row_map (InputRowMap),
                  col_map (InputColMap),
		  compressed (false),
		  graph (row_map.Comm().NumProc() > 1 ?
			 (std::auto_ptr<Epetra_FECrsGraph>
			  (new Epetra_FECrsGraph(Copy, row_map, 
						 (int*)const_cast<unsigned int*>
						 (&(n_entries_per_row[row_map.MinMyGID()])),
						 false)))
			 :
			 (std::auto_ptr<Epetra_FECrsGraph>
			  (new Epetra_FECrsGraph(Copy, row_map, col_map, 
						 (int*)const_cast<unsigned int*>
						 (&(n_entries_per_row[row_map.MinMyGID()])),
						 false)))
			 )
  {}

  SparsityPattern::SparsityPattern (const unsigned int m,
				    const unsigned int n,
				    const unsigned int n_entries_per_row)
		  :
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
                  row_map (m, 0, Epetra_MpiComm(MPI_COMM_SELF)),
                  col_map (n, 0, Epetra_MpiComm(MPI_COMM_SELF)),
#else
                  row_map (m, 0, Epetra_SerialComm()),
                  col_map (n, 0, Epetra_SerialComm()),
#endif
		  compressed (false),
		  graph (std::auto_ptr<Epetra_FECrsGraph>
			 (new Epetra_FECrsGraph(Copy, row_map, col_map,
						int(n_entries_per_row), false)))
  {}

  SparsityPattern::SparsityPattern (const unsigned int               m,
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
		  compressed (false),
		  graph (std::auto_ptr<Epetra_FECrsGraph>
			 (new Epetra_FECrsGraph(Copy, row_map, col_map, 
			  (int*)const_cast<unsigned int*>(&(n_entries_per_row[0])), 
						false)))
  {}

				   // Copy function is currently not working
				   // because the Trilinos Epetra_FECrsGraph
				   // does not implement a constructor from
				   // another graph.
  /*
  SparsityPattern::SparsityPattern (const SparsityPattern &InputSP)
  		  :
                  Subscriptor(),
		  row_map (InputSP.row_map),
 		  col_map (InputSP.col_map),
  		  compressed (false),
  		  graph (std::auto_ptr<Epetra_FECrsGraph>
  			  (new Epetra_FECrsGraph(*InputSP.graph)))
  {}
  */



  SparsityPattern::~SparsityPattern ()
  {}



  void 
  SparsityPattern::reinit (const Epetra_Map   &input_map,
			   const unsigned int  n_entries_per_row)
  {
    reinit (input_map, input_map, n_entries_per_row);
  }


  void 
  SparsityPattern::reinit (const unsigned int  m,
			   const unsigned int  n,
			   const unsigned int  n_entries_per_row)
  {
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
    Epetra_MpiComm    trilinos_communicator (MPI_COMM_SELF);
#else
    Epetra_SerialComm trilinos_communicator;
#endif

    const Epetra_Map rows (m, 0, trilinos_communicator);
    const Epetra_Map columns (n, 0, trilinos_communicator);

    reinit (rows, columns, n_entries_per_row);
  }


  void 
  SparsityPattern::reinit (const Epetra_Map   &input_row_map,
			   const Epetra_Map   &input_col_map,
			   const unsigned int  n_entries_per_row)
  {
    graph.reset();

    row_map = input_row_map;
    col_map = input_col_map;

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
    if (row_map.Comm().NumProc() > 1)
      graph = std::auto_ptr<Epetra_FECrsGraph>
	(new Epetra_FECrsGraph(Copy, row_map, n_entries_per_row, false));
    else
      graph = std::auto_ptr<Epetra_FECrsGraph>
	(new Epetra_FECrsGraph(Copy, row_map, col_map, n_entries_per_row, false));
  }



  void 
  SparsityPattern::reinit (const Epetra_Map   &input_map,
			   const std::vector<unsigned int> &n_entries_per_row)
  {
    reinit (input_map, input_map, n_entries_per_row);
  }



  void 
  SparsityPattern::reinit (const unsigned int  m,
			   const unsigned int  n,
			   const std::vector<unsigned int> &n_entries_per_row)
  {
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
    Epetra_MpiComm    trilinos_communicator (MPI_COMM_SELF);
#else
    Epetra_SerialComm trilinos_communicator;
#endif

    const Epetra_Map rows (m, 0, trilinos_communicator);
    const Epetra_Map columns (n, 0, trilinos_communicator);

    reinit (rows, columns, n_entries_per_row);
  }



  void 
  SparsityPattern::reinit (const Epetra_Map   &input_row_map,
			   const Epetra_Map   &input_col_map,
			   const std::vector<unsigned int> &n_entries_per_row)
  {
    graph.reset();

    Assert (n_entries_per_row.size() == 
	      static_cast<unsigned int>(input_row_map.NumGlobalElements()),
	    ExcDimensionMismatch (n_entries_per_row.size(),
				  input_row_map.NumGlobalElements()));
    row_map = input_row_map;
    col_map = input_col_map;

    if (row_map.Comm().NumProc() > 1)
      graph = std::auto_ptr<Epetra_FECrsGraph>
	(new Epetra_FECrsGraph(Copy, row_map, 
			       n_entries_per_row[input_row_map.MinMyGID()], 
			       false));
    else
      graph = std::auto_ptr<Epetra_FECrsGraph>
	(new Epetra_FECrsGraph(Copy, row_map, col_map,
			       n_entries_per_row[input_row_map.MinMyGID()], 
			       false));      
  }



  template <typename SparsityType>
  void 
  SparsityPattern::reinit (const Epetra_Map   &input_map,
			   const SparsityType &sp)
  {
    reinit (input_map, input_map, sp);
  }



  template <typename SparsityType>
  void 
  SparsityPattern::reinit (const Epetra_Map   &input_row_map,
			   const Epetra_Map   &input_col_map,
			   const SparsityType &sp)
  {
    graph.reset();

    Assert (sp.n_rows() == 
	      static_cast<unsigned int>(input_row_map.NumGlobalElements()),
	    ExcDimensionMismatch (sp.n_rows(),
				  input_row_map.NumGlobalElements()));
    Assert (sp.n_cols() == 
	      static_cast<unsigned int>(input_col_map.NumGlobalElements()),
	    ExcDimensionMismatch (sp.n_cols(),
				  input_col_map.NumGlobalElements()));

    row_map = input_row_map;
    col_map = input_col_map;

    const unsigned int n_rows = sp.n_rows();

    std::vector<int> n_entries_per_row(n_rows);

    for (unsigned int row=0; row<n_rows; ++row)
      n_entries_per_row[row] = sp.row_length(row);

    if (row_map.Comm().NumProc() > 1)
      graph = std::auto_ptr<Epetra_FECrsGraph>
	(new Epetra_FECrsGraph(Copy, row_map, 
			       n_entries_per_row[input_row_map.MinMyGID()], 
			       false));
    else
      graph = std::auto_ptr<Epetra_FECrsGraph>
	(new Epetra_FECrsGraph(Copy, row_map, col_map,
			       n_entries_per_row[input_row_map.MinMyGID()], 
			       false));      

    Assert (graph->NumGlobalRows() == (int)sp.n_rows(),
    	    ExcDimensionMismatch (graph->NumGlobalRows(),
    				  sp.n_rows()));


    std::vector<int>   row_indices;
    
    for (unsigned int row=0; row<n_rows; ++row)
      if (row_map.MyGID(row))
	{
	  const int row_length = sp.row_length(row);
	  row_indices.resize (row_length, -1);

	  for (int col=0; col < row_length; ++col)
	    row_indices[col] = sp.column_number (row, col);

	  graph->Epetra_CrsGraph::InsertGlobalIndices (row, row_length, 
						       &row_indices[0]);
	}

    compress();
  }



  template<>
  void 
  SparsityPattern::reinit (const Epetra_Map   &input_row_map,
			   const Epetra_Map   &input_col_map,
			   const CompressedSetSparsityPattern &sp)
  {
    graph.reset();

    Assert (sp.n_rows() == 
	      static_cast<unsigned int>(input_row_map.NumGlobalElements()),
	    ExcDimensionMismatch (sp.n_rows(),
				  input_row_map.NumGlobalElements()));
    Assert (sp.n_cols() ==
	      static_cast<unsigned int>(input_col_map.NumGlobalElements()),
	    ExcDimensionMismatch (sp.n_cols(),
				  input_col_map.NumGlobalElements()));

    row_map = input_row_map;
    col_map = input_col_map;

    const unsigned int n_rows = sp.n_rows();

    std::vector<int> n_entries_per_row(n_rows);

    for (unsigned int row=0; row<n_rows; ++row)
      n_entries_per_row[row] = sp.row_length(row);

    if (row_map.Comm().NumProc() > 1)
      graph = std::auto_ptr<Epetra_FECrsGraph>
	(new Epetra_FECrsGraph(Copy, row_map, 
			       n_entries_per_row[input_row_map.MinMyGID()], 
			       false));
    else
      graph = std::auto_ptr<Epetra_FECrsGraph>
	(new Epetra_FECrsGraph(Copy, row_map, col_map,
			       n_entries_per_row[input_row_map.MinMyGID()], 
			       false));      

    Assert (graph->NumGlobalRows() == (int)sp.n_rows(),
    	    ExcDimensionMismatch (graph->NumGlobalRows(),
    				  sp.n_rows()));


    std::vector<int>   row_indices;
    
    for (unsigned int row=0; row<n_rows; ++row)
      if (row_map.MyGID(row))
	{
	  const int row_length = sp.row_length(row);
	  row_indices.resize (row_length, -1);

	  CompressedSetSparsityPattern::row_iterator col_num = 
	    sp.row_begin (row);

	  for (unsigned int col = 0; 
	       col_num != sp.row_end (row); 
	       ++col_num, ++col)
	    row_indices[col] = *col_num;

	  graph->Epetra_CrsGraph::InsertGlobalIndices (row, row_length, 
						       &row_indices[0]);
	}

    compress();
  }



  /*  void
  SparsityPattern::copy_from (const SparsityPattern &sp)
  {
    graph.reset();
    row_map = sp.row_map;
    col_map = sp.col_map;

    graph = std::auto_ptr<Epetra_FECrsGraph> (new Epetra_FECrsGraph(*sp.graph));
  }
  */


  template <typename SparsityType>
  void 
  SparsityPattern::copy_from (const SparsityType &sp)
  {
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
    Epetra_MpiComm    trilinos_communicator (MPI_COMM_SELF);
#else
    Epetra_SerialComm trilinos_communicator;
#endif

    const Epetra_Map rows (sp.n_rows(), 0, trilinos_communicator);
    const Epetra_Map columns (sp.n_cols(), 0, trilinos_communicator);

    reinit (rows, columns, sp);
  }



  void
  SparsityPattern::clear ()
  {
				  // When we clear the matrix, reset
				  // the pointer and generate an
				  // empty matrix.
    graph.reset();
 
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
    row_map = Epetra_Map (0, 0, Epetra_MpiComm(MPI_COMM_SELF));
#else
    row_map = Epetra_Map (0, 0, Epetra_SerialComm());
#endif

    col_map = row_map;

    graph = std::auto_ptr<Epetra_FECrsGraph> 
	      (new Epetra_FECrsGraph(View, row_map, 0));

    graph->FillComplete();

    compressed = true;
  }



  void
  SparsityPattern::compress ()
  {
    int ierr;
    ierr = graph->GlobalAssemble (col_map, row_map, true);
    
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    ierr = graph->OptimizeStorage ();
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    compressed = true;
  }



  bool
  SparsityPattern::exists (const unsigned int i,
			   const unsigned int j) const
  {
				      // Extract local indices in
				      // the matrix.
    int trilinos_i = graph->LRID(i), trilinos_j = graph->LRID(j);

				      // If the data is not on the
				      // present processor, we throw
				      // an exception. This is on of
				      // the two tiny differences to
				      // the el(i,j) call, which does
				      // not throw any assertions.
    if (trilinos_i == -1)
      {
	return false;
      }
    else
      {
				      // Check whether the matrix 
				      // already is transformed to
				      // local indices.
	if (graph->Filled() == false)
	  {
	    int nnz_present = graph->NumGlobalIndices(i);
	    int nnz_extracted;
	    int *col_indices;

				      // Generate the view and make
				      // sure that we have not generated
				      // an error.
	    int ierr = graph->ExtractGlobalRowView(trilinos_i, nnz_extracted,
						   col_indices);
	    Assert (ierr==0, ExcTrilinosError(ierr));
	    Assert (nnz_present == nnz_extracted,
		    ExcDimensionMismatch(nnz_present, nnz_extracted));

				      // Search the index
	    int* el_find = std::find(col_indices, col_indices + nnz_present,
				     trilinos_j);

	    int local_col_index = (int)(el_find - col_indices);

	    if (local_col_index == nnz_present)
	      return false;
	  }
	else
	  {
				      // Prepare pointers for extraction
				      // of a view of the row.
	    int nnz_present = graph->NumGlobalIndices(i);
	    int nnz_extracted;
	    int *col_indices;

				      // Generate the view and make
				      // sure that we have not generated
				      // an error.
	    int ierr = graph->ExtractMyRowView(trilinos_i, nnz_extracted,
					       col_indices);
	    Assert (ierr==0, ExcTrilinosError(ierr));

	    Assert (nnz_present == nnz_extracted,
		    ExcDimensionMismatch(nnz_present, nnz_extracted));

				      // Search the index
	    int* el_find = std::find(col_indices, col_indices + nnz_present,
				     trilinos_j);

	    int local_col_index = (int)(el_find - col_indices);

	    if (local_col_index == nnz_present)
	      return false;
	  }
      }

    return true;
  }



  unsigned int
  SparsityPattern::bandwidth () const
  {
    unsigned int local_b=0;
    int global_b=0;
    for (unsigned int i=0; i<local_size(); ++i)
      {
	int * indices;
	int num_entries;
	graph->ExtractMyRowView(i, num_entries, indices);
	for (unsigned int j=0; j<(unsigned int)num_entries; ++j)
	  {
	    if (static_cast<unsigned int>(std::abs(static_cast<int>(i-indices[j]))) > local_b)
	      local_b = std::abs(static_cast<signed int>(i-indices[j]));
	  }
      }
    graph->Comm().MaxAll((int *)&local_b, &global_b, 1);
    return static_cast<unsigned int>(global_b);
  }



  unsigned int
  SparsityPattern::n_rows () const
  {
    const int n_rows = graph -> NumGlobalRows();
    return n_rows;
  }



  unsigned int
  SparsityPattern::n_cols () const
  {
    int n_cols;
    if (graph->Filled() == true)
      n_cols = graph -> NumGlobalCols();
    else
      n_cols = col_map.NumGlobalElements();

    return n_cols;
  }



  unsigned int
  SparsityPattern::local_size () const
  {
    int n_rows = graph -> NumMyRows();

    return n_rows;
  }



  std::pair<unsigned int, unsigned int>
  SparsityPattern::local_range () const
  {
    unsigned int begin, end;
    begin = graph -> RowMap().MinMyGID();
    end = graph -> RowMap().MaxMyGID()+1;
    
    return std::make_pair (begin, end);
  }



  unsigned int
  SparsityPattern::n_nonzero_elements () const
  {
    int nnz = graph->NumGlobalEntries();

    return static_cast<unsigned int>(nnz);
  }



  unsigned int
  SparsityPattern::max_entries_per_row () const
  {
    int nnz = graph->MaxNumIndices();

    return static_cast<unsigned int>(nnz);
  }



  unsigned int
  SparsityPattern::row_length (const unsigned int row) const
  {
    Assert (row < n_rows(), ExcInternalError());

				  // get a representation of the
				  // present row
    int ncols = -1;
    int local_row = graph->LRID(row);

				  // on the processor who owns this
				  // row, we'll have a non-negative
				  // value.
    if (local_row >= 0)
      ncols = graph->NumMyIndices (local_row);

    return static_cast<unsigned int>(ncols);
  }



  void
  SparsityPattern::write_ascii ()
  {
    Assert (false, ExcNotImplemented());
  }



				  // As of now, no particularly neat
				  // ouput is generated in case of
				  // multiple processors.
  void
  SparsityPattern::print (std::ostream &out) const
  {
    int * indices;
    int num_entries;
  
    for (int i=0; i<graph->NumMyRows(); ++i)
      {
	graph->ExtractMyRowView (i, num_entries, indices);
	for (int j=0; j<num_entries; ++j)
	  out << "(" << i << "," << indices[graph->GRID(j)] << ") " 
	      << std::endl;
      }
  
    AssertThrow (out, ExcIO());
  }



  void
  SparsityPattern::print_gnuplot (std::ostream &out) const
  { 
    Assert (graph->Filled() == true, ExcInternalError());
    for (unsigned int row=0; row<local_size(); ++row)
      {
	signed int * indices;
	int num_entries;
	graph->ExtractMyRowView (row, num_entries, indices);

	for (unsigned int j=0; j<(unsigned int)num_entries; ++j)
                                         // while matrix entries are usually
                                         // written (i,j), with i vertical and
                                         // j horizontal, gnuplot output is
                                         // x-y, that is we have to exchange
                                         // the order of output
	  out << indices[graph->GRID(j)] << " " << -static_cast<signed int>(row) 
	      << std::endl;
      }      

    AssertThrow (out, ExcIO());
  }




  // explicit instantiations
  //
  template void
  SparsityPattern::copy_from (const dealii::SparsityPattern &);
  template void
  SparsityPattern::copy_from (const dealii::CompressedSparsityPattern &);
  template void
  SparsityPattern::copy_from (const dealii::CompressedSetSparsityPattern &);
  template void
  SparsityPattern::copy_from (const dealii::CompressedSimpleSparsityPattern &);


  template void
  SparsityPattern::reinit (const Epetra_Map &,
			   const dealii::SparsityPattern &);
  template void
  SparsityPattern::reinit (const Epetra_Map &,
			   const dealii::CompressedSparsityPattern &);
  template void
  SparsityPattern::reinit (const Epetra_Map &,
			   const dealii::CompressedSetSparsityPattern &);
  template void
  SparsityPattern::reinit (const Epetra_Map &,
			   const dealii::CompressedSimpleSparsityPattern &);


  template void
  SparsityPattern::reinit (const Epetra_Map &,
			   const Epetra_Map &,
			   const dealii::SparsityPattern &);
  template void
  SparsityPattern::reinit (const Epetra_Map &,
			   const Epetra_Map &,
			   const dealii::CompressedSparsityPattern &);
  template void
  SparsityPattern::reinit (const Epetra_Map &,
			   const Epetra_Map &,
			   const dealii::CompressedSimpleSparsityPattern &);

}

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_USE_TRILINOS
