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

#ifdef DEAL_II_USE_TRILINOS

#  include <base/utilities.h>
#  include <lac/sparsity_pattern.h>
#  include <lac/compressed_sparsity_pattern.h>
#  include <lac/compressed_set_sparsity_pattern.h>
#  include <lac/compressed_simple_sparsity_pattern.h>

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
		  compressed (true)
  {
    column_space_map = std::auto_ptr<Epetra_Map>
      (new Epetra_Map (0, 0, Utilities::Trilinos::comm_self()));
    graph = std::auto_ptr<Epetra_FECrsGraph>
      (new Epetra_FECrsGraph(View, *column_space_map, *column_space_map, 0));
    graph->FillComplete();
  }


  SparsityPattern::SparsityPattern (const Epetra_Map  &input_map,
				    const unsigned int n_entries_per_row)
  {
    reinit (input_map, input_map, n_entries_per_row);
  }



  SparsityPattern::SparsityPattern (const Epetra_Map                &input_map,
				    const std::vector<unsigned int> &n_entries_per_row)
  {
    reinit (input_map, input_map, n_entries_per_row);
  }



  SparsityPattern::SparsityPattern (const Epetra_Map  &input_row_map,
				    const Epetra_Map  &input_col_map,
				    const unsigned int n_entries_per_row)
  {
    reinit (input_row_map, input_col_map, n_entries_per_row);
  }



  SparsityPattern::SparsityPattern (const Epetra_Map                &input_row_map,
				    const Epetra_Map                &input_col_map,
				    const std::vector<unsigned int> &n_entries_per_row)
  {
    reinit (input_row_map, input_col_map, n_entries_per_row);
  }



  SparsityPattern::SparsityPattern (const unsigned int m,
				    const unsigned int n,
				    const unsigned int n_entries_per_row)
  {
    reinit (m, n, n_entries_per_row);
  }



  SparsityPattern::SparsityPattern (const unsigned int               m,
				    const unsigned int               n,
				    const std::vector<unsigned int> &n_entries_per_row)
  {
    reinit (m, n, n_entries_per_row);
  }


				   // Copy function only works if the
				   // sparsity pattern is empty.
  SparsityPattern::SparsityPattern (const SparsityPattern &input_sparsity)
  		  :
                  Subscriptor(),
                  column_space_map (std::auto_ptr<Epetra_Map>
				    (new Epetra_Map(0, 0, Utilities::Trilinos::comm_self()))),
                  compressed (false),
                  graph (std::auto_ptr<Epetra_FECrsGraph>
			 (new Epetra_FECrsGraph(View, *column_space_map,
						*column_space_map, 0)))
  {
    Assert (input_sparsity.n_rows() == 0,
	    ExcMessage ("Copy constructor only works for empty sparsity patterns."));
  }



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
    const Epetra_Map rows (m, 0, Utilities::Trilinos::comm_self());
    const Epetra_Map columns (n, 0, Utilities::Trilinos::comm_self());

    reinit (rows, columns, n_entries_per_row);
  }


  void
  SparsityPattern::reinit (const Epetra_Map   &input_row_map,
			   const Epetra_Map   &input_col_map,
			   const unsigned int  n_entries_per_row)
  {
    column_space_map = std::auto_ptr<Epetra_Map> (new Epetra_Map (input_col_map));
    graph.reset();
    compressed = false;

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
    if (input_row_map.Comm().NumProc() > 1)
      graph = std::auto_ptr<Epetra_FECrsGraph>
	(new Epetra_FECrsGraph(Copy, input_row_map, n_entries_per_row, false));
    else
      graph = std::auto_ptr<Epetra_FECrsGraph>
	(new Epetra_FECrsGraph(Copy, input_row_map, input_col_map,
			       n_entries_per_row, false));
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
    const Epetra_Map rows (m, 0, Utilities::Trilinos::comm_self());
    const Epetra_Map columns (n, 0, Utilities::Trilinos::comm_self());

    reinit (rows, columns, n_entries_per_row);
  }



  void
  SparsityPattern::reinit (const Epetra_Map   &input_row_map,
			   const Epetra_Map   &input_col_map,
			   const std::vector<unsigned int> &n_entries_per_row)
  {
    Assert (n_entries_per_row.size() ==
	      static_cast<unsigned int>(input_row_map.NumGlobalElements()),
	    ExcDimensionMismatch (n_entries_per_row.size(),
				  input_row_map.NumGlobalElements()));

    column_space_map = std::auto_ptr<Epetra_Map> (new Epetra_Map (input_col_map));
    graph.reset();
    compressed = false;

    if (input_row_map.Comm().NumProc() > 1)
      graph = std::auto_ptr<Epetra_FECrsGraph>
	(new Epetra_FECrsGraph(Copy, input_row_map,
			       n_entries_per_row[input_row_map.MinMyGID()],
			       false));
    else
      graph = std::auto_ptr<Epetra_FECrsGraph>
	(new Epetra_FECrsGraph(Copy, input_row_map, input_col_map,
			       n_entries_per_row[input_row_map.MinMyGID()],
			       false));
  }



  template <typename SparsityType>
  void
  SparsityPattern::reinit (const Epetra_Map   &input_map,
			   const SparsityType &sp,
			   const bool          exchange_data)
  {
    reinit (input_map, input_map, sp, exchange_data);
  }



  template <typename SparsityType>
  void
  SparsityPattern::reinit (const Epetra_Map   &input_row_map,
			   const Epetra_Map   &input_col_map,
			   const SparsityType &sp,
			   const bool          exchange_data)
  {
    Assert (sp.n_rows() ==
	      static_cast<unsigned int>(input_row_map.NumGlobalElements()),
	    ExcDimensionMismatch (sp.n_rows(),
				  input_row_map.NumGlobalElements()));
    Assert (sp.n_cols() ==
	      static_cast<unsigned int>(input_col_map.NumGlobalElements()),
	    ExcDimensionMismatch (sp.n_cols(),
				  input_col_map.NumGlobalElements()));
    Assert (exchange_data == false, ExcNotImplemented());

    column_space_map = std::auto_ptr<Epetra_Map> (new Epetra_Map (input_col_map));
    graph.reset();
    compressed = false;

    Assert (input_row_map.LinearMap() == true,
	    ExcMessage ("This function is not efficient if the map is not contiguous."));

    std::vector<int> n_entries_per_row(input_row_map.MaxMyGID()-
				       input_row_map.MinMyGID() + 1);

    for (unsigned int row=input_row_map.MinMyGID();
	 row<static_cast<unsigned int>(input_row_map.MaxMyGID()+1);
	 ++row)
      n_entries_per_row[row-input_row_map.MinMyGID()] = sp.row_length(row);

    if (input_row_map.Comm().NumProc() > 1)
      graph = std::auto_ptr<Epetra_FECrsGraph>
	(new Epetra_FECrsGraph(Copy, input_row_map,
			       n_entries_per_row[0],
			       false));
    else
      graph = std::auto_ptr<Epetra_FECrsGraph>
	(new Epetra_FECrsGraph(Copy, input_row_map, input_col_map,
			       n_entries_per_row[0],
			       false));

    Assert (graph->NumGlobalRows() == (int)sp.n_rows(),
    	    ExcDimensionMismatch (graph->NumGlobalRows(),
    				  sp.n_rows()));

    std::vector<int>   row_indices;
    const unsigned int n_rows = sp.n_rows();

    for (unsigned int row=0; row<n_rows; ++row)
      if ( input_row_map.MyGID(row) )
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



  template <>
  void
  SparsityPattern::reinit (const Epetra_Map   &input_row_map,
			   const Epetra_Map   &input_col_map,
			   const CompressedSimpleSparsityPattern &sp,
			   const bool          exchange_data)
  {
    Assert (sp.n_rows() ==
	      static_cast<unsigned int>(input_row_map.NumGlobalElements()),
	    ExcDimensionMismatch (sp.n_rows(),
				  input_row_map.NumGlobalElements()));
    Assert (sp.n_cols() ==
	      static_cast<unsigned int>(input_col_map.NumGlobalElements()),
	    ExcDimensionMismatch (sp.n_cols(),
				  input_col_map.NumGlobalElements()));

    column_space_map = std::auto_ptr<Epetra_Map> (new Epetra_Map (input_col_map));
    graph.reset();
    compressed = false;

    Assert (input_row_map.LinearMap() == true,
	    ExcMessage ("This function is not efficient if the map is not contiguous."));

    std::vector<int> n_entries_per_row(input_row_map.MaxMyGID()-
				       input_row_map.MinMyGID() + 1);

    for (unsigned int row=input_row_map.MinMyGID();
	 row<static_cast<unsigned int>(input_row_map.MaxMyGID()+1);
	 ++row)
      n_entries_per_row[row-input_row_map.MinMyGID()] = sp.row_length(row);

    if (input_row_map.Comm().NumProc() > 1)
      graph = std::auto_ptr<Epetra_FECrsGraph>
	(new Epetra_FECrsGraph(Copy, input_row_map,
			       n_entries_per_row[0],
			       false));
    else
      graph = std::auto_ptr<Epetra_FECrsGraph>
	(new Epetra_FECrsGraph(Copy, input_row_map, input_col_map,
			       n_entries_per_row[0],
			       false));

    Assert (graph->NumGlobalRows() == (int)sp.n_rows(),
    	    ExcDimensionMismatch (graph->NumGlobalRows(),
    				  sp.n_rows()));

    const unsigned int n_rows = sp.n_rows();
    std::vector<int>   row_indices;

				// Include possibility to exchange data
				// since CompressedSimpleSparsityPattern is
				// able to do so
    for (unsigned int row=0; row<n_rows; ++row)
      if (input_row_map.MyGID(row) )
	{
	  const int row_length = sp.row_length(row);
	  row_indices.resize (row_length, -1);

	  for (int col=0; col < row_length; ++col)
	    row_indices[col] = sp.column_number (row, col);

	  graph->Epetra_CrsGraph::InsertGlobalIndices (row, row_length,
						       &row_indices[0]);
	}
      else if ( exchange_data && sp.row_index_set().is_element(row) )
	{
	  const int row_length = sp.row_length(row);
	  row_indices.resize (row_length, -1);

	  for (int col=0; col < row_length; ++col)
	    row_indices[col] = sp.column_number (row, col);

	  graph->InsertGlobalIndices (1, (int*)&row, row_length, &row_indices[0]);
	}

    compress();
  }



  template<>
  void
  SparsityPattern::reinit (const Epetra_Map   &input_row_map,
			   const Epetra_Map   &input_col_map,
			   const CompressedSetSparsityPattern &sp,
			   const bool          exchange_data)
  {
    Assert (exchange_data == false, ExcNotImplemented());
    Assert (sp.n_rows() ==
	      static_cast<unsigned int>(input_row_map.NumGlobalElements()),
	    ExcDimensionMismatch (sp.n_rows(),
				  input_row_map.NumGlobalElements()));
    Assert (sp.n_cols() ==
	      static_cast<unsigned int>(input_col_map.NumGlobalElements()),
	    ExcDimensionMismatch (sp.n_cols(),
				  input_col_map.NumGlobalElements()));

    column_space_map = std::auto_ptr<Epetra_Map> (new Epetra_Map (input_col_map));
    graph.reset();
    compressed = false;

    Assert (input_row_map.LinearMap() == true,
	    ExcMessage ("This function is not efficient if the map is not contiguous."));

    std::vector<int> n_entries_per_row(input_row_map.MaxMyGID()-
				       input_row_map.MinMyGID() + 1);
    for (unsigned int row=input_row_map.MinMyGID();
	 row<static_cast<unsigned int>(input_row_map.MaxMyGID()+1);
	 ++row)
      {
	n_entries_per_row[row-input_row_map.MinMyGID()] = sp.row_length(row);
      }


    if (input_row_map.Comm().NumProc() > 1)
      graph = std::auto_ptr<Epetra_FECrsGraph>
	(new Epetra_FECrsGraph(Copy, input_row_map,
			       n_entries_per_row[0],
			       false));
    else
      graph = std::auto_ptr<Epetra_FECrsGraph>
	(new Epetra_FECrsGraph(Copy, input_row_map, input_col_map,
			       n_entries_per_row[0],
			       false));

    Assert (graph->NumGlobalRows() == (int)sp.n_rows(),
    	    ExcDimensionMismatch (graph->NumGlobalRows(),
    				  sp.n_rows()));


    const unsigned int n_rows = sp.n_rows();
    std::vector<int>   row_indices;

    for (unsigned int row=0; row<n_rows; ++row)
     if (exchange_data || input_row_map.MyGID(row))
	{
	  const int row_length = sp.row_length(row);
	  if (row_length == 0)
	    continue;

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



  void
  SparsityPattern::operator = (const SparsityPattern &)
  {
    Assert (false, ExcNotImplemented());
  }



  template <typename SparsityType>
  void
  SparsityPattern::copy_from (const SparsityType &sp)
  {
    const Epetra_Map rows (sp.n_rows(), 0, Utilities::Trilinos::comm_self());
    const Epetra_Map columns (sp.n_cols(), 0, Utilities::Trilinos::comm_self());

    reinit (rows, columns, sp);
  }



  void
  SparsityPattern::clear ()
  {
				  // When we clear the matrix, reset
				  // the pointer and generate an
				  // empty sparsity pattern.
    column_space_map.reset();
    graph.reset();

     column_space_map = std::auto_ptr<Epetra_Map>
       (new Epetra_Map (0, 0, Utilities::Trilinos::comm_self()));
    graph = std::auto_ptr<Epetra_FECrsGraph>
      (new Epetra_FECrsGraph(View, *column_space_map, *column_space_map, 0));
    graph->FillComplete();

    compressed = true;
  }



  void
  SparsityPattern::compress ()
  {
    int ierr;
    Assert (&* column_space_map != 0, ExcInternalError());
    ierr = graph->GlobalAssemble (*column_space_map,
				  static_cast<const Epetra_Map&>(graph->RangeMap()),
				  true);

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
      n_cols = column_space_map->NumGlobalElements();

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
  SparsityPattern::print (std::ostream &out,
			  const bool    write_extended_trilinos_info) const
  {
    if (write_extended_trilinos_info)
      out << *graph;
    else
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
			   const dealii::SparsityPattern &,
			   bool);
  template void
  SparsityPattern::reinit (const Epetra_Map &,
			   const dealii::CompressedSparsityPattern &,
			   bool);
  template void
  SparsityPattern::reinit (const Epetra_Map &,
			   const dealii::CompressedSetSparsityPattern &,
			   bool);
  template void
  SparsityPattern::reinit (const Epetra_Map &,
			   const dealii::CompressedSimpleSparsityPattern &,
			   bool);


  template void
  SparsityPattern::reinit (const Epetra_Map &,
			   const Epetra_Map &,
			   const dealii::SparsityPattern &,
			   bool);
  template void
  SparsityPattern::reinit (const Epetra_Map &,
			   const Epetra_Map &,
			   const dealii::CompressedSparsityPattern &,
			   bool);

}

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_USE_TRILINOS
