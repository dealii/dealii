//----------------------------  filtered_matrix.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2001, 2002 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  filtered_matrix.cc  ---------------------------


#include <lac/filtered_matrix.templates.h>


template <>
void
FilteredMatrix<SparseMatrix<double>,Vector<double> >::
get_column_entries (const unsigned int           index,
		    std::vector<IndexValuePair> &column_entries,
		    const bool                   matrix_is_symmetric) const
{
				   // depending on whether the matrix
				   // can be assumed symmetric or not,
				   // either use a fast or a slow
				   // algorithm
  if (matrix_is_symmetric == true)
				     // ok, matrix is symmetric. we
				     // may determine the matrix
				     // entries in this column by
				     // looking at the matrix entries
				     // in this row which is
				     // significantly faster since we
				     // can traverse them linearly and
				     // do not have to check each row
				     // for the possible existence of
				     // a matrix entry
    {
      const unsigned int *
	col_nums   = &(matrix->get_sparsity_pattern().get_column_numbers()
		       [matrix->get_sparsity_pattern().get_rowstart_indices()[index]]);
      const unsigned int
	row_length = matrix->get_sparsity_pattern().row_length(index);

      for (unsigned int i=0; i<row_length; ++i)
	{
	  const unsigned int c = *(col_nums+i);

					   // if not diagonal entry,
					   // add to list
	  if (c != index)
	    column_entries.push_back (std::make_pair(c, (*matrix)(c,index)));
	};
    }
  else
    {
				       // otherwise check each row for
				       // occurrence of an entry in
				       // this column
      for (unsigned int row=0; row<n(); ++row)
	if (row != index)
	  {
	    const unsigned int
	      global_index = matrix->get_sparsity_pattern()(row,index);
	    if (global_index != SparsityPattern::invalid_entry)
	      column_entries.push_back (std::make_pair(row,
						       (*matrix)(row,index)));
	  };
    };
};



template <>
void
FilteredMatrix<BlockSparseMatrix<double>,BlockVector<double> >::
get_column_entries (const unsigned int           /*index*/,
		    std::vector<IndexValuePair> &/*column_entries*/,
		    const bool                   /*matrix_is_symmetric*/) const
{
				   // presently not implemented, but
				   // should be fairly simple to do
  Assert (false, ExcNotImplemented());
};




template <>
void
FilteredMatrix<SparseMatrix<double>,Vector<double> >::
allocate_tmp_vector () 
{
  tmp_mutex.acquire ();
  tmp_vector.reinit (matrix->n(), true);
  tmp_mutex.release ();
};



template <>
void
FilteredMatrix<SparseMatrix<float>,Vector<float> >::
allocate_tmp_vector () 
{
  tmp_mutex.acquire ();
  tmp_vector.reinit (matrix->n(), true);
  tmp_mutex.release ();
};



template <>
void
FilteredMatrix<BlockSparseMatrix<double>,BlockVector<double> >::
allocate_tmp_vector () 
{
  std::vector<unsigned int> block_sizes (matrix->n_block_rows());
  for (unsigned int i=0; i<block_sizes.size(); ++i)
    block_sizes[i] = matrix->block(i,i).n();
  
  tmp_mutex.acquire ();
  tmp_vector.reinit (block_sizes, true);
  tmp_mutex.release ();
};



template <>
void
FilteredMatrix<BlockSparseMatrix<float>,BlockVector<float> >::
allocate_tmp_vector () 
{
  std::vector<unsigned int> block_sizes (matrix->n_block_rows());
  for (unsigned int i=0; i<block_sizes.size(); ++i)
    block_sizes[i] = matrix->block(i,i).n();
  
  tmp_mutex.acquire ();
  tmp_vector.reinit (block_sizes, true);
  tmp_mutex.release ();
};


template class FilteredMatrix<SparseMatrix<double>,Vector<double> >;
template class FilteredMatrix<BlockSparseMatrix<double>,BlockVector<double> >;
