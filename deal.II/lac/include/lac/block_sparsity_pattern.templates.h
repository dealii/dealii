//----------------------------  block_sparsity_pattern.templates.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  block_sparsity_pattern.templates.h  ---------------------------
#ifndef __deal2__block_sparsity_pattern_templates_h
#define __deal2__block_sparsity_pattern_templates_h


#include <lac/block_sparsity_pattern.h>


template <int rows, int columns>
BlockSparsityPattern<rows,columns>::BlockSparsityPattern () 
{
  Assert (rows>0,    ExcInvalidSize (rows));
  Assert (columns>0, ExcInvalidSize (columns));
};

  

template <int rows, int columns>
BlockSparsityPattern<rows,columns> &
BlockSparsityPattern<rows,columns>::operator = (const BlockSparsityPattern<rows,columns> &bsp)
{
				   // copy objects
  for (unsigned int i=0; i<rows; ++i)
    for (unsigned int j=0; j<columns; ++j)
      sub_objects[i][j] = bsp.sub_objects[i][j];
				   // update index objects
  collect_sizes ();

  return *this;
};

  


template <int rows, int columns>
void
BlockSparsityPattern<rows,columns>::collect_sizes ()
{
  vector<unsigned int> row_sizes (rows);
  vector<unsigned int> col_sizes (columns);

				   // first find out the row sizes
				   // from the first block column
  for (unsigned int r=0; r<rows; ++r)
    row_sizes[r] = sub_objects[r][0].n_rows();
				   // then check that the following
				   // block columns have the same
				   // sizes
  for (unsigned int c=1; c<columns; ++c)
    for (unsigned int r=0; r<rows; ++r)
      Assert (row_sizes[r] == sub_objects[r][c].n_rows(),
	      ExcIncompatibleRowNumbers (r,0,r,c));

				   // finally initialize the row
				   // indices with this array
  row_indices.reinit (row_sizes);
  
  
				   // then do the same with the columns
  for (unsigned int c=0; c<columns; ++c)
    col_sizes[c] = sub_objects[0][c].n_cols();
  for (unsigned int r=1; r<rows; ++r)
    for (unsigned int c=0; c<columns; ++c)
      Assert (col_sizes[c] == sub_objects[r][c].n_cols(),
	      ExcIncompatibleRowNumbers (0,c,r,c));

				   // finally initialize the row
				   // indices with this array
  column_indices.reinit (col_sizes);
};



template <int rows, int columns>
void
BlockSparsityPattern<rows,columns>::compress ()
{
  for (unsigned int i=0; i<rows; ++i)
    for (unsigned int j=0; j<columns; ++j)
      sub_objects[i][j].compress ();
};



template <int rows, int columns>
bool
BlockSparsityPattern<rows,columns>::empty () const
{
  for (unsigned int i=0; i<rows; ++i)
    for (unsigned int j=0; j<columns; ++j)
      if (sub_objects[i][j].empty () == false)
	return false;
  return true;
};



template <int rows, int columns>
unsigned int
BlockSparsityPattern<rows,columns>::max_entries_per_row () const
{
  unsigned int max_entries = 0;
  for (unsigned int block_row=0; block_row<rows; ++block_row)
    {
      unsigned int this_row = 0;
      for (unsigned int c=0; c<columns; ++c)
	this_row += sub_objects[block_row][c].max_entries_per_row ();

      if (this_row > max_entries)
	max_entries = this_row;
    };
  return max_entries;
};


template <int rows, int columns>
unsigned int
BlockSparsityPattern<rows,columns>::n_rows () const
{
				   // only count in first column, since
				   // all rows should be equivalent
  unsigned int count = 0;
  for (unsigned int r=0; r<rows; ++r)
    count += sub_objects[r][0].n_rows();
  return count;
};



template <int rows, int columns>
unsigned int
BlockSparsityPattern<rows,columns>::n_cols () const
{
				   // only count in first row, since
				   // all rows should be equivalent
  unsigned int count = 0;
  for (unsigned int c=0; c<columns; ++c)
    count += sub_objects[0][c].n_cols();
  return count;
};





template <int rows, int columns>
unsigned int
BlockSparsityPattern<rows,columns>::n_nonzero_elements () const
{
  unsigned int count = 0;
  for (unsigned int i=0; i<rows; ++i)
    for (unsigned int j=0; j<columns; ++j)
      count += sub_objects[i][j].n_nonzero_elements ();
  return count;
};



template <int rows, int columns>
bool
BlockSparsityPattern<rows,columns>::is_compressed () const
{
  for (unsigned int i=0; i<rows; ++i)
    for (unsigned int j=0; j<columns; ++j)
      if (sub_objects[i][j].is_compressed () == false)
	return false;
  return true;
};



#endif // ifdef block_sparsity_pattern_templates_h
