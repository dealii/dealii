//----------------------------  block_block_sparsity_pattern.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  block_block_sparsity_pattern.cc  ---------------------------

#include <base/memory_consumption.h>
#include <lac/block_sparsity_pattern.h>



BlockSparsityPattern::BlockSparsityPattern ()
		:
		rows (0),
		columns (0)
{};



BlockSparsityPattern::BlockSparsityPattern (const unsigned int r,
					    const unsigned int c)
		:
		rows (0),
		columns (0)
{
  reinit (r,c);
};



BlockSparsityPattern::~BlockSparsityPattern ()
{
				   // clear all memory
  reinit (0,0);
};



void
BlockSparsityPattern::reinit (const unsigned int r,
			      const unsigned int c)
{
				   // delete previous content
  for (unsigned int i=0; i<rows; ++i)
    for (unsigned int j=0; j<columns; ++j)
      {
	SparsityPattern * sp = sub_objects[i][j];
	sub_objects[i][j] = 0;
	delete sp;
      };
  sub_objects.clear ();

				   // set new sizes
  rows = r;
  columns = c;
  sub_objects.resize (rows, vector<SmartPointer<SparsityPattern> > (columns));

				   // allocate new objects
  for (unsigned int i=0; i<rows; ++i)
    for (unsigned int j=0; j<columns; ++j)
      sub_objects[i][j] = new SparsityPattern;
};



BlockSparsityPattern &
BlockSparsityPattern::operator = (const BlockSparsityPattern &bsp)
{
  Assert (rows == bsp.rows, ExcIncompatibleSizes(rows, bsp.rows));
  Assert (columns == bsp.columns, ExcIncompatibleSizes(columns, bsp.columns));
				   // copy objects
  for (unsigned int i=0; i<rows; ++i)
    for (unsigned int j=0; j<columns; ++j)
      *sub_objects[i][j] = *bsp.sub_objects[i][j];
				   // update index objects
  collect_sizes ();

  return *this;
};

  


void
BlockSparsityPattern::collect_sizes ()
{
  vector<unsigned int> row_sizes (rows);
  vector<unsigned int> col_sizes (columns);

				   // first find out the row sizes
				   // from the first block column
  for (unsigned int r=0; r<rows; ++r)
    row_sizes[r] = sub_objects[r][0]->n_rows();
				   // then check that the following
				   // block columns have the same
				   // sizes
  for (unsigned int c=1; c<columns; ++c)
    for (unsigned int r=0; r<rows; ++r)
      Assert (row_sizes[r] == sub_objects[r][c]->n_rows(),
	      ExcIncompatibleRowNumbers (r,0,r,c));

				   // finally initialize the row
				   // indices with this array
  row_indices.reinit (row_sizes);
  
  
				   // then do the same with the columns
  for (unsigned int c=0; c<columns; ++c)
    col_sizes[c] = sub_objects[0][c]->n_cols();
  for (unsigned int r=1; r<rows; ++r)
    for (unsigned int c=0; c<columns; ++c)
      Assert (col_sizes[c] == sub_objects[r][c]->n_cols(),
	      ExcIncompatibleRowNumbers (0,c,r,c));

				   // finally initialize the row
				   // indices with this array
  column_indices.reinit (col_sizes);
};



void
BlockSparsityPattern::compress ()
{
  for (unsigned int i=0; i<rows; ++i)
    for (unsigned int j=0; j<columns; ++j)
      sub_objects[i][j]->compress ();
};



bool
BlockSparsityPattern::empty () const
{
  for (unsigned int i=0; i<rows; ++i)
    for (unsigned int j=0; j<columns; ++j)
      if (sub_objects[i][j]->empty () == false)
	return false;
  return true;
};



unsigned int
BlockSparsityPattern::max_entries_per_row () const
{
  unsigned int max_entries = 0;
  for (unsigned int block_row=0; block_row<rows; ++block_row)
    {
      unsigned int this_row = 0;
      for (unsigned int c=0; c<columns; ++c)
	this_row += sub_objects[block_row][c]->max_entries_per_row ();

      if (this_row > max_entries)
	max_entries = this_row;
    };
  return max_entries;
};


unsigned int
BlockSparsityPattern::n_rows () const
{
				   // only count in first column, since
				   // all rows should be equivalent
  unsigned int count = 0;
  for (unsigned int r=0; r<rows; ++r)
    count += sub_objects[r][0]->n_rows();
  return count;
};



unsigned int
BlockSparsityPattern::n_cols () const
{
				   // only count in first row, since
				   // all rows should be equivalent
  unsigned int count = 0;
  for (unsigned int c=0; c<columns; ++c)
    count += sub_objects[0][c]->n_cols();
  return count;
};





unsigned int
BlockSparsityPattern::n_nonzero_elements () const
{
  unsigned int count = 0;
  for (unsigned int i=0; i<rows; ++i)
    for (unsigned int j=0; j<columns; ++j)
      count += sub_objects[i][j]->n_nonzero_elements ();
  return count;
};



bool
BlockSparsityPattern::is_compressed () const
{
  for (unsigned int i=0; i<rows; ++i)
    for (unsigned int j=0; j<columns; ++j)
      if (sub_objects[i][j]->is_compressed () == false)
	return false;
  return true;
};


unsigned int
BlockSparsityPattern::memory_consumption () const
{
  unsigned int mem = 0;
  mem += (MemoryConsumption::memory_consumption (rows) +
	  MemoryConsumption::memory_consumption (columns) +
	  MemoryConsumption::memory_consumption (sub_objects) +
	  MemoryConsumption::memory_consumption (row_indices) +
	  MemoryConsumption::memory_consumption (column_indices));
  for (unsigned int r=0; r<rows; ++r)
    for (unsigned int c=0; c<columns; ++c)
      mem += MemoryConsumption::memory_consumption (*sub_objects[r][c]);
  
  return mem;
};
