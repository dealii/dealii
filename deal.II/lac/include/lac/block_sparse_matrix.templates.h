//-------------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-------------------------------------------------------------------------------
#ifndef __deal2__block_sparse_matrix_templates_h
#define __deal2__block_sparse_matrix_templates_h


#include <base/config.h>
#include <base/memory_consumption.h>
#include <lac/block_sparse_matrix.h>



template <typename number>
BlockSparseMatrix<number>::BlockSparseMatrix ()
                :
		rows (0),
		columns (0),
		sparsity_pattern (0)
{}



template <typename number>
BlockSparseMatrix<number>::
BlockSparseMatrix (const BlockSparsityPattern &sparsity)
                :
		rows (0),
		columns (0)
{
  reinit (sparsity);
}



template <typename number>
BlockSparseMatrix<number>::~BlockSparseMatrix ()
{
  				   // delete previous content of
				   // the subobjects array
  for (unsigned int r=0; r<rows; ++r)
    for (unsigned int c=0; c<columns; ++c)
      {
	BlockType *p = sub_objects[r][c];
	sub_objects[r][c] = 0;
	delete p;
      }
}



template <typename number>
BlockSparseMatrix<number> &
BlockSparseMatrix<number>::
operator = (const BlockSparseMatrix<number> &m) 
{
  Assert (rows == m.rows, ExcDimensionMismatch(rows, m.rows));
  Assert (columns == m.columns, ExcDimensionMismatch(columns, m.columns));
				   // this operator does not do
				   // anything except than checking
				   // whether the base objects want to
				   // do something
  for (unsigned int r=0; r<rows; ++r)
    for (unsigned int c=0; c<columns; ++c)
      block(r,c) = m.block(r,c);

  return *this;
}

 

template <typename number>
void
BlockSparseMatrix<number>::reinit ()
{
  for (unsigned int r=0; r<rows; ++r)
    for (unsigned int c=0; c<columns; ++c)
      block(r,c).reinit ();
}



template <typename number>
void
BlockSparseMatrix<number>::
reinit (const BlockSparsityPattern &sparsity)
{
				   // first delete previous content of
				   // the subobjects array
  for (unsigned int r=0; r<rows; ++r)
    for (unsigned int c=0; c<columns; ++c)
      {
	BlockType *p = sub_objects[r][c];
	sub_objects[r][c] = 0;
	delete p;
      }
  
  sub_objects.clear ();

				   // then associate new sparsity
				   // pattern and resize
  sparsity_pattern = &sparsity;
  rows = sparsity.n_block_rows();
  columns = sparsity.n_block_cols();
  sub_objects.reinit (rows, columns);

				   // and reinitialize the blocks
  for (unsigned int r=0; r<rows; ++r)
    for (unsigned int c=0; c<columns; ++c)
      {
        BlockType *p = new SparseMatrix<number>();
        p->reinit (sparsity.block(r,c));
	sub_objects[r][c] = p;
      }
}



template <typename number>
void
BlockSparseMatrix<number>::clear () 
{
  sparsity_pattern = 0;
  for (unsigned int r=0; r<rows; ++r)
    for (unsigned int c=0; c<columns; ++c)
      block(r,c).clear ();
}



template <typename number>
bool
BlockSparseMatrix<number>::empty () const
{
  for (unsigned int r=0; r<rows; ++r)
    for (unsigned int c=0; c<columns; ++c)
      if (block(r,c).empty () == false)
	return false;

  return true;
}




template <typename number>
unsigned int
BlockSparseMatrix<number>::n_nonzero_elements () const
{
  return sparsity_pattern->n_nonzero_elements ();
}



template <typename number>
unsigned int
BlockSparseMatrix<number>::n_actually_nonzero_elements () const
{
  unsigned int count = 0;
  for (unsigned int i=0; i<rows; ++i)
    for (unsigned int j=0; j<columns; ++j)
      count += sub_objects[i][j]->n_actually_nonzero_elements ();

  return count;
}



template <typename number>
const BlockSparsityPattern &
BlockSparseMatrix<number>::get_sparsity_pattern () const
{
  return *sparsity_pattern;
}



template <typename number>
void 
BlockSparseMatrix<number>::
print_formatted (std::ostream       &out,
                 const unsigned int  precision,
                 const bool          scientific,
                 const unsigned int  width,
                 const char         *zero_string,
                 const double        denominator) const
{
  for (unsigned int r=0;r<rows;++r)
    for (unsigned int c=0;c<columns;++c)
      {
        out << "Component (" << r << "," << c << ")" << std::endl;
        block(r,c).print_formatted (out, precision, scientific,
                                    width, zero_string, denominator);
      }
}



template <typename number>
unsigned int
BlockSparseMatrix<number>::memory_consumption () const
{
  unsigned int mem = sizeof(*this);
  mem += MemoryConsumption::memory_consumption (sub_objects);
  for (unsigned int r=0; r<rows; ++r)
    for (unsigned int c=0; c<columns; ++c)
      mem += MemoryConsumption::memory_consumption(*sub_objects[r][c]);

  return mem;
}



#endif // ifdef block_sparse_matrix_templates_h
