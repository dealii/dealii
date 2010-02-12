//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2008, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__block_sparse_matrix_templates_h
#define __deal2__block_sparse_matrix_templates_h


#include <base/config.h>
#include <base/memory_consumption.h>
#include <lac/block_sparse_matrix.h>

DEAL_II_NAMESPACE_OPEN


template <typename number>
BlockSparseMatrix<number>::BlockSparseMatrix ()
{}



template <typename number>
BlockSparseMatrix<number>::
BlockSparseMatrix (const BlockSparsityPattern &sparsity)
{
  reinit (sparsity);
}



template <typename number>
BlockSparseMatrix<number>::~BlockSparseMatrix ()
{
  				   // delete previous content of
				   // the subobjects array
  clear ();
  sparsity_pattern = 0;
}



template <typename number>
BlockSparseMatrix<number> &
BlockSparseMatrix<number>::
operator = (const BlockSparseMatrix<number> &m) 
{
  AssertDimension (this->row_block_indices, m.row_block_indices);
  AssertDimension (this->column_block_indices, m.column_block_indices);
  
                                   // this operator does not do
				   // anything except than checking
				   // whether the base objects want to
				   // do something
  for (unsigned int r=0; r<this->n_block_rows(); ++r)
    for (unsigned int c=0; c<this->n_block_cols(); ++c)
      this->block(r,c) = m.block(r,c);

  return *this;
}



template <typename number>
void
BlockSparseMatrix<number>::clear () 
{
  BlockMatrixBase<SparseMatrix<number> >::clear();
  sparsity_pattern = 0;
}



template <typename number>
void
BlockSparseMatrix<number>::
reinit (const BlockSparsityPattern &sparsity)
{
				   // first delete previous content of
				   // the subobjects array and delete
				   // the table completely
  clear ();

				   // then associate new sparsity
				   // pattern and resize
  sparsity_pattern = &sparsity;
  
  this->row_block_indices    = sparsity.row_indices;
  this->column_block_indices = sparsity.column_indices;

  this->sub_objects.reinit (sparsity.n_block_rows(),
                            sparsity.n_block_cols());

				   // and reinitialize the blocks
  for (unsigned int r=0; r<this->n_block_rows(); ++r)
    for (unsigned int c=0; c<this->n_block_cols(); ++c)
      {
        BlockType *p = new SparseMatrix<number>();
        p->reinit (sparsity.block(r,c));
	this->sub_objects[r][c] = p;
      }
}



template <typename number>
bool
BlockSparseMatrix<number>::empty () const
{
  for (unsigned int r=0; r<this->n_block_rows(); ++r)
    for (unsigned int c=0; c<this->n_block_cols(); ++c)
      if (this->block(r,c).empty () == false)
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
BlockSparseMatrix<number>::n_actually_nonzero_elements (const double threshold) const
{
  unsigned int count = 0;
  for (unsigned int i=0; i<this->n_block_rows(); ++i)
    for (unsigned int j=0; j<this->n_block_cols(); ++j)
      count += this->sub_objects[i][j]->n_actually_nonzero_elements (threshold);

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
  for (unsigned int r=0;r<this->n_block_rows();++r)
    for (unsigned int c=0;c<this->n_block_cols();++c)
      {
        out << "Component (" << r << "," << c << ")" << std::endl;
        this->block(r,c).print_formatted (out, precision, scientific,
                                          width, zero_string, denominator);
      }
}



template <typename number>
unsigned int
BlockSparseMatrix<number>::memory_consumption () const
{
  unsigned int mem = sizeof(*this);
  mem += MemoryConsumption::memory_consumption (this->sub_objects);
  for (unsigned int r=0; r<this->n_block_rows(); ++r)
    for (unsigned int c=0; c<this->n_block_cols(); ++c)
      mem += MemoryConsumption::memory_consumption(*this->sub_objects[r][c]);

  return mem;
}



DEAL_II_NAMESPACE_CLOSE

#endif // ifdef block_sparse_matrix_templates_h
