//----------------------------  block_sparse_matrix.templates.h  ---------------------------
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
//----------------------------  block_sparse_matrix.templates.h  ---------------------------
#ifndef __deal2__block_sparse_matrix_templates_h
#define __deal2__block_sparse_matrix_templates_h


#include <lac/block_sparse_matrix.h>



template <typename number, int  rows, int columns>
BlockSparseMatrix<number,rows,columns>::BlockSparseMatrix () :
		sparsity_pattern (0)
{};



template <typename number, int  rows, int columns>
BlockSparseMatrix<number,rows,columns>::
BlockSparseMatrix (const BlockSparsityPattern<rows,columns> &sparsity)
{
  reinit (sparsity);
};



template <typename number, int  rows, int columns>
BlockSparseMatrix<number,rows,columns> &
BlockSparseMatrix<number,rows,columns>::
operator = (const BlockSparseMatrix<number,rows,columns> &m) 
{
				   // this operator does not do
				   // anything except than checking
				   // whether the base objects want to
				   // do something
  for (unsigned int r=0; r<rows; ++r)
    for (unsigned int c=0; c<columns; ++c)
      block(r,c) = m.block(r,c);
  return *this;
};

 

template <typename number, int  rows, int columns>
void
BlockSparseMatrix<number,rows,columns>::reinit ()
{
  for (unsigned int r=0; r<rows; ++r)
    for (unsigned int c=0; c<columns; ++c)
      block(r,c).reinit ();
};



template <typename number, int  rows, int columns>
void
BlockSparseMatrix<number,rows,columns>::
reinit (const BlockSparsityPattern<rows,columns> &sparsity)
{
  sparsity_pattern = &sparsity;
  
  for (unsigned int r=0; r<rows; ++r)
    for (unsigned int c=0; c<columns; ++c)
      block(r,c).reinit (sparsity.block(r,c));
};


template <typename number, int  rows, int columns>
void
BlockSparseMatrix<number,rows,columns>::clear () 
{
  sparsity_pattern = 0;
  for (unsigned int r=0; r<rows; ++r)
    for (unsigned int c=0; c<columns; ++c)
      block(r,c).clear ();
};



template <typename number, int  rows, int columns>
bool
BlockSparseMatrix<number,rows,columns>::empty () const
{
  for (unsigned int r=0; r<rows; ++r)
    for (unsigned int c=0; c<columns; ++c)
      if (block(r,c).empty () == false)
	return false;
  return true;
};




template <typename number, int  rows, int columns>
unsigned int
BlockSparseMatrix<number,rows,columns>::n_nonzero_elements () const
{
  return sparsity_pattern->n_nonzero_elements ();
};



template <typename number, int  rows, int columns>
const BlockSparsityPattern<rows,columns> &
BlockSparseMatrix<number,rows,columns>::get_sparsity_pattern () const
{
  return *sparsity_pattern;
};



#endif // ifdef block_sparse_matrix_templates_h
