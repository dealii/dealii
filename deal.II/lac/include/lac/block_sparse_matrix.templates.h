//----------------------------  block_sparse_matrix.templates.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  block_sparse_matrix.templates.h  ---------------------------
#ifndef __deal2__block_sparse_matrix_templates_h
#define __deal2__block_sparse_matrix_templates_h


#include <base/memory_consumption.h>
#include <lac/block_sparse_matrix.h>



template <typename number>
BlockSparseMatrix<number>::BlockSparseMatrix () :
		rows (0),
		columns (0),
		sparsity_pattern (0)
{};



template <typename number>
BlockSparseMatrix<number>::
BlockSparseMatrix (const BlockSparsityPattern &sparsity) :
		rows (0),
		columns (0)
{
  reinit (sparsity);
};



template <typename number>
BlockSparseMatrix<number>::~BlockSparseMatrix ()
{
  				   // delete previous content of
				   // the subobjects array
  for (unsigned int r=0; r<rows; ++r)
    for (unsigned int c=0; c<columns; ++c)
      {
	SparseMatrix<number> *p = sub_objects[r][c];
	sub_objects[r][c] = 0;
	delete p;
      };
};



template <typename number>
BlockSparseMatrix<number> &
BlockSparseMatrix<number>::
operator = (const BlockSparseMatrix<number> &m) 
{
  Assert (rows == m.rows, ExcIncompatibleObjects());
  Assert (columns == m.columns, ExcIncompatibleObjects());
				   // this operator does not do
				   // anything except than checking
				   // whether the base objects want to
				   // do something
  for (unsigned int r=0; r<rows; ++r)
    for (unsigned int c=0; c<columns; ++c)
      block(r,c) = m.block(r,c);
  return *this;
};

 

template <typename number>
void
BlockSparseMatrix<number>::reinit ()
{
  for (unsigned int r=0; r<rows; ++r)
    for (unsigned int c=0; c<columns; ++c)
      block(r,c).reinit ();
};



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
	SparseMatrix<number> *p = sub_objects[r][c];
	sub_objects[r][c] = 0;
	delete p;
      };
  sub_objects.clear ();

				   // then associate new sparsity
				   // pattern and resize
  sparsity_pattern = &sparsity;
  rows = sparsity.n_block_rows();
  columns = sparsity.n_block_cols();
  sub_objects = std::vector<std::vector<SmartPointer<SparseMatrix<number> > > >
		(rows, std::vector<SmartPointer<SparseMatrix<number> > > (columns, 0));

				   // and reinitialize the blocks
  for (unsigned int r=0; r<rows; ++r)
    for (unsigned int c=0; c<columns; ++c)
      {
	sub_objects[r][c] = new SparseMatrix<number>();
	block(r,c).reinit (sparsity.block(r,c));
      };
};



template <typename number>
void
BlockSparseMatrix<number>::clear () 
{
  sparsity_pattern = 0;
  for (unsigned int r=0; r<rows; ++r)
    for (unsigned int c=0; c<columns; ++c)
      block(r,c).clear ();
};



template <typename number>
bool
BlockSparseMatrix<number>::empty () const
{
  for (unsigned int r=0; r<rows; ++r)
    for (unsigned int c=0; c<columns; ++c)
      if (block(r,c).empty () == false)
	return false;
  return true;
};




template <typename number>
unsigned int
BlockSparseMatrix<number>::n_nonzero_elements () const
{
  return sparsity_pattern->n_nonzero_elements ();
};



template <typename number>
const BlockSparsityPattern &
BlockSparseMatrix<number>::get_sparsity_pattern () const
{
  return *sparsity_pattern;
};



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
};



#endif // ifdef block_sparse_matrix_templates_h
