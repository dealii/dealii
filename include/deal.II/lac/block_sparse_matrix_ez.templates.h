// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef __deal2__block_sparse_matrix_ez_templates_h
#define __deal2__block_sparse_matrix_ez_templates_h


#include <deal.II/base/config.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/lac/block_sparse_matrix_ez.h>

DEAL_II_NAMESPACE_OPEN


template <typename number>
BlockSparseMatrixEZ<number>::BlockSparseMatrixEZ ()
{}



template <typename number>
BlockSparseMatrixEZ<number>::
BlockSparseMatrixEZ (const unsigned int rows,
                     const unsigned int cols)
  :
  row_indices (rows, 0),
  column_indices (cols, 0)
{}



//  template <typename number>
//  BlockSparseMatrixEZ<number>::~BlockSparseMatrixEZ ()
//  {
//                                 // delete previous content of
//                                 // the subobjects array
//    clear ();
//  };



template <typename number>
BlockSparseMatrixEZ<number> &
BlockSparseMatrixEZ<number>::
operator = (const BlockSparseMatrixEZ<number> &m)
{
  Assert (n_block_rows() == m.n_block_rows(),
          ExcDimensionMismatch(n_block_rows(), m.n_block_rows()));
  Assert (n_block_cols() == m.n_block_cols(),
          ExcDimensionMismatch(n_block_cols(), m.n_block_cols()));
  // this operator does not do
  // anything except than checking
  // whether the base objects want to
  // do something
  for (unsigned int r=0; r<n_block_rows(); ++r)
    for (unsigned int c=0; c<n_block_cols(); ++c)
      block(r,c) = m.block(r,c);
  return *this;
}



template <typename number>
BlockSparseMatrixEZ<number> &
BlockSparseMatrixEZ<number>::operator = (const double d)
{
  Assert (d==0, ExcScalarAssignmentOnlyForZeroValue());

  for (unsigned int r=0; r<n_block_rows(); ++r)
    for (unsigned int c=0; c<n_block_cols(); ++c)
      block(r,c) = 0;

  return *this;
}



template <typename number>
BlockSparseMatrixEZ<number>::BlockSparseMatrixEZ (
  const BlockSparseMatrixEZ<number> &m)
  :
  Subscriptor (m),
  row_indices(m.row_indices),
  column_indices(m.column_indices),
  blocks(m.blocks)
{}



template <typename number>
void
BlockSparseMatrixEZ<number>::reinit (const unsigned int rows,
                                     const unsigned int cols)
{
  row_indices.reinit(rows, 0);
  column_indices.reinit(cols, 0);
  blocks.reinit(rows, cols);
}



template <typename number>
void
BlockSparseMatrixEZ<number>::clear ()
{
  row_indices.reinit(0, 0);
  column_indices.reinit(0, 0);
  blocks.reinit(0, 0);
}



template <typename number>
bool
BlockSparseMatrixEZ<number>::empty () const
{
  for (unsigned int r=0; r<n_block_rows(); ++r)
    for (unsigned int c=0; c<n_block_cols(); ++c)
      if (block(r,c).empty () == false)
        return false;
  return true;
}



template <typename number>
void
BlockSparseMatrixEZ<number>::collect_sizes ()
{
  const unsigned int rows = n_block_rows();
  const unsigned int columns = n_block_cols();
  std::vector<size_type> row_sizes (rows);
  std::vector<size_type> col_sizes (columns);

  // first find out the row sizes
  // from the first block column
  for (unsigned int r=0; r<rows; ++r)
    row_sizes[r] = blocks[r][0].m();
  // then check that the following
  // block columns have the same
  // sizes
  for (unsigned int c=1; c<columns; ++c)
    for (unsigned int r=0; r<rows; ++r)
      Assert (row_sizes[r] == blocks[r][c].m(),
              ExcDimensionMismatch (row_sizes[r], blocks[r][c].m()));

  // finally initialize the row
  // indices with this array
  row_indices.reinit (row_sizes);


  // then do the same with the columns
  for (unsigned int c=0; c<columns; ++c)
    col_sizes[c] = blocks[0][c].n();
  for (unsigned int r=1; r<rows; ++r)
    for (unsigned int c=0; c<columns; ++c)
      Assert (col_sizes[c] == blocks[r][c].n(),
              ExcDimensionMismatch (col_sizes[c], blocks[r][c].n()));

  // finally initialize the row
  // indices with this array
  column_indices.reinit (col_sizes);
}




//  template <typename number>
//  unsigned int
//  BlockSparseMatrixEZ<number>::n_nonzero_elements () const
//  {
//    return sparsity_pattern->n_nonzero_elements ();
//  };



//  template <typename number>
//  unsigned int
//  BlockSparseMatrixEZ<number>::n_actually_nonzero_elements () const
//  {
//    unsigned int count = 0;
//    for (unsigned int i=0; i<rows; ++i)
//      for (unsigned int j=0; j<columns; ++j)
//        count += sub_objects[i][j]->n_actually_nonzero_elements ();
//    return count;
//  };



//  template <typename number>
//  std::size_t
//  BlockSparseMatrixEZ<number>::memory_consumption () const
//  {
//    std::size_t mem = sizeof(*this);
//    mem += MemoryConsumption::memory_consumption (sub_objects);
//    for (unsigned int r=0; r<rows; ++r)
//      for (unsigned int c=0; c<columns; ++c)
//        mem += MemoryConsumption::memory_consumption(*sub_objects[r][c]);
//    return mem;
//  };



DEAL_II_NAMESPACE_CLOSE

#endif // ifdef block_sparse_matrix_templates_h
