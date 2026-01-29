// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2002 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_block_sparse_matrix_ez_templates_h
#define dealii_block_sparse_matrix_ez_templates_h


#include <deal.II/base/config.h>

#include <deal.II/base/memory_consumption.h>

#include <deal.II/lac/block_sparse_matrix_ez.h>

DEAL_II_NAMESPACE_OPEN

template <typename number>
BlockSparseMatrixEZ<number>::BlockSparseMatrixEZ(const unsigned int rows,
                                                 const unsigned int cols)
{
  reinit(rows, cols);
}



template <typename number>
BlockSparseMatrixEZ<number>::BlockSparseMatrixEZ(
  const BlockSparseMatrixEZ<number> &m)
{
  Assert(
    m.empty(),
    ExcMessage(
      "This operator can only be called if the provided right "
      "hand side is a block matrix with empty blocks. This operator can not be "
      "used to copy a non-empty matrix."));
  this->reinit(m.n_block_rows(), m.n_block_cols());
}



template <typename number>
BlockSparseMatrixEZ<number> &
BlockSparseMatrixEZ<number>::operator=(const BlockSparseMatrixEZ<number> &m)
{
  Assert(this->n_block_rows() == m.n_block_rows(),
         ExcDimensionMismatch(this->n_block_rows(), m.n_block_rows()));
  Assert(this->n_block_cols() == m.n_block_cols(),
         ExcDimensionMismatch(this->n_block_cols(), m.n_block_cols()));
  // this operator does not do
  // anything except than checking
  // whether the base objects want to
  // do something
  for (unsigned int r = 0; r < this->n_block_rows(); ++r)
    for (unsigned int c = 0; c < this->n_block_cols(); ++c)
      *(this->sub_objects(r, c)) = m.block(r, c);
  return *this;
}



template <typename number>
BlockSparseMatrixEZ<number> &
BlockSparseMatrixEZ<number>::operator=(const double d)
{
  (void)d;
  Assert(d == 0, ExcScalarAssignmentOnlyForZeroValue());

  for (unsigned int r = 0; r < this->n_block_rows(); ++r)
    for (unsigned int c = 0; c < this->n_block_cols(); ++c)
      *(this->sub_objects(r, c)) = 0;

  return *this;
}



template <typename number>
void
BlockSparseMatrixEZ<number>::reinit(const unsigned int br,
                                    const unsigned int bc)
{
  clear(); // calls BlockMatrixBase::clear(), deletes old blocks
  this->sub_objects.reinit(br, bc);

  for (unsigned int r = 0; r < br; ++r)
    for (unsigned int c = 0; c < bc; ++c)
      this->sub_objects[r][c] = new SparseMatrixEZ<number>();

  this->row_block_indices    = BlockIndices(br, 0);
  this->column_block_indices = BlockIndices(bc, 0);
}



template <typename number>
void
BlockSparseMatrixEZ<number>::clear()
{
  BlockMatrixBase<SparseMatrixEZ<number>>::clear();
}



template <typename number>
bool
BlockSparseMatrixEZ<number>::empty() const
{
  for (unsigned int r = 0; r < this->n_block_rows(); ++r)
    for (unsigned int c = 0; c < this->n_block_cols(); ++c)
      if (this->sub_objects(r, c)->empty() == false)
        return false;
  return true;
}



template <typename number>
void
BlockSparseMatrixEZ<number>::collect_sizes()
{
  const unsigned int     rows    = this->n_block_rows();
  const unsigned int     columns = this->n_block_cols();
  std::vector<size_type> row_sizes(rows);
  std::vector<size_type> col_sizes(columns);

  // first find out the row sizes
  // from the first block column
  for (unsigned int r = 0; r < rows; ++r)
    row_sizes[r] = this->sub_objects[r][0]->m();
  // then check that the following
  // block columns have the same
  // sizes
  for (unsigned int c = 1; c < columns; ++c)
    for (unsigned int r = 0; r < rows; ++r)
      Assert(row_sizes[r] == this->sub_objects[r][c]->m(),
             ExcDimensionMismatch(row_sizes[r], this->sub_objects[r][c]->m()));

  // finally initialize the row
  // indices with this array
  this->row_block_indices.reinit(row_sizes);


  // then do the same with the columns
  for (unsigned int c = 0; c < columns; ++c)
    col_sizes[c] = this->sub_objects[0][c]->n();
  for (unsigned int r = 1; r < rows; ++r)
    for (unsigned int c = 0; c < columns; ++c)
      Assert(col_sizes[c] == this->sub_objects[r][c]->n(),
             ExcDimensionMismatch(col_sizes[c], this->sub_objects[r][c]->n()));

  // finally initialize the row
  // indices with this array
  this->column_block_indices.reinit(col_sizes);
}



DEAL_II_NAMESPACE_CLOSE

#endif // ifdef block_sparse_matrix_templates_h
