// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2000 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_block_sparse_matrix_templates_h
#define dealii_block_sparse_matrix_templates_h


#include <deal.II/base/config.h>

#include <deal.II/base/memory_consumption.h>

#include <deal.II/lac/block_sparse_matrix.h>

DEAL_II_NAMESPACE_OPEN

template <typename number>
BlockSparseMatrix<number>::BlockSparseMatrix(
  const BlockSparsityPattern &sparsity)
{
  // virtual functions called in constructors and destructors never use the
  // override in a derived class
  // for clarity be explicit on which function is called
  BlockSparseMatrix<number>::reinit(sparsity);
}



template <typename number>
BlockSparseMatrix<number>::~BlockSparseMatrix()
{
  // delete previous content of
  // the subobjects array
  try
    {
      clear();
    }
  catch (...)
    {}
  sparsity_pattern = nullptr;
}



template <typename number>
BlockSparseMatrix<number> &
BlockSparseMatrix<number>::operator=(const BlockSparseMatrix<number> &m)
{
  Assert(this->row_block_indices == m.row_block_indices,
         ExcBlockDimensionMismatch());
  Assert(this->column_block_indices == m.column_block_indices,
         ExcBlockDimensionMismatch());

  // this operator does not do
  // anything except than checking
  // whether the base objects want to
  // do something
  for (size_type r = 0; r < this->n_block_rows(); ++r)
    for (size_type c = 0; c < this->n_block_cols(); ++c)
      this->block(r, c) = m.block(r, c);

  return *this;
}



template <typename number>
void
BlockSparseMatrix<number>::clear()
{
  BlockMatrixBase<SparseMatrix<number>>::clear();
  sparsity_pattern = nullptr;
}



template <typename number>
void
BlockSparseMatrix<number>::reinit(const BlockSparsityPattern &sparsity)
{
  // first delete previous content of
  // the subobjects array and delete
  // the table completely
  clear();

  // then associate new sparsity
  // pattern and resize
  sparsity_pattern = &sparsity;

  this->row_block_indices    = sparsity.row_indices;
  this->column_block_indices = sparsity.column_indices;

  this->sub_objects.reinit(sparsity.n_block_rows(), sparsity.n_block_cols());

  // and reinitialize the blocks
  for (size_type r = 0; r < this->n_block_rows(); ++r)
    for (size_type c = 0; c < this->n_block_cols(); ++c)
      {
        BlockType *p = new SparseMatrix<number>();
        p->reinit(sparsity.block(r, c));
        this->sub_objects[r][c] = p;
      }
}



template <typename number>
bool
BlockSparseMatrix<number>::empty() const
{
  for (size_type r = 0; r < this->n_block_rows(); ++r)
    for (size_type c = 0; c < this->n_block_cols(); ++c)
      if (this->block(r, c).empty() == false)
        return false;

  return true;
}



template <typename number>
typename BlockSparseMatrix<number>::size_type
BlockSparseMatrix<number>::get_row_length(const size_type row) const
{
  return sparsity_pattern->row_length(row);
}



template <typename number>
typename BlockSparseMatrix<number>::size_type
BlockSparseMatrix<number>::n_nonzero_elements() const
{
  return sparsity_pattern->n_nonzero_elements();
}



template <typename number>
typename BlockSparseMatrix<number>::size_type
BlockSparseMatrix<number>::n_actually_nonzero_elements(
  const double threshold) const
{
  size_type count = 0;
  for (size_type i = 0; i < this->n_block_rows(); ++i)
    for (size_type j = 0; j < this->n_block_cols(); ++j)
      count += this->sub_objects[i][j]->n_actually_nonzero_elements(threshold);

  return count;
}



template <typename number>
const BlockSparsityPattern &
BlockSparseMatrix<number>::get_sparsity_pattern() const
{
  return *sparsity_pattern;
}



template <typename number>
void
BlockSparseMatrix<number>::print_formatted(std::ostream      &out,
                                           const unsigned int precision,
                                           const bool         scientific,
                                           const unsigned int width,
                                           const char        *zero_string,
                                           const double       denominator,
                                           const char        *separator) const
{
  for (size_type r = 0; r < this->n_block_rows(); ++r)
    for (size_type c = 0; c < this->n_block_cols(); ++c)
      {
        out << "Component (" << r << ',' << c << ')' << std::endl;
        this->block(r, c).print_formatted(out,
                                          precision,
                                          scientific,
                                          width,
                                          zero_string,
                                          denominator,
                                          separator);
      }
}



template <typename number>
std::size_t
BlockSparseMatrix<number>::memory_consumption() const
{
  std::size_t mem = sizeof(*this);
  mem += MemoryConsumption::memory_consumption(this->sub_objects);
  for (size_type r = 0; r < this->n_block_rows(); ++r)
    for (size_type c = 0; c < this->n_block_cols(); ++c)
      mem += MemoryConsumption::memory_consumption(*this->sub_objects[r][c]);

  return mem;
}



DEAL_II_NAMESPACE_CLOSE

#endif // ifdef block_sparse_matrix_templates_h
