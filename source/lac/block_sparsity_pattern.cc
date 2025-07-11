// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2000 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/memory_consumption.h>

#include <deal.II/lac/block_sparsity_pattern.h>

DEAL_II_NAMESPACE_OPEN


template <typename SparsityPatternType>
BlockSparsityPatternBase<SparsityPatternType>::BlockSparsityPatternBase()
  : block_rows(0)
  , block_columns(0)
{}



template <typename SparsityPatternType>
BlockSparsityPatternBase<SparsityPatternType>::BlockSparsityPatternBase(
  const size_type block_rows,
  const size_type block_columns)
  : BlockSparsityPatternBase()
{
  reinit(block_rows, block_columns);
}



template <typename SparsityPatternType>
BlockSparsityPatternBase<SparsityPatternType>::BlockSparsityPatternBase(
  const BlockSparsityPatternBase &s)
  : BlockSparsityPatternBase()
{
  Assert(s.n_block_rows() == 0 && s.n_block_cols() == 0,
         ExcMessage(
           "This constructor can only be called if the provided argument "
           "is the sparsity pattern for an empty matrix. This constructor can "
           "not be used to copy-construct a non-empty sparsity pattern."));
}



template <typename SparsityPatternType>
void
BlockSparsityPatternBase<SparsityPatternType>::reinit(
  const size_type new_block_rows,
  const size_type new_block_columns)
{
  sub_objects.reinit(0, 0);

  block_rows    = new_block_rows;
  block_columns = new_block_columns;

  sub_objects.reinit(block_rows, block_columns);
  for (size_type i = 0; i < n_block_rows(); ++i)
    for (size_type j = 0; j < n_block_cols(); ++j)
      sub_objects[i][j] = std::make_unique<SparsityPatternType>();
}


template <typename SparsityPatternType>
BlockSparsityPatternBase<SparsityPatternType> &
BlockSparsityPatternBase<SparsityPatternType>::operator=(
  const BlockSparsityPatternBase<SparsityPatternType> &bsp)
{
  AssertDimension(n_block_rows(), bsp.n_block_rows());
  AssertDimension(n_block_cols(), bsp.n_block_cols());
  // copy objects
  for (size_type i = 0; i < n_block_rows(); ++i)
    for (size_type j = 0; j < n_block_cols(); ++j)
      *sub_objects[i][j] = *bsp.sub_objects[i][j];
  // update index objects
  collect_sizes();

  return *this;
}



template <typename SparsityPatternType>
typename BlockSparsityPatternBase<SparsityPatternType>::size_type
BlockSparsityPatternBase<SparsityPatternType>::compute_n_rows() const
{
  // only count in first column, since
  // all rows should be equivalent
  size_type count = 0;
  for (size_type r = 0; r < n_block_rows(); ++r)
    count += sub_objects[r][0]->n_rows();
  return count;
}



template <typename SparsityPatternType>
typename BlockSparsityPatternBase<SparsityPatternType>::size_type
BlockSparsityPatternBase<SparsityPatternType>::compute_n_cols() const
{
  // only count in first row, since
  // all rows should be equivalent
  size_type count = 0;
  for (size_type c = 0; c < n_block_cols(); ++c)
    count += sub_objects[0][c]->n_cols();
  return count;
}



template <typename SparsityPatternType>
void
BlockSparsityPatternBase<SparsityPatternType>::collect_sizes()
{
  SparsityPatternBase::resize(compute_n_rows(), compute_n_cols());

  std::vector<size_type> row_sizes(n_block_rows());
  std::vector<size_type> col_sizes(n_block_cols());

  // first find out the row sizes
  // from the first block column
  for (size_type r = 0; r < n_block_rows(); ++r)
    row_sizes[r] = sub_objects[r][0]->n_rows();
  // then check that the following
  // block columns have the same
  // sizes
  for (size_type c = 1; c < n_block_cols(); ++c)
    for (size_type r = 0; r < n_block_rows(); ++r)
      Assert(row_sizes[r] == sub_objects[r][c]->n_rows(),
             ExcIncompatibleRowNumbers(r, 0, r, c));

  // finally initialize the row
  // indices with this array
  row_indices.reinit(row_sizes);


  // then do the same with the columns
  for (size_type c = 0; c < n_block_cols(); ++c)
    col_sizes[c] = sub_objects[0][c]->n_cols();
  for (size_type r = 1; r < n_block_rows(); ++r)
    for (size_type c = 0; c < n_block_cols(); ++c)
      Assert(col_sizes[c] == sub_objects[r][c]->n_cols(),
             ExcIncompatibleRowNumbers(0, c, r, c));

  // finally initialize the row
  // indices with this array
  column_indices.reinit(col_sizes);

  // Resize scratch arrays
  block_column_indices.resize(n_block_cols());
  counter_within_block.resize(n_block_cols());
}



template <typename SparsityPatternType>
void
BlockSparsityPatternBase<SparsityPatternType>::compress()
{
  for (size_type i = 0; i < n_block_rows(); ++i)
    for (size_type j = 0; j < n_block_cols(); ++j)
      sub_objects[i][j]->compress();
}



template <typename SparsityPatternType>
bool
BlockSparsityPatternBase<SparsityPatternType>::empty() const
{
  for (size_type i = 0; i < n_block_rows(); ++i)
    for (size_type j = 0; j < n_block_cols(); ++j)
      if (sub_objects[i][j]->empty() == false)
        return false;
  return true;
}



template <typename SparsityPatternType>
typename BlockSparsityPatternBase<SparsityPatternType>::size_type
BlockSparsityPatternBase<SparsityPatternType>::max_entries_per_row() const
{
  size_type max_entries = 0;
  for (size_type block_row = 0; block_row < n_block_rows(); ++block_row)
    {
      size_type this_row = 0;
      for (size_type c = 0; c < n_block_cols(); ++c)
        this_row += sub_objects[block_row][c]->max_entries_per_row();

      if (this_row > max_entries)
        max_entries = this_row;
    }
  return max_entries;
}



template <typename SparsityPatternType>
typename BlockSparsityPatternBase<SparsityPatternType>::size_type
BlockSparsityPatternBase<SparsityPatternType>::n_nonzero_elements() const
{
  size_type count = 0;
  for (size_type i = 0; i < n_block_rows(); ++i)
    for (size_type j = 0; j < n_block_cols(); ++j)
      count += sub_objects[i][j]->n_nonzero_elements();
  return count;
}



template <typename SparsityPatternType>
void
BlockSparsityPatternBase<SparsityPatternType>::print(std::ostream &out) const
{
  size_type k = 0;
  for (size_type ib = 0; ib < n_block_rows(); ++ib)
    {
      for (size_type i = 0; i < block(ib, 0).n_rows(); ++i)
        {
          out << '[' << i + k;
          size_type l = 0;
          for (size_type jb = 0; jb < n_block_cols(); ++jb)
            {
              const SparsityPatternType &b = block(ib, jb);
              for (size_type j = 0; j < b.n_cols(); ++j)
                if (b.exists(i, j))
                  out << ',' << l + j;
              l += b.n_cols();
            }
          out << ']' << std::endl;
        }
      k += block(ib, 0).n_rows();
    }
}


#ifndef DOXYGEN
template <>
void
BlockSparsityPatternBase<DynamicSparsityPattern>::print(std::ostream &out) const
{
  size_type k = 0;
  for (size_type ib = 0; ib < n_block_rows(); ++ib)
    {
      for (size_type i = 0; i < block(ib, 0).n_rows(); ++i)
        {
          out << '[' << i + k;
          size_type l = 0;
          for (size_type jb = 0; jb < n_block_cols(); ++jb)
            {
              const DynamicSparsityPattern &b = block(ib, jb);
              if (b.row_index_set().size() == 0 ||
                  b.row_index_set().is_element(i))
                for (size_type j = 0; j < b.n_cols(); ++j)
                  if (b.exists(i, j))
                    out << ',' << l + j;
              l += b.n_cols();
            }
          out << ']' << std::endl;
        }
      k += block(ib, 0).n_rows();
    }
}
#endif



template <typename SparsityPatternType>
void
BlockSparsityPatternBase<SparsityPatternType>::print_gnuplot(
  std::ostream &out) const
{
  size_type k = 0;
  for (size_type ib = 0; ib < n_block_rows(); ++ib)
    {
      for (size_type i = 0; i < block(ib, 0).n_rows(); ++i)
        {
          size_type l = 0;
          for (size_type jb = 0; jb < n_block_cols(); ++jb)
            {
              const SparsityPatternType &b = block(ib, jb);
              for (size_type j = 0; j < b.n_cols(); ++j)
                if (b.exists(i, j))
                  out << l + j << " " << -static_cast<signed int>(i + k)
                      << std::endl;
              l += b.n_cols();
            }
        }
      k += block(ib, 0).n_rows();
    }
}



template <typename SparsityPatternType>
void
BlockSparsityPatternBase<SparsityPatternType>::print_svg(
  std::ostream &out) const
{
  const unsigned int m = this->n_rows();
  const unsigned int n = this->n_cols();
  out
    << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" viewBox=\"0 0 "
    << n + 2 << " " << m + 2
    << " \">\n"
       "<style type=\"text/css\" >\n"
       "     <![CDATA[\n"
       "      rect.pixel {\n"
       "          fill:   #ff0000;\n"
       "      }\n"
       "    ]]>\n"
       "  </style>\n\n"
       "   <rect width=\""
    << n + 2 << "\" height=\"" << m + 2
    << "\" fill=\"rgb(128, 128, 128)\"/>\n"
       "   <rect x=\"1\" y=\"1\" width=\""
    << n + 0.1 << "\" height=\"" << m + 0.1
    << "\" fill=\"rgb(255, 255, 255)\"/>\n\n";

  for (unsigned int block_i = 0; block_i < n_block_rows(); ++block_i)
    for (unsigned int block_j = 0; block_j < n_block_cols(); ++block_j)
      for (const auto &entry : block(block_i, block_j))
        {
          out << "  <rect class=\"pixel\" x=\""
              << column_indices.local_to_global(block_j, entry.column()) + 1
              << "\" y=\""
              << row_indices.local_to_global(block_i, entry.row()) + 1
              << "\" width=\".9\" height=\".9\"/>\n";
        }

  out << "</svg>" << std::endl;
}



template <typename SparsityPatternType>
std::size_t
BlockSparsityPatternBase<SparsityPatternType>::memory_consumption() const
{
  std::size_t mem = 0;
  mem += (MemoryConsumption::memory_consumption(n_block_rows()) +
          MemoryConsumption::memory_consumption(n_block_cols()) +
          MemoryConsumption::memory_consumption(sub_objects) +
          MemoryConsumption::memory_consumption(row_indices) +
          MemoryConsumption::memory_consumption(column_indices));
  for (size_type r = 0; r < n_block_rows(); ++r)
    for (size_type c = 0; c < n_block_cols(); ++c)
      mem += MemoryConsumption::memory_consumption(*sub_objects[r][c]);

  return mem;
}



BlockSparsityPattern::BlockSparsityPattern(const size_type n_rows,
                                           const size_type n_columns)
  : BlockSparsityPatternBase<SparsityPattern>(n_rows, n_columns)
{}


void
BlockSparsityPattern::reinit(
  const BlockIndices                           &rows,
  const BlockIndices                           &cols,
  const std::vector<std::vector<unsigned int>> &row_lengths)
{
  AssertDimension(row_lengths.size(), cols.size());

  this->reinit(rows.size(), cols.size());
  for (size_type j = 0; j < cols.size(); ++j)
    for (size_type i = 0; i < rows.size(); ++i)
      {
        const size_type start  = rows.local_to_global(i, 0);
        const size_type length = rows.block_size(i);

        if (row_lengths[j].size() == 1)
          block(i, j).reinit(rows.block_size(i),
                             cols.block_size(j),
                             row_lengths[j][0]);
        else
          {
            Assert(row_lengths[j].begin() + start + length <=
                     row_lengths[j].end(),
                   ExcInternalError());
            ArrayView<const unsigned int> block_rows(row_lengths[j].data() +
                                                       start,
                                                     length);
            block(i, j).reinit(rows.block_size(i),
                               cols.block_size(j),
                               block_rows);
          }
      }
  this->collect_sizes();
  Assert(this->row_indices == rows, ExcInternalError());
  Assert(this->column_indices == cols, ExcInternalError());
}


bool
BlockSparsityPattern::is_compressed() const
{
  for (size_type i = 0; i < n_block_rows(); ++i)
    for (size_type j = 0; j < n_block_cols(); ++j)
      if (sub_objects[i][j]->is_compressed() == false)
        return false;
  return true;
}



void
BlockSparsityPattern::copy_from(const BlockDynamicSparsityPattern &dsp)
{
  // delete old content, set block
  // sizes anew
  reinit(dsp.n_block_rows(), dsp.n_block_cols());

  // copy over blocks
  for (size_type i = 0; i < n_block_rows(); ++i)
    for (size_type j = 0; j < n_block_cols(); ++j)
      block(i, j).copy_from(dsp.block(i, j));

  // and finally enquire their new
  // sizes
  collect_sizes();
}



BlockDynamicSparsityPattern::BlockDynamicSparsityPattern(
  const size_type n_rows,
  const size_type n_columns)
  : BlockSparsityPatternBase<DynamicSparsityPattern>(n_rows, n_columns)
{}



BlockDynamicSparsityPattern::BlockDynamicSparsityPattern(
  const std::vector<size_type> &row_indices,
  const std::vector<size_type> &col_indices)
  : BlockSparsityPatternBase<DynamicSparsityPattern>(row_indices.size(),
                                                     col_indices.size())
{
  for (size_type i = 0; i < row_indices.size(); ++i)
    for (size_type j = 0; j < col_indices.size(); ++j)
      this->block(i, j).reinit(row_indices[i], col_indices[j]);
  this->collect_sizes();
}



BlockDynamicSparsityPattern::BlockDynamicSparsityPattern(
  const std::vector<IndexSet> &partitioning)
  : BlockSparsityPatternBase<DynamicSparsityPattern>(partitioning.size(),
                                                     partitioning.size())
{
  for (size_type i = 0; i < partitioning.size(); ++i)
    for (size_type j = 0; j < partitioning.size(); ++j)
      this->block(i, j).reinit(partitioning[i].size(),
                               partitioning[j].size(),
                               partitioning[i]);
  this->collect_sizes();
}



BlockDynamicSparsityPattern::BlockDynamicSparsityPattern(
  const BlockIndices &row_indices,
  const BlockIndices &col_indices)
{
  reinit(row_indices, col_indices);
}



void
BlockDynamicSparsityPattern::reinit(
  const std::vector<size_type> &row_block_sizes,
  const std::vector<size_type> &col_block_sizes)
{
  BlockSparsityPatternBase<DynamicSparsityPattern>::reinit(
    row_block_sizes.size(), col_block_sizes.size());
  for (size_type i = 0; i < row_block_sizes.size(); ++i)
    for (size_type j = 0; j < col_block_sizes.size(); ++j)
      this->block(i, j).reinit(row_block_sizes[i], col_block_sizes[j]);
  this->collect_sizes();
}



void
BlockDynamicSparsityPattern::reinit(const std::vector<IndexSet> &partitioning)
{
  BlockSparsityPatternBase<DynamicSparsityPattern>::reinit(partitioning.size(),
                                                           partitioning.size());
  for (size_type i = 0; i < partitioning.size(); ++i)
    for (size_type j = 0; j < partitioning.size(); ++j)
      this->block(i, j).reinit(partitioning[i].size(),
                               partitioning[j].size(),
                               partitioning[i]);
  this->collect_sizes();
}



void
BlockDynamicSparsityPattern::reinit(const BlockIndices &row_indices,
                                    const BlockIndices &col_indices)
{
  BlockSparsityPatternBase<DynamicSparsityPattern>::reinit(row_indices.size(),
                                                           col_indices.size());
  for (size_type i = 0; i < row_indices.size(); ++i)
    for (size_type j = 0; j < col_indices.size(); ++j)
      this->block(i, j).reinit(row_indices.block_size(i),
                               col_indices.block_size(j));
  this->collect_sizes();
}


#ifdef DEAL_II_WITH_TRILINOS
namespace TrilinosWrappers
{
  BlockSparsityPattern::BlockSparsityPattern(const size_type n_rows,
                                             const size_type n_columns)
    : dealii::BlockSparsityPatternBase<SparsityPattern>(n_rows, n_columns)
  {}



  BlockSparsityPattern::BlockSparsityPattern(
    const std::vector<size_type> &row_indices,
    const std::vector<size_type> &col_indices)
    : BlockSparsityPatternBase<SparsityPattern>(row_indices.size(),
                                                col_indices.size())
  {
    for (size_type i = 0; i < row_indices.size(); ++i)
      for (size_type j = 0; j < col_indices.size(); ++j)
        this->block(i, j).reinit(row_indices[i], col_indices[j]);
    this->collect_sizes();
  }



  BlockSparsityPattern::BlockSparsityPattern(
    const std::vector<IndexSet> &parallel_partitioning,
    const MPI_Comm               communicator)
    : BlockSparsityPatternBase<SparsityPattern>(parallel_partitioning.size(),
                                                parallel_partitioning.size())
  {
    for (size_type i = 0; i < parallel_partitioning.size(); ++i)
      for (size_type j = 0; j < parallel_partitioning.size(); ++j)
        this->block(i, j).reinit(parallel_partitioning[i],
                                 parallel_partitioning[j],
                                 communicator);
    this->collect_sizes();
  }



  BlockSparsityPattern::BlockSparsityPattern(
    const std::vector<IndexSet> &row_parallel_partitioning,
    const std::vector<IndexSet> &col_parallel_partitioning,
    const std::vector<IndexSet> &writable_rows,
    const MPI_Comm               communicator)
    : BlockSparsityPatternBase<SparsityPattern>(
        row_parallel_partitioning.size(),
        col_parallel_partitioning.size())
  {
    for (size_type i = 0; i < row_parallel_partitioning.size(); ++i)
      for (size_type j = 0; j < col_parallel_partitioning.size(); ++j)
        this->block(i, j).reinit(row_parallel_partitioning[i],
                                 col_parallel_partitioning[j],
                                 writable_rows[i],
                                 communicator);
    this->collect_sizes();
  }



  void
  BlockSparsityPattern::reinit(const std::vector<size_type> &row_block_sizes,
                               const std::vector<size_type> &col_block_sizes)
  {
    dealii::BlockSparsityPatternBase<SparsityPattern>::reinit(
      row_block_sizes.size(), col_block_sizes.size());
    for (size_type i = 0; i < row_block_sizes.size(); ++i)
      for (size_type j = 0; j < col_block_sizes.size(); ++j)
        this->block(i, j).reinit(row_block_sizes[i], col_block_sizes[j]);
    this->collect_sizes();
  }



  void
  BlockSparsityPattern::reinit(
    const std::vector<IndexSet> &parallel_partitioning,
    const MPI_Comm               communicator)
  {
    dealii::BlockSparsityPatternBase<SparsityPattern>::reinit(
      parallel_partitioning.size(), parallel_partitioning.size());
    for (size_type i = 0; i < parallel_partitioning.size(); ++i)
      for (size_type j = 0; j < parallel_partitioning.size(); ++j)
        this->block(i, j).reinit(parallel_partitioning[i],
                                 parallel_partitioning[j],
                                 communicator);
    this->collect_sizes();
  }



  void
  BlockSparsityPattern::reinit(
    const std::vector<IndexSet> &row_parallel_partitioning,
    const std::vector<IndexSet> &col_parallel_partitioning,
    const MPI_Comm               communicator)
  {
    dealii::BlockSparsityPatternBase<SparsityPattern>::reinit(
      row_parallel_partitioning.size(), col_parallel_partitioning.size());
    for (size_type i = 0; i < row_parallel_partitioning.size(); ++i)
      for (size_type j = 0; j < col_parallel_partitioning.size(); ++j)
        this->block(i, j).reinit(row_parallel_partitioning[i],
                                 col_parallel_partitioning[j],
                                 communicator);
    this->collect_sizes();
  }



  void
  BlockSparsityPattern::reinit(
    const std::vector<IndexSet> &row_parallel_partitioning,
    const std::vector<IndexSet> &col_parallel_partitioning,
    const std::vector<IndexSet> &writable_rows,
    const MPI_Comm               communicator)
  {
    AssertDimension(writable_rows.size(), row_parallel_partitioning.size());
    dealii::BlockSparsityPatternBase<SparsityPattern>::reinit(
      row_parallel_partitioning.size(), col_parallel_partitioning.size());
    for (size_type i = 0; i < row_parallel_partitioning.size(); ++i)
      for (size_type j = 0; j < col_parallel_partitioning.size(); ++j)
        this->block(i, j).reinit(row_parallel_partitioning[i],
                                 col_parallel_partitioning[j],
                                 writable_rows[i],
                                 communicator);
    this->collect_sizes();
  }

} // namespace TrilinosWrappers

#endif

template class BlockSparsityPatternBase<SparsityPattern>;
template class BlockSparsityPatternBase<DynamicSparsityPattern>;
#ifdef DEAL_II_WITH_TRILINOS
template class BlockSparsityPatternBase<TrilinosWrappers::SparsityPattern>;
#endif

DEAL_II_NAMESPACE_CLOSE
