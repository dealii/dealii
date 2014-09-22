// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2013 by the deal.II authors
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


#include <deal.II/base/vector_slice.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/lac/block_sparsity_pattern.h>

DEAL_II_NAMESPACE_OPEN


template <class SparsityPatternBase>
BlockSparsityPatternBase<SparsityPatternBase>::BlockSparsityPatternBase ()
  :
  rows (0),
  columns (0)
{}



template <class SparsityPatternBase>
BlockSparsityPatternBase<SparsityPatternBase>::
BlockSparsityPatternBase (const size_type n_block_rows,
                          const size_type n_block_columns)
  :
  rows (0),
  columns (0)
{
  reinit (n_block_rows, n_block_columns);
}



template <class SparsityPatternBase>
BlockSparsityPatternBase<SparsityPatternBase>::
BlockSparsityPatternBase (const BlockSparsityPatternBase &s)
  :
  Subscriptor ()
{
  Assert(s.rows==0, ExcInvalidConstructorCall());
  Assert(s.columns==0, ExcInvalidConstructorCall());

  rows = 0;
  columns=0;
}



template <class SparsityPatternBase>
BlockSparsityPatternBase<SparsityPatternBase>::~BlockSparsityPatternBase ()
{
  // clear all memory
  reinit (0,0);
}



template <class SparsityPatternBase>
void
BlockSparsityPatternBase<SparsityPatternBase>::
reinit (const size_type n_block_rows,
        const size_type n_block_columns)
{
  // delete previous content and
  // clean the sub_objects array
  // completely
  for (size_type i=0; i<rows; ++i)
    for (size_type j=0; j<columns; ++j)
      {
        SparsityPatternBase *sp = sub_objects[i][j];
        sub_objects[i][j] = 0;
        delete sp;
      };
  sub_objects.reinit (0,0);

  // then set new sizes
  rows = n_block_rows;
  columns = n_block_columns;
  sub_objects.reinit (rows, columns);

  // allocate new objects
  for (size_type i=0; i<rows; ++i)
    for (size_type j=0; j<columns; ++j)
      {
        SparsityPatternBase *p = new SparsityPatternBase;
        sub_objects[i][j] = p;
      }
}


template <class SparsityPatternBase>
BlockSparsityPatternBase<SparsityPatternBase> &
BlockSparsityPatternBase<SparsityPatternBase>::
operator = (const BlockSparsityPatternBase<SparsityPatternBase> &bsp)
{
  Assert (rows == bsp.rows, ExcDimensionMismatch(rows, bsp.rows));
  Assert (columns == bsp.columns, ExcDimensionMismatch(columns, bsp.columns));
  // copy objects
  for (size_type i=0; i<rows; ++i)
    for (size_type j=0; j<columns; ++j)
      *sub_objects[i][j] = *bsp.sub_objects[i][j];
  // update index objects
  collect_sizes ();

  return *this;
}



template <class SparsityPatternBase>
void
BlockSparsityPatternBase<SparsityPatternBase>::collect_sizes ()
{
  std::vector<size_type > row_sizes (rows);
  std::vector<size_type > col_sizes (columns);

  // first find out the row sizes
  // from the first block column
  for (size_type r=0; r<rows; ++r)
    row_sizes[r] = sub_objects[r][0]->n_rows();
  // then check that the following
  // block columns have the same
  // sizes
  for (size_type c=1; c<columns; ++c)
    for (size_type r=0; r<rows; ++r)
      Assert (row_sizes[r] == sub_objects[r][c]->n_rows(),
              ExcIncompatibleRowNumbers (r,0,r,c));

  // finally initialize the row
  // indices with this array
  row_indices.reinit (row_sizes);


  // then do the same with the columns
  for (size_type c=0; c<columns; ++c)
    col_sizes[c] = sub_objects[0][c]->n_cols();
  for (size_type r=1; r<rows; ++r)
    for (size_type c=0; c<columns; ++c)
      Assert (col_sizes[c] == sub_objects[r][c]->n_cols(),
              ExcIncompatibleRowNumbers (0,c,r,c));

  // finally initialize the row
  // indices with this array
  column_indices.reinit (col_sizes);
}



template <class SparsityPatternBase>
void
BlockSparsityPatternBase<SparsityPatternBase>::compress ()
{
  for (size_type i=0; i<rows; ++i)
    for (size_type j=0; j<columns; ++j)
      sub_objects[i][j]->compress ();
}



template <class SparsityPatternBase>
bool
BlockSparsityPatternBase<SparsityPatternBase>::empty () const
{
  for (size_type i=0; i<rows; ++i)
    for (size_type j=0; j<columns; ++j)
      if (sub_objects[i][j]->empty () == false)
        return false;
  return true;
}



template <class SparsityPatternBase>
typename BlockSparsityPatternBase<SparsityPatternBase>::size_type
BlockSparsityPatternBase<SparsityPatternBase>::max_entries_per_row () const
{
  size_type max_entries = 0;
  for (size_type block_row=0; block_row<rows; ++block_row)
    {
      size_type this_row = 0;
      for (size_type c=0; c<columns; ++c)
        this_row += sub_objects[block_row][c]->max_entries_per_row ();

      if (this_row > max_entries)
        max_entries = this_row;
    };
  return max_entries;
}



template <class SparsityPatternBase>
typename BlockSparsityPatternBase<SparsityPatternBase>::size_type
BlockSparsityPatternBase<SparsityPatternBase>::n_rows () const
{
  // only count in first column, since
  // all rows should be equivalent
  size_type count = 0;
  for (size_type r=0; r<rows; ++r)
    count += sub_objects[r][0]->n_rows();
  return count;
}



template <class SparsityPatternBase>
typename BlockSparsityPatternBase<SparsityPatternBase>::size_type
BlockSparsityPatternBase<SparsityPatternBase>::n_cols () const
{
  // only count in first row, since
  // all rows should be equivalent
  size_type count = 0;
  for (size_type c=0; c<columns; ++c)
    count += sub_objects[0][c]->n_cols();
  return count;
}



template <class SparsityPatternBase>
typename BlockSparsityPatternBase<SparsityPatternBase>::size_type
BlockSparsityPatternBase<SparsityPatternBase>::n_nonzero_elements () const
{
  size_type count = 0;
  for (size_type i=0; i<rows; ++i)
    for (size_type j=0; j<columns; ++j)
      count += sub_objects[i][j]->n_nonzero_elements ();
  return count;
}



template <class SparsityPatternBase>
void
BlockSparsityPatternBase<SparsityPatternBase>::print(std::ostream &out) const
{
  size_type k=0;
  for (size_type ib=0; ib<n_block_rows(); ++ib)
    {
      for (size_type i=0; i<block(ib,0).n_rows(); ++i)
        {
          out << '[' << i+k;
          size_type l=0;
          for (size_type jb=0; jb<n_block_cols(); ++jb)
            {
              const SparsityPatternBase &b = block(ib,jb);
              for (size_type j=0; j<b.n_cols(); ++j)
                if (b.exists(i,j))
                  out << ',' << l+j;
              l += b.n_cols();
            }
          out << ']' << std::endl;
        }
      k += block(ib,0).n_rows();
    }
}


template <>
void
BlockSparsityPatternBase<CompressedSimpleSparsityPattern>::print(std::ostream &out) const
{
  size_type k=0;
  for (size_type ib=0; ib<n_block_rows(); ++ib)
    {
      for (size_type i=0; i<block(ib,0).n_rows(); ++i)
        {
          out << '[' << i+k;
          size_type l=0;
          for (size_type jb=0; jb<n_block_cols(); ++jb)
            {
              const CompressedSimpleSparsityPattern &b = block(ib,jb);
              if (b.row_index_set().size()==0 || b.row_index_set().is_element(i))
                for (size_type j=0; j<b.n_cols(); ++j)
                  if (b.exists(i,j))
                    out << ',' << l+j;
              l += b.n_cols();
            }
          out << ']' << std::endl;
        }
      k += block(ib,0).n_rows();
    }
}


template <class SparsityPatternBase>
void
BlockSparsityPatternBase<SparsityPatternBase>::print_gnuplot(std::ostream &out) const
{
  size_type k=0;
  for (size_type ib=0; ib<n_block_rows(); ++ib)
    {
      for (size_type i=0; i<block(ib,0).n_rows(); ++i)
        {
          size_type l=0;
          for (size_type jb=0; jb<n_block_cols(); ++jb)
            {
              const SparsityPatternBase &b = block(ib,jb);
              for (size_type j=0; j<b.n_cols(); ++j)
                if (b.exists(i,j))
                  out << l+j << " " << -static_cast<signed int>(i+k) << std::endl;
              l += b.n_cols();
            }
        }
      k += block(ib,0).n_rows();
    }
}




// Remark: The following explicit instantiations needed to be moved to
// this place here to work around a problem with gcc3.3 on Apple MacOSX.
// The reason is that some of the functions instantiated here are used
// further down; if they are not explicitly instantiated here, then the
// compiler will do an implicit instantiation and give it internal linkage
// (despite the later explicit instantiation that should make sure it
// gets external linkage). To make sure the functions have external
// linkage, we need to place the explicit instantiation before the first
// use.
//
// For more information, see http://gcc.gnu.org/bugzilla/show_bug.cgi?id=24331
// +++++++++++++


template class BlockSparsityPatternBase<SparsityPattern>;
template class BlockSparsityPatternBase<CompressedSparsityPattern>;
template class BlockSparsityPatternBase<CompressedSimpleSparsityPattern>;
template class BlockSparsityPatternBase<CompressedSetSparsityPattern>;
#ifdef DEAL_II_WITH_TRILINOS
template class BlockSparsityPatternBase<TrilinosWrappers::SparsityPattern>;
#endif



BlockSparsityPattern::BlockSparsityPattern ()
{}



BlockSparsityPattern::BlockSparsityPattern (const size_type n_rows,
                                            const size_type n_columns)
  :
  BlockSparsityPatternBase<SparsityPattern>(n_rows,
                                            n_columns)
{}


void
BlockSparsityPattern::reinit(
  const BlockIndices &rows,
  const BlockIndices &cols,
  const std::vector<std::vector<unsigned int> > &row_lengths)
{
  AssertDimension (row_lengths.size(), cols.size());

  this->reinit(rows.size(), cols.size());
  for (size_type j=0; j<cols.size(); ++j)
    for (size_type i=0; i<rows.size(); ++i)
      {
        const size_type start = rows.local_to_global(i, 0);
        const size_type length = rows.block_size(i);

        if (row_lengths[j].size()==1)
          block(i,j).reinit(rows.block_size(i),
                            cols.block_size(j), row_lengths[j][0]);
        else
          {
            VectorSlice<const std::vector<unsigned int> >
            block_rows(row_lengths[j], start, length);
            block(i,j).reinit(rows.block_size(i),
                              cols.block_size(j),
                              block_rows);
          }
      }
  this->collect_sizes();
  Assert (this->row_indices == rows, ExcInternalError());
  Assert (this->column_indices == cols, ExcInternalError());
}


bool
BlockSparsityPattern::is_compressed () const
{
  for (size_type i=0; i<rows; ++i)
    for (size_type j=0; j<columns; ++j)
      if (sub_objects[i][j]->is_compressed () == false)
        return false;
  return true;
}


std::size_t
BlockSparsityPattern::memory_consumption () const
{
  std::size_t mem = 0;
  mem += (MemoryConsumption::memory_consumption (rows) +
          MemoryConsumption::memory_consumption (columns) +
          MemoryConsumption::memory_consumption (sub_objects) +
          MemoryConsumption::memory_consumption (row_indices) +
          MemoryConsumption::memory_consumption (column_indices));
  for (size_type r=0; r<rows; ++r)
    for (size_type c=0; c<columns; ++c)
      mem += MemoryConsumption::memory_consumption (*sub_objects[r][c]);

  return mem;
}



void
BlockSparsityPattern::copy_from  (const BlockCompressedSparsityPattern &csp)
{
  // delete old content, set block
  // sizes anew
  reinit (csp.n_block_rows(), csp.n_block_cols());

  // copy over blocks
  for (size_type i=0; i<n_block_rows(); ++i)
    for (size_type j=0; j<n_block_cols(); ++j)
      block(i,j).copy_from (csp.block(i,j));

  // and finally enquire their new
  // sizes
  collect_sizes();
}

void
BlockSparsityPattern::copy_from  (const BlockCompressedSimpleSparsityPattern &csp)
{
  // delete old content, set block
  // sizes anew
  reinit (csp.n_block_rows(), csp.n_block_cols());

  // copy over blocks
  for (size_type i=0; i<rows; ++i)
    for (size_type j=0; j<rows; ++j)
      block(i,j).copy_from (csp.block(i,j));

  // and finally enquire their new
  // sizes
  collect_sizes();
}

void
BlockSparsityPattern::copy_from  (const BlockCompressedSetSparsityPattern &csp)
{
  // delete old content, set block
  // sizes anew
  reinit (csp.n_block_rows(), csp.n_block_cols());

  // copy over blocks
  for (size_type i=0; i<rows; ++i)
    for (size_type j=0; j<rows; ++j)
      block(i,j).copy_from (csp.block(i,j));

  // and finally enquire their new
  // sizes
  collect_sizes();
}



BlockCompressedSparsityPattern::BlockCompressedSparsityPattern ()
{}



BlockCompressedSparsityPattern::
BlockCompressedSparsityPattern (
  const size_type n_rows,
  const size_type n_columns)
  :
  BlockSparsityPatternBase<CompressedSparsityPattern>(n_rows,
                                                      n_columns)
{}


BlockCompressedSparsityPattern::
BlockCompressedSparsityPattern (
  const std::vector<size_type> &row_indices,
  const std::vector<size_type> &col_indices)
{
  reinit(row_indices, col_indices);
}


BlockCompressedSparsityPattern::
BlockCompressedSparsityPattern (
  const BlockIndices &row_indices,
  const BlockIndices &col_indices)
{
  reinit(row_indices, col_indices);
}


void
BlockCompressedSparsityPattern::reinit (
  const std::vector<size_type> &row_block_sizes,
  const std::vector<size_type> &col_block_sizes)
{
  BlockSparsityPatternBase<CompressedSparsityPattern>::reinit(row_block_sizes.size(), col_block_sizes.size());
  for (size_type i=0; i<row_block_sizes.size(); ++i)
    for (size_type j=0; j<col_block_sizes.size(); ++j)
      this->block(i,j).reinit(row_block_sizes[i],col_block_sizes[j]);
  this->collect_sizes();
}



void
BlockCompressedSparsityPattern::reinit (
  const BlockIndices &row_indices,
  const BlockIndices &col_indices)
{
  BlockSparsityPatternBase<CompressedSparsityPattern>::reinit(row_indices.size(),
                                                              col_indices.size());
  for (size_type i=0; i<row_indices.size(); ++i)
    for (size_type j=0; j<col_indices.size(); ++j)
      this->block(i,j).reinit(row_indices.block_size(i),
                              col_indices.block_size(j));
  this->collect_sizes();
}



BlockCompressedSetSparsityPattern::BlockCompressedSetSparsityPattern ()
{}



BlockCompressedSetSparsityPattern::
BlockCompressedSetSparsityPattern (
  const size_type n_rows,
  const size_type n_columns)
  :
  BlockSparsityPatternBase<CompressedSetSparsityPattern>(n_rows,
                                                         n_columns)
{}


BlockCompressedSetSparsityPattern::
BlockCompressedSetSparsityPattern (
  const std::vector<size_type> &row_indices,
  const std::vector<size_type> &col_indices)
{
  reinit(row_indices, col_indices);
}


BlockCompressedSetSparsityPattern::
BlockCompressedSetSparsityPattern (
  const BlockIndices &row_indices,
  const BlockIndices &col_indices)
{
  reinit(row_indices, col_indices);
}


void
BlockCompressedSetSparsityPattern::reinit (
  const std::vector<size_type> &row_block_sizes,
  const std::vector<size_type> &col_block_sizes)
{
  BlockSparsityPatternBase<CompressedSetSparsityPattern>::reinit(row_block_sizes.size(), col_block_sizes.size());
  for (size_type i=0; i<row_block_sizes.size(); ++i)
    for (size_type j=0; j<col_block_sizes.size(); ++j)
      this->block(i,j).reinit(row_block_sizes[i],col_block_sizes[j]);
  this->collect_sizes();
}



void
BlockCompressedSetSparsityPattern::reinit (
  const BlockIndices &row_indices,
  const BlockIndices &col_indices)
{
  BlockSparsityPatternBase<CompressedSetSparsityPattern>::reinit(row_indices.size(),
      col_indices.size());
  for (size_type i=0; i<row_indices.size(); ++i)
    for (size_type j=0; j<col_indices.size(); ++j)
      this->block(i,j).reinit(row_indices.block_size(i),
                              col_indices.block_size(j));
  this->collect_sizes();
}



BlockCompressedSimpleSparsityPattern::BlockCompressedSimpleSparsityPattern ()
{}



BlockCompressedSimpleSparsityPattern::
BlockCompressedSimpleSparsityPattern (const size_type n_rows,
                                      const size_type n_columns)
  :
  BlockSparsityPatternBase<CompressedSimpleSparsityPattern>(n_rows,
                                                            n_columns)
{}



BlockCompressedSimpleSparsityPattern::
BlockCompressedSimpleSparsityPattern (const std::vector<size_type> &row_indices,
                                      const std::vector<size_type> &col_indices)
  :
  BlockSparsityPatternBase<CompressedSimpleSparsityPattern>(row_indices.size(),
                                                            col_indices.size())
{
  for (size_type i=0; i<row_indices.size(); ++i)
    for (size_type j=0; j<col_indices.size(); ++j)
      this->block(i,j).reinit(row_indices[i],col_indices[j]);
  this->collect_sizes();
}


BlockCompressedSimpleSparsityPattern::
BlockCompressedSimpleSparsityPattern (const std::vector<IndexSet> &partitioning)
  :
  BlockSparsityPatternBase<CompressedSimpleSparsityPattern>(partitioning.size(),
                                                            partitioning.size())
{
  for (size_type i=0; i<partitioning.size(); ++i)
    for (size_type j=0; j<partitioning.size(); ++j)
      this->block(i,j).reinit(partitioning[i].size(),
                              partitioning[j].size(),
                              partitioning[i]);
  this->collect_sizes();
}



void
BlockCompressedSimpleSparsityPattern::reinit (
  const std::vector<size_type> &row_block_sizes,
  const std::vector<size_type> &col_block_sizes)
{
  BlockSparsityPatternBase<CompressedSimpleSparsityPattern>::
  reinit(row_block_sizes.size(), col_block_sizes.size());
  for (size_type i=0; i<row_block_sizes.size(); ++i)
    for (size_type j=0; j<col_block_sizes.size(); ++j)
      this->block(i,j).reinit(row_block_sizes[i],col_block_sizes[j]);
  this->collect_sizes();
}

void
BlockCompressedSimpleSparsityPattern::reinit (
  const std::vector< IndexSet > &partitioning)
{
  BlockSparsityPatternBase<CompressedSimpleSparsityPattern>::
  reinit(partitioning.size(), partitioning.size());
  for (size_type i=0; i<partitioning.size(); ++i)
    for (size_type j=0; j<partitioning.size(); ++j)
      this->block(i,j).reinit(partitioning[i].size(),
                              partitioning[j].size(),
                              partitioning[i]);
  this->collect_sizes();
}


#ifdef DEAL_II_WITH_TRILINOS
namespace TrilinosWrappers
{

  BlockSparsityPattern::BlockSparsityPattern ()
  {}



  BlockSparsityPattern::
  BlockSparsityPattern (const size_type n_rows,
                        const size_type n_columns)
    :
    dealii::BlockSparsityPatternBase<SparsityPattern>(n_rows,
                                                      n_columns)
  {}



  BlockSparsityPattern::
  BlockSparsityPattern (const std::vector<size_type> &row_indices,
                        const std::vector<size_type> &col_indices)
    :
    BlockSparsityPatternBase<SparsityPattern>(row_indices.size(),
                                              col_indices.size())
  {
    for (size_type i=0; i<row_indices.size(); ++i)
      for (size_type j=0; j<col_indices.size(); ++j)
        this->block(i,j).reinit(row_indices[i],col_indices[j]);
    this->collect_sizes();
  }



  BlockSparsityPattern::
  BlockSparsityPattern (const std::vector<Epetra_Map> &parallel_partitioning)
    :
    BlockSparsityPatternBase<SparsityPattern>
    (parallel_partitioning.size(),
     parallel_partitioning.size())
  {
    for (size_type i=0; i<parallel_partitioning.size(); ++i)
      for (size_type j=0; j<parallel_partitioning.size(); ++j)
        this->block(i,j).reinit(parallel_partitioning[i],
                                parallel_partitioning[j]);
    this->collect_sizes();
  }



  BlockSparsityPattern::
  BlockSparsityPattern (const std::vector<IndexSet> &parallel_partitioning,
                        const MPI_Comm              &communicator)
    :
    BlockSparsityPatternBase<SparsityPattern>
    (parallel_partitioning.size(),
     parallel_partitioning.size())
  {
    for (size_type i=0; i<parallel_partitioning.size(); ++i)
      for (size_type j=0; j<parallel_partitioning.size(); ++j)
        this->block(i,j).reinit(parallel_partitioning[i],
                                parallel_partitioning[j],
                                communicator);
    this->collect_sizes();
  }



  BlockSparsityPattern::
  BlockSparsityPattern (const std::vector<IndexSet> &row_parallel_partitioning,
                        const std::vector<IndexSet> &col_parallel_partitioning,
                        const std::vector<IndexSet> &writable_rows,
                        const MPI_Comm              &communicator)
    :
    BlockSparsityPatternBase<SparsityPattern>
    (row_parallel_partitioning.size(),
     col_parallel_partitioning.size())
  {
    for (size_type i=0; i<row_parallel_partitioning.size(); ++i)
      for (size_type j=0; j<col_parallel_partitioning.size(); ++j)
        this->block(i,j).reinit(row_parallel_partitioning[i],
                                col_parallel_partitioning[j],
                                writable_rows[i],
                                communicator);
    this->collect_sizes();
  }



  void
  BlockSparsityPattern::reinit (const std::vector<size_type> &row_block_sizes,
                                const std::vector<size_type> &col_block_sizes)
  {
    dealii::BlockSparsityPatternBase<SparsityPattern>::
    reinit(row_block_sizes.size(), col_block_sizes.size());
    for (size_type i=0; i<row_block_sizes.size(); ++i)
      for (size_type j=0; j<col_block_sizes.size(); ++j)
        this->block(i,j).reinit(row_block_sizes[i],col_block_sizes[j]);
    this->collect_sizes();
  }



  void
  BlockSparsityPattern::reinit (const std::vector<Epetra_Map> &parallel_partitioning)
  {
    dealii::BlockSparsityPatternBase<SparsityPattern>::
    reinit(parallel_partitioning.size(),
           parallel_partitioning.size());
    for (size_type i=0; i<parallel_partitioning.size(); ++i)
      for (size_type j=0; j<parallel_partitioning.size(); ++j)
        this->block(i,j).reinit(parallel_partitioning[i],
                                parallel_partitioning[j]);
    this->collect_sizes();
  }



  void
  BlockSparsityPattern::reinit (const std::vector<IndexSet> &parallel_partitioning,
                                const MPI_Comm &communicator)
  {
    dealii::BlockSparsityPatternBase<SparsityPattern>::
    reinit(parallel_partitioning.size(),
           parallel_partitioning.size());
    for (size_type i=0; i<parallel_partitioning.size(); ++i)
      for (size_type j=0; j<parallel_partitioning.size(); ++j)
        this->block(i,j).reinit(parallel_partitioning[i],
                                parallel_partitioning[j],
                                communicator);
    this->collect_sizes();
  }



  void
  BlockSparsityPattern::reinit (const std::vector<IndexSet> &row_parallel_partitioning,
                                const std::vector<IndexSet> &col_parallel_partitioning,
                                const MPI_Comm &communicator)
  {
    dealii::BlockSparsityPatternBase<SparsityPattern>::
    reinit(row_parallel_partitioning.size(),
           col_parallel_partitioning.size());
    for (size_type i=0; i<row_parallel_partitioning.size(); ++i)
      for (size_type j=0; j<col_parallel_partitioning.size(); ++j)
        this->block(i,j).reinit(row_parallel_partitioning[i],
                                col_parallel_partitioning[j],
                                communicator);
    this->collect_sizes();
  }



  void
  BlockSparsityPattern::reinit (const std::vector<IndexSet> &row_parallel_partitioning,
                                const std::vector<IndexSet> &col_parallel_partitioning,
                                const std::vector<IndexSet> &writable_rows,
                                const MPI_Comm &communicator)
  {
    AssertDimension(writable_rows.size(), row_parallel_partitioning.size());
    dealii::BlockSparsityPatternBase<SparsityPattern>::
    reinit(row_parallel_partitioning.size(),
           col_parallel_partitioning.size());
    for (size_type i=0; i<row_parallel_partitioning.size(); ++i)
      for (size_type j=0; j<col_parallel_partitioning.size(); ++j)
        this->block(i,j).reinit(row_parallel_partitioning[i],
                                col_parallel_partitioning[j],
                                writable_rows[i],
                                communicator);
    this->collect_sizes();
  }

}

#endif

// Remark: The explicit instantiations for "BlockSparsityPatternBase" were moved
// to the top of this source file. The reason is a slightly buggy version
// of the Apple gcc v.3.3.
// For more information, see http://gcc.gnu.org/bugzilla/show_bug.cgi?id=24331

DEAL_II_NAMESPACE_CLOSE
