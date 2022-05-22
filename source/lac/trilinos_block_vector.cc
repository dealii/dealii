// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#include <deal.II/lac/trilinos_parallel_block_vector.h>

#ifdef DEAL_II_WITH_TRILINOS

#  include <deal.II/lac/trilinos_block_sparse_matrix.h>
#  include <deal.II/lac/trilinos_index_access.h>


DEAL_II_NAMESPACE_OPEN

namespace TrilinosWrappers
{
  namespace MPI
  {
    BlockVector &
    BlockVector::operator=(const value_type s)
    {
      BaseClass::operator=(s);
      return *this;
    }



    BlockVector &
    BlockVector::operator=(const BlockVector &v)
    {
      // we only allow assignment to vectors with the same number of blocks
      // or to an empty BlockVector
      Assert(n_blocks() == 0 || n_blocks() == v.n_blocks(),
             ExcDimensionMismatch(n_blocks(), v.n_blocks()));

      if (this->n_blocks() != v.n_blocks())
        reinit(v.n_blocks());

      for (size_type i = 0; i < this->n_blocks(); ++i)
        this->components[i] = v.block(i);

      collect_sizes();

      return *this;
    }



    BlockVector &
    BlockVector::operator=(BlockVector &&v) noexcept
    {
      swap(v);
      return *this;
    }



    void
    BlockVector::reinit(const std::vector<IndexSet> &parallel_partitioning,
                        const MPI_Comm &             communicator,
                        const bool                   omit_zeroing_entries)
    {
      const size_type        no_blocks = parallel_partitioning.size();
      std::vector<size_type> block_sizes(no_blocks);

      for (size_type i = 0; i < no_blocks; ++i)
        {
          block_sizes[i] = parallel_partitioning[i].size();
        }

      this->block_indices.reinit(block_sizes);
      if (components.size() != n_blocks())
        components.resize(n_blocks());

      for (size_type i = 0; i < n_blocks(); ++i)
        components[i].reinit(parallel_partitioning[i],
                             communicator,
                             omit_zeroing_entries);

      collect_sizes();
    }



    void
    BlockVector::reinit(const std::vector<IndexSet> &parallel_partitioning,
                        const std::vector<IndexSet> &ghost_values,
                        const MPI_Comm &             communicator,
                        const bool                   vector_writable)
    {
      const size_type        no_blocks = parallel_partitioning.size();
      std::vector<size_type> block_sizes(no_blocks);

      for (size_type i = 0; i < no_blocks; ++i)
        {
          block_sizes[i] = parallel_partitioning[i].size();
        }

      this->block_indices.reinit(block_sizes);
      if (components.size() != n_blocks())
        components.resize(n_blocks());

      for (size_type i = 0; i < n_blocks(); ++i)
        components[i].reinit(parallel_partitioning[i],
                             ghost_values[i],
                             communicator,
                             vector_writable);

      collect_sizes();
    }



    void
    BlockVector::reinit(const BlockVector &v, const bool omit_zeroing_entries)
    {
      block_indices = v.get_block_indices();
      if (components.size() != n_blocks())
        components.resize(n_blocks());

      for (size_type i = 0; i < n_blocks(); ++i)
        components[i].reinit(v.block(i), omit_zeroing_entries, false);

      collect_sizes();
    }



    void
    BlockVector::reinit(const size_type num_blocks)
    {
      std::vector<size_type> block_sizes(num_blocks, 0);
      this->block_indices.reinit(block_sizes);
      if (this->components.size() != this->n_blocks())
        this->components.resize(this->n_blocks());

      for (size_type i = 0; i < this->n_blocks(); ++i)
        components[i].clear();

      collect_sizes();
    }



    void
    BlockVector::import_nonlocal_data_for_fe(
      const TrilinosWrappers::BlockSparseMatrix &m,
      const BlockVector &                        v)
    {
      Assert(m.n_block_rows() == v.n_blocks(),
             ExcDimensionMismatch(m.n_block_rows(), v.n_blocks()));
      Assert(m.n_block_cols() == v.n_blocks(),
             ExcDimensionMismatch(m.n_block_cols(), v.n_blocks()));

      if (v.n_blocks() != n_blocks())
        {
          block_indices = v.get_block_indices();
          components.resize(v.n_blocks());
        }

      for (size_type i = 0; i < this->n_blocks(); ++i)
        components[i].import_nonlocal_data_for_fe(m.block(i, i), v.block(i));

      collect_sizes();
    }



    void
    BlockVector::print(std::ostream &     out,
                       const unsigned int precision,
                       const bool         scientific,
                       const bool         across) const
    {
      for (size_type i = 0; i < this->n_blocks(); ++i)
        {
          if (across)
            out << 'C' << i << ':';
          else
            out << "Component " << i << std::endl;
          this->components[i].print(out, precision, scientific, across);
        }
    }

  } /* end of namespace MPI */
} /* end of namespace TrilinosWrappers */


DEAL_II_NAMESPACE_CLOSE

#endif
