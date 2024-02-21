// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2008 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

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
      Assert(this->n_blocks() == 0 || this->n_blocks() == v.n_blocks(),
             ExcDimensionMismatch(this->n_blocks(), v.n_blocks()));

      if (this->n_blocks() != v.n_blocks())
        this->block_indices = v.block_indices;

      this->components.resize(this->n_blocks());
      for (unsigned int i = 0; i < this->n_blocks(); ++i)
        this->components[i] = v.components[i];

      this->collect_sizes();

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
                        const MPI_Comm               communicator,
                        const bool                   omit_zeroing_entries)
    {
      // update the number of blocks
      this->block_indices.reinit(parallel_partitioning.size(), 0);

      // initialize each block
      this->components.resize(this->n_blocks());
      for (unsigned int i = 0; i < this->n_blocks(); ++i)
        this->components[i].reinit(parallel_partitioning[i],
                                   communicator,
                                   omit_zeroing_entries);

      // update block_indices content
      this->collect_sizes();
    }



    void
    BlockVector::reinit(const std::vector<IndexSet> &parallel_partitioning,
                        const std::vector<IndexSet> &ghost_values,
                        const MPI_Comm               communicator,
                        const bool                   vector_writable)
    {
      AssertDimension(parallel_partitioning.size(), ghost_values.size());

      // update the number of blocks
      this->block_indices.reinit(parallel_partitioning.size(), 0);

      // initialize each block
      this->components.resize(this->n_blocks());
      for (unsigned int i = 0; i < this->n_blocks(); ++i)
        this->components[i].reinit(parallel_partitioning[i],
                                   ghost_values[i],
                                   communicator,
                                   vector_writable);

      // update block_indices content
      this->collect_sizes();
    }



    void
    BlockVector::reinit(
      const std::vector<std::shared_ptr<const Utilities::MPI::Partitioner>>
                &partitioners,
      const bool make_ghosted,
      const bool vector_writable)
    {
      // update the number of blocks
      this->block_indices.reinit(partitioners.size(), 0);

      // initialize each block
      this->components.resize(this->n_blocks());
      for (unsigned int i = 0; i < this->n_blocks(); ++i)
        this->components[i].reinit(partitioners[i],
                                   make_ghosted,
                                   vector_writable);

      // update block_indices content
      this->collect_sizes();
    }



    void
    BlockVector::reinit(const BlockVector &v, const bool omit_zeroing_entries)
    {
      if (this->n_blocks() != v.n_blocks())
        this->block_indices = v.get_block_indices();

      this->components.resize(this->n_blocks());
      for (unsigned int i = 0; i < this->n_blocks(); ++i)
        this->components[i].reinit(v.components[i],
                                   omit_zeroing_entries,
                                   false);

      this->collect_sizes();
    }



    void
    BlockVector::reinit(const size_type num_blocks)
    {
      this->block_indices.reinit(num_blocks, 0);

      this->components.resize(this->n_blocks());
      for (auto &c : this->components)
        c.clear();
    }



    void
    BlockVector::import_nonlocal_data_for_fe(
      const TrilinosWrappers::BlockSparseMatrix &m,
      const BlockVector                         &v)
    {
      AssertDimension(m.n_block_rows(), v.n_blocks());
      AssertDimension(m.n_block_cols(), v.n_blocks());

      if (this->n_blocks() != v.n_blocks())
        this->block_indices = v.get_block_indices();

      this->components.resize(this->n_blocks());
      for (unsigned int i = 0; i < this->n_blocks(); ++i)
        this->components[i].import_nonlocal_data_for_fe(m.block(i, i),
                                                        v.block(i));

      this->collect_sizes();
    }



    void
    BlockVector::print(std::ostream      &out,
                       const unsigned int precision,
                       const bool         scientific,
                       const bool         across) const
    {
      for (unsigned int i = 0; i < this->n_blocks(); ++i)
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
