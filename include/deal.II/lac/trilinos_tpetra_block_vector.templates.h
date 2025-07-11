// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2024 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_trilinos_tpetra_block_vector_templates_h
#define dealii_trilinos_tpetra_block_vector_templates_h

#include <deal.II/base/config.h>

#include <deal.II/lac/trilinos_tpetra_block_vector.h>

#ifdef DEAL_II_TRILINOS_WITH_TPETRA

#  include <deal.II/base/index_set.h>
#  include <deal.II/base/trilinos_utilities.h>

DEAL_II_NAMESPACE_OPEN

namespace LinearAlgebra
{
  namespace TpetraWrappers
  {
    template <typename Number, typename MemorySpace>
    BlockVector<Number, MemorySpace>::BlockVector(
      const std::vector<IndexSet> &parallel_partitioning,
      const MPI_Comm               communicator)
    {
      reinit(parallel_partitioning, communicator, false);
    }



    template <typename Number, typename MemorySpace>
    BlockVector<Number, MemorySpace>::BlockVector(
      const std::vector<IndexSet> &parallel_partitioning,
      const std::vector<IndexSet> &ghost_values,
      const MPI_Comm               communicator,
      const bool                   vector_writable)
    {
      reinit(parallel_partitioning,
             ghost_values,
             communicator,
             vector_writable);
    }



    template <typename Number, typename MemorySpace>
    BlockVector<Number, MemorySpace>::BlockVector(
      const BlockVector<Number, MemorySpace> &v)
      : dealii::BlockVectorBase<TpetraWrappers::Vector<Number, MemorySpace>>()
    {
      this->block_indices = v.block_indices;

      this->components.resize(this->n_blocks());
      for (unsigned int i = 0; i < this->n_blocks(); ++i)
        this->components[i] = v.components[i];

      this->collect_sizes();
    }



    template <typename Number, typename MemorySpace>
    BlockVector<Number, MemorySpace>::BlockVector(const size_type num_blocks)
      : dealii::BlockVectorBase<TpetraWrappers::Vector<Number, MemorySpace>>()
    {
      reinit(num_blocks);
    }



    template <typename Number, typename MemorySpace>
    BlockVector<Number, MemorySpace> &
    BlockVector<Number, MemorySpace>::operator=(const Number s)
    {
      BaseClass::operator=(s);
      return *this;
    }



    template <typename Number, typename MemorySpace>
    BlockVector<Number, MemorySpace> &
    BlockVector<Number, MemorySpace>::operator=(const BlockVector &v)
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



    template <typename Number, typename MemorySpace>
    void
    BlockVector<Number, MemorySpace>::reinit(
      const std::vector<IndexSet> &parallel_partitioning,
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



    template <typename Number, typename MemorySpace>
    void
    BlockVector<Number, MemorySpace>::reinit(
      const std::vector<IndexSet> &parallel_partitioning,
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



    template <typename Number, typename MemorySpace>
    void
    BlockVector<Number, MemorySpace>::reinit(
      const BlockVector<Number, MemorySpace> &v,
      const bool                              omit_zeroing_entries)
    {
      if (this->n_blocks() != v.n_blocks())
        this->block_indices = v.get_block_indices();

      this->components.resize(this->n_blocks());
      for (unsigned int i = 0; i < this->n_blocks(); ++i)
        this->components[i].reinit(v.components[i], omit_zeroing_entries);

      this->collect_sizes();
    }



    template <typename Number, typename MemorySpace>
    void
    BlockVector<Number, MemorySpace>::reinit(const size_type num_blocks)
    {
      this->block_indices.reinit(num_blocks, 0);

      this->components.resize(this->n_blocks());
      for (auto &c : this->components)
        c.clear();
    }



    template <typename Number, typename MemorySpace>
    bool
    BlockVector<Number, MemorySpace>::has_ghost_elements() const
    {
      bool ghosted = this->block(0).has_ghost_elements();
      if constexpr (running_in_debug_mode())
        {
          for (unsigned int i = 0; i < this->n_blocks(); ++i)
            Assert(this->block(i).has_ghost_elements() == ghosted,
                   ExcInternalError());
        }
      return ghosted;
    }



    template <typename Number, typename MemorySpace>
    void
    BlockVector<Number, MemorySpace>::print(std::ostream      &out,
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

  } // namespace TpetraWrappers
} // namespace LinearAlgebra

DEAL_II_NAMESPACE_CLOSE

#else

// Make sure the scripts that create the C++20 module input files have
// something to latch on if the preprocessor #ifdef above would
// otherwise lead to an empty content of the file.
DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_TRILINOS_WITH_TPETRA

#endif // dealii_trilinos_tpetra_block_vector_templates_h
