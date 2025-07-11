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

#ifndef dealii_tpetra_trilinos_block_sparse_matrix_templates_h
#define dealii_tpetra_trilinos_block_sparse_matrix_templates_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_TRILINOS_WITH_TPETRA

#  include <deal.II/lac/trilinos_tpetra_block_sparse_matrix.h>

DEAL_II_NAMESPACE_OPEN

namespace LinearAlgebra
{

  namespace TpetraWrappers
  {

    template <typename Number, typename MemorySpace>
    BlockSparseMatrix<Number, MemorySpace>::~BlockSparseMatrix()
    {
      // delete previous content of
      // the subobjects array
      try
        {
          clear();
        }
      catch (...)
        {}
    }



    template <typename Number, typename MemorySpace>
    BlockSparseMatrix<Number, MemorySpace> &
    BlockSparseMatrix<Number, MemorySpace>::operator=(const Number d)
    {
      Assert(d == 0, ExcScalarAssignmentOnlyForZeroValue());

      for (size_type r = 0; r < this->n_block_rows(); ++r)
        for (size_type c = 0; c < this->n_block_cols(); ++c)
          this->block(r, c) = d;

      return *this;
    }



    template <typename Number, typename MemorySpace>
    void
    BlockSparseMatrix<Number, MemorySpace>::reinit(
      const size_type n_block_rows,
      const size_type n_block_columns)
    {
      // first delete previous content of
      // the subobjects array
      clear();

      // then resize. set sizes of blocks to
      // zero. user will later have to call
      // collect_sizes for this
      this->sub_objects.reinit(n_block_rows, n_block_columns);
      this->row_block_indices.reinit(n_block_rows, 0);
      this->column_block_indices.reinit(n_block_columns, 0);

      // and reinitialize the blocks
      for (size_type r = 0; r < this->n_block_rows(); ++r)
        for (size_type c = 0; c < this->n_block_cols(); ++c)
          {
            BlockType *p = new BlockType();

            Assert(this->sub_objects[r][c] == nullptr, ExcInternalError());
            this->sub_objects[r][c] = p;
          }
    }



    template <typename Number, typename MemorySpace>
    template <typename BlockSparsityPatternType>
    void
    BlockSparseMatrix<Number, MemorySpace>::reinit(
      const std::vector<IndexSet>    &parallel_partitioning,
      const BlockSparsityPatternType &block_sparsity_pattern,
      const MPI_Comm                  communicator,
      const bool                      exchange_data)
    {
      if constexpr (running_in_debug_mode())
        {
          std::vector<typename TpetraTypes::MapType<MemorySpace>> tpetra_maps;
          for (size_type i = 0; i < block_sparsity_pattern.n_block_rows(); ++i)
            tpetra_maps.push_back(
              parallel_partitioning[i]
                .template make_tpetra_map<
                  typename TpetraTypes::NodeType<MemorySpace>>(communicator,
                                                               false));

          Assert(tpetra_maps.size() == block_sparsity_pattern.n_block_rows(),
                 ExcDimensionMismatch(tpetra_maps.size(),
                                      block_sparsity_pattern.n_block_rows()));
          Assert(tpetra_maps.size() == block_sparsity_pattern.n_block_cols(),
                 ExcDimensionMismatch(tpetra_maps.size(),
                                      block_sparsity_pattern.n_block_cols()));

          const size_type n_block_rows = tpetra_maps.size();
          Assert(n_block_rows == block_sparsity_pattern.n_block_rows(),
                 ExcDimensionMismatch(n_block_rows,
                                      block_sparsity_pattern.n_block_rows()));
          Assert(n_block_rows == block_sparsity_pattern.n_block_cols(),
                 ExcDimensionMismatch(n_block_rows,
                                      block_sparsity_pattern.n_block_cols()));
        }


      // Call the other basic reinit function, ...
      reinit(block_sparsity_pattern.n_block_rows(),
             block_sparsity_pattern.n_block_cols());

      // ... set the correct sizes, ...
      this->row_block_indices    = block_sparsity_pattern.get_row_indices();
      this->column_block_indices = block_sparsity_pattern.get_column_indices();

      // ... and then assign the correct
      // data to the blocks.
      for (size_type r = 0; r < this->n_block_rows(); ++r)
        for (size_type c = 0; c < this->n_block_cols(); ++c)
          {
            this->sub_objects[r][c]->reinit(parallel_partitioning[r],
                                            parallel_partitioning[c],
                                            block_sparsity_pattern.block(r, c),
                                            communicator,
                                            exchange_data);
          }
    }



    template <typename Number, typename MemorySpace>
    template <typename BlockSparsityPatternType>
    void
    BlockSparseMatrix<Number, MemorySpace>::reinit(
      const BlockSparsityPatternType &block_sparsity_pattern)
    {
      std::vector<IndexSet> parallel_partitioning;
      for (size_type i = 0; i < block_sparsity_pattern.n_block_rows(); ++i)
        parallel_partitioning.emplace_back(
          complete_index_set(block_sparsity_pattern.block(i, 0).n_rows()));

      reinit(parallel_partitioning, block_sparsity_pattern);
    }



    template <typename Number, typename MemorySpace>
    void
    BlockSparseMatrix<Number, MemorySpace>::reinit(
      const std::vector<IndexSet>               &parallel_partitioning,
      const ::dealii::BlockSparseMatrix<double> &dealii_block_sparse_matrix,
      const MPI_Comm                             communicator,
      const double                               drop_tolerance)
    {
      const size_type n_block_rows = parallel_partitioning.size();

      Assert(n_block_rows == dealii_block_sparse_matrix.n_block_rows(),
             ExcDimensionMismatch(n_block_rows,
                                  dealii_block_sparse_matrix.n_block_rows()));
      Assert(n_block_rows == dealii_block_sparse_matrix.n_block_cols(),
             ExcDimensionMismatch(n_block_rows,
                                  dealii_block_sparse_matrix.n_block_cols()));

      // Call the other basic reinit function ...
      reinit(n_block_rows, n_block_rows);

      // ... and then assign the correct
      // data to the blocks.
      for (size_type r = 0; r < this->n_block_rows(); ++r)
        for (size_type c = 0; c < this->n_block_cols(); ++c)
          {
            this->sub_objects[r][c]->reinit(parallel_partitioning[r],
                                            parallel_partitioning[c],
                                            dealii_block_sparse_matrix.block(r,
                                                                             c),
                                            communicator,
                                            drop_tolerance);
          }

      collect_sizes();
    }



    template <typename Number, typename MemorySpace>
    void
    BlockSparseMatrix<Number, MemorySpace>::reinit(
      const ::dealii::BlockSparseMatrix<double> &dealii_block_sparse_matrix,
      const double                               drop_tolerance)
    {
      Assert(dealii_block_sparse_matrix.n_block_rows() ==
               dealii_block_sparse_matrix.n_block_cols(),
             ExcDimensionMismatch(dealii_block_sparse_matrix.n_block_rows(),
                                  dealii_block_sparse_matrix.n_block_cols()));
      Assert(dealii_block_sparse_matrix.m() == dealii_block_sparse_matrix.n(),
             ExcDimensionMismatch(dealii_block_sparse_matrix.m(),
                                  dealii_block_sparse_matrix.n()));

      std::vector<IndexSet> parallel_partitioning;
      for (size_type i = 0; i < dealii_block_sparse_matrix.n_block_rows(); ++i)
        parallel_partitioning.emplace_back(
          complete_index_set(dealii_block_sparse_matrix.block(i, 0).m()));

      reinit(parallel_partitioning,
             dealii_block_sparse_matrix,
             MPI_COMM_SELF,
             drop_tolerance);
    }



    template <typename Number, typename MemorySpace>
    inline bool
    BlockSparseMatrix<Number, MemorySpace>::is_compressed() const
    {
      bool compressed = true;
      for (size_type row = 0; row < this->n_block_rows(); ++row)
        for (size_type col = 0; col < this->n_block_cols(); ++col)
          if (this->block(row, col).is_compressed() == false)
            {
              compressed = false;
              break;
            }

      return compressed;
    }



    template <typename Number, typename MemorySpace>
    void
    BlockSparseMatrix<Number, MemorySpace>::collect_sizes()
    {
      // simply forward to the (non-public) function of the base class
      BaseClass::collect_sizes();
    }



    template <typename Number, typename MemorySpace>
    std::uint64_t
    BlockSparseMatrix<Number, MemorySpace>::n_nonzero_elements() const
    {
      std::uint64_t n_nonzero = 0;
      for (size_type rows = 0; rows < this->n_block_rows(); ++rows)
        for (size_type cols = 0; cols < this->n_block_cols(); ++cols)
          n_nonzero += this->block(rows, cols).n_nonzero_elements();

      return n_nonzero;
    }



    template <typename Number, typename MemorySpace>
    MPI_Comm
    BlockSparseMatrix<Number, MemorySpace>::get_mpi_communicator() const
    {
      Assert(this->n_block_cols() != 0, ExcNotInitialized());
      Assert(this->n_block_rows() != 0, ExcNotInitialized());
      return this->sub_objects[0][0]->get_mpi_communicator();
    }



    template <typename Number, typename MemorySpace>
    inline std::vector<IndexSet>
    BlockSparseMatrix<Number, MemorySpace>::locally_owned_domain_indices() const
    {
      Assert(this->n_block_cols() != 0, ExcNotInitialized());
      Assert(this->n_block_rows() != 0, ExcNotInitialized());

      std::vector<IndexSet> domain_indices;
      for (size_type c = 0; c < this->n_block_cols(); ++c)
        domain_indices.push_back(
          this->sub_objects[0][c]->locally_owned_domain_indices());

      return domain_indices;
    }



    template <typename Number, typename MemorySpace>
    inline std::vector<IndexSet>
    BlockSparseMatrix<Number, MemorySpace>::locally_owned_range_indices() const
    {
      Assert(this->n_block_cols() != 0, ExcNotInitialized());
      Assert(this->n_block_rows() != 0, ExcNotInitialized());

      std::vector<IndexSet> range_indices;
      for (size_type r = 0; r < this->n_block_rows(); ++r)
        range_indices.push_back(
          this->sub_objects[r][0]->locally_owned_range_indices());

      return range_indices;
    }



    template <typename Number, typename MemorySpace>
    template <typename VectorType1, typename VectorType2, typename VectorType3>
    Number
    BlockSparseMatrix<Number, MemorySpace>::residual(VectorType1       &dst,
                                                     const VectorType2 &x,
                                                     const VectorType3 &b) const
    {
      vmult(dst, x);
      dst -= b;
      dst *= -1.;

      return dst.l2_norm();
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

#endif // dealii_tpetra_trilinos_block_sparse_matrix_templates_h
