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

#ifndef dealii_trilinos_tpetra_sparse_matrix_templates_h
#define dealii_trilinos_tpetra_sparse_matrix_templates_h

#include <deal.II/base/config.h>

#include <deal.II/base/types.h>

#include <deal.II/lac/trilinos_tpetra_types.h>

#ifdef DEAL_II_TRILINOS_WITH_TPETRA

#  include <deal.II/lac/dynamic_sparsity_pattern.h>
#  include <deal.II/lac/full_matrix.h>
#  include <deal.II/lac/trilinos_tpetra_sparse_matrix.h>
#  include <deal.II/lac/trilinos_tpetra_sparsity_pattern.h>

DEAL_II_NAMESPACE_OPEN

namespace LinearAlgebra
{

  namespace TpetraWrappers
  {
    namespace internal
    {
      template <typename Number, typename MemorySpace>
      void
      apply(const SparseMatrix<Number, MemorySpace> &M,
            const Vector<Number, MemorySpace>       &src,
            Vector<Number, MemorySpace>             &dst,
            Teuchos::ETransp                         mode = Teuchos::NO_TRANS,
            Number alpha = Teuchos::ScalarTraits<Number>::one(),
            Number beta  = Teuchos::ScalarTraits<Number>::zero())
      {
        Assert(
          &src != &dst,
          (typename SparseMatrix<double,
                                 MemorySpace>::ExcSourceEqualsDestination()));
        Assert(M.trilinos_matrix().isFillComplete(),
               (typename SparseMatrix<double,
                                      MemorySpace>::ExcMatrixNotCompressed()));

        if (mode == Teuchos::NO_TRANS)
          {
            Assert(src.trilinos_vector().getMap()->isSameAs(
                     *M.trilinos_matrix().getDomainMap()),
                   (typename SparseMatrix<double,
                                          MemorySpace>::ExcColMapMismatch()));
            Assert(
              dst.trilinos_vector().getMap()->isSameAs(
                *M.trilinos_matrix().getRangeMap()),
              (typename SparseMatrix<double,
                                     MemorySpace>::ExcDomainMapMismatch()));
          }
        else
          {
            Assert(dst.trilinos_vector().getMap()->isSameAs(
                     *M.trilinos_matrix().getDomainMap()),
                   (typename SparseMatrix<double,
                                          MemorySpace>::ExcColMapMismatch()));
            Assert(
              src.trilinos_vector().getMap()->isSameAs(
                *M.trilinos_matrix().getRangeMap()),
              (typename SparseMatrix<double,
                                     MemorySpace>::ExcDomainMapMismatch()));
          }

        M.trilinos_matrix().apply(
          src.trilinos_vector(), dst.trilinos_vector(), mode, alpha, beta);
      }



      template <typename Number, typename MemorySpace>
      void
      apply(const SparseMatrix<Number, MemorySpace> &M,
            const dealii::Vector<Number>            &src,
            dealii::Vector<Number>                  &dst,
            Teuchos::ETransp                         mode = Teuchos::NO_TRANS,
            Number alpha = Teuchos::ScalarTraits<Number>::one(),
            Number beta  = Teuchos::ScalarTraits<Number>::zero())
      {
        Assert(
          &src != &dst,
          (typename SparseMatrix<double,
                                 MemorySpace>::ExcSourceEqualsDestination()));
        Assert(M.trilinos_matrix().isFillComplete(),
               (typename SparseMatrix<double,
                                      MemorySpace>::ExcMatrixNotCompressed()));

        // get the size of the input vectors:
        const size_type dst_local_size = dst.end() - dst.begin();
        const size_type src_local_size = src.end() - src.begin();

        // For the dst vector:
        Kokkos::View<Number **, Kokkos::LayoutLeft, Kokkos::HostSpace>
          kokkos_view_dst(dst.begin(), dst_local_size, 1);

        // get a Kokkos::DualView
        auto mirror_view_dst = Kokkos::create_mirror_view_and_copy(
          typename MemorySpace::kokkos_space{}, kokkos_view_dst);
        typename TpetraTypes::VectorType<Number, MemorySpace>::dual_view_type
          kokkos_dual_view_dst(mirror_view_dst, kokkos_view_dst);

        // create the Tpetra::Vector
        typename TpetraTypes::VectorType<Number, MemorySpace> tpetra_dst(
          M.trilinos_matrix().getRangeMap(), kokkos_dual_view_dst);

        // For the src vector:
        // create a Kokkos::View from the src vector
        Kokkos::View<Number **, Kokkos::LayoutLeft, Kokkos::HostSpace>
          kokkos_view_src(const_cast<Number *>(src.begin()), src_local_size, 1);

        // get a Kokkos::DualView
        auto mirror_view_src = Kokkos::create_mirror_view_and_copy(
          typename MemorySpace::kokkos_space{}, kokkos_view_src);
        typename TpetraTypes::VectorType<Number, MemorySpace>::dual_view_type
          kokkos_dual_view_src(mirror_view_src, kokkos_view_src);

        // create the Tpetra::Vector
        typename TpetraTypes::VectorType<Number, MemorySpace> tpetra_src(
          M.trilinos_matrix().getDomainMap(), kokkos_dual_view_src);

        M.trilinos_matrix().apply(tpetra_src, tpetra_dst, mode, alpha, beta);
      }



      template <typename Number, typename MemorySpace, typename VectorType>
      void
      apply(SparseMatrix<Number, MemorySpace> &,
            const VectorType &,
            VectorType &,
            Teuchos::ETransp,
            Number,
            Number)
      {
        DEAL_II_NOT_IMPLEMENTED();
      }
    } // namespace internal



    // reinit_matrix():
    namespace SparseMatrixImpl
    {
      using size_type = dealii::types::signed_global_dof_index;

      template <typename Number,
                typename MemorySpace,
                typename SparsityPatternType>
      void
      reinit_matrix(
        const IndexSet            &row_parallel_partitioning,
        const IndexSet            &column_parallel_partitioning,
        const SparsityPatternType &sparsity_pattern,
        const bool                 exchange_data,
        const MPI_Comm             communicator,
        Teuchos::RCP<TpetraTypes::MapType<MemorySpace>> &column_space_map,
        Teuchos::RCP<TpetraTypes::MatrixType<Number, MemorySpace>> &matrix)
      {
        // release memory before reallocation
        matrix.reset();

        // Get the Tpetra::Maps
        Teuchos::RCP<TpetraTypes::MapType<MemorySpace>> row_space_map =
          row_parallel_partitioning
            .template make_tpetra_map_rcp<TpetraTypes::NodeType<MemorySpace>>(
              communicator, false);

        column_space_map =
          column_parallel_partitioning
            .template make_tpetra_map_rcp<TpetraTypes::NodeType<MemorySpace>>(
              communicator, false);

        if (column_space_map->getComm()->getRank() == 0)
          {
            AssertDimension(sparsity_pattern.n_rows(),
                            row_parallel_partitioning.size());
            AssertDimension(sparsity_pattern.n_cols(),
                            column_parallel_partitioning.size());
          }

        // if we want to exchange data, build a usual Trilinos sparsity pattern
        // and let that handle the exchange. otherwise, manually create a
        // CrsGraph, which consumes considerably less memory because it can set
        // correct number of indices right from the start
        if (exchange_data)
          {
            SparsityPattern<MemorySpace> trilinos_sparsity;
            trilinos_sparsity.reinit(row_parallel_partitioning,
                                     column_parallel_partitioning,
                                     sparsity_pattern,
                                     communicator,
                                     exchange_data);
            matrix = Utilities::Trilinos::internal::make_rcp<
              TpetraTypes::MatrixType<Number, MemorySpace>>(
              trilinos_sparsity.trilinos_sparsity_pattern());

            return;
          }

        // compute the number of entries per row
        const size_type first_row = row_space_map->getMinGlobalIndex();
        const size_type last_row  = row_space_map->getMaxGlobalIndex() + 1;

        Teuchos::Array<size_t> n_entries_per_row(last_row - first_row);
        for (size_type row = first_row; row < last_row; ++row)
          n_entries_per_row[row - first_row] = sparsity_pattern.row_length(row);

        // The deal.II notion of a 'sparsity pattern' corresponds to the
        // Tpetra concept of a 'graph'. Hence, we generate a graph by copying
        // the sparsity pattern into it, and then build up the matrix from the
        // graph. This is considerable faster than directly filling elements
        // into the matrix. Moreover, it consumes less memory, since the
        // internal reordering is done on ints only, and we can leave the
        // doubles aside.
        Teuchos::RCP<TpetraTypes::GraphType<MemorySpace>> graph =
          Utilities::Trilinos::internal::make_rcp<
            TpetraTypes::GraphType<MemorySpace>>(row_space_map,
                                                 n_entries_per_row);

        // This functions assumes that the sparsity pattern sits on all
        // processors (completely). The parallel version uses a Tpetra graph
        // that is already distributed.

        // now insert the indices
        Teuchos::Array<TrilinosWrappers::types::int_type> row_indices;

        for (size_type global_row = first_row; global_row < last_row;
             ++global_row)
          {
            const int row_length = sparsity_pattern.row_length(global_row);
            if (row_length == 0)
              continue;

            row_indices.resize(row_length, -1);
            for (size_type col = 0; col < row_length; ++col)
              row_indices[col] =
                sparsity_pattern.column_number(global_row, col);

            AssertIndexRange(global_row, row_space_map->getGlobalNumElements());
            graph->insertGlobalIndices(global_row, row_indices);
          }

        // Eventually, optimize the graph structure (sort indices, make memory
        // contiguous, etc.). note that the documentation of the function indeed
        // states that we first need to provide the column (domain) map and then
        // the row (range) map
        graph->fillComplete(column_space_map, row_space_map);

        // check whether we got the number of columns right.
        AssertDimension(sparsity_pattern.n_cols(), graph->getGlobalNumCols());

        // And now finally generate the matrix.
        matrix = Utilities::Trilinos::internal::make_rcp<
          TpetraTypes::MatrixType<Number, MemorySpace>>(graph);
      }



      template <typename Number, typename MemorySpace>
      void
      reinit_matrix(
        const IndexSet               &row_parallel_partitioning,
        const IndexSet               &column_parallel_partitioning,
        const DynamicSparsityPattern &sparsity_pattern,
        const bool                    exchange_data,
        const MPI_Comm                communicator,
        Teuchos::RCP<TpetraTypes::MapType<MemorySpace>> &column_space_map,
        Teuchos::RCP<TpetraTypes::MatrixType<Number, MemorySpace>> &matrix)
      {
        // release memory before reallocation
        matrix.reset();

        // Get the Tpetra::Maps
        Teuchos::RCP<TpetraTypes::MapType<MemorySpace>> row_space_map =
          row_parallel_partitioning
            .template make_tpetra_map_rcp<TpetraTypes::NodeType<MemorySpace>>(
              communicator, false);

        column_space_map =
          column_parallel_partitioning
            .template make_tpetra_map_rcp<TpetraTypes::NodeType<MemorySpace>>(
              communicator, false);

        if (column_space_map->getComm()->getRank() == 0)
          {
            AssertDimension(sparsity_pattern.n_rows(),
                            row_parallel_partitioning.size());
            AssertDimension(sparsity_pattern.n_cols(),
                            column_parallel_partitioning.size());
          }

        // if we want to exchange data, build a usual Trilinos sparsity pattern
        // and let that handle the exchange. otherwise, manually create a
        // CrsGraph, which consumes considerably less memory because it can set
        // correct number of indices right from the start
        if (exchange_data)
          {
            SparsityPattern<MemorySpace> trilinos_sparsity;
            trilinos_sparsity.reinit(row_parallel_partitioning,
                                     column_parallel_partitioning,
                                     sparsity_pattern,
                                     communicator,
                                     exchange_data);
            matrix = Utilities::Trilinos::internal::make_rcp<
              TpetraTypes::MatrixType<Number, MemorySpace>>(
              trilinos_sparsity.trilinos_sparsity_pattern());

            return;
          }

        IndexSet relevant_rows(sparsity_pattern.row_index_set());
        // serial case
        if (relevant_rows.size() == 0)
          {
            relevant_rows.set_size(row_space_map->getGlobalNumElements());
            relevant_rows.add_range(0, row_space_map->getGlobalNumElements());
          }
        relevant_rows.compress();


        std::vector<TrilinosWrappers::types::int_type> ghost_rows;
        Teuchos::Array<size_t>                         n_entries_per_row(
#  if DEAL_II_TRILINOS_VERSION_GTE(14, 0, 0)
          row_space_map->getLocalNumElements());
#  else
          row_space_map->getNodeNumElements());
#  endif
        {
          size_type own = 0;
          for (const auto global_row : relevant_rows)
            {
              if (row_space_map->isNodeGlobalElement(global_row))
                n_entries_per_row[own++] =
                  sparsity_pattern.row_length(global_row);
            }
        }

        // The deal.II notation of a Sparsity pattern corresponds to the Tpetra
        // concept of a Graph. Hence, we generate a graph by copying the
        // sparsity pattern into it, and then build up the matrix from the
        // graph. This is considerable faster than directly filling elements
        // into the matrix. Moreover, it consumes less memory, since the
        // internal reordering is done on ints only, and we can leave the
        // doubles aside.
        Teuchos::RCP<TpetraTypes::GraphType<MemorySpace>> graph =
          Utilities::Trilinos::internal::make_rcp<
            TpetraTypes::GraphType<MemorySpace>>(row_space_map,
                                                 n_entries_per_row);

        // This functions assumes that the sparsity pattern sits on all
        // processors (completely). The parallel version uses a Tpetra graph
        // that is already distributed.

        // now insert the indices
        Teuchos::Array<TrilinosWrappers::types::int_type> row_indices;

        for (const auto global_row : relevant_rows)
          {
            const int row_length = sparsity_pattern.row_length(global_row);
            if (row_length == 0)
              continue;

            row_indices.resize(row_length, -1);
            for (size_type col = 0; col < row_length; ++col)
              row_indices[col] =
                sparsity_pattern.column_number(global_row, col);

            AssertIndexRange(global_row, row_space_map->getGlobalNumElements());
            graph->insertGlobalIndices(global_row, row_indices);
          }

        // Eventually, optimize the graph structure (sort indices, make memory
        // contiguous, etc.). note that the documentation of the function indeed
        // states that we first need to provide the column (domain) map and then
        // the row (range) map
        graph->fillComplete(column_space_map, row_space_map);

        // check whether we got the number of columns right.
        AssertDimension(sparsity_pattern.n_cols(), graph->getGlobalNumCols());

        // And now finally generate the matrix.
        matrix = Utilities::Trilinos::internal::make_rcp<
          TpetraTypes::MatrixType<Number, MemorySpace>>(graph);
      }
    } // namespace SparseMatrixImpl



    // Constructors and initialization:

    // The constructor is actually the only point where we have to check
    // whether we build a serial or a parallel Trilinos matrix.
    // Actually, it does not even matter how many threads there are, but
    // only if we use an MPI compiler or a standard compiler. So, even one
    // thread on a configuration with MPI will still get a parallel interface.
    template <typename Number, typename MemorySpace>
    SparseMatrix<Number, MemorySpace>::SparseMatrix()
      : column_space_map(Utilities::Trilinos::internal::make_rcp<
                         TpetraTypes::MapType<MemorySpace>>(
          0,
          0,
          Utilities::Trilinos::tpetra_comm_self()))
    {
      // Prepare the graph
      Teuchos::RCP<TpetraTypes::GraphType<MemorySpace>> graph =
        Utilities::Trilinos::internal::make_rcp<
          TpetraTypes::GraphType<MemorySpace>>(column_space_map,
                                               column_space_map,
                                               0);
      graph->fillComplete();

      // Create the matrix from the graph
      matrix = Utilities::Trilinos::internal::make_rcp<
        TpetraTypes::MatrixType<Number, MemorySpace>>(graph);

      compressed = false;
    }



    template <typename Number, typename MemorySpace>
    SparseMatrix<Number, MemorySpace>::SparseMatrix(
      const SparsityPattern<MemorySpace> &sparsity_pattern)
      : matrix(Utilities::Trilinos::internal::make_rcp<
               TpetraTypes::MatrixType<Number, MemorySpace>>(
          sparsity_pattern.trilinos_sparsity_pattern()))
    {
      column_space_map =
        Teuchos::rcp_const_cast<TpetraTypes::MapType<MemorySpace>>(
          sparsity_pattern.domain_partitioner());
      compressed = false;
      compress(VectorOperation::add);
    }



    template <typename Number, typename MemorySpace>
    SparseMatrix<Number, MemorySpace>::SparseMatrix(
      const size_type    m,
      const size_type    n,
      const unsigned int n_max_entries_per_row)
      : column_space_map(Utilities::Trilinos::internal::make_rcp<
                         TpetraTypes::MapType<MemorySpace>>(
          static_cast<dealii::types::signed_global_dof_index>(n),
          0,
          Utilities::Trilinos::tpetra_comm_self()))
      , matrix(Utilities::Trilinos::internal::make_rcp<
               TpetraTypes::MatrixType<Number, MemorySpace>>(
          Utilities::Trilinos::internal::make_rcp<
            TpetraTypes::MapType<MemorySpace>>(
            static_cast<dealii::types::signed_global_dof_index>(m),
            0,
            Utilities::Trilinos::tpetra_comm_self()),
          column_space_map,
          n_max_entries_per_row))
      , compressed(false)
    {}



    template <typename Number, typename MemorySpace>
    SparseMatrix<Number, MemorySpace>::SparseMatrix(
      const size_type                  m,
      const size_type                  n,
      const std::vector<unsigned int> &n_entries_per_row)
      : column_space_map(Utilities::Trilinos::internal::make_rcp<
                         TpetraTypes::MapType<MemorySpace>>(
          n,
          0,
          Utilities::Trilinos::tpetra_comm_self()))
      , compressed(false)
    {
      std::vector<size_t> entries_per_row_size_type(n_entries_per_row.begin(),
                                                    n_entries_per_row.end());
      matrix = Utilities::Trilinos::internal::make_rcp<
        TpetraTypes::MatrixType<Number, MemorySpace>>(
        Utilities::Trilinos::internal::make_rcp<
          TpetraTypes::MapType<MemorySpace>>(
          m, 0, Utilities::Trilinos::tpetra_comm_self()),
        column_space_map,
        Teuchos::ArrayView<size_t>{entries_per_row_size_type});
    }



    template <typename Number, typename MemorySpace>
    SparseMatrix<Number, MemorySpace>::SparseMatrix(
      SparseMatrix<Number, MemorySpace> &&other) noexcept
      : column_space_map(std::move(other.column_space_map))
      , matrix(std::move(other.matrix))
      , compressed(std::move(other.compressed))
    {
      other.compressed = false;
    }



    template <typename Number, typename MemorySpace>
    SparseMatrix<Number, MemorySpace> &
    SparseMatrix<Number, MemorySpace>::operator=(
      SparseMatrix<Number, MemorySpace> &&other) noexcept
    {
      column_space_map = std::move(other.column_space_map);
      matrix           = std::move(other.matrix);
      compressed       = std::move(other.compressed);

      return *this;
    }



    template <typename Number, typename MemorySpace>
    template <typename SparsityPatternType>
    void
    SparseMatrix<Number, MemorySpace>::reinit(
      const SparsityPatternType &sparsity_pattern)
    {
      SparseMatrixImpl::reinit_matrix<Number, MemorySpace, SparsityPatternType>(
        complete_index_set(sparsity_pattern.n_rows()),
        complete_index_set(sparsity_pattern.n_cols()),
        sparsity_pattern,
        false,
        MPI_COMM_SELF,
        column_space_map,
        matrix);

      compressed = false;
      compress(VectorOperation::add);
    }



    template <typename Number, typename MemorySpace>
    void
    SparseMatrix<Number, MemorySpace>::reinit(
      const SparsityPattern<MemorySpace> &sparsity_pattern)
    {
      column_space_map.reset();
      matrix.reset();

      // reinit with a (distributed) Trilinos sparsity pattern.
      column_space_map =
        Teuchos::rcp_const_cast<TpetraTypes::MapType<MemorySpace>>(
          sparsity_pattern.domain_partitioner());
      matrix = Utilities::Trilinos::internal::make_rcp<
        TpetraTypes::MatrixType<Number, MemorySpace>>(
        sparsity_pattern.trilinos_sparsity_pattern());

      compressed = false;
      compress(VectorOperation::add);
    }



    // Constructors and initialization using an IndexSet description:

    template <typename Number, typename MemorySpace>
    SparseMatrix<Number, MemorySpace>::SparseMatrix(
      const IndexSet    &parallel_partitioning,
      const MPI_Comm     communicator,
      const unsigned int n_max_entries_per_row)
      : column_space_map(
          parallel_partitioning.template make_tpetra_map_rcp<
            TpetraTypes::NodeType<MemorySpace>>(communicator, false))
      , matrix(Utilities::Trilinos::internal::make_rcp<
               TpetraTypes::MatrixType<Number, MemorySpace>>(
          column_space_map,
          n_max_entries_per_row))
      , compressed(false)
    {}



    template <typename Number, typename MemorySpace>
    SparseMatrix<Number, MemorySpace>::SparseMatrix(
      const IndexSet                  &parallel_partitioning,
      const MPI_Comm                   communicator,
      const std::vector<unsigned int> &n_entries_per_row)
      : column_space_map(
          parallel_partitioning.template make_tpetra_map_rcp<
            TpetraTypes::NodeType<MemorySpace>>(communicator, false))
      , compressed(false)
    {
      Teuchos::Array<size_t> n_entries_per_row_array(n_entries_per_row.begin(),
                                                     n_entries_per_row.end());
      matrix = Utilities::Trilinos::internal::make_rcp<
        TpetraTypes::MatrixType<Number, MemorySpace>>(column_space_map,
                                                      n_entries_per_row_array);
    }



    template <typename Number, typename MemorySpace>
    SparseMatrix<Number, MemorySpace>::SparseMatrix(
      const IndexSet &row_parallel_partitioning,
      const IndexSet &col_parallel_partitioning,
      const MPI_Comm  communicator,
      const size_type n_max_entries_per_row)
      : column_space_map(
          col_parallel_partitioning.template make_tpetra_map_rcp<
            TpetraTypes::NodeType<MemorySpace>>(communicator, false))
      , matrix(Utilities::Trilinos::internal::make_rcp<
               TpetraTypes::MatrixType<Number, MemorySpace>>(
          row_parallel_partitioning.template make_tpetra_map_rcp<
            TpetraTypes::NodeType<MemorySpace>>(communicator, false),
          n_max_entries_per_row))
      , compressed(false)
    {}



    template <typename Number, typename MemorySpace>
    SparseMatrix<Number, MemorySpace>::SparseMatrix(
      const IndexSet                  &row_parallel_partitioning,
      const IndexSet                  &col_parallel_partitioning,
      const MPI_Comm                   communicator,
      const std::vector<unsigned int> &n_entries_per_row)
      : column_space_map(
          col_parallel_partitioning.template make_tpetra_map_rcp<
            TpetraTypes::NodeType<MemorySpace>>(communicator, false))
      , compressed(false)
    {
      Teuchos::Array<size_t> n_entries_per_row_array(n_entries_per_row.begin(),
                                                     n_entries_per_row.end());

      matrix = Utilities::Trilinos::internal::make_rcp<
        TpetraTypes::MatrixType<Number, MemorySpace>>(
        row_parallel_partitioning
          .template make_tpetra_map_rcp<TpetraTypes::NodeType<MemorySpace>>(
            communicator, false),
        n_entries_per_row_array);
    }



    template <typename Number, typename MemorySpace>
    template <typename SparsityPatternType>
    std::enable_if_t<
      !std::is_same_v<SparsityPatternType, dealii::SparseMatrix<double>>>
    SparseMatrix<Number, MemorySpace>::reinit(
      const IndexSet            &parallel_partitioning,
      const SparsityPatternType &sparsity_pattern,
      const MPI_Comm             communicator,
      const bool                 exchange_data)
    {
      reinit(parallel_partitioning,
             parallel_partitioning,
             sparsity_pattern,
             communicator,
             exchange_data);
    }



    template <typename Number, typename MemorySpace>
    template <typename SparsityPatternType>
    std::enable_if_t<
      !std::is_same_v<SparsityPatternType, dealii::SparseMatrix<double>>>
    SparseMatrix<Number, MemorySpace>::reinit(
      const IndexSet            &row_parallel_partitioning,
      const IndexSet            &col_parallel_partitioning,
      const SparsityPatternType &sparsity_pattern,
      const MPI_Comm             communicator,
      const bool                 exchange_data)
    {
      SparseMatrixImpl::reinit_matrix<Number, MemorySpace, SparsityPatternType>(
        row_parallel_partitioning,
        col_parallel_partitioning,
        sparsity_pattern,
        exchange_data,
        communicator,
        column_space_map,
        matrix);

      compressed = false;
      compress(VectorOperation::add);
    }



    template <typename Number, typename MemorySpace>
    void
    SparseMatrix<Number, MemorySpace>::reinit(
      const IndexSet                     &row_parallel_partitioning,
      const IndexSet                     &col_parallel_partitioning,
      const dealii::SparseMatrix<Number> &dealii_sparse_matrix,
      const MPI_Comm                      communicator,
      const double                        drop_tolerance,
      const bool                          copy_values,
      const dealii::SparsityPattern      *use_this_sparsity)
    {
      const size_type n_rows = dealii_sparse_matrix.m();
      AssertDimension(row_parallel_partitioning.size(), n_rows);
      AssertDimension(col_parallel_partitioning.size(),
                      dealii_sparse_matrix.n());

      const dealii::SparsityPattern &sparsity_pattern =
        (use_this_sparsity != nullptr) ?
          *use_this_sparsity :
          dealii_sparse_matrix.get_sparsity_pattern();

      if (matrix.is_null() || m() != n_rows ||
          n_nonzero_elements() != sparsity_pattern.n_nonzero_elements() ||
          copy_values)
        if (use_this_sparsity == nullptr)
          reinit(row_parallel_partitioning,
                 col_parallel_partitioning,
                 sparsity_pattern,
                 communicator,
                 false);

      // in case we do not copy values, we are done
      if (copy_values == false)
        return;

        // fill the values: go through all rows of the
        // matrix, and then all columns. since the sparsity patterns of the
        // input matrix and the specified sparsity pattern might be different,
        // need to go through the row for both these sparsity structures
        // simultaneously in order to really set the correct values.
#  if DEAL_II_TRILINOS_VERSION_GTE(14, 0, 0)
      size_type maximum_row_length = matrix->getLocalMaxNumRowEntries();
#  else
      size_type maximum_row_length = matrix->getNodeMaxNumRowEntries();
#  endif
      std::vector<size_type> row_indices(maximum_row_length);
      std::vector<Number>    values(maximum_row_length);

      for (size_type row = 0; row < n_rows; ++row)
        // see if the row is locally stored on this processor
        if (row_parallel_partitioning.is_element(row) == true)
          {
            dealii::SparsityPattern::iterator select_index =
              sparsity_pattern.begin(row);
            typename dealii::SparseMatrix<Number>::const_iterator it =
              dealii_sparse_matrix.begin(row);
            size_type col = 0;
            if (sparsity_pattern.n_rows() == sparsity_pattern.n_cols())
              {
                // optimized diagonal
                AssertDimension(it->column(), row);
                if (std::fabs(it->value()) > drop_tolerance)
                  {
                    values[col]        = it->value();
                    row_indices[col++] = it->column();
                  }
                ++select_index;
                ++it;
              }

            while (it != dealii_sparse_matrix.end(row) &&
                   select_index != sparsity_pattern.end(row))
              {
                while (select_index->column() < it->column() &&
                       select_index != sparsity_pattern.end(row))
                  ++select_index;
                while (it->column() < select_index->column() &&
                       it != dealii_sparse_matrix.end(row))
                  ++it;

                if (it == dealii_sparse_matrix.end(row))
                  break;
                if (std::fabs(it->value()) > drop_tolerance)
                  {
                    values[col]        = it->value();
                    row_indices[col++] = it->column();
                  }
                ++select_index;
                ++it;
              }
            set(row,
                col,
                reinterpret_cast<size_type *>(row_indices.data()),
                values.data(),
                false);
          }

      compress(VectorOperation::insert);
    }



    // Information on the matrix

    template <typename Number, typename MemorySpace>
    inline unsigned int
    SparseMatrix<Number, MemorySpace>::local_size() const
    {
#  if DEAL_II_TRILINOS_VERSION_GTE(14, 0, 0)
      return matrix->getLocalNumRows();
#  else
      return matrix->getNodeNumRows();
#  endif
    }



    template <typename Number, typename MemorySpace>
    inline std::pair<typename SparseMatrix<Number, MemorySpace>::size_type,
                     typename SparseMatrix<Number, MemorySpace>::size_type>
    SparseMatrix<Number, MemorySpace>::local_range() const
    {
      size_type begin, end;
      begin = matrix->getRowMap()->getMinLocalIndex();
      end   = matrix->getRowMap()->getMaxLocalIndex() + 1;

      return std::make_pair(begin, end);
    }



    template <typename Number, typename MemorySpace>
    inline size_t
    SparseMatrix<Number, MemorySpace>::n_nonzero_elements() const
    {
      return matrix->getGlobalNumEntries();
    }



    template <typename Number, typename MemorySpace>
    MPI_Comm
    SparseMatrix<Number, MemorySpace>::get_mpi_communicator() const
    {
      return Utilities::Trilinos::teuchos_comm_to_mpi_comm(matrix->getComm());
    }



    // Modifying entries

    template <typename Number, typename MemorySpace>
    SparseMatrix<Number, MemorySpace> &
    SparseMatrix<Number, MemorySpace>::operator=(const double d)
    {
      Assert(d == 0, ExcScalarAssignmentOnlyForZeroValue());

      if (compressed)
        {
          matrix->resumeFill();
          compressed = false;
        }

      // As checked above, we are only allowed to use d==0.0, so pass
      // a constant zero (instead of a run-time value 'd' that *happens* to
      // have a zero value) to the underlying class in hopes that the compiler
      // can optimize this somehow.
      matrix->setAllToScalar(/*d=*/0.0);

      return *this;
    }



    template <typename Number, typename MemorySpace>
    SparseMatrix<Number, MemorySpace> &
    SparseMatrix<Number, MemorySpace>::operator*=(const Number a)
    {
      if (compressed)
        {
          matrix->resumeFill();
          compressed = false;
        }

      matrix->scale(a);

      matrix->fillComplete();
      compressed = true;

      return *this;
    }



    template <typename Number, typename MemorySpace>
    SparseMatrix<Number, MemorySpace> &
    SparseMatrix<Number, MemorySpace>::operator/=(const Number a)
    {
      Assert(a != 0, ExcDivideByZero());

      if (compressed)
        {
          matrix->resumeFill();
          compressed = false;
        }

      const Number factor = 1.0 / a;
      matrix->scale(factor);

      matrix->fillComplete();
      compressed = true;

      return *this;
    }



    template <typename Number, typename MemorySpace>
    Number
    SparseMatrix<Number, MemorySpace>::matrix_norm_square(
      const Vector<Number, MemorySpace> &v) const
    {
      AssertDimension(m(), v.size());
      Assert(matrix->isFillComplete(), ExcMatrixNotCompressed());
      Assert(matrix->getRowMap()->isSameAs(*matrix->getDomainMap()),
             ExcNotQuadratic());

      Vector<Number, MemorySpace> temp_vector;
      temp_vector.reinit(v, true);

      vmult(temp_vector, v);
      return temp_vector * v;
    }



    template <typename Number, typename MemorySpace>
    Number
    SparseMatrix<Number, MemorySpace>::matrix_scalar_product(
      const Vector<Number, MemorySpace> &u,
      const Vector<Number, MemorySpace> &v) const
    {
      AssertDimension(m(), u.size());
      AssertDimension(m(), v.size());
      Assert(matrix->isFillComplete(), ExcMatrixNotCompressed());
      Assert(matrix->getRowMap()->isSameAs(*matrix->getDomainMap()),
             ExcNotQuadratic());

      Vector<Number, MemorySpace> temp_vector;
      temp_vector.reinit(v, true);

      vmult(temp_vector, v);
      return u * temp_vector;
    }



    template <typename Number, typename MemorySpace>
    Number
    SparseMatrix<Number, MemorySpace>::frobenius_norm() const
    {
      Assert(matrix->isFillComplete(), ExcMatrixNotCompressed());
      return matrix->getFrobeniusNorm();
    }



    template <typename Number, typename MemorySpace>
    void
    SparseMatrix<Number, MemorySpace>::set(const size_type  row,
                                           const size_type  n_cols,
                                           const size_type *col_indices,
                                           const Number    *values,
                                           const bool       elide_zero_values)
    {
      AssertIndexRange(row, this->m());
      const types::signed_global_dof_index *col_index_ptr;
      const Number                         *col_value_ptr;
      const types::signed_global_dof_index  trilinos_row = row;
      types::signed_global_dof_index        n_columns;

      boost::container::small_vector<Number, 200> local_value_array(
        elide_zero_values ? n_cols : 0);
      boost::container::small_vector<types::signed_global_dof_index, 200>
        local_index_array(elide_zero_values ? n_cols : 0);

      // If we don't elide zeros, the pointers are already available... need to
      // cast to non-const pointers as that is the format taken by Trilinos (but
      // we will not modify const data)
      if (elide_zero_values == false)
        {
          col_index_ptr =
            reinterpret_cast<const dealii::types::signed_global_dof_index *>(
              col_indices);
          col_value_ptr = values;
          n_columns     = n_cols;
        }
      else
        {
          // Otherwise, extract nonzero values in each row and get the
          // respective indices.
          col_index_ptr = local_index_array.data();
          col_value_ptr = local_value_array.data();

          n_columns = 0;
          for (size_type j = 0; j < n_cols; ++j)
            {
              const double value = values[j];
              AssertIsFinite(value);
              if (value != 0)
                {
                  local_index_array[n_columns] = col_indices[j];
                  local_value_array[n_columns] = value;
                  ++n_columns;
                }
            }

          AssertIndexRange(n_columns, n_cols + 1);
        }

      // We distinguish between two cases: the first one is when the matrix is
      // not filled (i.e., it is possible to add new elements to the sparsity
      // pattern), and the second one is when the pattern is already fixed. In
      // the former case, we add the possibility to insert new values, and in
      // the second we just replace data.

      // If the matrix is marked as compressed, we need to
      // call resumeFill() first.
      if (compressed || matrix->isFillComplete())
        {
          matrix->resumeFill();
          compressed = false;
        }

      if (!matrix->isStaticGraph())
        matrix->insertGlobalValues(trilinos_row,
                                   n_columns,
                                   col_value_ptr,
                                   col_index_ptr);
      else
        matrix->replaceGlobalValues(trilinos_row,
                                    n_columns,
                                    col_value_ptr,
                                    col_index_ptr);
    }



    template <typename Number, typename MemorySpace>
    void
    SparseMatrix<Number, MemorySpace>::set(
      const std::vector<size_type> &indices,
      const FullMatrix<Number>     &values,
      const bool                    elide_zero_values)
    {
      Assert(indices.size() == values.m(),
             ExcDimensionMismatch(indices.size(), values.m()));
      Assert(values.m() == values.n(), ExcNotQuadratic());

      for (size_type i = 0; i < indices.size(); ++i)
        set(indices[i],
            indices.size(),
            indices.data(),
            &values(i, 0),
            elide_zero_values);
    }



    template <typename Number, typename MemorySpace>
    void
    SparseMatrix<Number, MemorySpace>::set(
      const std::vector<size_type> &row_indices,
      const std::vector<size_type> &col_indices,
      const FullMatrix<Number>     &full_matrix,
      const bool                    elide_zero_values)
    {
      Assert(row_indices.size() == full_matrix.m(),
             ExcDimensionMismatch(row_indices.size(), full_matrix.m()));
      Assert(col_indices.size() == full_matrix.n(),
             ExcDimensionMismatch(col_indices.size(), full_matrix.n()));

      for (size_type i = 0; i < row_indices.size(); ++i)
        set(row_indices[i],
            col_indices.size(),
            col_indices.data(),
            &full_matrix(i, 0),
            elide_zero_values);
    }



    template <typename Number, typename MemorySpace>
    void
    SparseMatrix<Number, MemorySpace>::set(
      const size_type               row,
      const std::vector<size_type> &col_indices,
      const std::vector<Number>    &values,
      const bool                    elide_zero_values)
    {
      Assert(col_indices.size() == values.size(),
             ExcDimensionMismatch(col_indices.size(), values.size()));

      set(row,
          col_indices.size(),
          col_indices.data(),
          &values[0],
          elide_zero_values);
    }



    template <typename Number, typename MemorySpace>
    void
    SparseMatrix<Number, MemorySpace>::add(
      const size_type  row,
      const size_type  n_cols,
      const size_type *col_indices,
      const Number    *values,
      const bool       elide_zero_values,
      const bool /*col_indices_are_sorted*/)
    {
      AssertIndexRange(row, this->m());
      for (size_t n = 0; n < n_cols; ++n)
        AssertIndexRange(col_indices[n], this->n());

      // If the matrix is marked as compressed, we need to
      // call resumeFill() first.
      if (compressed || matrix->isFillComplete())
        {
          matrix->resumeFill();
          compressed = false;
        }

      // count zero entries;
      const size_t n_zero_entries =
        (elide_zero_values ? std::count(values, values + n_cols, Number(0)) :
                             0);

      // Exit early if there is nothing to do
      if (n_zero_entries == n_cols)
        return;

      // Convert the input into Teuchos::Array
      Teuchos::Array<types::signed_global_dof_index> col_indices_array(
        n_cols - n_zero_entries);
      Teuchos::Array<Number> values_array(n_cols - n_zero_entries);
      if (elide_zero_values)
        {
          size_t n_columns = 0;
          for (size_t i = 0; i < n_cols; ++i)
            {
              // skip all zero entries, while filling the
              if (values[i] != 0)
                {
                  AssertIsFinite(values[i]);
                  AssertIndexRange(n_columns, n_zero_entries);
                  col_indices_array[n_columns] = col_indices[i];
                  values_array[n_columns]      = values[i];
                  ++n_columns;
                }
            }
        }
      else
        for (size_t i = 0; i < n_cols; ++i)
          {
            AssertIsFinite(values[i]);
            col_indices_array[i] = col_indices[i];
            values_array[i]      = values[i];
          }

      // Sum the values into the global matrix.
      matrix->sumIntoGlobalValues(row, col_indices_array, values_array);
    }



    template <typename Number, typename MemorySpace>
    void
    SparseMatrix<Number, MemorySpace>::add(
      const Number                             factor,
      const SparseMatrix<Number, MemorySpace> &source)
    {
      AssertDimension(source.m(), m());
      AssertDimension(source.n(), n());
      AssertDimension(source.local_range().first, local_range().first);
      AssertDimension(source.local_range().second, local_range().second);
      Assert(matrix->getRowMap()->isSameAs(*source.matrix->getRowMap()),
             ExcMessage(
               "Can only add matrices with same distribution of rows."));
      Assert(matrix->isFillComplete() && source.matrix->isFillComplete(),
             ExcMessage("Addition of matrices is only allowed if the matrices "
                        "are filled, i.e., if compress() has been called."));

      // Tpetra does not have a function that adds a matrix in-place to the
      // we're currently storing. But it (inefficiently) has one that returns a
      // new matrix with the result. As a consequence, replace the current one
      // by the one that is returned.
      //
      // Inconveniently, however, the Tpetra::CrsMatrix::add() function returns
      // a RCP<Tpetra::RowMatrix>, though it guarantees that the stored object
      // is actually a Tpetra::CrsMatrix for our combination of arguments. So,
      // we have to do some casting around to assign things back to 'matrix':
      auto result = matrix->add(factor,
                                *source.matrix,
                                1.0,
                                matrix->getDomainMap(),
                                matrix->getRangeMap(),
                                Teuchos::null);
      // needs a local typedef or Assert thinks it gets more arguments than it
      // does...
      using MatrixType = TpetraTypes::MatrixType<Number, MemorySpace>;
      Assert(dynamic_cast<MatrixType *>(result.get()), ExcInternalError());
      matrix.reset(dynamic_cast<TpetraTypes::MatrixType<Number, MemorySpace> *>(
        result.release().get()));
    }



    template <typename Number, typename MemorySpace>
    void
    SparseMatrix<Number, MemorySpace>::clear_row(const size_type row,
                                                 const Number    new_diag_value)
    {
      clear_rows(ArrayView<const size_type>(row), new_diag_value);
    }



    template <typename Number, typename MemorySpace>
    void
    SparseMatrix<Number, MemorySpace>::clear_rows(
      const ArrayView<const size_type> &rows,
      const Number                      new_diag_value)
    {
      // If the matrix is marked as compressed, we need to
      // call resumeFill() first.
      if (compressed || matrix->isFillComplete())
        {
          matrix->resumeFill();
          compressed = false;
        }

      std::vector<int>    col_indices_vector;
      std::vector<Number> values_vector;

      for (size_type row : rows)
        {
          // Only do this on the rows owned locally on this processor.
          int local_row = matrix->getRowMap()->getLocalElement(row);
          if (local_row != Teuchos::OrdinalTraits<int>::invalid())
            {
              size_t nnz = matrix->getNumEntriesInLocalRow(local_row);
              col_indices_vector.resize(nnz);
              values_vector.resize(nnz);

              typename TpetraTypes::MatrixType<Number, MemorySpace>::
                nonconst_local_inds_host_view_type col_indices(
                  col_indices_vector.data(), nnz);
              typename TpetraTypes::MatrixType<Number, MemorySpace>::
                nonconst_values_host_view_type values(values_vector.data(),
                                                      nnz);

              matrix->getLocalRowCopy(local_row, col_indices, values, nnz);

              const size_t diag_index = std::find(col_indices_vector.begin(),
                                                  col_indices_vector.end(),
                                                  local_row) -
                                        col_indices_vector.begin();

              for (size_t j = 0; j < nnz; ++j)
                if (diag_index != j || new_diag_value == 0)
                  values[j] = 0.;

              if (diag_index != nnz)
                values[diag_index] = new_diag_value;

              [[maybe_unused]] int n_replacements =
                matrix->replaceLocalValues(local_row, col_indices, values);
              AssertDimension(n_replacements, nnz);
            }
        }
      compress(VectorOperation::insert);
    }



    template <typename Number, typename MemorySpace>
    void
    SparseMatrix<Number, MemorySpace>::copy_from(
      const SparseMatrix<Number, MemorySpace> &source)
    {
      if (this == &source)
        return;

      // release memory before reallocation
      matrix.reset();
      column_space_map.reset();

      // TODO:
      // If the source and the target matrix have the same structure, we do
      // not need to perform a deep copy.

      // Perform a deep copy
      matrix = Utilities::Trilinos::internal::make_rcp<
        TpetraTypes::MatrixType<Number, MemorySpace>>(*source.matrix,
                                                      Teuchos::Copy);

      column_space_map =
        Teuchos::rcp_const_cast<TpetraTypes::MapType<MemorySpace>>(
          matrix->getColMap());
      compressed = source.compressed;
    }



    template <typename Number, typename MemorySpace>
    void
    SparseMatrix<Number, MemorySpace>::clear()
    {
      // When we clear the matrix, reset
      // the pointer and generate an
      // empty matrix.
      column_space_map = Utilities::Trilinos::internal::make_rcp<
        TpetraTypes::MapType<MemorySpace>>(
        0, 0, Utilities::Trilinos::tpetra_comm_self());

      // Prepare the graph
      Teuchos::RCP<TpetraTypes::GraphType<MemorySpace>> graph =
        Utilities::Trilinos::internal::make_rcp<
          TpetraTypes::GraphType<MemorySpace>>(column_space_map,
                                               column_space_map,
                                               0);
      graph->fillComplete();

      // Create the matrix from the graph
      matrix = Utilities::Trilinos::internal::make_rcp<
        TpetraTypes::MatrixType<Number, MemorySpace>>(graph);

      compressed = true;
    }



    // Multiplications
    template <typename Number, typename MemorySpace>
    template <typename InputVectorType>
    void
    SparseMatrix<Number, MemorySpace>::vmult(InputVectorType       &dst,
                                             const InputVectorType &src) const
    {
      internal::apply(*this, src, dst);
    }



    template <typename Number, typename MemorySpace>
    template <typename InputVectorType>
    void
    SparseMatrix<Number, MemorySpace>::Tvmult(InputVectorType       &dst,
                                              const InputVectorType &src) const
    {
      internal::apply(*this, src, dst, Teuchos::TRANS);
    }



    template <typename Number, typename MemorySpace>
    template <typename InputVectorType>
    void
    SparseMatrix<Number, MemorySpace>::vmult_add(
      InputVectorType       &dst,
      const InputVectorType &src) const
    {
      internal::apply(*this,
                      src,
                      dst,
                      Teuchos::NO_TRANS,
                      Teuchos::ScalarTraits<Number>::one(),
                      Teuchos::ScalarTraits<Number>::one());
    }



    template <typename Number, typename MemorySpace>
    template <typename InputVectorType>
    void
    SparseMatrix<Number, MemorySpace>::Tvmult_add(
      InputVectorType       &dst,
      const InputVectorType &src) const
    {
      internal::apply(*this,
                      src,
                      dst,
                      Teuchos::TRANS,
                      Teuchos::ScalarTraits<Number>::one(),
                      Teuchos::ScalarTraits<Number>::one());
    }



    template <typename Number, typename MemorySpace>
    void
    SparseMatrix<Number, MemorySpace>::print(
      std::ostream &out,
      const bool    print_detailed_trilinos_information) const
    {
      if (print_detailed_trilinos_information)
        {
          auto teuchos_out = Teuchos::getFancyOStream(Teuchos::rcpFromRef(out));
          matrix->describe(*teuchos_out, Teuchos::VERB_EXTREME);
        }
      else
        {
          typename TpetraTypes::MatrixType<Number,
                                           MemorySpace>::values_host_view_type
            values;
          typename TpetraTypes::MatrixType<Number, MemorySpace>::
            local_inds_host_view_type indices;

#  if DEAL_II_TRILINOS_VERSION_GTE(14, 0, 0)
          for (size_t i = 0; i < matrix->getLocalNumRows(); ++i)
#  else
          for (size_t i = 0; i < matrix->getNodeNumRows(); ++i)
#  endif
            {
              matrix->getLocalRowView(i, indices, values);

              for (size_type j = 0; j < static_cast<size_type>(indices.size());
                   ++j)
                out << "(" << matrix->getRowMap()->getGlobalElement(i) << ","
                    << matrix->getColMap()->getGlobalElement(indices[j]) << ") "
                    << values[j] << std::endl;
            }
        }

      AssertThrow(out.fail() == false, ExcIO());
    }



    template <typename Number, typename MemorySpace>
    void
    SparseMatrix<Number, MemorySpace>::compress(VectorOperation::values)
    {
      if (!compressed)
        {
          matrix->fillComplete(column_space_map, matrix->getRowMap());
          compressed = true;
        }
    }



    template <typename Number, typename MemorySpace>
    void
    SparseMatrix<Number, MemorySpace>::resume_fill()
    {
      if (compressed)
        {
          matrix->resumeFill();
          compressed = false;
        }
    }


    template <typename Number, typename MemorySpace>
    Number
    SparseMatrix<Number, MemorySpace>::element(const size_type i,
                                               const size_type j,
                                               const bool      no_error) const
    {
      // Extract local indices in the matrix.
      const int trilinos_i = matrix->getRowMap()->getLocalElement(i);
      const int trilinos_j = matrix->getColMap()->getLocalElement(j);
      Number    value      = 0.;

      if (trilinos_i == Teuchos::OrdinalTraits<int>::invalid() ||
          trilinos_j == Teuchos::OrdinalTraits<int>::invalid())
        {
          if (no_error)
            return {};
          Assert(false,
                 ExcAccessToNonLocalElement(
                   i, j, local_range().first, local_range().second - 1));
        }
      else
        {
          Assert(matrix->isFillComplete(), ExcMatrixNotCompressed());

          // Prepare pointers for extraction of a view of the row.
          size_t nnz_present = matrix->getNumEntriesInLocalRow(trilinos_i);

          typename TpetraTypes::MatrixType<Number, MemorySpace>::
            nonconst_local_inds_host_view_type col_indices("indices",
                                                           nnz_present);
          typename TpetraTypes::MatrixType<Number, MemorySpace>::
            nonconst_values_host_view_type values("values", nnz_present);

          matrix->getLocalRowCopy(trilinos_i, col_indices, values, nnz_present);

          // Search the index where we look for the value, and then finally get
          // it.
          int local_col_index = 0;
          for (; local_col_index < static_cast<int>(nnz_present);
               ++local_col_index)
            {
              if (col_indices[local_col_index] == trilinos_j)
                break;
            }

          if (local_col_index == static_cast<int>(nnz_present))
            {
              if (no_error)
                return {};
              Assert(false, ExcInvalidIndex(i, j));
            }
          else
            value = values[local_col_index];
        }

      return value;
    }



    template <typename Number, typename MemorySpace>
    Number
    SparseMatrix<Number, MemorySpace>::operator()(const size_type i,
                                                  const size_type j) const
    {
      return element(i, j, /* no_error */ false);
    }



    template <typename Number, typename MemorySpace>
    Number
    SparseMatrix<Number, MemorySpace>::el(const size_type i,
                                          const size_type j) const
    {
      return element(i, j, /* no_error */ true);
    }



    template <typename Number, typename MemorySpace>
    Number
    SparseMatrix<Number, MemorySpace>::diag_element(const size_type i) const
    {
      Assert(m() == n(), ExcNotQuadratic());

      if constexpr (running_in_debug_mode())
        {
          // use operator() in debug mode because it checks if this is a valid
          // element (in parallel)
          return operator()(i, i);
        }
      else
        {
          // Trilinos doesn't seem to have a more efficient way to access the
          // diagonal than by just using the standard el(i,j) function.
          return el(i, i);
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

#endif // dealii_trilinos_tpetra_sparse_matrix_templates_h
