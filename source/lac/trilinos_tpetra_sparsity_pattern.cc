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

#include <deal.II/base/config.h>

#ifdef DEAL_II_TRILINOS_WITH_TPETRA

#  include <deal.II/base/mpi.h>
#  include <deal.II/base/trilinos_utilities.h>

#  include <deal.II/lac/dynamic_sparsity_pattern.h>
#  include <deal.II/lac/sparsity_pattern.h>
#  include <deal.II/lac/trilinos_index_access.h>
#  include <deal.II/lac/trilinos_tpetra_sparsity_pattern.h>

#  include <limits>

DEAL_II_NAMESPACE_OPEN

namespace LinearAlgebra
{

  namespace TpetraWrappers
  {
    namespace SparsityPatternIterators
    {
      template <typename MemorySpace>
      void
      Accessor<MemorySpace>::visit_present_row()
      {
        // if we are asked to visit the past-the-end line, then simply
        // release all our caches and go on with life
        if (static_cast<std::size_t>(this->a_row) == sparsity_pattern->n_rows())
          {
            colnum_cache.reset();
            return;
          }

        // otherwise first flush Trilinos caches if necessary
        if (!sparsity_pattern->is_compressed())
          sparsity_pattern->compress();

        colnum_cache =
          std::make_shared<std::vector<dealii::types::signed_global_dof_index>>(
            sparsity_pattern->row_length(this->a_row));

        if (colnum_cache->size() > 0)
          {
            // get a representation of the present row
            std::size_t ncols;
            typename TpetraTypes::GraphType<
              MemorySpace>::nonconst_global_inds_host_view_type
              column_indices_view(colnum_cache->data(), colnum_cache->size());

            sparsity_pattern->graph->getGlobalRowCopy(this->a_row,
                                                      column_indices_view,
                                                      ncols);
            AssertThrow(ncols == colnum_cache->size(), ExcInternalError());
          }
      }
    } // namespace SparsityPatternIterators


    // The constructor is actually the only point where we have to check whether
    // we build a serial or a parallel Trilinos matrix. Actually, it does not
    // even matter how many threads there are, but only if we use an MPI
    // compiler or a standard compiler. So, even one thread on a configuration
    // with MPI will still get a parallel interface.
    template <typename MemorySpace>
    SparsityPattern<MemorySpace>::SparsityPattern()
    {
      column_space_map = Utilities::Trilinos::internal::make_rcp<
        TpetraTypes::MapType<MemorySpace>>(
        TrilinosWrappers::types::int_type(0),
        TrilinosWrappers::types::int_type(0),
        Utilities::Trilinos::tpetra_comm_self());
      graph = Utilities::Trilinos::internal::make_rcp<
        TpetraTypes::GraphType<MemorySpace>>(column_space_map,
                                             column_space_map,
                                             0);
      graph->fillComplete();
    }



    template <typename MemorySpace>
    SparsityPattern<MemorySpace>::SparsityPattern(
      const size_type m,
      const size_type n,
      const size_type n_entries_per_row)
    {
      reinit(m, n, n_entries_per_row);
    }



    template <typename MemorySpace>
    SparsityPattern<MemorySpace>::SparsityPattern(
      const size_type               m,
      const size_type               n,
      const std::vector<size_type> &n_entries_per_row)
    {
      reinit(m, n, n_entries_per_row);
    }



    template <typename MemorySpace>
    SparsityPattern<MemorySpace>::SparsityPattern(
      SparsityPattern<MemorySpace> &&other) noexcept
      : SparsityPatternBase(std::move(other))
      , column_space_map(std::move(other.column_space_map))
      , graph(std::move(other.graph))
      , nonlocal_graph(std::move(other.nonlocal_graph))
    {}



    // Copy function only works if the sparsity pattern is empty.
    template <typename MemorySpace>
    SparsityPattern<MemorySpace>::SparsityPattern(
      const SparsityPattern<MemorySpace> &input_sparsity)
      : SparsityPatternBase(input_sparsity)
      , column_space_map(Utilities::Trilinos::internal::make_rcp<
                         TpetraTypes::MapType<MemorySpace>>(
          0,
          0,
          Utilities::Trilinos::tpetra_comm_self()))
      , graph(Utilities::Trilinos::internal::make_rcp<
              TpetraTypes::GraphType<MemorySpace>>(column_space_map,
                                                   column_space_map,
                                                   0))
    {
      (void)input_sparsity;
      Assert(input_sparsity.n_rows() == 0,
             ExcMessage(
               "Copy constructor only works for empty sparsity patterns."));
    }



    template <typename MemorySpace>
    SparsityPattern<MemorySpace>::SparsityPattern(
      const IndexSet &parallel_partitioning,
      const MPI_Comm  communicator,
      const size_type n_entries_per_row)
    {
      reinit(parallel_partitioning,
             parallel_partitioning,
             communicator,
             n_entries_per_row);
    }



    template <typename MemorySpace>
    SparsityPattern<MemorySpace>::SparsityPattern(
      const IndexSet               &parallel_partitioning,
      const MPI_Comm                communicator,
      const std::vector<size_type> &n_entries_per_row)
    {
      reinit(parallel_partitioning,
             parallel_partitioning,
             communicator,
             n_entries_per_row);
    }



    template <typename MemorySpace>
    SparsityPattern<MemorySpace>::SparsityPattern(
      const IndexSet &row_parallel_partitioning,
      const IndexSet &col_parallel_partitioning,
      const MPI_Comm  communicator,
      const size_type n_entries_per_row)
    {
      reinit(row_parallel_partitioning,
             col_parallel_partitioning,
             communicator,
             n_entries_per_row);
    }



    template <typename MemorySpace>
    SparsityPattern<MemorySpace>::SparsityPattern(
      const IndexSet               &row_parallel_partitioning,
      const IndexSet               &col_parallel_partitioning,
      const MPI_Comm                communicator,
      const std::vector<size_type> &n_entries_per_row)
    {
      reinit(row_parallel_partitioning,
             col_parallel_partitioning,
             communicator,
             n_entries_per_row);
    }



    template <typename MemorySpace>
    SparsityPattern<MemorySpace>::SparsityPattern(
      const IndexSet &row_parallel_partitioning,
      const IndexSet &col_parallel_partitioning,
      const IndexSet &writable_rows,
      const MPI_Comm  communicator,
      const size_type n_max_entries_per_row)
    {
      reinit(row_parallel_partitioning,
             col_parallel_partitioning,
             writable_rows,
             communicator,
             n_max_entries_per_row);
    }



    template <typename MemorySpace>
    void
    SparsityPattern<MemorySpace>::reinit(const size_type m,
                                         const size_type n,
                                         const size_type n_entries_per_row)
    {
      reinit(complete_index_set(m),
             complete_index_set(n),
             MPI_COMM_SELF,
             n_entries_per_row);
    }



    template <typename MemorySpace>
    void
    SparsityPattern<MemorySpace>::reinit(
      const size_type               m,
      const size_type               n,
      const std::vector<size_type> &n_entries_per_row)
    {
      reinit(complete_index_set(m),
             complete_index_set(n),
             MPI_COMM_SELF,
             n_entries_per_row);
    }



    namespace SparsityPatternImpl
    {
      template <typename MemorySpace>
      using size_type = typename SparsityPattern<MemorySpace>::size_type;

      template <typename MemorySpace>
      void
      reinit_sp(
        const Teuchos::RCP<TpetraTypes::MapType<MemorySpace>> &row_map,
        const Teuchos::RCP<TpetraTypes::MapType<MemorySpace>> &col_map,
        const size_type<MemorySpace>                       n_entries_per_row,
        Teuchos::RCP<TpetraTypes::MapType<MemorySpace>>   &column_space_map,
        Teuchos::RCP<TpetraTypes::GraphType<MemorySpace>> &graph,
        Teuchos::RCP<TpetraTypes::GraphType<MemorySpace>> &nonlocal_graph)
      {
        Assert(row_map->isOneToOne(),
               ExcMessage("Row map must be 1-to-1, i.e., no overlap between "
                          "the maps of different processors."));
        Assert(col_map->isOneToOne(),
               ExcMessage("Column map must be 1-to-1, i.e., no overlap between "
                          "the maps of different processors."));

        nonlocal_graph.reset();
        graph.reset();
        column_space_map = col_map;

        // for more than one processor, need to specify only row map first and
        // let the matrix entries decide about the column map (which says which
        // columns are present in the matrix, not to be confused with the
        // col_map that tells how the domain dofs of the matrix will be
        // distributed). for only one processor, we can directly assign the
        // columns as well. If we use a recent Trilinos version, we can also
        // require building a non-local graph which gives us thread-safe
        // initialization.
        graph = Utilities::Trilinos::internal::make_rcp<
          TpetraTypes::GraphType<MemorySpace>>(row_map,
                                               row_map,
                                               n_entries_per_row);
      }



      template <typename MemorySpace>
      void
      reinit_sp(
        const Teuchos::RCP<TpetraTypes::MapType<MemorySpace>> &row_map,
        const Teuchos::RCP<TpetraTypes::MapType<MemorySpace>> &col_map,
        const std::vector<size_type<MemorySpace>>         &n_entries_per_row,
        Teuchos::RCP<TpetraTypes::MapType<MemorySpace>>   &column_space_map,
        Teuchos::RCP<TpetraTypes::GraphType<MemorySpace>> &graph,
        Teuchos::RCP<TpetraTypes::GraphType<MemorySpace>> &nonlocal_graph)
      {
        Assert(row_map->isOneToOne(),
               ExcMessage("Row map must be 1-to-1, i.e., no overlap between "
                          "the maps of different processors."));
        Assert(col_map->isOneToOne(),
               ExcMessage("Column map must be 1-to-1, i.e., no overlap between "
                          "the maps of different processors."));

        // release memory before reallocation
        nonlocal_graph.reset();
        graph.reset();
        AssertDimension(n_entries_per_row.size(),
                        row_map->getGlobalNumElements());

        column_space_map = col_map;

        // Translate the vector of row lengths into one that only stores
        // those entries that related to the locally stored rows of the matrix:
        Kokkos::DualView<size_t *, typename MemorySpace::kokkos_space>
          local_entries_per_row("local_entries_per_row",
                                row_map->getMaxGlobalIndex() -
                                  row_map->getMinGlobalIndex());

        auto local_entries_per_row_host =
          local_entries_per_row
            .template view<Kokkos::DefaultHostExecutionSpace>();

        std::uint64_t total_size = 0;
        for (unsigned int i = 0; i < local_entries_per_row.extent(0); ++i)
          {
            local_entries_per_row_host(i) =
              n_entries_per_row[row_map->getMinGlobalIndex() + i];
            total_size += local_entries_per_row_host[i];
          }
        local_entries_per_row
          .template modify<Kokkos::DefaultHostExecutionSpace>();
        local_entries_per_row
          .template sync<typename MemorySpace::kokkos_space>();

        AssertThrow(
          total_size < static_cast<std::uint64_t>(
                         std::numeric_limits<
                           dealii::types::signed_global_dof_index>::max()),
          ExcMessage(
            "You are requesting to store more elements than global ordinal type allows."));

        graph = Utilities::Trilinos::internal::make_rcp<
          TpetraTypes::GraphType<MemorySpace>>(row_map,
                                               col_map,
                                               local_entries_per_row);
      }



      template <typename SparsityPatternType, typename MemorySpace>
      void
      reinit_sp(
        const Teuchos::RCP<TpetraTypes::MapType<MemorySpace>> &row_map,
        const Teuchos::RCP<TpetraTypes::MapType<MemorySpace>> &col_map,
        const SparsityPatternType                             &sp,
        [[maybe_unused]] const bool                            exchange_data,
        Teuchos::RCP<TpetraTypes::MapType<MemorySpace>>       &column_space_map,
        Teuchos::RCP<TpetraTypes::GraphType<MemorySpace>>     &graph,
        Teuchos::RCP<TpetraTypes::GraphType<MemorySpace>>     &nonlocal_graph)
      {
        nonlocal_graph.reset();
        graph.reset();

        AssertDimension(sp.n_rows(), row_map->getGlobalNumElements());
        AssertDimension(sp.n_cols(), col_map->getGlobalNumElements());

        column_space_map = Utilities::Trilinos::internal::make_rcp<
          TpetraTypes::MapType<MemorySpace>>(*col_map);

        Assert(row_map->isContiguous() == true,
               ExcMessage(
                 "This function only works if the row map is contiguous."));

        const size_type<MemorySpace> first_row = row_map->getMinGlobalIndex(),
                                     last_row =
                                       row_map->getMaxGlobalIndex() + 1;
        Teuchos::Array<size_t> n_entries_per_row(last_row - first_row);

        for (size_type<MemorySpace> row = first_row; row < last_row; ++row)
          n_entries_per_row[row - first_row] = sp.row_length(row);

        AssertThrow(
          std::accumulate(n_entries_per_row.begin(),
                          n_entries_per_row.end(),
                          std::uint64_t(0)) <
            static_cast<std::uint64_t>(std::numeric_limits<int>::max()),
          ExcMessage(
            "The TrilinosWrappers use Tpetra internally, and "
            "Trilinos/Tpetra was compiled with 'local ordinate = int'. "
            "Therefore, 'signed int' is used to represent local indices, "
            "and only 2,147,483,647 nonzero matrix entries can be stored "
            "on a single process, but you are requesting more than "
            "that. Either use more MPI processes or recompile Trilinos "
            "with 'local ordinate = long long' "));

        if (row_map->getComm()->getSize() > 1)
          graph = Utilities::Trilinos::internal::make_rcp<
            TpetraTypes::GraphType<MemorySpace>>(row_map, n_entries_per_row);
        else
          graph = Utilities::Trilinos::internal::make_rcp<
            TpetraTypes::GraphType<MemorySpace>>(row_map,
                                                 col_map,
                                                 n_entries_per_row);

        AssertDimension(sp.n_rows(), graph->getGlobalNumRows());
        AssertDimension(sp.n_cols(), graph->getGlobalNumEntries());

        std::vector<TrilinosWrappers::types::int_type> row_indices;

        for (size_type<MemorySpace> row = first_row; row < last_row; ++row)
          {
            const TrilinosWrappers::types::int_type row_length =
              sp.row_length(row);
            if (row_length == 0)
              continue;

            row_indices.resize(row_length, -1);
            {
              typename SparsityPatternType::iterator p = sp.begin(row);
              // avoid incrementing p over the end of the current row because
              // it is slow for DynamicSparsityPattern in parallel
              for (int col = 0; col < row_length;)
                {
                  row_indices[col++] = p->column();
                  if (col < row_length)
                    ++p;
                }
            }
            graph->insertGlobalIndices(row, row_length, row_indices.data());
          }

        graph->globalAssemble();
      }
    } // namespace SparsityPatternImpl


    template <typename MemorySpace>
    void
    SparsityPattern<MemorySpace>::reinit(const IndexSet &parallel_partitioning,
                                         const MPI_Comm  communicator,
                                         const size_type n_entries_per_row)
    {
      SparsityPatternBase::resize(parallel_partitioning.size(),
                                  parallel_partitioning.size());
      Teuchos::RCP<TpetraTypes::MapType<MemorySpace>> map =
        parallel_partitioning
          .template make_tpetra_map_rcp<TpetraTypes::NodeType<MemorySpace>>(
            communicator, false);
      SparsityPatternImpl::reinit_sp<MemorySpace>(
        map, map, n_entries_per_row, column_space_map, graph, nonlocal_graph);
    }



    template <typename MemorySpace>
    void
    SparsityPattern<MemorySpace>::reinit(
      const IndexSet               &parallel_partitioning,
      const MPI_Comm                communicator,
      const std::vector<size_type> &n_entries_per_row)
    {
      SparsityPatternBase::resize(parallel_partitioning.size(),
                                  parallel_partitioning.size());
      Teuchos::RCP<TpetraTypes::MapType<MemorySpace>> map =
        parallel_partitioning
          .template make_tpetra_map_rcp<TpetraTypes::NodeType<MemorySpace>>(
            communicator, false);
      SparsityPatternImpl::reinit_sp<MemorySpace>(
        map, map, n_entries_per_row, column_space_map, graph, nonlocal_graph);
    }



    template <typename MemorySpace>
    void
    SparsityPattern<MemorySpace>::reinit(
      const IndexSet &row_parallel_partitioning,
      const IndexSet &col_parallel_partitioning,
      const MPI_Comm  communicator,
      const size_type n_entries_per_row)
    {
      SparsityPatternBase::resize(row_parallel_partitioning.size(),
                                  col_parallel_partitioning.size());
      Teuchos::RCP<TpetraTypes::MapType<MemorySpace>> row_map =
        row_parallel_partitioning
          .template make_tpetra_map_rcp<TpetraTypes::NodeType<MemorySpace>>(
            communicator, false);
      Teuchos::RCP<TpetraTypes::MapType<MemorySpace>> col_map =
        col_parallel_partitioning
          .template make_tpetra_map_rcp<TpetraTypes::NodeType<MemorySpace>>(
            communicator, false);
      SparsityPatternImpl::reinit_sp<MemorySpace>(row_map,
                                                  col_map,
                                                  n_entries_per_row,
                                                  column_space_map,
                                                  graph,
                                                  nonlocal_graph);
    }



    template <typename MemorySpace>
    void
    SparsityPattern<MemorySpace>::reinit(
      const IndexSet               &row_parallel_partitioning,
      const IndexSet               &col_parallel_partitioning,
      const MPI_Comm                communicator,
      const std::vector<size_type> &n_entries_per_row)
    {
      SparsityPatternBase::resize(row_parallel_partitioning.size(),
                                  col_parallel_partitioning.size());
      Teuchos::RCP<TpetraTypes::MapType<MemorySpace>> row_map =
        row_parallel_partitioning
          .template make_tpetra_map_rcp<TpetraTypes::NodeType<MemorySpace>>(
            communicator, false);
      Teuchos::RCP<TpetraTypes::MapType<MemorySpace>> col_map =
        col_parallel_partitioning
          .template make_tpetra_map_rcp<TpetraTypes::NodeType<MemorySpace>>(
            communicator, false);
      SparsityPatternImpl::reinit_sp<MemorySpace>(row_map,
                                                  col_map,
                                                  n_entries_per_row,
                                                  column_space_map,
                                                  graph,
                                                  nonlocal_graph);
    }



    template <typename MemorySpace>
    void
    SparsityPattern<MemorySpace>::reinit(
      const IndexSet &row_parallel_partitioning,
      const IndexSet &col_parallel_partitioning,
      const IndexSet &writable_rows,
      const MPI_Comm  communicator,
      const size_type n_entries_per_row)
    {
      SparsityPatternBase::resize(row_parallel_partitioning.size(),
                                  col_parallel_partitioning.size());
      Teuchos::RCP<TpetraTypes::MapType<MemorySpace>> row_map =
        row_parallel_partitioning
          .template make_tpetra_map_rcp<TpetraTypes::NodeType<MemorySpace>>(
            communicator, false);
      Teuchos::RCP<TpetraTypes::MapType<MemorySpace>> col_map =
        col_parallel_partitioning
          .template make_tpetra_map_rcp<TpetraTypes::NodeType<MemorySpace>>(
            communicator, false);
      SparsityPatternImpl::reinit_sp<MemorySpace>(row_map,
                                                  col_map,
                                                  n_entries_per_row,
                                                  column_space_map,
                                                  graph,
                                                  nonlocal_graph);

      IndexSet nonlocal_partitioner = writable_rows;
      AssertDimension(nonlocal_partitioner.size(),
                      row_parallel_partitioning.size());
      if constexpr (running_in_debug_mode())
        {
          {
            IndexSet tmp = writable_rows & row_parallel_partitioning;
            Assert(tmp == row_parallel_partitioning,
                   ExcMessage(
                     "The set of writable rows passed to this method does not "
                     "contain the locally owned rows, which is not allowed."));
          }
        }
      nonlocal_partitioner.subtract_set(row_parallel_partitioning);
      if (Utilities::MPI::n_mpi_processes(communicator) > 1)
        {
          Teuchos::RCP<TpetraTypes::MapType<MemorySpace>> nonlocal_map =
            nonlocal_partitioner
              .template make_tpetra_map_rcp<TpetraTypes::NodeType<MemorySpace>>(
                communicator, true);
          nonlocal_graph = Utilities::Trilinos::internal::make_rcp<
            TpetraTypes::GraphType<MemorySpace>>(nonlocal_map, col_map, 0);
        }
      else
        Assert(nonlocal_partitioner.n_elements() == 0, ExcInternalError());
    }



    template <typename MemorySpace>
    template <typename SparsityPatternType>
    void
    SparsityPattern<MemorySpace>::reinit(
      const IndexSet            &row_parallel_partitioning,
      const IndexSet            &col_parallel_partitioning,
      const SparsityPatternType &nontrilinos_sparsity_pattern,
      const MPI_Comm             communicator,
      const bool                 exchange_data)
    {
      SparsityPatternBase::resize(row_parallel_partitioning.size(),
                                  col_parallel_partitioning.size());
      Teuchos::RCP<TpetraTypes::MapType<MemorySpace>> row_map =
        row_parallel_partitioning
          .template make_tpetra_map_rcp<TpetraTypes::NodeType<MemorySpace>>(
            communicator, false);
      Teuchos::RCP<TpetraTypes::MapType<MemorySpace>> col_map =
        col_parallel_partitioning
          .template make_tpetra_map_rcp<TpetraTypes::NodeType<MemorySpace>>(
            communicator, false);
      SparsityPatternImpl::reinit_sp<SparsityPatternType, MemorySpace>(
        row_map,
        col_map,
        nontrilinos_sparsity_pattern,
        exchange_data,
        column_space_map,
        graph,
        nonlocal_graph);
    }



    template <typename MemorySpace>
    template <typename SparsityPatternType>
    void
    SparsityPattern<MemorySpace>::reinit(
      const IndexSet            &parallel_partitioning,
      const SparsityPatternType &nontrilinos_sparsity_pattern,
      const MPI_Comm             communicator,
      const bool                 exchange_data)
    {
      AssertDimension(nontrilinos_sparsity_pattern.n_rows(),
                      parallel_partitioning.size());
      AssertDimension(nontrilinos_sparsity_pattern.n_cols(),
                      parallel_partitioning.size());
      SparsityPatternBase::resize(parallel_partitioning.size(),
                                  parallel_partitioning.size());
      Teuchos::RCP<TpetraTypes::MapType<MemorySpace>> map =
        parallel_partitioning
          .template make_tpetra_map_rcp<TpetraTypes::NodeType<MemorySpace>>(
            communicator, false);
      SparsityPatternImpl::reinit_sp<SparsityPatternType, MemorySpace>(
        map,
        map,
        nontrilinos_sparsity_pattern,
        exchange_data,
        column_space_map,
        graph,
        nonlocal_graph);
    }



    template <typename MemorySpace>
    SparsityPattern<MemorySpace> &
    SparsityPattern<MemorySpace>::operator=(
      const SparsityPattern<MemorySpace> &)
    {
      DEAL_II_NOT_IMPLEMENTED();
      return *this;
    }



    template <typename MemorySpace>
    void
    SparsityPattern<MemorySpace>::copy_from(
      const SparsityPattern<MemorySpace> &sp)
    {
      SparsityPatternBase::resize(sp.n_rows(), sp.n_cols());
      column_space_map = Utilities::Trilinos::internal::make_rcp<
        TpetraTypes::MapType<MemorySpace>>(*sp.column_space_map);
      graph = Utilities::Trilinos::internal::make_rcp<
        TpetraTypes::GraphType<MemorySpace>>(*sp.graph);

      if (sp.nonlocal_graph.get() != nullptr)
        nonlocal_graph = Utilities::Trilinos::internal::make_rcp<
          TpetraTypes::GraphType<MemorySpace>>(*sp.nonlocal_graph);
      else
        nonlocal_graph.reset();
    }



    template <typename MemorySpace>
    template <typename SparsityPatternType>
    void
    SparsityPattern<MemorySpace>::copy_from(const SparsityPatternType &sp)
    {
      SparsityPatternBase::resize(sp.n_rows(), sp.n_cols());
      Teuchos::RCP<TpetraTypes::MapType<MemorySpace>> rows =
        Utilities::Trilinos::internal::make_rcp<
          TpetraTypes::MapType<MemorySpace>>(
          TrilinosWrappers::types::int_type(sp.n_rows()),
          0,
          Utilities::Trilinos::tpetra_comm_self());
      Teuchos::RCP<TpetraTypes::MapType<MemorySpace>> columns =
        Utilities::Trilinos::internal::make_rcp<
          TpetraTypes::MapType<MemorySpace>>(
          TrilinosWrappers::types::int_type(sp.n_cols()),
          0,
          Utilities::Trilinos::tpetra_comm_self());

      SparsityPatternImpl::reinit_sp<SparsityPatternType, MemorySpace>(
        rows, columns, sp, false, column_space_map, graph, nonlocal_graph);
    }



    template <typename MemorySpace>
    void
    SparsityPattern<MemorySpace>::clear()
    {
      SparsityPatternBase::resize(0, 0);
      // When we clear the matrix, reset
      // the pointer and generate an
      // empty sparsity pattern.
      column_space_map = Utilities::Trilinos::internal::make_rcp<
        TpetraTypes::MapType<MemorySpace>>(
        TrilinosWrappers::types::int_type(0),
        TrilinosWrappers::types::int_type(0),
        Utilities::Trilinos::tpetra_comm_self());
      graph = Utilities::Trilinos::internal::make_rcp<
        TpetraTypes::GraphType<MemorySpace>>(column_space_map,
                                             column_space_map,
                                             0);
      graph->fillComplete();

      nonlocal_graph.reset();
    }



    template <typename MemorySpace>
    void
    SparsityPattern<MemorySpace>::compress()
    {
      Assert(column_space_map.get(), ExcInternalError());
      if (nonlocal_graph.get() != nullptr)
        {
#  if DEAL_II_TRILINOS_VERSION_GTE(14, 0, 0)
          if (nonlocal_graph->getRowMap()->getLocalNumElements() > 0 &&
              column_space_map->getGlobalNumElements() > 0)
#  else
          if (nonlocal_graph->getRowMap()->getNodeNumElements() > 0 &&
              column_space_map->getGlobalNumElements() > 0)
#  endif
            {
              // Insert dummy element at (row, column) that corresponds to row 0
              // in local index counting.
              TrilinosWrappers::types::int_type row =
                nonlocal_graph->getRowMap()->getGlobalElement(0);
              TrilinosWrappers::types::int_type column = 0;

              // in case we have a square sparsity pattern, add the entry on the
              // diagonal
              if (column_space_map->getGlobalNumElements() ==
                  graph->getRangeMap()->getGlobalNumElements())
                column = row;
                // if not, take a column index that we have ourselves since we
                // know for sure it is there (and it will not create spurious
                // messages to many ranks like putting index 0 on many
                // processors)
#  if DEAL_II_TRILINOS_VERSION_GTE(14, 0, 0)
              else if (column_space_map->getLocalNumElements() > 0)
#  else
              else if (column_space_map->getNodeNumElements() > 0)
#  endif
                column = column_space_map->getGlobalElement(0);
              nonlocal_graph->insertGlobalIndices(row, 1, &column);
            }
#  if DEAL_II_TRILINOS_VERSION_GTE(14, 0, 0)
          Assert(nonlocal_graph->getRowMap()->getLocalNumElements() == 0 ||
                   column_space_map->getGlobalNumElements() == 0,
                 ExcInternalError());
#  else
          Assert(nonlocal_graph->getRowMap()->getNodeNumElements() == 0 ||
                   column_space_map->getGlobalNumElements() == 0,
                 ExcInternalError());
#  endif

          nonlocal_graph->fillComplete(column_space_map, graph->getRowMap());
        }
      graph->fillComplete(column_space_map, graph->getRowMap());

      // Check consistency between the sizes set at the beginning and what
      // Trilinos stores:
      using namespace deal_II_exceptions::internals;
      AssertDimension(n_rows(), graph->getGlobalNumRows());
      AssertDimension(n_cols(), graph->getGlobalNumCols());
    }



    template <typename MemorySpace>
    bool
    SparsityPattern<MemorySpace>::row_is_stored_locally(const size_type i) const
    {
      return graph->getRowMap()->getLocalElement(i) !=
             Teuchos::OrdinalTraits<int>::invalid();
    }



    template <typename MemorySpace>
    bool
    SparsityPattern<MemorySpace>::exists(const size_type i,
                                         const size_type j) const
    {
      if (!row_is_stored_locally(i))
        return false;

      // Extract local indices in  the matrix.
      const auto trilinos_i = graph->getRowMap()->getLocalElement(i);
      const auto trilinos_j = graph->getColMap()->getLocalElement(j);

      typename TpetraTypes::GraphType<MemorySpace>::local_inds_host_view_type
        col_indices;

      // Generate the view.
      graph->getLocalRowView(trilinos_i, col_indices);

      // Search the index
      const size_type local_col_index =
        std::find(col_indices.data(),
                  col_indices.data() + col_indices.size(),
                  trilinos_j) -
        col_indices.data();

      return static_cast<std::size_t>(local_col_index) != col_indices.size();
    }



    template <typename MemorySpace>
    typename SparsityPattern<MemorySpace>::size_type
    SparsityPattern<MemorySpace>::bandwidth() const
    {
      size_type local_b = 0;
      for (int i = 0; i < static_cast<int>(local_size()); ++i)
        {
#  if DEAL_II_TRILINOS_VERSION_GTE(13, 2, 0)
          typename TpetraTypes::GraphType<
            MemorySpace>::local_inds_host_view_type indices;
#  else
          Teuchos::ArrayView<const int> indices;
#  endif

          graph->getLocalRowView(i, indices);
          const auto num_entries = indices.size();
          for (unsigned int j = 0; j < static_cast<unsigned int>(num_entries);
               ++j)
            {
              if (static_cast<size_type>(std::abs(i - indices[j])) > local_b)
                local_b = std::abs(i - indices[j]);
            }
        }

      TrilinosWrappers::types::int_type global_b =
        Utilities::MPI::max(local_b,
                            Utilities::Trilinos::teuchos_comm_to_mpi_comm(
                              graph->getComm()));
      return static_cast<size_type>(global_b);
    }



    template <typename MemorySpace>
    unsigned int
    SparsityPattern<MemorySpace>::local_size() const
    {
#  if DEAL_II_TRILINOS_VERSION_GTE(14, 0, 0)
      return graph->getLocalNumRows();
#  else
      return graph->getNodeNumRows();
#  endif
    }



    template <typename MemorySpace>
    std::pair<typename SparsityPattern<MemorySpace>::size_type,
              typename SparsityPattern<MemorySpace>::size_type>
    SparsityPattern<MemorySpace>::local_range() const
    {
      const size_type begin = graph->getRowMap()->getMinGlobalIndex();
      const size_type end   = graph->getRowMap()->getMaxGlobalIndex() + 1;

      return {begin, end};
    }



    template <typename MemorySpace>
    std::uint64_t
    SparsityPattern<MemorySpace>::n_nonzero_elements() const
    {
      return graph->getGlobalNumEntries();
    }



    template <typename MemorySpace>
    unsigned int
    SparsityPattern<MemorySpace>::max_entries_per_row() const
    {
#  if DEAL_II_TRILINOS_VERSION_GTE(14, 0, 0)
      return graph->getLocalMaxNumRowEntries();
#  else
      return graph->getNodeMaxNumRowEntries();
#  endif
    }



    template <typename MemorySpace>
    typename SparsityPattern<MemorySpace>::size_type
    SparsityPattern<MemorySpace>::row_length(const size_type row) const
    {
      Assert(row < (size_type)n_rows(), ExcInternalError());

      // Get a representation of the where the present row is located on
      // the current processor
      TrilinosWrappers::types::int_type local_row =
        graph->getRowMap()->getLocalElement(row);

      // On the processor who owns this row, we'll have a non-negative
      // value for `local_row` and can ask for the length of the row.
      if (local_row >= 0)
        return graph->getNumEntriesInLocalRow(local_row);
      else
        return static_cast<size_type>(-1);
    }



    template <typename MemorySpace>
    void
    SparsityPattern<MemorySpace>::add_row_entries(
      const dealii::types::global_dof_index                  &row,
      const ArrayView<const dealii::types::global_dof_index> &columns,
      const bool indices_are_sorted)
    {
      add_entries(row, columns.begin(), columns.end(), indices_are_sorted);
    }



    template <typename MemorySpace>
    Teuchos::RCP<const typename TpetraTypes::MapType<MemorySpace>>
    SparsityPattern<MemorySpace>::domain_partitioner() const
    {
      return graph->getDomainMap();
    }



    template <typename MemorySpace>
    Teuchos::RCP<const typename TpetraTypes::MapType<MemorySpace>>
    SparsityPattern<MemorySpace>::range_partitioner() const
    {
      return graph->getRangeMap();
    }



    template <typename MemorySpace>
    MPI_Comm
    SparsityPattern<MemorySpace>::get_mpi_communicator() const
    {
      return Utilities::Trilinos::teuchos_comm_to_mpi_comm(
        graph->getRangeMap()->getComm());
    }



    template <typename MemorySpace>
    Teuchos::RCP<const Teuchos::Comm<int>>
    SparsityPattern<MemorySpace>::get_teuchos_mpi_communicator() const
    {
      return graph->getRangeMap()->getComm();
    }



    // As of now, no particularly neat
    // output is generated in case of
    // multiple processors.
    template <typename MemorySpace>
    void
    SparsityPattern<MemorySpace>::print(
      std::ostream &out,
      const bool    write_extended_trilinos_info) const
    {
      if (write_extended_trilinos_info)
        out << *graph;
      else
        {
#  if DEAL_II_TRILINOS_VERSION_GTE(14, 0, 0)
          for (unsigned int i = 0; i < graph->getLocalNumRows(); ++i)
#  else
          for (unsigned int i = 0; i < graph->getNodeNumRows(); ++i)
#  endif
            {
#  if DEAL_II_TRILINOS_VERSION_GTE(13, 2, 0)
              typename TpetraTypes::GraphType<
                MemorySpace>::local_inds_host_view_type indices;
#  else
              Teuchos::ArrayView<const int> indices;
#  endif
              graph->getLocalRowView(i, indices);
              int num_entries = indices.size();
              for (int j = 0; j < num_entries; ++j)
                out << "(" << graph->getRowMap()->getGlobalElement(i) << ","
                    << graph->getColMap()->getGlobalElement(indices[j]) << ") "
                    << std::endl;
            }
        }

      AssertThrow(out.fail() == false, ExcIO());
    }



    template <typename MemorySpace>
    void
    SparsityPattern<MemorySpace>::print_gnuplot(std::ostream &out) const
    {
      Assert(graph->isFillComplete() == true, ExcInternalError());

      for (unsigned int row = 0; row < local_size(); ++row)
        {
#  if DEAL_II_TRILINOS_VERSION_GTE(13, 2, 0)
          typename TpetraTypes::GraphType<
            MemorySpace>::local_inds_host_view_type indices;
#  else
          Teuchos::ArrayView<const int> indices;
#  endif
          graph->getLocalRowView(row, indices);
          int num_entries = indices.size();

          Assert(num_entries >= 0, ExcInternalError());
          // avoid sign comparison warning
          const dealii::types::signed_global_dof_index num_entries_ =
            num_entries;
          for (dealii::types::signed_global_dof_index j = 0; j < num_entries_;
               ++j)
            // while matrix entries are usually
            // written (i,j), with i vertical and
            // j horizontal, gnuplot output is
            // x-y, that is we have to exchange
            // the order of output
            out << static_cast<int>(
                     graph->getColMap()->getGlobalElement(indices[j]))
                << " "
                << -static_cast<int>(graph->getRowMap()->getGlobalElement(row))
                << std::endl;
        }

      AssertThrow(out.fail() == false, ExcIO());
    }

    // TODO: Implement!
    template <typename MemorySpace>
    std::size_t
    SparsityPattern<MemorySpace>::memory_consumption() const
    {
      DEAL_II_NOT_IMPLEMENTED();
      return 0;
    }


#  ifndef DOXYGEN
    // explicit instantiations
    template class SparsityPattern<dealii::MemorySpace::Host>;

    template void
    SparsityPattern<dealii::MemorySpace::Host>::copy_from(
      const dealii::SparsityPattern &);
    template void
    SparsityPattern<dealii::MemorySpace::Host>::copy_from(
      const dealii::DynamicSparsityPattern &);

    template void
    SparsityPattern<dealii::MemorySpace::Host>::reinit(
      const IndexSet &,
      const dealii::SparsityPattern &,
      const MPI_Comm,
      bool);
    template void
    SparsityPattern<dealii::MemorySpace::Host>::reinit(
      const IndexSet &,
      const dealii::DynamicSparsityPattern &,
      const MPI_Comm,
      bool);


    template void
    SparsityPattern<dealii::MemorySpace::Host>::reinit(
      const IndexSet &,
      const IndexSet &,
      const dealii::SparsityPattern &,
      const MPI_Comm,
      bool);
    template void
    SparsityPattern<dealii::MemorySpace::Host>::reinit(
      const IndexSet &,
      const IndexSet &,
      const dealii::DynamicSparsityPattern &,
      const MPI_Comm,
      bool);


    template class SparsityPattern<dealii::MemorySpace::Default>;

    template void
    SparsityPattern<dealii::MemorySpace::Default>::copy_from(
      const dealii::SparsityPattern &);
    template void
    SparsityPattern<dealii::MemorySpace::Default>::copy_from(
      const dealii::DynamicSparsityPattern &);

    template void
    SparsityPattern<dealii::MemorySpace::Default>::reinit(
      const IndexSet &,
      const dealii::SparsityPattern &,
      const MPI_Comm,
      bool);
    template void
    SparsityPattern<dealii::MemorySpace::Default>::reinit(
      const IndexSet &,
      const dealii::DynamicSparsityPattern &,
      const MPI_Comm,
      bool);


    template void
    SparsityPattern<dealii::MemorySpace::Default>::reinit(
      const IndexSet &,
      const IndexSet &,
      const dealii::SparsityPattern &,
      const MPI_Comm,
      bool);
    template void
    SparsityPattern<dealii::MemorySpace::Default>::reinit(
      const IndexSet &,
      const IndexSet &,
      const dealii::DynamicSparsityPattern &,
      const MPI_Comm,
      bool);

#  endif

  } // namespace TpetraWrappers

} // namespace LinearAlgebra

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_TRILINOS_WITH_TPETRA
