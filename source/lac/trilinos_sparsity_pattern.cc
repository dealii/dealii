// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2023 by the deal.II authors
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

#include <deal.II/lac/trilinos_index_access.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>

#ifdef DEAL_II_WITH_TRILINOS

#  include <deal.II/base/mpi.h>
#  include <deal.II/base/trilinos_utilities.h>

#  include <deal.II/lac/dynamic_sparsity_pattern.h>
#  include <deal.II/lac/sparsity_pattern.h>

#  include <Epetra_Export.h>

#  include <limits>

DEAL_II_NAMESPACE_OPEN

namespace TrilinosWrappers
{
  namespace SparsityPatternIterators
  {
    void
    Accessor::visit_present_row()
    {
      // if we are asked to visit the past-the-end line, then simply
      // release all our caches and go on with life
      if (this->a_row == sparsity_pattern->n_rows())
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
          Tpetra::CrsGraph<int, dealii::types::signed_global_dof_index>::
            nonconst_global_inds_host_view_type column_indices_view(
              colnum_cache->data(), colnum_cache->size());
          sparsity_pattern->graph->getGlobalRowCopy(this->a_row,
                                                    column_indices_view,
                                                    ncols);
          AssertThrow(ncols == colnum_cache->size(), ExcInternalError());
        }
    }
  } // namespace SparsityPatternIterators


  // The constructor is actually the only point where we have to check whether
  // we build a serial or a parallel Trilinos matrix. Actually, it does not even
  // matter how many threads there are, but only if we use an MPI compiler or a
  // standard compiler. So, even one thread on a configuration with MPI will
  // still get a parallel interface.
  SparsityPattern::SparsityPattern()
  {
    column_space_map =
      Teuchos::rcp(new Tpetra::Map<int, dealii::types::signed_global_dof_index>(
        TrilinosWrappers::types::int_type(0),
        TrilinosWrappers::types::int_type(0),
        Utilities::Trilinos::tpetra_comm_self()));
    graph = Teuchos::rcp(
      new Tpetra::FECrsGraph<int, dealii::types::signed_global_dof_index>(
        column_space_map, column_space_map, 0));
    graph->fillComplete();
  }



  SparsityPattern::SparsityPattern(const size_type m,
                                   const size_type n,
                                   const size_type n_entries_per_row)
  {
    reinit(m, n, n_entries_per_row);
  }



  SparsityPattern::SparsityPattern(
    const size_type               m,
    const size_type               n,
    const std::vector<size_type> &n_entries_per_row)
  {
    reinit(m, n, n_entries_per_row);
  }



  SparsityPattern::SparsityPattern(SparsityPattern &&other) noexcept
    : SparsityPatternBase(std::move(other))
    , column_space_map(std::move(other.column_space_map))
    , graph(std::move(other.graph))
    , nonlocal_graph(std::move(other.nonlocal_graph))
  {}



  // Copy function only works if the sparsity pattern is empty.
  SparsityPattern::SparsityPattern(const SparsityPattern &input_sparsity)
    : SparsityPatternBase(input_sparsity)
    , column_space_map(
        new Tpetra::Map<int, dealii::types::signed_global_dof_index>(
          0,
          0,
          Utilities::Trilinos::tpetra_comm_self()))
    , graph(new Tpetra::FECrsGraph<int, dealii::types::signed_global_dof_index>(
        column_space_map,
        column_space_map,
        0))
  {
    (void)input_sparsity;
    Assert(input_sparsity.n_rows() == 0,
           ExcMessage(
             "Copy constructor only works for empty sparsity patterns."));
  }



  SparsityPattern::SparsityPattern(const IndexSet &parallel_partitioning,
                                   const MPI_Comm  communicator,
                                   const size_type n_entries_per_row)
  {
    reinit(parallel_partitioning,
           parallel_partitioning,
           communicator,
           n_entries_per_row);
  }



  SparsityPattern::SparsityPattern(
    const IndexSet &              parallel_partitioning,
    const MPI_Comm                communicator,
    const std::vector<size_type> &n_entries_per_row)
  {
    reinit(parallel_partitioning,
           parallel_partitioning,
           communicator,
           n_entries_per_row);
  }



  SparsityPattern::SparsityPattern(const IndexSet &row_parallel_partitioning,
                                   const IndexSet &col_parallel_partitioning,
                                   const MPI_Comm  communicator,
                                   const size_type n_entries_per_row)
  {
    reinit(row_parallel_partitioning,
           col_parallel_partitioning,
           communicator,
           n_entries_per_row);
  }



  SparsityPattern::SparsityPattern(
    const IndexSet &              row_parallel_partitioning,
    const IndexSet &              col_parallel_partitioning,
    const MPI_Comm                communicator,
    const std::vector<size_type> &n_entries_per_row)
  {
    reinit(row_parallel_partitioning,
           col_parallel_partitioning,
           communicator,
           n_entries_per_row);
  }



  SparsityPattern::SparsityPattern(const IndexSet &row_parallel_partitioning,
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



  void
  SparsityPattern::reinit(const size_type m,
                          const size_type n,
                          const size_type n_entries_per_row)
  {
    reinit(complete_index_set(m),
           complete_index_set(n),
           MPI_COMM_SELF,
           n_entries_per_row);
  }



  void
  SparsityPattern::reinit(const size_type               m,
                          const size_type               n,
                          const std::vector<size_type> &n_entries_per_row)
  {
    reinit(complete_index_set(m),
           complete_index_set(n),
           MPI_COMM_SELF,
           n_entries_per_row);
  }



  namespace
  {
    using size_type = SparsityPattern::size_type;

    void
    reinit_sp(
      const Teuchos::RCP<
        Tpetra::Map<int, dealii::types::signed_global_dof_index>> &row_map,
      const Teuchos::RCP<
        Tpetra::Map<int, dealii::types::signed_global_dof_index>> &col_map,
      const size_type n_entries_per_row,
      Teuchos::RCP<Tpetra::Map<int, dealii::types::signed_global_dof_index>>
        &column_space_map,
      Teuchos::RCP<
        Tpetra::FECrsGraph<int, dealii::types::signed_global_dof_index>> &graph,
      Teuchos::RCP<
        Tpetra::CrsGraph<int, dealii::types::signed_global_dof_index>>
        &nonlocal_graph)
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
        graph = Teuchos::rcp(
          new Tpetra::FECrsGraph<int, dealii::types::signed_global_dof_index>(
            row_map, row_map, n_entries_per_row));
    }



    void
    reinit_sp(
      const Teuchos::RCP<
        Tpetra::Map<int, dealii::types::signed_global_dof_index>> &row_map,
      const Teuchos::RCP<
        Tpetra::Map<int, dealii::types::signed_global_dof_index>> &col_map,
      const std::vector<size_type> &n_entries_per_row,
      Teuchos::RCP<Tpetra::Map<int, dealii::types::signed_global_dof_index>>
        &column_space_map,
      Teuchos::RCP<
        Tpetra::FECrsGraph<int, dealii::types::signed_global_dof_index>> &graph,
      Teuchos::RCP<
        Tpetra::CrsGraph<int, dealii::types::signed_global_dof_index>>
        &nonlocal_graph)
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
      Kokkos::DualView<dealii::types::signed_global_dof_index*> local_entries_per_row("local_entries_per_row",
        row_map->getMaxGlobalIndex() - row_map->getMinGlobalIndex());
      auto local_entries_per_row_host = local_entries_per_row.view<Kokkos::DefaultHostExecutionSpace>();
      std::uint64_t total_size = 0;
      for (unsigned int i = 0; i < local_entries_per_row.extent(0); ++i) {
        local_entries_per_row_host(i) =
          n_entries_per_row[row_map->getMinGlobalIndex() + i];
        total_size += local_entries_per_row_host[i];
      }  
      local_entries_per_row.modify<Kokkos::DefaultHostExecutionSpace>();
      local_entries_per_row.sync<Kokkos::DefaultExecutionSpace>();

      AssertThrow(total_size
                                 <
                    static_cast<std::uint64_t>(std::numeric_limits< dealii::types::signed_global_dof_index>::max()),
                  ExcMessage(
                             "You are requesting to store more elements than global ordinal type allows."));

        graph = Teuchos::rcp(new
          Tpetra::FECrsGraph<int, dealii::types::signed_global_dof_index>(
          row_map, row_map, local_entries_per_row));
    }



    template <typename SparsityPatternType>
    void
    reinit_sp(
      const Teuchos::RCP<Tpetra::Map<int, dealii::types::signed_global_dof_index>> &row_map,
      const Teuchos::RCP<Tpetra::Map<int, dealii::types::signed_global_dof_index>> &col_map,
      const SparsityPatternType &                                     sp,
      const bool exchange_data,
      Teuchos::RCP<Tpetra::Map<int, dealii::types::signed_global_dof_index>>
        &column_space_map,
      Teuchos::RCP<
        Tpetra::FECrsGraph<int, dealii::types::signed_global_dof_index>> &graph,
      Teuchos::RCP<
        Tpetra::CrsGraph<int, dealii::types::signed_global_dof_index>>
        &nonlocal_graph)
    {
      nonlocal_graph.reset();
      graph.reset();

      AssertDimension(sp.n_rows(),
                      row_map->getGlobalNumElements());
      AssertDimension(sp.n_cols(),
                      col_map->getGlobalNumElements());

      column_space_map = Teuchos::rcp(new 
        Tpetra::Map<int, dealii::types::signed_global_dof_index>(col_map));

      Assert(row_map->isContiguous() == true,
             ExcMessage(
               "This function only works if the row map is contiguous."));

      const size_type first_row = row_map->getMinGlobalIndex(),
                      last_row  = row_map->getMaxGlobalIndex() + 1;
      std::vector<int> n_entries_per_row(last_row - first_row);

      // Trilinos wants the row length as an int this is hopefully never going
      // to be a problem.
      for (size_type row = first_row; row < last_row; ++row)
        n_entries_per_row[row - first_row] =
          static_cast<int>(sp.row_length(row));

      AssertThrow(std::accumulate(n_entries_per_row.begin(),
                                  n_entries_per_row.end(),
                                  std::uint64_t(0)) <
                    static_cast<std::uint64_t>(std::numeric_limits<int>::max()),
                  ExcMessage("The TrilinosWrappers use Epetra internally which "
                             "uses 'signed int' to represent local indices. "
                             "Therefore, only 2,147,483,647 nonzero matrix "
                             "entries can be stored on a single process, "
                             "but you are requesting more than that. "
                             "If possible, use more MPI processes."));

      if (row_map->getComm()->getSize() > 1)
        graph = Teuchos::rcp(new
          Tpetra::FECrsGraph<int, dealii::types::signed_global_dof_index>(
          Copy, row_map, n_entries_per_row.data(), false));
      else
        graph = Teuchos::rcp(new
          Tpetra::FECrsGraph<int, dealii::types::signed_global_dof_index>(
          Copy, row_map, col_map, n_entries_per_row.data(), false));

      AssertDimension(sp.n_rows(), graph->getGlobalNumEntries());

      std::vector<TrilinosWrappers::types::int_type> row_indices;

      // Include possibility to exchange data since DynamicSparsityPattern is
      // able to do so
      if (exchange_data == false)
        for (size_type row = first_row; row < last_row; ++row)
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
            graph
              ->Tpetra::CrsGraph<int, dealii::types::signed_global_dof_index>::
                insertGlobalIndices(row, row_length, row_indices.data());
          }
      else
        for (size_type row = 0; row < sp.n_rows(); ++row)
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
            const TrilinosWrappers::types::int_type trilinos_row = row;
            graph->insertGlobalIndices(trilinos_row,
                                       row_length,
                                       row_indices.data());
          }

      // TODO A dynamic_cast fails here, this is suspicious.
      // const auto &range_map =
      //  static_cast<const Tpetra::Map<int,
      //  dealii::types::signed_global_dof_index> &>(*graph->getRangeMap()); //
      //  NOLINT
      graph->globalAssemble();

      // graph->OptimizeStorage();
    }
  } // namespace



  void
  SparsityPattern::reinit(const IndexSet &parallel_partitioning,
                          const MPI_Comm  communicator,
                          const size_type n_entries_per_row)
  {
    SparsityPatternBase::resize(parallel_partitioning.size(),
                                parallel_partitioning.size());
    Teuchos::RCP<Tpetra::Map<int, dealii::types::signed_global_dof_index>> map =
      parallel_partitioning.make_tpetra_map(communicator, false);
    reinit_sp(
      map, map, n_entries_per_row, column_space_map, graph, nonlocal_graph);
  }



  void
  SparsityPattern::reinit(const IndexSet &              parallel_partitioning,
                          const MPI_Comm                communicator,
                          const std::vector<size_type> &n_entries_per_row)
  {
    SparsityPatternBase::resize(parallel_partitioning.size(),
                                parallel_partitioning.size());
    Teuchos::RCP<Tpetra::Map<int, dealii::types::signed_global_dof_index>> map =
      parallel_partitioning.make_tpetra_map(communicator, false);
    reinit_sp(
      map, map, n_entries_per_row, column_space_map, graph, nonlocal_graph);
  }



  void
  SparsityPattern::reinit(const IndexSet &row_parallel_partitioning,
                          const IndexSet &col_parallel_partitioning,
                          const MPI_Comm  communicator,
                          const size_type n_entries_per_row)
  {
    SparsityPatternBase::resize(row_parallel_partitioning.size(),
                                col_parallel_partitioning.size());
    Teuchos::RCP<Tpetra::Map<int, dealii::types::signed_global_dof_index>> row_map =
      row_parallel_partitioning.make_tpetra_map(communicator, false);
    Teuchos::RCP<Tpetra::Map<int, dealii::types::signed_global_dof_index>> col_map =
      col_parallel_partitioning.make_tpetra_map(communicator, false);
    reinit_sp(row_map,
              col_map,
              n_entries_per_row,
              column_space_map,
              graph,
              nonlocal_graph);
  }



  void
  SparsityPattern::reinit(const IndexSet &row_parallel_partitioning,
                          const IndexSet &col_parallel_partitioning,
                          const MPI_Comm  communicator,
                          const std::vector<size_type> &n_entries_per_row)
  {
    SparsityPatternBase::resize(row_parallel_partitioning.size(),
                                col_parallel_partitioning.size());
    Teuchos::RCP<Tpetra::Map<int, dealii::types::signed_global_dof_index>> row_map =
      row_parallel_partitioning.make_tpetra_map(communicator, false);
    Teuchos::RCP<Tpetra::Map<int, dealii::types::signed_global_dof_index>> col_map =
      col_parallel_partitioning.make_tpetra_map(communicator, false);
    reinit_sp(row_map,
              col_map,
              n_entries_per_row,
              column_space_map,
              graph,
              nonlocal_graph);
  }



  void
  SparsityPattern::reinit(const IndexSet &row_parallel_partitioning,
                          const IndexSet &col_parallel_partitioning,
                          const IndexSet &writable_rows,
                          const MPI_Comm  communicator,
                          const size_type n_entries_per_row)
  {
    SparsityPatternBase::resize(row_parallel_partitioning.size(),
                                col_parallel_partitioning.size());
    Teuchos::RCP<Tpetra::Map<int, dealii::types::signed_global_dof_index>> row_map =
      row_parallel_partitioning.make_tpetra_map(communicator, false);
    Teuchos::RCP<Tpetra::Map<int, dealii::types::signed_global_dof_index>> col_map =
      col_parallel_partitioning.make_tpetra_map(communicator, false);
    reinit_sp(row_map,
              col_map,
              n_entries_per_row,
              column_space_map,
              graph,
              nonlocal_graph);

    IndexSet nonlocal_partitioner = writable_rows;
    AssertDimension(nonlocal_partitioner.size(),
                    row_parallel_partitioning.size());
#  ifdef DEBUG
    {
      IndexSet tmp = writable_rows & row_parallel_partitioning;
      Assert(tmp == row_parallel_partitioning,
             ExcMessage(
               "The set of writable rows passed to this method does not "
               "contain the locally owned rows, which is not allowed."));
    }
#  endif
    nonlocal_partitioner.subtract_set(row_parallel_partitioning);
    if (Utilities::MPI::n_mpi_processes(communicator) > 1)
      {
        Teuchos::RCP<Tpetra::Map<int, dealii::types::signed_global_dof_index>> nonlocal_map =
          nonlocal_partitioner.make_tpetra_map(communicator, true);
        nonlocal_graph = Teuchos::rcp(new
          Tpetra::CrsGraph<int, dealii::types::signed_global_dof_index>(
          nonlocal_map, nonlocal_map, 0));
      }
    else
      Assert(nonlocal_partitioner.n_elements() == 0, ExcInternalError());
  }



  template <typename SparsityPatternType>
  void
  SparsityPattern::reinit(
    const IndexSet &           row_parallel_partitioning,
    const IndexSet &           col_parallel_partitioning,
    const SparsityPatternType &nontrilinos_sparsity_pattern,
    const MPI_Comm             communicator,
    const bool                 exchange_data)
  {
    SparsityPatternBase::resize(row_parallel_partitioning.size(),
                                col_parallel_partitioning.size());
    Teuchos::RCP<Tpetra::Map<int, dealii::types::signed_global_dof_index>> row_map =
      row_parallel_partitioning.make_tpetra_map(communicator, false);
    Teuchos::RCP<Tpetra::Map<int, dealii::types::signed_global_dof_index>> col_map =
      col_parallel_partitioning.make_tpetra_map(communicator, false);
    reinit_sp(row_map,
              col_map,
              nontrilinos_sparsity_pattern,
              exchange_data,
              column_space_map,
              graph,
              nonlocal_graph);
  }



  template <typename SparsityPatternType>
  void
  SparsityPattern::reinit(
    const IndexSet &           parallel_partitioning,
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
    Teuchos::RCP<Tpetra::Map<int, dealii::types::signed_global_dof_index>> map =
      parallel_partitioning.make_tpetra_map(communicator, false);
    reinit_sp(map,
              map,
              nontrilinos_sparsity_pattern,
              exchange_data,
              column_space_map,
              graph,
              nonlocal_graph);
  }



  SparsityPattern &
  SparsityPattern::operator=(const SparsityPattern &)
  {
    Assert(false, ExcNotImplemented());
    return *this;
  }



  void
  SparsityPattern::copy_from(const SparsityPattern &sp)
  {
    SparsityPatternBase::resize(sp.n_rows(), sp.n_cols());
    column_space_map = Teuchos::rcp(new
      Tpetra::Map<int, dealii::types::signed_global_dof_index>(
      *sp.column_space_map));
    graph = sp.graph;

    if (sp.nonlocal_graph.get() != nullptr)
      nonlocal_graph = Teuchos::rcp(new
        Tpetra::CrsGraph<int, dealii::types::signed_global_dof_index>(
        *sp.nonlocal_graph));
    else
      nonlocal_graph.reset();
  }



  template <typename SparsityPatternType>
  void
  SparsityPattern::copy_from(const SparsityPatternType &sp)
  {
    SparsityPatternBase::resize(sp.n_rows(), sp.n_cols());
    auto rows = Teuchos::rcp(new Tpetra::Map<int, dealii::types::signed_global_dof_index>(
      TrilinosWrappers::types::int_type(sp.n_rows()),
      0,
      Utilities::Trilinos::tpetra_comm_self()));
    auto columns = Teuchos::rcp(new Tpetra::Map<int, dealii::types::signed_global_dof_index>(
      TrilinosWrappers::types::int_type(sp.n_cols()),
      0,
      Utilities::Trilinos::tpetra_comm_self()));

    reinit_sp(
      rows, columns, sp, false, column_space_map, graph, nonlocal_graph);
  }



  void
  SparsityPattern::clear()
  {
    SparsityPatternBase::resize(0, 0);
    // When we clear the matrix, reset
    // the pointer and generate an
    // empty sparsity pattern.
    column_space_map = Teuchos::rcp(new
      Tpetra::Map<int, dealii::types::signed_global_dof_index>(
      TrilinosWrappers::types::int_type(0),
      TrilinosWrappers::types::int_type(0),
      Utilities::Trilinos::tpetra_comm_self()));
    graph = Teuchos::rcp(new 
      Tpetra::FECrsGraph<int, dealii::types::signed_global_dof_index>(
      column_space_map, column_space_map, 0));
    graph->fillComplete();

    nonlocal_graph.reset();
  }



  void
  SparsityPattern::compress()
  {
    Assert(column_space_map.get(), ExcInternalError());
    if (nonlocal_graph.get() != nullptr)
      {
        if (
            nonlocal_graph->getRowMap()->getLocalNumElements() > 0 &&
            column_space_map->getGlobalNumElements() > 0)
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
            // messages to many ranks like putting index 0 on many processors)
            else if (column_space_map->getLocalNumElements() > 0)
              column = column_space_map->getGlobalElement(0);
            nonlocal_graph->insertGlobalIndices(row, 1, &column);
          }
        Assert(nonlocal_graph->getRowMap()->getLocalNumElements() == 0 ||
                 column_space_map->getGlobalNumElements() == 0,
               ExcInternalError());

          nonlocal_graph->fillComplete(column_space_map, graph->getRangeMap());
        //ierr = nonlocal_graph->OptimizeStorage();
        //AssertThrow(ierr >= 0, ExcTrilinosError(ierr));
        Tpetra::Export<int, dealii::types::signed_global_dof_index> exporter(
          nonlocal_graph->getRowMap(), graph->getRowMap());
        //ierr = graph->export(*nonlocal_graph, exporter, Add);
        //AssertThrow(ierr == 0, ExcTrilinosError(ierr));
        graph->fillComplete(column_space_map, graph->getRangeMap());
      }
    else
      {
        graph->globalAssemble();
      }

    /*try
      {
        ierr = graph->OptimizeStorage();
      }
    catch (const int error_code)
      {
        AssertThrow(
          false,
          ExcMessage(
            "The Tpetra::CrsGraph<int, dealii::types::signed_global_dof_index>::OptimizeStorage() function "
            "has thrown an error with code " +
            std::to_string(error_code) +
            ". You will have to look up the exact meaning of this error "
            "in the Trilinos source code, but oftentimes, this function "
            "throwing an error indicates that you are trying to allocate "
            "more than 2,147,483,647 nonzero entries in the sparsity "
            "pattern on the local process; this will not work because "
            "Epetra indexes entries with a simple 'signed int'."));
      }*/

    // Check consistency between the sizes set at the beginning and what
    // Trilinos stores:
    using namespace deal_II_exceptions::internals;
    Assert(compare_for_equality(n_rows(), graph->getGlobalNumEntries()),
           ExcInternalError());
    Assert(compare_for_equality(n_cols(), graph->getGlobalNumEntries()),
           ExcInternalError());
  }



  bool
  SparsityPattern::row_is_stored_locally(const size_type i) const
  {
    return graph->getRowMap()->getLocalElement(i) != Teuchos::OrdinalTraits<int>::invalid();
  }



  bool
  SparsityPattern::exists(const size_type i, const size_type j) const
  {
    if (!row_is_stored_locally(i))
        return false;
        
// Extract local indices in  the matrix.
        auto trilinos_i =
              graph->getRowMap()->getLocalElement(i);

  auto trilinos_j =
              graph->getColMap()->getLocalElement(j);

        Tpetra::FECrsGraph<int, dealii::types::signed_global_dof_index>::local_inds_host_view_type col_indices;

            // Generate the view.
            graph->getLocalRowView(trilinos_i,
                                                   col_indices);
            // Search the index
            const std::ptrdiff_t local_col_index =
              std::find(col_indices.data(), col_indices.data() + col_indices.size(), trilinos_j) -
              col_indices.data();

            return local_col_index != col_indices.size();
  }



  SparsityPattern::size_type
  SparsityPattern::bandwidth() const
  {
    size_type                         local_b  = 0;
    TrilinosWrappers::types::int_type global_b = 0;
    for (int i = 0; i < static_cast<int>(local_size()); ++i)
      {
        Tpetra::FECrsGraph<int, dealii::types::signed_global_dof_index>::local_inds_host_view_type indices;
        graph->getLocalRowView(i, indices);
        int num_entries = indices.size();
        for (unsigned int j = 0; j < static_cast<unsigned int>(num_entries);
             ++j)
          {
            if (static_cast<size_type>(std::abs(i - indices[j])) > local_b)
              local_b = std::abs(i - indices[j]);
          }
      }
    graph->getComm()->MaxAll(reinterpret_cast<TrilinosWrappers::types::int_type *>(
                           &local_b),
                         &global_b,
                         1);
    return static_cast<size_type>(global_b);
  }



  unsigned int
  SparsityPattern::local_size() const
  {
    int n_rows = graph->getLocalNumRows();

    return n_rows;
  }



  std::pair<SparsityPattern::size_type, SparsityPattern::size_type>
  SparsityPattern::local_range() const
  {
    size_type begin, end;
    begin = graph->getRowMap()->getMinGlobalIndex();
    end   = graph->getRowMap()->getMaxGlobalIndex() + 1;

    return std::make_pair(begin, end);
  }



  std::uint64_t
  SparsityPattern::n_nonzero_elements() const
  {
    TrilinosWrappers::types::int64_type nnz = graph->getGlobalNumEntries();

    return static_cast<std::uint64_t>(nnz);
  }



  unsigned int
  SparsityPattern::max_entries_per_row() const
  {
    int nnz = graph->getLocalMaxNumRowEntries();

    return static_cast<unsigned int>(nnz);
  }



  SparsityPattern::size_type
  SparsityPattern::row_length(const size_type row) const
  {
    Assert(row < n_rows(), ExcInternalError());

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



  void
  SparsityPattern::add_row_entries(const size_type &                 row,
                                   const ArrayView<const size_type> &columns,
                                   const bool indices_are_sorted)
  {
    add_entries(row, columns.begin(), columns.end(), indices_are_sorted);
  }



  const Tpetra::Map<int, dealii::types::signed_global_dof_index> &
  SparsityPattern::domain_partitioner() const
  {
    return *graph->getDomainMap();
  }



  const Tpetra::Map<int, dealii::types::signed_global_dof_index> &
  SparsityPattern::range_partitioner() const
  {
    return *graph->getRangeMap();
  }



  MPI_Comm
  SparsityPattern::get_mpi_communicator() const
  {
    return graph->getRangeMap()->getComm();
  }



  void
  SparsityPattern::write_ascii()
  {
    Assert(false, ExcNotImplemented());
  }



  // As of now, no particularly neat
  // output is generated in case of
  // multiple processors.
  void
  SparsityPattern::print(std::ostream &out,
                         const bool    write_extended_trilinos_info) const
  {
    if (write_extended_trilinos_info)
      out << *graph;
    else
      {
        int *indices;

        for (int i = 0; i < graph->getLocalNumRows(); ++i)
          {
            Tpetra::FECrsGraph<int, dealii::types::signed_global_dof_index>::local_inds_host_view_type indices;
            graph->getLocalRowView(i, indices);
            int num_entries = indices.size();
            for (int j = 0; j < num_entries; ++j)
              out << "(" << graph->getRowMap()->getGlobalElement(i)
                  << ","
                  << graph->getColMap()->getGlobalElement(indices[j])
                  << ") " << std::endl;
          }
      }

    AssertThrow(out.fail() == false, ExcIO());
  }



  void
  SparsityPattern::print_gnuplot(std::ostream &out) const
  {
    Assert(graph->isFillComplete() == true, ExcInternalError());
    for (dealii::types::signed_global_dof_index row = 0; row < local_size();
         ++row)
      {
        Tpetra::FECrsGraph<int, dealii::types::signed_global_dof_index>::local_inds_host_view_type indices;
        graph->getLocalRowView(row, indices);
        int num_entries = indices.size();

        Assert(num_entries >= 0, ExcInternalError());
        // avoid sign comparison warning
        const dealii::types::signed_global_dof_index num_entries_ = num_entries;
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
              << -static_cast<int>(
                   graph->getRowMap()->getGlobalElement(row))
              << std::endl;
      }

    AssertThrow(out.fail() == false, ExcIO());
  }

  // TODO: Implement!
  std::size_t
  SparsityPattern::memory_consumption() const
  {
    Assert(false, ExcNotImplemented());
    return 0;
  }


#  ifndef DOXYGEN
  // explicit instantiations
  //
  template void
  SparsityPattern::copy_from(const dealii::SparsityPattern &);
  template void
  SparsityPattern::copy_from(const dealii::DynamicSparsityPattern &);

  template void
  SparsityPattern::reinit(const IndexSet &,
                          const dealii::SparsityPattern &,
                          const MPI_Comm,
                          bool);
  template void
  SparsityPattern::reinit(const IndexSet &,
                          const dealii::DynamicSparsityPattern &,
                          const MPI_Comm,
                          bool);


  template void
  SparsityPattern::reinit(const IndexSet &,
                          const IndexSet &,
                          const dealii::SparsityPattern &,
                          const MPI_Comm,
                          bool);
  template void
  SparsityPattern::reinit(const IndexSet &,
                          const IndexSet &,
                          const dealii::DynamicSparsityPattern &,
                          const MPI_Comm,
                          bool);
#  endif

} // namespace TrilinosWrappers

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_TRILINOS
