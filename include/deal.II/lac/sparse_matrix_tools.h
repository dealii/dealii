// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_sparse_matrix_tools_h
#define dealii_sparse_matrix_tools_h

#include <deal.II/base/config.h>

#include <deal.II/base/mpi_consensus_algorithms.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>

DEAL_II_NAMESPACE_OPEN

/**
 * A namespace to process sparse matrices.
 */
namespace SparseMatrixTools
{
  /**
   * Given a sparse matrix (@p system_matrix, @p sparsity_pattern),
   * construct a new sparse matrix (@p system_matrix_out, @p sparsity_pattern_out)
   * by restriction
   * @f[
   *  A_i = R_i A R_i^T,
   * @f]
   * where the Boolean matrix $R_i$ is defined by the entries of @p requested_is.
   *
   * The function can be called by multiple processes with different sets
   * of indices, allowing to assign each process a different $A_i$.
   *
   * Such a function is useful to implement Schwarz methods, where
   * operations of type
   * @f[
   *  u^{n} = u^{n-1} + \sum_{i} R_i^T A_i^{-1} R_i (f - A u^{n-1})
   * @f]
   * are performed to iteratively solve a system of type $Au=f$.
   *
   * @warning This is a collective call that needs to be executed by all
   * processes in the communicator of @p sparse_matrix_in.
   */
  template <typename SparseMatrixType,
            typename SparsityPatternType,
            typename SparseMatrixType2,
            typename SparsityPatternType2>
  void
  restrict_to_serial_sparse_matrix(const SparseMatrixType    &sparse_matrix_in,
                                   const SparsityPatternType &sparsity_pattern,
                                   const IndexSet            &requested_is,
                                   SparseMatrixType2         &system_matrix_out,
                                   SparsityPatternType2 &sparsity_pattern_out);

  /**
   * Similar to the above function, but taking two index sets
   * (@p index_set_0, @p index_set_1), allowing to block the matrix. This
   * is particularly useful, when dealing with vectors of type
   * parallel::distributed::Vector, where the vector is blocked according
   * to locally owned and ghost indices. As a consequence, the most
   * typical usecase will be to pass in the set of locally owned DoFs and set
   * of active or locally relevant DoFs.
   *
   * @warning This is a collective call that needs to be executed by all
   * processes in the communicator of @p sparse_matrix_in.
   */
  template <typename SparseMatrixType,
            typename SparsityPatternType,
            typename SparseMatrixType2,
            typename SparsityPatternType2>
  void
  restrict_to_serial_sparse_matrix(const SparseMatrixType    &sparse_matrix_in,
                                   const SparsityPatternType &sparsity_pattern,
                                   const IndexSet            &index_set_0,
                                   const IndexSet            &index_set_1,
                                   SparseMatrixType2         &system_matrix_out,
                                   SparsityPatternType2 &sparsity_pattern_out);

  /**
   * A restriction operation similar to the above one. However, the operation
   * is performed for each locally owned active cell individually and index sets
   * are given by their DoFs. The correct entries in the resulting vector can
   * accessed by CellAccessor::active_cell_index().
   *
   * @note In a certain sense, this is the reversion of the cell loop during
   * matrix assembly. However, doing this on a distributed matrix is not
   * trivial, since 1) rows might be owned by different processes and 2) degrees
   * of freedom might be constrained, resulting in "missing" entries in the
   * matrix.
   *
   * @warning This is a collective call that needs to be executed by all
   * processes in the communicator of @p sparse_matrix_in.
   */
  template <int dim,
            int spacedim,
            typename SparseMatrixType,
            typename SparsityPatternType,
            typename Number>
  void
  restrict_to_cells(const SparseMatrixType          &system_matrix,
                    const SparsityPatternType       &sparsity_pattern,
                    const DoFHandler<dim, spacedim> &dof_handler,
                    std::vector<FullMatrix<Number>> &blocks);

  /**
   * A restriction operation similar to the above one. However, the indices
   * of the blocks can be chosen arbitrarily. If the indices of cells are
   * given, the output is the same as of the above function. However, one
   * can also provide, e.g., indices that are also part of a halo around
   * a cell to implement element-block based overlapping Schwarz methods.
   *
   * If no indices are provided for a block, the resulting matrix of this
   * block is empty.
   *
   * @warning This is a collective call that needs to be executed by all
   * processes in the communicator of @p sparse_matrix_in.
   */
  template <typename SparseMatrixType,
            typename SparsityPatternType,
            typename Number>
  void
  restrict_to_full_matrices(
    const SparseMatrixType                                  &sparse_matrix_in,
    const SparsityPatternType                               &sparsity_pattern,
    const std::vector<std::vector<types::global_dof_index>> &indices,
    std::vector<FullMatrix<Number>>                         &blocks);


#ifndef DOXYGEN
  /*---------------------- Inline functions ---------------------------------*/

  namespace internal
  {
    template <typename T>
    using get_mpi_communicator_t =
      decltype(std::declval<const T>().get_mpi_communicator());

    template <typename T>
    constexpr bool has_get_mpi_communicator =
      dealii::internal::is_supported_operation<get_mpi_communicator_t, T>;

    template <typename T>
    using local_size_t = decltype(std::declval<const T>().local_size());

    template <typename T>
    constexpr bool has_local_size =
      dealii::internal::is_supported_operation<local_size_t, T>;

    template <typename SparseMatrixType,
              std::enable_if_t<has_get_mpi_communicator<SparseMatrixType>,
                               SparseMatrixType> * = nullptr>
    MPI_Comm
    get_mpi_communicator(const SparseMatrixType &sparse_matrix)
    {
      return sparse_matrix.get_mpi_communicator();
    }

    template <typename SparseMatrixType,
              std::enable_if_t<!has_get_mpi_communicator<SparseMatrixType>,
                               SparseMatrixType> * = nullptr>
    MPI_Comm
    get_mpi_communicator(const SparseMatrixType & /*sparse_matrix*/)
    {
      return MPI_COMM_SELF;
    }

    template <typename SparseMatrixType,
              std::enable_if_t<has_local_size<SparseMatrixType>,
                               SparseMatrixType> * = nullptr>
    unsigned int
    get_local_size(const SparseMatrixType &sparse_matrix)
    {
      return sparse_matrix.local_size();
    }

    template <typename SparseMatrixType,
              std::enable_if_t<!has_local_size<SparseMatrixType>,
                               SparseMatrixType> * = nullptr>
    unsigned int
    get_local_size(const SparseMatrixType &sparse_matrix)
    {
      AssertDimension(sparse_matrix.m(), sparse_matrix.n());

      return sparse_matrix.m();
    }

    // Helper function to extract for a distributed sparse matrix rows
    // potentially not owned by the current process.
    template <typename Number,
              typename SparseMatrixType,
              typename SparsityPatternType>
    std::vector<std::vector<std::pair<types::global_dof_index, Number>>>
    extract_remote_rows(const SparseMatrixType    &system_matrix,
                        const SparsityPatternType &sparsity_pattern,
                        const IndexSet            &locally_active_dofs,
                        const MPI_Comm             comm)
    {
      std::vector<unsigned int> dummy(locally_active_dofs.n_elements());

      const auto local_size = get_local_size(system_matrix);
      const auto [prefix_sum, total_sum] =
        Utilities::MPI::partial_and_total_sum(local_size, comm);
      IndexSet locally_owned_dofs(total_sum);
      locally_owned_dofs.add_range(prefix_sum, prefix_sum + local_size);

      using T1 = std::vector<
        std::pair<types::global_dof_index,
                  std::vector<std::pair<types::global_dof_index, Number>>>>;

      std::map<unsigned int, IndexSet> requesters;
      std::tie(std::ignore, requesters) =
        Utilities::MPI::compute_index_owner_and_requesters(locally_owned_dofs,
                                                           locally_active_dofs,
                                                           comm);

      std::vector<std::vector<std::pair<types::global_dof_index, Number>>>
        locally_relevant_matrix_entries(locally_active_dofs.n_elements());


      std::vector<unsigned int> ranks;
      ranks.reserve(requesters.size());

      for (const auto &i : requesters)
        ranks.push_back(i.first);

      std::vector<std::vector<unsigned int>> row_to_procs(
        locally_owned_dofs.n_elements());

      for (const auto &requester : requesters)
        for (const auto &index : requester.second)
          row_to_procs[locally_owned_dofs.index_within_set(index)].push_back(
            requester.first);

      std::map<unsigned int, T1> data;

      std::pair<types::global_dof_index,
                std::vector<std::pair<types::global_dof_index, Number>>>
        buffer;

      for (unsigned int i = 0; i < row_to_procs.size(); ++i)
        {
          if (row_to_procs[i].empty())
            continue;

          const auto row   = locally_owned_dofs.nth_index_in_set(i);
          auto       entry = system_matrix.begin(row);

          const unsigned int row_length = sparsity_pattern.row_length(row);

          buffer.first = row;
          buffer.second.resize(row_length);

          for (unsigned int j = 0; j < row_length; ++j, ++entry)
            buffer.second[j] = {entry->column(), entry->value()};

          for (const auto &proc : row_to_procs[i])
            data[proc].emplace_back(buffer);
        }

      dealii::Utilities::MPI::ConsensusAlgorithms::selector<T1>(
        ranks,
        [&](const unsigned int other_rank) { return data[other_rank]; },
        [&](const unsigned int &, const T1 &buffer_recv) {
          for (const auto &i : buffer_recv)
            {
              auto &dst =
                locally_relevant_matrix_entries[locally_active_dofs
                                                  .index_within_set(i.first)];
              dst = i.second;
              std::sort(dst.begin(),
                        dst.end(),
                        [](const auto &a, const auto &b) {
                          return a.first < b.first;
                        });
            }
        },
        comm);

      return locally_relevant_matrix_entries;
    }
  } // namespace internal



  template <typename SparseMatrixType,
            typename SparsityPatternType,
            typename SparseMatrixType2,
            typename SparsityPatternType2>
  void
  restrict_to_serial_sparse_matrix(const SparseMatrixType    &system_matrix,
                                   const SparsityPatternType &sparsity_pattern,
                                   const IndexSet            &index_set_0,
                                   const IndexSet            &index_set_1,
                                   SparseMatrixType2         &system_matrix_out,
                                   SparsityPatternType2 &sparsity_pattern_out)
  {
    Assert(index_set_1.size() == 0 || index_set_0.size() == index_set_1.size(),
           ExcInternalError());

    auto index_set_1_cleared = index_set_1;
    if (index_set_1.size() != 0)
      index_set_1_cleared.subtract_set(index_set_0);

    const auto index_within_set = [&index_set_0,
                                   &index_set_1_cleared](const auto n) {
      if (index_set_0.is_element(n))
        return index_set_0.index_within_set(n);
      else
        return index_set_0.n_elements() +
               index_set_1_cleared.index_within_set(n);
    };

    // 1) collect needed rows
    auto index_set_union = index_set_0;

    if (index_set_1.size() != 0)
      index_set_union.add_indices(index_set_1_cleared);

    // TODO: actually only communicate remote rows as in the case of
    // SparseMatrixTools::restrict_to_cells()
    const auto locally_relevant_matrix_entries =
      internal::extract_remote_rows<typename SparseMatrixType2::value_type>(
        system_matrix,
        sparsity_pattern,
        index_set_union,
        internal::get_mpi_communicator(system_matrix));


    // 2) create sparsity pattern
    const unsigned int n_rows = index_set_union.n_elements();
    const unsigned int n_cols = index_set_union.n_elements();
    const unsigned int entries_per_row =
      locally_relevant_matrix_entries.empty() ?
        0 :
        std::max_element(locally_relevant_matrix_entries.begin(),
                         locally_relevant_matrix_entries.end(),
                         [](const auto &a, const auto &b) {
                           return a.size() < b.size();
                         })
          ->size();

    sparsity_pattern_out.reinit(n_rows, n_cols, entries_per_row);

    std::vector<types::global_dof_index>                temp_indices;
    std::vector<typename SparseMatrixType2::value_type> temp_values;

    for (unsigned int row = 0; row < index_set_union.n_elements(); ++row)
      {
        const auto &global_row_entries = locally_relevant_matrix_entries[row];

        temp_indices.clear();
        temp_indices.reserve(global_row_entries.size());

        for (const auto &global_row_entry : global_row_entries)
          {
            const auto global_index = std::get<0>(global_row_entry);

            if (index_set_union.is_element(global_index))
              temp_indices.push_back(index_within_set(global_index));
          }

        sparsity_pattern_out.add_entries(
          index_within_set(index_set_union.nth_index_in_set(row)),
          temp_indices.begin(),
          temp_indices.end());
      }

    sparsity_pattern_out.compress();

    // 3) setup matrix
    system_matrix_out.reinit(sparsity_pattern_out);

    // 4) fill entries
    for (unsigned int row = 0; row < index_set_union.n_elements(); ++row)
      {
        const auto &global_row_entries = locally_relevant_matrix_entries[row];

        temp_indices.clear();
        temp_values.clear();

        temp_indices.reserve(global_row_entries.size());
        temp_values.reserve(global_row_entries.size());

        for (const auto &global_row_entry : global_row_entries)
          {
            const auto global_index = std::get<0>(global_row_entry);

            if (index_set_union.is_element(global_index))
              {
                temp_indices.push_back(index_within_set(global_index));
                temp_values.push_back(std::get<1>(global_row_entry));
              }
          }

        system_matrix_out.add(index_within_set(
                                index_set_union.nth_index_in_set(row)),
                              temp_indices,
                              temp_values);
      }

    system_matrix_out.compress(VectorOperation::add);
  }



  template <typename SparseMatrixType,
            typename SparsityPatternType,
            typename SparseMatrixType2,
            typename SparsityPatternType2>
  void
  restrict_to_serial_sparse_matrix(const SparseMatrixType    &system_matrix,
                                   const SparsityPatternType &sparsity_pattern,
                                   const IndexSet            &requested_is,
                                   SparseMatrixType2         &system_matrix_out,
                                   SparsityPatternType2 &sparsity_pattern_out)
  {
    restrict_to_serial_sparse_matrix(system_matrix,
                                     sparsity_pattern,
                                     requested_is,
                                     IndexSet(), // simply pass empty index set
                                     system_matrix_out,
                                     sparsity_pattern_out);
  }



  template <typename SparseMatrixType,
            typename SparsityPatternType,
            typename Number>
  void
  restrict_to_full_matrices(
    const SparseMatrixType                                  &system_matrix,
    const SparsityPatternType                               &sparsity_pattern,
    const std::vector<std::vector<types::global_dof_index>> &indices,
    std::vector<FullMatrix<Number>>                         &blocks)
  {
    // 0) determine which rows are locally owned and which ones are remote
    const auto local_size = internal::get_local_size(system_matrix);
    const auto prefix_sum = Utilities::MPI::partial_and_total_sum(
      local_size, internal::get_mpi_communicator(system_matrix));
    IndexSet locally_owned_dofs(std::get<1>(prefix_sum));
    locally_owned_dofs.add_range(std::get<0>(prefix_sum),
                                 std::get<0>(prefix_sum) + local_size);

    std::vector<dealii::types::global_dof_index> ghost_indices_vector;

    for (const auto &i : indices)
      ghost_indices_vector.insert(ghost_indices_vector.end(),
                                  i.begin(),
                                  i.end());

    std::sort(ghost_indices_vector.begin(), ghost_indices_vector.end());

    IndexSet locally_active_dofs(std::get<1>(prefix_sum));
    locally_active_dofs.add_indices(ghost_indices_vector.begin(),
                                    ghost_indices_vector.end());

    locally_active_dofs.subtract_set(locally_owned_dofs);

    // 1) collect remote rows of sparse matrix
    const auto locally_relevant_matrix_entries =
      internal::extract_remote_rows<Number>(system_matrix,
                                            sparsity_pattern,
                                            locally_active_dofs,
                                            internal::get_mpi_communicator(
                                              system_matrix));


    // 2) loop over all cells and "revert" assembly
    blocks.clear();
    blocks.resize(indices.size());

    for (unsigned int c = 0; c < indices.size(); ++c)
      {
        if (indices[c].empty())
          continue;

        const auto &local_dof_indices = indices[c];
        auto       &cell_matrix       = blocks[c];

        // allocate memory
        const unsigned int dofs_per_cell = indices[c].size();

        cell_matrix = FullMatrix<Number>(dofs_per_cell, dofs_per_cell);

        // loop over all entries of the restricted element matrix and
        // do different things if rows are locally owned or not and
        // if column entries of that row exist or not
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          for (unsigned int j = 0; j < dofs_per_cell; ++j)
            {
              if (locally_owned_dofs.is_element(
                    local_dof_indices[i])) // row is local
                {
                  if constexpr (std::is_same_v<SparseMatrixType,
                                               dealii::SparseMatrix<Number>>)
                    {
                      const types::global_dof_index ind =
                        system_matrix.get_sparsity_pattern()(
                          local_dof_indices[i], local_dof_indices[j]);

                      // If SparsityPattern::operator()` found the entry, then
                      // we can access the corresponding value without a
                      // second search in the sparse matrix, otherwise the
                      // matrix entry at that index is zero because it does
                      // not exist in the sparsity pattern
                      if (ind != SparsityPattern::invalid_entry)
                        {
                          const SparseMatrixIterators::Accessor<Number, true>
                            accessor(&system_matrix, ind);
                          cell_matrix(i, j) = accessor.value();
                        }
                      else
                        cell_matrix(i, j) = 0.0;
                    }
                  else
                    cell_matrix(i, j) =
                      sparsity_pattern.exists(local_dof_indices[i],
                                              local_dof_indices[j]) ?
                        system_matrix(local_dof_indices[i],
                                      local_dof_indices[j]) :
                        0.0;
                }
              else // row is ghost
                {
                  Assert(locally_active_dofs.is_element(local_dof_indices[i]),
                         ExcInternalError());

                  const auto &row_entries =
                    locally_relevant_matrix_entries[locally_active_dofs
                                                      .index_within_set(
                                                        local_dof_indices[i])];

                  const auto ptr =
                    std::lower_bound(row_entries.begin(),
                                     row_entries.end(),
                                     std::pair<types::global_dof_index, Number>{
                                       local_dof_indices[j], /*dummy*/ 0.0},
                                     [](const auto a, const auto b) {
                                       return a.first < b.first;
                                     });

                  if (ptr != row_entries.end() &&
                      local_dof_indices[j] == ptr->first)
                    cell_matrix(i, j) = ptr->second;
                  else
                    cell_matrix(i, j) = 0.0;
                }
            }
      }
  }



  template <int dim,
            int spacedim,
            typename SparseMatrixType,
            typename SparsityPatternType,
            typename Number>
  void
  restrict_to_cells(const SparseMatrixType          &system_matrix,
                    const SparsityPatternType       &sparsity_pattern,
                    const DoFHandler<dim, spacedim> &dof_handler,
                    std::vector<FullMatrix<Number>> &blocks)
  {
    std::vector<std::vector<types::global_dof_index>> all_dof_indices;
    all_dof_indices.resize(dof_handler.get_triangulation().n_active_cells());

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        if (cell->is_locally_owned() == false)
          continue;

        auto &local_dof_indices = all_dof_indices[cell->active_cell_index()];
        local_dof_indices.resize(cell->get_fe().n_dofs_per_cell());
        cell->get_dof_indices(local_dof_indices);
      }

    restrict_to_full_matrices(system_matrix,
                              sparsity_pattern,
                              all_dof_indices,
                              blocks);
  }
#endif

} // namespace SparseMatrixTools

DEAL_II_NAMESPACE_CLOSE

#endif
