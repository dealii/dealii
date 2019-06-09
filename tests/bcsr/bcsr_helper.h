// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2017 by the deal.II authors
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


// auxiliary header to setup a 1D problem with sparse vectors

#include <deal.II/base/index_set.h>
#include <deal.II/base/mpi.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/mapping_q1.h>

#include <deal.II/lac/block_csr_matrix.h>

#include <algorithm>
#include <vector>

using namespace dealii;

/**
 * Build sparsity pattern of BCSR based on input @p local_support and
 * row and column blocking.
 */
void
setup_sparsity_pattern(const std::vector<types::global_dof_index> &row_blocks,
                       const std::vector<types::global_dof_index> &col_blocks,
                       const std::vector<IndexSet> & local_support,
                       const types::global_dof_index owned_start,
                       DynamicSparsityPattern &      dsp)
{
  AssertDimension(std::accumulate(col_blocks.begin(),
                                  col_blocks.end(),
                                  types::global_dof_index(0)),
                  local_support.size());

  dsp.reinit(row_blocks.size(), col_blocks.size());

  // go through all row blocks and then through all cols block to figure out
  // sparsity
  types::global_dof_index row_start = 0;
  for (types::global_dof_index row = 0; row < row_blocks.size(); ++row)
    {
      IndexSet row_block(local_support[0].size());
      row_block.add_range(owned_start + row_start,
                          owned_start + row_start + row_blocks[row]);
      std::vector<types::global_dof_index> sparsity;
      types::global_dof_index              col_start = 0;
      for (types::global_dof_index col = 0; col < col_blocks.size(); ++col)
        {
          const types::global_dof_index col_end = col_start + col_blocks[col];
          for (types::global_dof_index s = col_start; s < col_end; ++s)
            {
              const IndexSet intersection = local_support[s] & row_block;
              if (!intersection.is_empty())
                {
                  // if intersection is non-empty, break loop over this
                  // column block and add to the sparsity
                  sparsity.push_back(col);
                  break;
                }
            }
          col_start += col_blocks[col];
        }
      dsp.add_entries(row,
                      sparsity.begin(),
                      sparsity.end(),
                      /*indices_are_unique_and_sorted*/ true);

      row_start += row_blocks[row];
    }
}


/**
 * Gather sparsity pattern based on @p owned_rows partitioning.
 */
void
gather(DynamicSparsityPattern &      local,
       const DynamicSparsityPattern &ghost,
       const IndexSet &              owned_rows,
       const MPI_Comm                comm)
{
  const unsigned int myid = dealii::Utilities::MPI::this_mpi_process(comm);

  const auto owned =
    dealii::Utilities::MPI::all_gather(comm, owned_rows.n_elements());

  std::vector<DynamicSparsityPattern::size_type> start_index(owned.size() + 1);
  start_index[0] = 0;
  for (DynamicSparsityPattern::size_type i = 0; i < owned.size(); ++i)
    start_index[i + 1] = start_index[i] + owned[i];

  const auto b = start_index[myid];
  const auto e = start_index[myid + 1];

  local.reinit(e - b, ghost.n_cols());

  using Cols = std::vector<types::global_dof_index>;

  // what we send to owning procs:
  std::map<unsigned int, std::vector<std::pair<types::global_dof_index, Cols>>>
    to_send;

  Cols cols;
  cols.reserve(ghost.max_entries_per_row());
  // 2. go through sparsity and populate what we need to send
  unsigned int rank     = 0;
  auto         rank_end = ++start_index.begin();
  for (types::global_dof_index r = 0; r < ghost.n_rows(); ++r)
    {
      // 2.1. adjust rank if needed:
      while (r >= *rank_end)
        {
          ++rank;
          ++rank_end;
        }

      // 2.2. collect columns
      cols.resize(0);
      for (auto it = ghost.begin(r); it != ghost.end(r); ++it)
        cols.push_back(it->column());

      // 2.3. add to local sparsity or to to-be-sent buffer
      if (rank == myid)
        {
          AssertDimension(*rank_end, e);
          Assert(r >= b && r < e, ExcInternalError());
          local.add_entries(r - b, cols.begin(), cols.end(), true);
        }
      else
        {
          to_send[rank].push_back({r, cols});
        }
    }

  // 3. exchange sparsity
  const auto received = dealii::Utilities::MPI::some_to_some(comm, to_send);

  // 4. add to our sparsity based on what we received
  for (auto &el : received)
    for (auto &sp : el.second)
      {
        cols = sp.second;
        Assert(sp.first >= b && sp.first < e, ExcInternalError());
        local.add_entries(sp.first - b, cols.begin(), cols.end(), true);
      }
}



void
get_owned_columns(IndexSet &                       owned_columns,
                  std::vector<unsigned int> &      blocks_local,
                  const std::vector<unsigned int> &blocks_global,
                  const MPI_Comm                   comm)
{
  const unsigned int myid   = dealii::Utilities::MPI::this_mpi_process(comm);
  const unsigned int n_proc = dealii::Utilities::MPI::n_mpi_processes(comm);

  owned_columns.clear();
  owned_columns.set_size(blocks_global.size());

  blocks_local.resize(0);

  const auto &M     = blocks_global.size();
  const auto  block = M / n_proc;
  const auto  start = block * myid;
  const auto  end   = myid == n_proc - 1 ? M : block * (myid + 1);

  owned_columns.add_range(start, end);
  blocks_local.reserve(end - start);
  for (unsigned int b = start; b < end; ++b)
    blocks_local.push_back(blocks_global[b]);
}


/**
 * Given @p owned_size locally owned rows and MPI communicator @p comm,
 * create an index set of locally owned rows on this MPI process.
 */
IndexSet
get_owned_dofs(const types::global_dof_index owned_size, const MPI_Comm comm)
{
  const unsigned int myid = dealii::Utilities::MPI::this_mpi_process(comm);
  const auto all_owned = dealii::Utilities::MPI::all_gather(comm, owned_size);

  const auto start      = std::accumulate(all_owned.begin(),
                                     all_owned.begin() + myid,
                                     types::global_dof_index(0));
  const auto total_size = std::accumulate(all_owned.begin(),
                                          all_owned.end(),
                                          types::global_dof_index(0));

  IndexSet owned_rows(total_size);
  owned_rows.add_range(start, start + owned_size);
  return owned_rows;
}



/**
 * A generalization of the classic MPI_Bcast function, that accepts arbitrary
 * data types T, as long as boost::serialize accepts T as an argument.
 */
template <typename T>
void
bcast(T &object, const MPI_Comm &comm, const unsigned int root)
{
#ifndef DEAL_II_WITH_MPI
  (void)comm;
  (void)root;
  (void)object;
  return;
#else
  AssertIndexRange(root, dealii::Utilities::MPI::n_mpi_processes(comm));
  const auto this_proc = dealii::Utilities::MPI::this_mpi_process(comm);

  // 1. serialize on root and let others know the size
  std::vector<char> buffer;
  int               n_local_data;
  if (this_proc == root)
    {
      buffer       = dealii::Utilities::pack(object);
      n_local_data = buffer.size();
    }

  // Broadcast to others the size
  int ierr = MPI_Bcast(&n_local_data, 1, MPI_INT, root, comm);
  AssertThrowMPI(ierr);

  // resize buffer on non-root processes
  if (this_proc != root)
    buffer.resize(n_local_data);

  // 2. Bcast the actual packed object
  ierr = MPI_Bcast(buffer.data(), n_local_data, MPI_CHAR, root, comm);
  AssertThrowMPI(ierr);

  // 3. unpack the received object on non-root:
  if (this_proc != root)
    object = Utilities::unpack<T>(buffer);
#endif
}



/**
 * Broadcast sparsity pattern
 */
void
bcast_sp(DynamicSparsityPattern &sp,
         const MPI_Comm          comm,
         const unsigned int      root)
{
  std::vector<std::vector<DynamicSparsityPattern::size_type>> all_cols;
  DynamicSparsityPattern::size_type                           n_cols;
  const auto this_proc = dealii::Utilities::MPI::this_mpi_process(comm);
  if (this_proc == root)
    {
      n_cols = sp.n_cols();
      std::vector<DynamicSparsityPattern::size_type> cols;
      cols.reserve(sp.max_entries_per_row());
      all_cols.resize(sp.n_rows());
      for (DynamicSparsityPattern::size_type r = 0; r < sp.n_rows(); ++r)
        {
          cols.resize(0);
          for (auto it = sp.begin(r); it != sp.end(r); ++it)
            cols.push_back(it->column());

          all_cols[r] = cols;
        }
    }

  bcast(all_cols, comm, root);
  bcast(n_cols, comm, root);

  if (this_proc != root)
    {
      sp.reinit(all_cols.size(), n_cols);
      for (DynamicSparsityPattern::size_type r = 0; r < sp.n_rows(); ++r)
        sp.add_entries(r, all_cols[r].begin(), all_cols[r].end(), true);
    }
}


/**
 * Given column partitioning @p column_partitioning (with just one range),
 * partition this range into blocks with smallest element being @p smallest .
 * Reminder from the division will be put into the last block.
 */
std::vector<unsigned int>
get_local_col_blocks(const IndexSet &   column_partitioning,
                     const unsigned int smallest)
{
  Assert(smallest > 0, ExcInternalError());
  AssertDimension(column_partitioning.n_intervals(), 1);
  const auto b = column_partitioning.nth_index_in_set(0);
  const auto e =
    column_partitioning.nth_index_in_set(column_partitioning.n_elements() - 1) +
    1;
  const auto                size     = e - b;
  const auto                n_blocks = size / smallest;
  std::vector<unsigned int> blocks(n_blocks, smallest);
  if (n_blocks * smallest != size)
    blocks.back() += size - n_blocks * smallest;

  return blocks;
}



/**
 * Given sparsity pattern @p global_dsp, extract its
 * subset based on @p owned_col_blocks.
 */
void
get_view(DynamicSparsityPattern &      local_dsp,
         const DynamicSparsityPattern &global_dsp,
         const IndexSet &              owned_col_blocks)
{
  AssertDimension(owned_col_blocks.n_intervals(), 1);
  const auto begin = owned_col_blocks.nth_index_in_set(0);
  const auto end =
    owned_col_blocks.nth_index_in_set(owned_col_blocks.n_elements() - 1) + 1;

  local_dsp.reinit(end - begin, global_dsp.n_cols());

  unsigned int                                   local_row = 0;
  std::vector<DynamicSparsityPattern::size_type> cols;
  cols.reserve(global_dsp.max_entries_per_row());
  for (DynamicSparsityPattern::size_type r = begin; r < end; ++r, ++local_row)
    {
      cols.resize(0);
      for (auto it = global_dsp.begin(r); it != global_dsp.end(r); ++it)
        cols.push_back(it->column());

      local_dsp.add_entries(local_row, cols.begin(), cols.end(), true);
    }
}



/**
 * Return index set of columns stored in the sparsity pattern @p dsp.
 */
IndexSet
columns(const DynamicSparsityPattern &dsp)
{
  std::set<DynamicSparsityPattern::size_type> cols;
  for (auto it = dsp.begin(); it != dsp.end(); ++it)
    cols.insert(it->column());

  IndexSet columns(dsp.n_cols());
  columns.add_indices(cols.begin(), cols.end());
  return columns;
}


/**
 * Given owned column blocks, setup global column blocks and index set
 * representing owned blocks.
 */
void
setup_column_blocks(IndexSet &                       owned_col_blocks,
                    std::vector<unsigned int> &      col_blocks,
                    const std::vector<unsigned int> &col_blocks_local,
                    const MPI_Comm                   mpi_communicator)
{
  const unsigned int myid =
    dealii::Utilities::MPI::this_mpi_process(mpi_communicator);
  col_blocks.resize(0);

  const auto col_blocks_gathered =
    Utilities::MPI::all_gather(mpi_communicator, col_blocks_local);
  for (auto v : col_blocks_gathered)
    col_blocks.insert(col_blocks.end(), v.begin(), v.end());

  unsigned int begin = 0;
  for (unsigned int i = 0; i < myid; ++i)
    begin += col_blocks_gathered[i].size();

  const unsigned int end = begin + col_blocks_gathered[myid].size();

  owned_col_blocks.set_size(col_blocks.size());
  owned_col_blocks.add_range(begin, end);
}



/**
 * Given global sparsity pattern of projected matrix @p global_dsp,
 * column blocks that are owned @p owned_col_blocks and local sparsity
 * pattern of `A` in `C=A^T B`, setup column partitioner and local
 * sparsity pattern for `C`.
 */
void
setup_col_partitioner_and_sparsity(
  std::shared_ptr<dealii::Utilities::MPI::Partitioner> &col_partitioner,
  DynamicSparsityPattern &                              local_dsp,
  const DynamicSparsityPattern &                        global_dsp,
  const IndexSet &                                      owned_col_blocks,
  const DynamicSparsityPattern &                        local_dsp_A,
  const MPI_Comm                                        mpi_communicator)
{
  get_view(local_dsp, global_dsp, owned_col_blocks);
  const IndexSet ghost_columns = columns(local_dsp_A);

  AssertThrow((ghost_columns & owned_col_blocks) == owned_col_blocks,
              ExcMessage(
                "ghost_columns is not a superset of owned_col_blocks "
                " on process " +
                std::to_string(
                  dealii::Utilities::MPI::this_mpi_process(mpi_communicator))));

  col_partitioner =
    std::make_shared<dealii::Utilities::MPI::Partitioner>(owned_col_blocks,
                                                          ghost_columns,
                                                          mpi_communicator);
}



template <typename Number>
void
init_bcsr(BlockCSRMatrix<Number> &A)
{
  const auto &rb = A.get_row_blocks();
  const auto &cb = A.get_col_blocks();

  const unsigned int row_mult = static_cast<unsigned int>(
    std::pow(10, std::ceil(std::log10(static_cast<double>(A.n())))));

  for (unsigned int r = 0; r < rb->size(); ++r)
    {
      const auto end       = A.end_local(r);
      const auto row_start = rb->block_start(r);
      const auto row_size  = rb->block_size(r);
      for (auto it = A.begin_local(r); it != end; ++it)
        {
          const auto c         = it->column();
          const auto col_start = cb->block_start(c);
          const auto col_size  = cb->block_size(c);

          for (unsigned int ii = 0; ii < row_size; ++ii)
            for (unsigned int jj = 0; jj < col_size; ++jj)
              *(it->data() + BlockCSRMatrix<Number>::local_index(
                               ii, jj, row_size, col_size)) =
                (row_start + ii + 1) * row_mult + (col_start + jj + 1);
        }
    }
}

template <typename NumberType>
void
init_bcsr(BlockCSRMatrix<NumberType> &A,
          const NumberType            row_mult,
          const NumberType            col_mult,
          const NumberType            shift)
{
  const auto  local_range = A.local_range();
  const auto &rb          = A.get_row_blocks();
  const auto &cb          = A.get_col_blocks();
  for (unsigned int r = 0; r < A.n_local_row_blocks(); ++r)
    {
      const auto end       = A.end_local(r);
      const auto row_start = local_range.first + rb->block_start(r);
      const auto row_size  = rb->block_size(r);
      for (auto it = A.begin_local(r); it != end; ++it)
        {
          const auto c         = it->column();
          const auto col_start = cb->block_start(c);
          const auto col_size  = cb->block_size(c);

          for (unsigned int ii = 0; ii < row_size; ++ii)
            for (unsigned int jj = 0; jj < col_size; ++jj)
              *(it->data() + BlockCSRMatrix<double>::local_index(
                               ii, jj, row_size, col_size)) =
                row_mult * (row_start + ii + 1) +
                col_mult * (col_start + jj + 1) + shift;
        }
    }
}

/**
 * A helper function to renumber based on coordinates of node in lexicographic
 * order.
 *
 * Return row blocking for locally owned DoFs so that nodes are grouped
 * according to x coordinate.
 */
template <int dim>
std::vector<unsigned int>
renumber_based_on_nodes(DoFHandler<dim> &dh)
{
  const IndexSet &                              owned = dh.locally_owned_dofs();
  std::map<types::global_dof_index, Point<dim>> support_points;
  DoFTools::map_dofs_to_support_points(StaticMappingQ1<dim>::mapping,
                                       dh,
                                       support_points);

  // rework map into a vector
  using PAIR = std::pair<Point<dim>, types::global_dof_index>;
  std::vector<PAIR> to_sort;
  to_sort.resize(dh.n_locally_owned_dofs());
  for (auto m : support_points)
    if (owned.is_element(m.first))
      to_sort[owned.index_within_set(m.first)] = {m.second, m.first};

  std::sort(to_sort.begin(),
            to_sort.end(),
            [](const PAIR &a, const PAIR &b) -> bool {
              // return a < b
              static const double eps = 1e-8;
              for (unsigned int d = 0; d < dim; ++d)
                if (std::abs(a.first[d] - b.first[d]) > eps)
                  return a.first[d] < b.first[d];

              // if points are the same according to the eps,
              // compare based on current DoFs
              return a.second < b.second;
            });

  std::vector<types::global_dof_index> new_numbers(dh.n_locally_owned_dofs());

  unsigned int              block_size = 0;
  std::vector<unsigned int> row_blocks;
  double                    x = to_sort[0].first[0];
  for (unsigned int ind = 0; ind < to_sort.size(); ++ind, ++block_size)
    {
      const auto old = owned.index_within_set(to_sort[ind].second);
      AssertIndexRange(old, new_numbers.size());
      new_numbers[old] = owned.nth_index_in_set(ind);
      // deallog << old << " -> " << ind << std::endl;
      const auto diff = std::fabs(x - to_sort[ind].first[0]);
      if (diff > 1e-8)
        {
          x = to_sort[ind].first[0];
          row_blocks.push_back(block_size);
          block_size = 0;
        }
    }

  if (block_size > 0)
    row_blocks.push_back(block_size);

  AssertDimension(std::accumulate(row_blocks.begin(), row_blocks.end(), 0),
                  to_sort.size());

  dh.renumber_dofs(new_numbers);

  return row_blocks;
}

/**
 * Setups data structures to test sparse vector in a 1D linear FEM / FD case.
 *
 * Centers of localization domains equidistributed in [0,N)
 * as follows Ci = {C0,C0+step,C0+2step,...}, each with radius L.
 *
 * @parameter local_support for each localized vector, contains its support
 * @parameter locally_owned_dofs  locally owned DoFs from decomposing [0,N) using @p mpi_communicator
 * @parameter column_partitioning partitioning of columns based on the locaction
 * of localization center.
 * @parameter mpi_communicator MPI communicator to be used
 * @parameter L sparsity radius
 * @parameter N number of DoFs
 * @parameter C0 first localization center
 * @parameter step increment/step between localization centers
 */
void
setup_1d_sparsity(std::vector<IndexSet> &global_support,
                  IndexSet &             locally_owned_dofs,
                  IndexSet &             locally_relevant_dofs,
                  IndexSet &             column_partitioning,
                  const MPI_Comm         mpi_communicator,
                  const unsigned int     L           = 22,
                  const unsigned int     N           = 100,
                  const unsigned int     C0          = 5,
                  const unsigned int     step        = 10,
                  const unsigned int     ghost_width = 3)
{
  const unsigned int n_proc = Utilities::MPI::n_mpi_processes(mpi_communicator);
  const unsigned int this_proc =
    Utilities::MPI::this_mpi_process(mpi_communicator);
  // setup locally owned dofs:
  locally_owned_dofs.clear();
  locally_owned_dofs.set_size(N);
  const auto start = this_proc * N / n_proc;
  const auto end   = std::min((this_proc + 1) * N / n_proc, N);
  locally_owned_dofs.add_range(start, end);

  locally_relevant_dofs = locally_owned_dofs;
  // setup ghost ranges:
  if (this_proc > 0)
    {
      const auto l_width = ghost_width;
      const auto g_start = start >= l_width ? start - l_width : 0;
      locally_relevant_dofs.add_range(g_start, start);
    }
  if (this_proc < n_proc - 1)
    {
      const auto r_width = ghost_width;
      const auto g_end   = end + r_width <= N ? end + r_width : N;
      locally_relevant_dofs.add_range(end, g_end);
    }

  // now setup local support:
  global_support.clear();
  std::vector<unsigned int> owned_rows;
  unsigned int              counter = 0;
  for (unsigned int center = C0; center < N; center += step)
    {
      // support of this center:
      IndexSet           sup(N);
      const unsigned int begin = center - std::min(center, L);
      const unsigned int end   = std::min(center + L + 1, N);
      sup.add_range(begin, end);

      // intersect with locally owned
      // sup = sup & locally_owned_dofs;

      global_support.push_back(sup);

      // finally decide who owns this center:
      if (locally_owned_dofs.is_element(center))
        owned_rows.push_back(counter);

      counter++;
    }

  Assert(counter == global_support.size(), ExcInternalError());
  column_partitioning.clear();
  column_partitioning.set_size(global_support.size());
  column_partitioning.add_indices(owned_rows.begin(), owned_rows.end());
}



/**
 * Given local row blocks @p row_blocks, global column blocks @p col_blocks,
 * local support of vectors @p local_support and locally owned dofs @p locally_owned_dofs,
 * create a local sparsity of Block Compressed Sparse Row matrix.
 */
void
setup_bcsr_sparsity(const std::vector<types::global_dof_index> &row_blocks,
                    const std::vector<types::global_dof_index> &col_blocks,
                    const std::vector<IndexSet> &               local_support,
                    const IndexSet &        locally_owned_dofs,
                    DynamicSparsityPattern &dsp)
{
  const types::global_dof_index n_cols =
    std::accumulate(col_blocks.begin(),
                    col_blocks.end(),
                    types::global_dof_index(0));
  (void)n_cols;
  Assert(n_cols == local_support.size(), ExcInternalError());

  dsp.reinit(row_blocks.size(), col_blocks.size());

  // go through all row blocks and then through all cols block to figure out
  // sparsity
  types::global_dof_index row_start = locally_owned_dofs.nth_index_in_set(0);
  for (types::global_dof_index row = 0; row < row_blocks.size(); ++row)
    {
      IndexSet row_block(locally_owned_dofs.size());
      row_block.add_range(row_start, row_start + row_blocks[row]);
      row_block = (row_block & locally_owned_dofs);
      std::vector<types::global_dof_index> sparsity;
      types::global_dof_index              col_start = 0;
      for (types::global_dof_index col = 0; col < col_blocks.size(); ++col)
        {
          const types::global_dof_index col_end = col_start + col_blocks[col];
          for (types::global_dof_index s = col_start; s < col_end; ++s)
            {
              const IndexSet intersection = local_support[s] & row_block;
              if (!intersection.is_empty())
                {
                  // if intersection is non-empty, break loop over this
                  // column block and add to the sparsity
                  sparsity.push_back(col);
                  break;
                }
            }
          col_start += col_blocks[col];
        }
      dsp.add_entries(row,
                      sparsity.begin(),
                      sparsity.end(),
                      /*indices_are_unique_and_sorted*/ true);

      row_start += row_blocks[row];
    }
}



/**
 * Take local sparsity pattern @p local on each process in @p mpi_communicator
 * and gather them into the global one @p global
 */
void
gather_sparsity(DynamicSparsityPattern &      global,
                const DynamicSparsityPattern &local,
                const MPI_Comm                mpi_communicator)
{
  // both DSP and SP won't work through all_gather since they don't have
  // operator= for non-empty objects. So save local sparsity in some reasonable
  // format and add it manually
  std::vector<DynamicSparsityPattern::size_type> cols;
  std::vector<std::pair<DynamicSparsityPattern::size_type,
                        std::vector<DynamicSparsityPattern::size_type>>>
    sp_send;
  cols.reserve(local.max_entries_per_row());
  for (DynamicSparsityPattern::size_type r = 0; r < local.n_rows(); ++r)
    {
      cols.resize(0);
      for (auto it = local.begin(r); it != local.end(r); ++it)
        cols.push_back(it->column());

      if (cols.size() > 0)
        sp_send.push_back({r, cols});
    }

  global.reinit(local.n_rows(), local.n_cols());
  auto sparsities = Utilities::MPI::all_gather(mpi_communicator, sp_send);
  for (auto &sp : sparsities)
    for (auto &rows : sp)
      global.add_entries(rows.first,
                         rows.second.begin(),
                         rows.second.end(),
                         true);
}
