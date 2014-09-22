// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2014 by the deal.II authors
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

#include <deal.II/lac/trilinos_sparse_matrix.h>

#ifdef DEAL_II_WITH_TRILINOS

#  include <deal.II/base/utilities.h>
#  include <deal.II/lac/sparse_matrix.h>
#  include <deal.II/lac/trilinos_sparsity_pattern.h>
#  include <deal.II/lac/sparsity_pattern.h>
#  include <deal.II/lac/compressed_sparsity_pattern.h>
#  include <deal.II/lac/compressed_set_sparsity_pattern.h>
#  include <deal.II/lac/compressed_simple_sparsity_pattern.h>
#  include <deal.II/lac/sparsity_tools.h>

#  include <Epetra_Export.h>
#  include <ml_epetra_utils.h>
#  include <ml_struct.h>
#  include <Teuchos_RCP.hpp>

DEAL_II_NAMESPACE_OPEN

namespace TrilinosWrappers
{
  namespace
  {
#ifndef DEAL_II_WITH_64BIT_INDICES
    // define a helper function that queries the size of an Epetra_Map object
    // by calling either the 32- or 64-bit function necessary, and returns the
    // result in the correct data type so that we can use it in calling other
    // Epetra member functions that are overloaded by index type
    int n_global_elements (const Epetra_BlockMap &map)
    {
      return map.NumGlobalElements();
    }

    int min_my_gid(const Epetra_BlockMap &map)
    {
      return map.MinMyGID();
    }

    int max_my_gid(const Epetra_BlockMap &map)
    {
      return map.MaxMyGID();
    }

    int n_global_cols(const Epetra_CrsGraph &graph)
    {
      return graph.NumGlobalCols();
    }

    int global_column_index(const Epetra_CrsMatrix &matrix, int i)
    {
      return matrix.GCID(i);
    }

    int global_row_index(const Epetra_CrsMatrix &matrix, int i)
    {
      return matrix.GRID(i);
    }
#else
    // define a helper function that queries the size of an Epetra_Map object
    // by calling either the 32- or 64-bit function necessary, and returns the
    // result in the correct data type so that we can use it in calling other
    // Epetra member functions that are overloaded by index type
    long long int n_global_elements (const Epetra_BlockMap &map)
    {
      return map.NumGlobalElements64();
    }

    long long int min_my_gid(const Epetra_BlockMap &map)
    {
      return map.MinMyGID64();
    }

    long long int max_my_gid(const Epetra_BlockMap &map)
    {
      return map.MaxMyGID64();
    }

    long long int n_global_cols(const Epetra_CrsGraph &graph)
    {
      return graph.NumGlobalCols64();
    }

    long long int global_column_index(const Epetra_CrsMatrix &matrix, int i)
    {
      return matrix.GCID64(i);
    }

    long long int global_row_index(const Epetra_CrsMatrix &matrix, int i)
    {
      return matrix.GRID64(i);
    }
#endif
  }


  namespace SparseMatrixIterators
  {
    void
    AccessorBase::visit_present_row ()
    {
      // if we are asked to visit the
      // past-the-end line, then simply
      // release all our caches and go on
      // with life
      if (this->a_row == matrix->m())
        {
          colnum_cache.reset ();
          value_cache.reset ();

          return;
        }

      // otherwise first flush Trilinos caches
      matrix->compress ();

      // get a representation of the present
      // row
      int ncols;
      TrilinosWrappers::types::int_type colnums = matrix->n();
      if (value_cache.get() == 0)
        {
          value_cache.reset (new std::vector<TrilinosScalar> (matrix->n()));
          colnum_cache.reset (new std::vector<size_type> (matrix->n()));
        }
      else
        {
          value_cache->resize (matrix->n());
          colnum_cache->resize (matrix->n());
        }

      int ierr = matrix->trilinos_matrix().
                 ExtractGlobalRowCopy((TrilinosWrappers::types::int_type)this->a_row,
                                      colnums,
                                      ncols, &((*value_cache)[0]),
                                      reinterpret_cast<TrilinosWrappers::types::int_type *>(&((*colnum_cache)[0])));
      value_cache->resize (ncols);
      colnum_cache->resize (ncols);
      AssertThrow (ierr == 0, ExcTrilinosError(ierr));

      // copy it into our caches if the
      // line isn't empty. if it is, then
      // we've done something wrong, since
      // we shouldn't have initialized an
      // iterator for an empty line (what
      // would it point to?)
    }
  }


  // The constructor is actually the
  // only point where we have to check
  // whether we build a serial or a
  // parallel Trilinos matrix.
  // Actually, it does not even matter
  // how many threads there are, but
  // only if we use an MPI compiler or
  // a standard compiler. So, even one
  // thread on a configuration with
  // MPI will still get a parallel
  // interface.
  SparseMatrix::SparseMatrix ()
    :
    column_space_map (new Epetra_Map (0, 0,
                                      Utilities::Trilinos::comm_self())),
    matrix (new Epetra_FECrsMatrix(View, *column_space_map,
                                   *column_space_map, 0)),
    last_action (Zero),
    compressed (true)
  {
    matrix->FillComplete();
  }



  SparseMatrix::SparseMatrix (const Epetra_Map  &input_map,
                              const size_type n_max_entries_per_row)
    :
    column_space_map (new Epetra_Map (input_map)),
    matrix (new Epetra_FECrsMatrix(Copy, *column_space_map,
                                   TrilinosWrappers::types::int_type(n_max_entries_per_row), false)),
    last_action (Zero),
    compressed (false)
  {}



  SparseMatrix::SparseMatrix (const Epetra_Map                &input_map,
                              const std::vector<unsigned int> &n_entries_per_row)
    :
    column_space_map (new Epetra_Map (input_map)),
    matrix (new Epetra_FECrsMatrix
            (Copy, *column_space_map,
             (int *)const_cast<unsigned int *>(&(n_entries_per_row[0])),
             false)),
    last_action (Zero),
    compressed (false)
  {}



  SparseMatrix::SparseMatrix (const Epetra_Map  &input_row_map,
                              const Epetra_Map  &input_col_map,
                              const size_type n_max_entries_per_row)
    :
    column_space_map (new Epetra_Map (input_col_map)),
    matrix (new Epetra_FECrsMatrix(Copy, input_row_map,
                                   TrilinosWrappers::types::int_type(n_max_entries_per_row), false)),
    last_action (Zero),
    compressed (false)
  {}



  SparseMatrix::SparseMatrix (const Epetra_Map                &input_row_map,
                              const Epetra_Map                &input_col_map,
                              const std::vector<unsigned int> &n_entries_per_row)
    :
    column_space_map (new Epetra_Map (input_col_map)),
    matrix (new Epetra_FECrsMatrix(Copy, input_row_map,
                                   (int *)const_cast<unsigned int *>(&(n_entries_per_row[0])),
                                   false)),
    last_action (Zero),
    compressed (false)
  {}



  SparseMatrix::SparseMatrix (const size_type m,
                              const size_type n,
                              const unsigned int n_max_entries_per_row)
    :
    column_space_map (new Epetra_Map (static_cast<TrilinosWrappers::types::int_type>(n), 0,
                                      Utilities::Trilinos::comm_self())),

    // on one processor only, we know how the
    // columns of the matrix will be
    // distributed (everything on one
    // processor), so we can hand in this
    // information to the constructor. we
    // can't do so in parallel, where the
    // information from columns is only
    // available when entries have been added
    matrix (new Epetra_FECrsMatrix(Copy,
                                   Epetra_Map (static_cast<TrilinosWrappers::types::int_type>(m), 0,
                                               Utilities::Trilinos::comm_self()),
                                   *column_space_map,
                                   n_max_entries_per_row,
                                   false)),
    last_action (Zero),
    compressed (false)
  {}



  SparseMatrix::SparseMatrix (const size_type                  m,
                              const size_type                  n,
                              const std::vector<unsigned int> &n_entries_per_row)
    :
    column_space_map (new Epetra_Map (static_cast<TrilinosWrappers::types::int_type>(n), 0,
                                      Utilities::Trilinos::comm_self())),
    matrix (new Epetra_FECrsMatrix(Copy,
                                   Epetra_Map (static_cast<TrilinosWrappers::types::int_type>(m), 0,
                                               Utilities::Trilinos::comm_self()),
                                   *column_space_map,
                                   (int *)const_cast<unsigned int *>(&(n_entries_per_row[0])),
                                   false)),
    last_action (Zero),
    compressed (false)
  {}



  SparseMatrix::SparseMatrix (const IndexSet     &parallel_partitioning,
                              const MPI_Comm     &communicator,
                              const unsigned int  n_max_entries_per_row)
    :
    column_space_map (new Epetra_Map(parallel_partitioning.
                                     make_trilinos_map(communicator, false))),
    matrix (new Epetra_FECrsMatrix(Copy,
                                   *column_space_map,
                                   n_max_entries_per_row,
                                   false)),
    last_action (Zero),
    compressed (false)
  {}



  SparseMatrix::SparseMatrix (const IndexSet                  &parallel_partitioning,
                              const MPI_Comm                  &communicator,
                              const std::vector<unsigned int> &n_entries_per_row)
    :
    column_space_map (new Epetra_Map(parallel_partitioning.
                                     make_trilinos_map(communicator, false))),
    matrix (new Epetra_FECrsMatrix(Copy,
                                   *column_space_map,
                                   (int *)const_cast<unsigned int *>(&(n_entries_per_row[0])),
                                   false)),
    last_action (Zero),
    compressed (false)
  {}



  SparseMatrix::SparseMatrix (const IndexSet  &row_parallel_partitioning,
                              const IndexSet  &col_parallel_partitioning,
                              const MPI_Comm  &communicator,
                              const size_type  n_max_entries_per_row)
    :
    column_space_map (new Epetra_Map(col_parallel_partitioning.
                                     make_trilinos_map(communicator, false))),
    matrix (new Epetra_FECrsMatrix(Copy,
                                   row_parallel_partitioning.
                                   make_trilinos_map(communicator, false),
                                   n_max_entries_per_row,
                                   false)),
    last_action (Zero),
    compressed (false)
  {}



  SparseMatrix::SparseMatrix (const IndexSet                  &row_parallel_partitioning,
                              const IndexSet                  &col_parallel_partitioning,
                              const MPI_Comm                  &communicator,
                              const std::vector<unsigned int> &n_entries_per_row)
    :
    column_space_map (new Epetra_Map(col_parallel_partitioning.
                                     make_trilinos_map(communicator, false))),
    matrix (new Epetra_FECrsMatrix(Copy,
                                   row_parallel_partitioning.
                                   make_trilinos_map(communicator, false),
                                   (int *)const_cast<unsigned int *>(&(n_entries_per_row[0])),
                                   false)),
    last_action (Zero),
    compressed (false)
  {}



  SparseMatrix::SparseMatrix (const SparsityPattern &sparsity_pattern)
    :
    column_space_map (new Epetra_Map (sparsity_pattern.domain_partitioner())),
    matrix (new Epetra_FECrsMatrix(Copy,
                                   sparsity_pattern.trilinos_sparsity_pattern(),
                                   false)),
    last_action (Zero),
    compressed (true)
  {
    Assert(sparsity_pattern.trilinos_sparsity_pattern().Filled() == true,
           ExcMessage("The Trilinos sparsity pattern has not been compressed."));
    compress();
  }



  SparseMatrix::SparseMatrix (const SparseMatrix &input_matrix)
    :
    Subscriptor(),
    column_space_map (new Epetra_Map (input_matrix.domain_partitioner())),
    matrix (new Epetra_FECrsMatrix(*input_matrix.matrix)),
    last_action (Zero),
    compressed (true)
  {}



  SparseMatrix::~SparseMatrix ()
  {}



  void
  SparseMatrix::copy_from (const SparseMatrix &m)
  {
    nonlocal_matrix.reset();
    nonlocal_matrix_exporter.reset();

    // check whether we need to update the
    // partitioner or can just copy the data:
    // in case we have the same distribution,
    // we can just copy the data.
    if (local_range() == m.local_range())
      *matrix = *m.matrix;
    else
      {
        column_space_map.reset (new Epetra_Map (m.domain_partitioner()));

        // release memory before reallocation
        matrix.reset ();
        matrix.reset (new Epetra_FECrsMatrix(*m.matrix));
      }

    if (m.nonlocal_matrix.get() != 0)
      nonlocal_matrix.reset(new Epetra_CrsMatrix(Copy, m.nonlocal_matrix->Graph()));

    compress();
  }



  template <typename SparsityType>
  void
  SparseMatrix::reinit (const SparsityType &sparsity_pattern)
  {
    const Epetra_Map rows (static_cast<TrilinosWrappers::types::int_type>(sparsity_pattern.n_rows()),
                           0,
                           Utilities::Trilinos::comm_self());
    const Epetra_Map columns (static_cast<TrilinosWrappers::types::int_type>(sparsity_pattern.n_cols()),
                              0,
                              Utilities::Trilinos::comm_self());

    reinit (rows, columns, sparsity_pattern);
  }



  template <typename SparsityType>
  void
  SparseMatrix::reinit (const Epetra_Map    &input_map,
                        const SparsityType  &sparsity_pattern,
                        const bool           exchange_data)
  {
    reinit (input_map, input_map, sparsity_pattern, exchange_data);
  }



  namespace internal
  {
    namespace
    {
      // distinguish between compressed sparsity types that define row_begin()
      // and SparsityPattern that uses begin() as iterator type
      template <typename Sparsity>
      void copy_row (const Sparsity        &csp,
                     const size_type        row,
                     std::vector<TrilinosWrappers::types::int_type> &row_indices)
      {
        typename Sparsity::row_iterator col_num = csp.row_begin (row);
        for (size_type col=0; col_num != csp.row_end (row); ++col_num, ++col)
          row_indices[col] = *col_num;
      }

      void copy_row (const dealii::SparsityPattern &csp,
                     const size_type                row,
                     std::vector<TrilinosWrappers::types::int_type> &row_indices)
      {
        dealii::SparsityPattern::iterator col_num = csp.begin (row);
        for (size_type col=0; col_num != csp.end (row); ++col_num, ++col)
          row_indices[col] = col_num->column();
      }
    }
  }



  template <typename SparsityType>
  void
  SparseMatrix::reinit (const Epetra_Map    &input_row_map,
                        const Epetra_Map    &input_col_map,
                        const SparsityType  &sparsity_pattern,
                        const bool           exchange_data)
  {
    // release memory before reallocation
    matrix.reset();
    nonlocal_matrix.reset();
    nonlocal_matrix_exporter.reset();

    // if we want to exchange data, build a usual Trilinos sparsity pattern
    // and let that handle the exchange. otherwise, manually create a
    // CrsGraph, which consumes considerably less memory because it can set
    // correct number of indices right from the start
    if (exchange_data)
      {
        SparsityPattern trilinos_sparsity;
        trilinos_sparsity.reinit (input_row_map, input_col_map,
                                  sparsity_pattern, exchange_data);
        reinit (trilinos_sparsity);

        return;
      }

    Assert (exchange_data == false, ExcNotImplemented());
    if (input_row_map.Comm().MyPID() == 0)
      {
        AssertDimension (sparsity_pattern.n_rows(),
                         static_cast<size_type>(n_global_elements(input_row_map)));
        AssertDimension (sparsity_pattern.n_cols(),
                         static_cast<size_type>(n_global_elements(input_col_map)));
      }

    column_space_map.reset (new Epetra_Map (input_col_map));

    const size_type first_row = min_my_gid(input_row_map),
                    last_row = max_my_gid(input_row_map)+1;
    std::vector<int> n_entries_per_row(last_row-first_row);

    for (size_type row=first_row; row<last_row; ++row)
      n_entries_per_row[row-first_row] = sparsity_pattern.row_length(row);

    // The deal.II notation of a Sparsity pattern corresponds to the Epetra
    // concept of a Graph. Hence, we generate a graph by copying the sparsity
    // pattern into it, and then build up the matrix from the graph. This is
    // considerable faster than directly filling elements into the
    // matrix. Moreover, it consumes less memory, since the internal
    // reordering is done on ints only, and we can leave the doubles aside.

    // for more than one processor, need to specify only row map first and let
    // the matrix entries decide about the column map (which says which
    // columns are present in the matrix, not to be confused with the col_map
    // that tells how the domain dofs of the matrix will be distributed). for
    // only one processor, we can directly assign the columns as well. Compare
    // this with bug # 4123 in the Sandia Bugzilla.
    std_cxx11::shared_ptr<Epetra_CrsGraph> graph;
    if (input_row_map.Comm().NumProc() > 1)
      graph.reset (new Epetra_CrsGraph (Copy, input_row_map,
                                        &n_entries_per_row[0], true));
    else
      graph.reset (new Epetra_CrsGraph (Copy, input_row_map, input_col_map,
                                        &n_entries_per_row[0], true));

    // This functions assumes that the sparsity pattern sits on all processors
    // (completely). The parallel version uses an Epetra graph that is already
    // distributed.

    // now insert the indices
    std::vector<TrilinosWrappers::types::int_type>   row_indices;

    for (size_type row=first_row; row<last_row; ++row)
      {
        const int row_length = sparsity_pattern.row_length(row);
        if (row_length == 0)
          continue;

        row_indices.resize (row_length, -1);
        internal::copy_row(sparsity_pattern, row, row_indices);
        graph->Epetra_CrsGraph::InsertGlobalIndices (row, row_length,
                                                     &row_indices[0]);
      }

    // Eventually, optimize the graph structure (sort indices, make memory
    // contiguous, etc). note that the documentation of the function indeed
    // states that we first need to provide the column (domain) map and then
    // the row (range) map
    graph->FillComplete(input_col_map, input_row_map);
    graph->OptimizeStorage();

    // check whether we got the number of columns right.
    AssertDimension (sparsity_pattern.n_cols(),static_cast<size_type>(
                       n_global_cols(*graph)));

    // And now finally generate the matrix.
    matrix.reset (new Epetra_FECrsMatrix(Copy, *graph, false));
    last_action = Zero;

    // In the end, the matrix needs to be compressed in order to be really
    // ready.
    compress();
  }



  // specialization for CompressedSimpleSparsityPattern which can provide us
  // with more information about the non-locally owned rows
  template <>
  void
  SparseMatrix::reinit (const Epetra_Map    &input_row_map,
                        const Epetra_Map    &input_col_map,
                        const CompressedSimpleSparsityPattern &sparsity_pattern,
                        const bool           exchange_data)
  {
    matrix.reset();
    nonlocal_matrix.reset();
    nonlocal_matrix_exporter.reset();

    AssertDimension (sparsity_pattern.n_rows(),
                     static_cast<size_type>(n_global_elements(input_row_map)));
    AssertDimension (sparsity_pattern.n_cols(),
                     static_cast<size_type>(n_global_elements(input_col_map)));

    column_space_map.reset (new Epetra_Map (input_col_map));

    IndexSet relevant_rows (sparsity_pattern.row_index_set());
    // serial case
    if (relevant_rows.size() == 0)
      {
        relevant_rows.set_size(n_global_elements(input_row_map));
        relevant_rows.add_range(0, n_global_elements(input_row_map));
      }
    relevant_rows.compress();
    Assert(relevant_rows.n_elements() >= static_cast<unsigned int>(input_row_map.NumMyElements()),
           ExcMessage("Locally relevant rows of sparsity pattern must contain "
                      "all locally owned rows"));

    // check whether the relevant rows correspond to exactly the same map as
    // the owned rows. In that case, do not create the nonlocal graph and fill
    // the columns by demand
    bool have_ghost_rows = false;
    {
      std::vector<dealii::types::global_dof_index> indices;
      relevant_rows.fill_index_vector(indices);
      Epetra_Map relevant_map (TrilinosWrappers::types::int_type(-1),
                               TrilinosWrappers::types::int_type(relevant_rows.n_elements()),
                               (indices.empty() ? 0 :
                                reinterpret_cast<TrilinosWrappers::types::int_type *>(&indices[0])),
                               0, input_row_map.Comm());
      if (relevant_map.SameAs(input_row_map))
        have_ghost_rows = false;
      else
        have_ghost_rows = true;
    }

    const unsigned int n_rows = relevant_rows.n_elements();
    std::vector<TrilinosWrappers::types::int_type> ghost_rows;
    std::vector<int> n_entries_per_row(input_row_map.NumMyElements());
    std::vector<int> n_entries_per_ghost_row;
    for (unsigned int i=0, own=0; i<n_rows; ++i)
      {
        const TrilinosWrappers::types::int_type global_row =
          relevant_rows.nth_index_in_set(i);
        if (input_row_map.MyGID(global_row))
          n_entries_per_row[own++] = sparsity_pattern.row_length(global_row);
        else if (sparsity_pattern.row_length(global_row) > 0)
          {
            ghost_rows.push_back(global_row);
            n_entries_per_ghost_row.push_back(sparsity_pattern.row_length(global_row));
          }
      }

    // make sure all processors create an off-processor matrix with at least
    // one entry
    if (have_ghost_rows == true && ghost_rows.empty() == true)
      {
        ghost_rows.push_back(0);
        n_entries_per_ghost_row.push_back(1);
      }

    Epetra_Map off_processor_map(-1, ghost_rows.size(), &ghost_rows[0],
                                 0, input_row_map.Comm());

    std_cxx11::shared_ptr<Epetra_CrsGraph> graph, nonlocal_graph;
    if (input_row_map.Comm().NumProc() > 1)
      {
        graph.reset (new Epetra_CrsGraph (Copy, input_row_map,
                                          &n_entries_per_row[0],
                                          exchange_data ? false : true));
        if (have_ghost_rows == true)
          nonlocal_graph.reset (new Epetra_CrsGraph (Copy, off_processor_map,
                                                     &n_entries_per_ghost_row[0],
                                                     true));
      }
    else
      graph.reset (new Epetra_CrsGraph (Copy, input_row_map, input_col_map,
                                        &n_entries_per_row[0], true));

    // now insert the indices, select between the right matrix
    std::vector<TrilinosWrappers::types::int_type> row_indices;

    for (unsigned int i=0; i<n_rows; ++i)
      {
        const TrilinosWrappers::types::int_type global_row =
          relevant_rows.nth_index_in_set(i);
        const int row_length = sparsity_pattern.row_length(global_row);
        if (row_length == 0)
          continue;

        row_indices.resize (row_length, -1);
        internal::copy_row(sparsity_pattern, global_row, row_indices);

        if (input_row_map.MyGID(global_row))
          graph->InsertGlobalIndices (global_row, row_length, &row_indices[0]);
        else
          {
            Assert(nonlocal_graph.get() != 0, ExcInternalError());
            nonlocal_graph->InsertGlobalIndices (global_row, row_length,
                                                 &row_indices[0]);
          }
      }

    // finalize nonlocal graph and create nonlocal matrix
    if (nonlocal_graph.get() != 0)
      {
        if (nonlocal_graph->IndicesAreGlobal() == false &&
            nonlocal_graph->RowMap().NumMyElements() > 0)
          {
            // insert dummy element
            TrilinosWrappers::types::int_type row =
              nonlocal_graph->RowMap().MyGID(TrilinosWrappers::types::int_type(0));
            nonlocal_graph->InsertGlobalIndices(row, 1, &row);
          }
        Assert(nonlocal_graph->IndicesAreGlobal() == true,
               ExcInternalError());
        nonlocal_graph->FillComplete(input_col_map, input_row_map);
        nonlocal_graph->OptimizeStorage();

        // insert data from nonlocal graph into the final sparsity pattern
        if (exchange_data)
          {
            Epetra_Export exporter(nonlocal_graph->RowMap(), input_row_map);
            int ierr = graph->Export(*nonlocal_graph, exporter, Add);
            Assert (ierr==0, ExcTrilinosError(ierr));
          }

        nonlocal_matrix.reset (new Epetra_CrsMatrix(Copy, *nonlocal_graph));
      }

    graph->FillComplete(input_col_map, input_row_map);
    graph->OptimizeStorage();

    AssertDimension (sparsity_pattern.n_cols(),static_cast<size_type>(
                       n_global_cols(*graph)));

    matrix.reset (new Epetra_FECrsMatrix(Copy, *graph, false));
    last_action = Zero;

    // In the end, the matrix needs to be compressed in order to be really
    // ready.
    compress();
  }




  void
  SparseMatrix::reinit (const SparsityPattern &sparsity_pattern)
  {
    matrix.reset ();
    nonlocal_matrix_exporter.reset();

    // reinit with a (parallel) Trilinos sparsity pattern.
    column_space_map.reset (new Epetra_Map
                            (sparsity_pattern.domain_partitioner()));
    matrix.reset (new Epetra_FECrsMatrix
                  (Copy, sparsity_pattern.trilinos_sparsity_pattern(), false));

    if (sparsity_pattern.nonlocal_graph.get() != 0)
      nonlocal_matrix.reset (new Epetra_CrsMatrix(Copy, *sparsity_pattern.nonlocal_graph));
    else
      nonlocal_matrix.reset ();

    compress();
    last_action = Zero;
  }



  void
  SparseMatrix::reinit (const SparseMatrix &sparse_matrix)
  {
    column_space_map.reset (new Epetra_Map (sparse_matrix.domain_partitioner()));
    matrix.reset ();
    nonlocal_matrix_exporter.reset();
    matrix.reset (new Epetra_FECrsMatrix
                  (Copy, sparse_matrix.trilinos_sparsity_pattern(), false));
    if (sparse_matrix.nonlocal_matrix != 0)
      nonlocal_matrix.reset (new Epetra_CrsMatrix
                             (Copy, sparse_matrix.nonlocal_matrix->Graph()));
    else
      nonlocal_matrix.reset();

    compress();
  }



  template <typename number>
  void
  SparseMatrix::reinit (const ::dealii::SparseMatrix<number> &dealii_sparse_matrix,
                        const double                          drop_tolerance,
                        const bool                            copy_values,
                        const ::dealii::SparsityPattern      *use_this_sparsity)
  {
    const Epetra_Map rows (static_cast<TrilinosWrappers::types::int_type>(dealii_sparse_matrix.m()),
                           0,
                           Utilities::Trilinos::comm_self());
    const Epetra_Map columns (static_cast<TrilinosWrappers::types::int_type>(dealii_sparse_matrix.n()),
                              0,
                              Utilities::Trilinos::comm_self());
    reinit (rows, columns, dealii_sparse_matrix, drop_tolerance,
            copy_values, use_this_sparsity);
  }



  template <typename number>
  void
  SparseMatrix::reinit (const Epetra_Map                     &input_map,
                        const ::dealii::SparseMatrix<number> &dealii_sparse_matrix,
                        const double                          drop_tolerance,
                        const bool                            copy_values,
                        const ::dealii::SparsityPattern      *use_this_sparsity)
  {
    reinit (input_map, input_map, dealii_sparse_matrix, drop_tolerance,
            copy_values, use_this_sparsity);
  }



  template <typename number>
  void
  SparseMatrix::reinit (const Epetra_Map                     &input_row_map,
                        const Epetra_Map                     &input_col_map,
                        const ::dealii::SparseMatrix<number> &dealii_sparse_matrix,
                        const double                          drop_tolerance,
                        const bool                            copy_values,
                        const ::dealii::SparsityPattern      *use_this_sparsity)
  {
    if (copy_values == false)
      {
        // in case we do not copy values, just
        // call the other function.
        if (use_this_sparsity == 0)
          reinit (input_row_map, input_col_map,
                  dealii_sparse_matrix.get_sparsity_pattern());
        else
          reinit (input_row_map, input_col_map,
                  *use_this_sparsity);
        return;
      }

    const size_type n_rows = dealii_sparse_matrix.m();

    Assert (static_cast<size_type>(n_global_elements(input_row_map)) == n_rows,
            ExcDimensionMismatch (n_global_elements(input_row_map),
                                  n_rows));
    Assert (n_global_elements(input_row_map) == (TrilinosWrappers::types::int_type)n_rows,
            ExcDimensionMismatch (n_global_elements(input_row_map),
                                  n_rows));
    Assert (n_global_elements(input_col_map) == (TrilinosWrappers::types::int_type)dealii_sparse_matrix.n(),
            ExcDimensionMismatch (n_global_elements(input_col_map),
                                  dealii_sparse_matrix.n()));

    const ::dealii::SparsityPattern &sparsity_pattern =
      (use_this_sparsity!=0)? *use_this_sparsity :
      dealii_sparse_matrix.get_sparsity_pattern();

    if (matrix.get() == 0 ||
        m() != n_rows ||
        n_nonzero_elements() != sparsity_pattern.n_nonzero_elements())
      {
        SparsityPattern trilinos_sparsity;
        trilinos_sparsity.reinit (input_row_map, input_col_map, sparsity_pattern);
        reinit (trilinos_sparsity);
      }

    // fill the values. the same as above: go through all rows of the matrix,
    // and then all columns. since the sparsity patterns of the input matrix
    // and the specified sparsity pattern might be different, need to go
    // through the row for both these sparsity structures simultaneously in
    // order to really set the correct values.
    size_type maximum_row_length = matrix->MaxNumEntries();
    std::vector<size_type> row_indices (maximum_row_length);
    std::vector<TrilinosScalar> values (maximum_row_length);

    for (size_type row=0; row<n_rows; ++row)
      // see if the row is locally stored on this processor
      if (input_row_map.MyGID(static_cast<TrilinosWrappers::types::int_type>(row)) == true)
        {
          ::dealii::SparsityPattern::iterator select_index =
            sparsity_pattern.begin(row);
          typename ::dealii::SparseMatrix<number>::const_iterator it =
            dealii_sparse_matrix.begin(row);
          size_type col = 0;
          if (sparsity_pattern.n_rows() == sparsity_pattern.n_cols())
            {
              // optimized diagonal
              AssertDimension(it->column(), row);
              if (std::fabs(it->value()) > drop_tolerance)
                {
                  values[col] = it->value();
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
                  values[col] = it->value();
                  row_indices[col++] = it->column();
                }
              ++select_index;
              ++it;
            }
          set (row, col, reinterpret_cast<size_type *>(&row_indices[0]),
               &values[0], false);
        }

    compress();
  }



  void
  SparseMatrix::reinit (const Epetra_CrsMatrix &input_matrix,
                        const bool              copy_values)
  {
    Assert (input_matrix.Filled()==true,
            ExcMessage("Input CrsMatrix has not called FillComplete()!"));

    column_space_map.reset (new Epetra_Map (input_matrix.DomainMap()));

    const Epetra_CrsGraph *graph = &input_matrix.Graph();

    nonlocal_matrix.reset();
    nonlocal_matrix_exporter.reset();
    matrix.reset ();
    matrix.reset (new Epetra_FECrsMatrix(Copy, *graph, false));

    matrix->FillComplete (*column_space_map, input_matrix.RangeMap(), true);

    if (copy_values == true)
      {
        // point to the first data entry in the two
        // matrices and copy the content
        const TrilinosScalar *in_values = input_matrix[0];
        TrilinosScalar *values = (*matrix)[0];
        const size_type my_nonzeros = input_matrix.NumMyNonzeros();
        std::memcpy (&values[0], &in_values[0],
                     my_nonzeros*sizeof (TrilinosScalar));
      }

    compress();
  }



  void
  SparseMatrix::compress (::dealii::VectorOperation::values operation)
  {

    Epetra_CombineMode mode = last_action;
    if (last_action == Zero)
      {
        if ((operation==::dealii::VectorOperation::add) ||
            (operation==::dealii::VectorOperation::unknown))
          mode = Add;
        else if (operation==::dealii::VectorOperation::insert)
          mode = Insert;
      }
    else
      {
        Assert(
          ((last_action == Add) && (operation!=::dealii::VectorOperation::insert))
          ||
          ((last_action == Insert) && (operation!=::dealii::VectorOperation::add)),
          ExcMessage("operation and argument to compress() do not match"));
      }

    // flush buffers
    int ierr;
    if (nonlocal_matrix.get() != 0 && mode == Add)
      {
        // do only export in case of an add() operation, otherwise the owning
        // processor must have set the correct entry
        nonlocal_matrix->FillComplete(*column_space_map, matrix->RowMap());
        if (nonlocal_matrix_exporter.get() == 0)
          nonlocal_matrix_exporter.reset
          (new Epetra_Export(nonlocal_matrix->RowMap(), matrix->RowMap()));
        ierr = matrix->Export(*nonlocal_matrix, *nonlocal_matrix_exporter, mode);
        AssertThrow(ierr == 0, ExcTrilinosError(ierr));
        ierr = matrix->FillComplete(*column_space_map, matrix->RowMap());
        nonlocal_matrix->PutScalar(0);
      }
    else
      ierr = matrix->GlobalAssemble (*column_space_map, matrix->RowMap(),
                                     true, mode);

    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    ierr = matrix->OptimizeStorage ();
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    last_action = Zero;

    compressed = true;
  }



  void
  SparseMatrix::clear ()
  {
    // When we clear the matrix, reset
    // the pointer and generate an
    // empty matrix.
    column_space_map.reset (new Epetra_Map (0, 0,
                                            Utilities::Trilinos::comm_self()));
    matrix.reset (new Epetra_FECrsMatrix(View, *column_space_map, 0));
    nonlocal_matrix.reset();
    nonlocal_matrix_exporter.reset();

    matrix->FillComplete();

    compressed = true;
  }



  void
  SparseMatrix::clear_row (const size_type      row,
                           const TrilinosScalar new_diag_value)
  {
    Assert (matrix->Filled()==true, ExcMatrixNotCompressed());

    // Only do this on the rows owned
    // locally on this processor.
    int local_row =
      matrix->LRID(static_cast<TrilinosWrappers::types::int_type>(row));
    if (local_row >= 0)
      {
        TrilinosScalar *values;
        int *col_indices;
        int num_entries;
        const int ierr = matrix->ExtractMyRowView(local_row, num_entries,
                                                  values, col_indices);

        Assert (ierr == 0,
                ExcTrilinosError(ierr));

        int *diag_find = std::find(col_indices,col_indices+num_entries,local_row);
        int diag_index = (int)(diag_find - col_indices);

        for (TrilinosWrappers::types::int_type j=0; j<num_entries; ++j)
          if (diag_index != j || new_diag_value == 0)
            values[j] = 0.;

        if (diag_find && std::fabs(values[diag_index]) == 0.0 &&
            new_diag_value != 0.0)
          values[diag_index] = new_diag_value;
      }
  }



  void
  SparseMatrix::clear_rows (const std::vector<size_type> &rows,
                            const TrilinosScalar          new_diag_value)
  {
    compress();
    for (size_type row=0; row<rows.size(); ++row)
      clear_row(rows[row], new_diag_value);

    // This function needs to be called
    // on all processors. We change some
    // data, so we need to flush the
    // buffers to make sure that the
    // right data is used.
    compress();
  }



  TrilinosScalar
  SparseMatrix::operator() (const size_type i,
                            const size_type j) const
  {
    // Extract local indices in
    // the matrix.
    int trilinos_i = matrix->LRID(static_cast<TrilinosWrappers::types::int_type>(i)),
        trilinos_j = matrix->LCID(static_cast<TrilinosWrappers::types::int_type>(j));
    TrilinosScalar value = 0.;

    // If the data is not on the
    // present processor, we throw
    // an exception. This is one of
    // the two tiny differences to
    // the el(i,j) call, which does
    // not throw any assertions.
    if (trilinos_i == -1)
      {
        Assert (false, ExcAccessToNonLocalElement(i, j, local_range().first,
                                                  local_range().second));
      }
    else
      {
        // Check whether the matrix has
        // already been transformed to local
        // indices.
        Assert (matrix->Filled(), ExcMatrixNotCompressed());

        // Prepare pointers for extraction
        // of a view of the row.
        int nnz_present = matrix->NumMyEntries(trilinos_i);
        int nnz_extracted;
        int *col_indices;
        TrilinosScalar *values;

        // Generate the view and make
        // sure that we have not generated
        // an error.
        // TODO Check that col_indices are int and not long long
        int ierr = matrix->ExtractMyRowView(trilinos_i, nnz_extracted,
                                            values, col_indices);
        Assert (ierr==0, ExcTrilinosError(ierr));

        Assert (nnz_present == nnz_extracted,
                ExcDimensionMismatch(nnz_present, nnz_extracted));

        // Search the index where we
        // look for the value, and then
        // finally get it.

        int *el_find = std::find(col_indices, col_indices + nnz_present, trilinos_j);

        int local_col_index = (int)(el_find - col_indices);

        // This is actually the only
        // difference to the el(i,j)
        // function, which means that
        // we throw an exception in
        // this case instead of just
        // returning zero for an
        // element that is not present
        // in the sparsity pattern.
        if (local_col_index == nnz_present)
          {
            Assert (false, ExcInvalidIndex (i,j));
          }
        else
          value = values[local_col_index];
      }

    return value;
  }



  TrilinosScalar
  SparseMatrix::el (const size_type i,
                    const size_type j) const
  {
    // Extract local indices in
    // the matrix.
    int trilinos_i = matrix->LRID(static_cast<TrilinosWrappers::types::int_type>(i)),
        trilinos_j = matrix->LCID(static_cast<TrilinosWrappers::types::int_type>(j));
    TrilinosScalar value = 0.;

    // If the data is not on the
    // present processor, we can't
    // continue. Just print out zero
    // as discussed in the
    // documentation of this
    // function. if you want error
    // checking, use operator().
    if ((trilinos_i == -1 ) || (trilinos_j == -1))
      return 0.;
    else
      {
        // Check whether the matrix
        // already is transformed to
        // local indices.
        Assert (matrix->Filled(), ExcMatrixNotCompressed());

        // Prepare pointers for extraction
        // of a view of the row.
        int nnz_present = matrix->NumMyEntries(trilinos_i);
        int nnz_extracted;
        int *col_indices;
        TrilinosScalar *values;

        // Generate the view and make
        // sure that we have not generated
        // an error.
        int ierr = matrix->ExtractMyRowView(trilinos_i, nnz_extracted,
                                            values, col_indices);
        Assert (ierr==0, ExcTrilinosError(ierr));

        Assert (nnz_present == nnz_extracted,
                ExcDimensionMismatch(nnz_present, nnz_extracted));

        // Search the index where we
        // look for the value, and then
        // finally get it.
        int *el_find = std::find(col_indices, col_indices + nnz_present, trilinos_j);

        int local_col_index = (int)(el_find - col_indices);


        // This is actually the only
        // difference to the () function
        // querying (i,j), where we throw an
        // exception instead of just
        // returning zero for an element
        // that is not present in the
        // sparsity pattern.
        if (local_col_index == nnz_present)
          value = 0;
        else
          value = values[local_col_index];
      }

    return value;
  }



  TrilinosScalar
  SparseMatrix::diag_element (const size_type i) const
  {
    Assert (m() == n(), ExcNotQuadratic());

#ifdef DEBUG
    // use operator() in debug mode because
    // it checks if this is a valid element
    // (in parallel)
    return operator()(i,i);
#else
    // Trilinos doesn't seem to have a
    // more efficient way to access the
    // diagonal than by just using the
    // standard el(i,j) function.
    return el(i,i);
#endif
  }



  unsigned int
  SparseMatrix::row_length (const size_type row) const
  {
    Assert (row < m(), ExcInternalError());

    // get a representation of the
    // present row
    int ncols = -1;
    int local_row = matrix->LRID(static_cast<TrilinosWrappers::types::int_type>(row));

    // on the processor who owns this
    // row, we'll have a non-negative
    // value.
    if (local_row >= 0)
      {
        int ierr = matrix->NumMyRowEntries (local_row, ncols);
        AssertThrow (ierr == 0, ExcTrilinosError(ierr));
      }

    return ncols;
  }



  namespace internals
  {
    typedef dealii::types::global_dof_index size_type;

    void perform_mmult (const SparseMatrix &inputleft,
                        const SparseMatrix &inputright,
                        SparseMatrix       &result,
                        const VectorBase   &V,
                        const bool          transpose_left)
    {
#ifdef DEAL_II_WITH_64BIT_INDICES
      Assert(false,ExcNotImplemented())
#endif
      const bool use_vector = (V.size() == inputright.m() ? true : false);
      if (transpose_left == false)
        {
          Assert (inputleft.n() == inputright.m(),
                  ExcDimensionMismatch(inputleft.n(), inputright.m()));
          Assert (inputleft.domain_partitioner().SameAs(inputright.range_partitioner()),
                  ExcMessage ("Parallel partitioning of A and B does not fit."));
        }
      else
        {
          Assert (inputleft.m() == inputright.m(),
                  ExcDimensionMismatch(inputleft.m(), inputright.m()));
          Assert (inputleft.range_partitioner().SameAs(inputright.range_partitioner()),
                  ExcMessage ("Parallel partitioning of A and B does not fit."));
        }

      result.clear();

      // create a suitable operator B: in case
      // we do not use a vector, all we need to
      // do is to set the pointer. Otherwise,
      // we insert the data from B, but
      // multiply each row with the respective
      // vector element.
      Teuchos::RCP<Epetra_CrsMatrix> mod_B;
      if (use_vector == false)
        {
          mod_B = Teuchos::rcp(const_cast<Epetra_CrsMatrix *>
                               (&inputright.trilinos_matrix()),
                               false);
        }
      else
        {
          mod_B = Teuchos::rcp(new Epetra_CrsMatrix
                               (Copy, inputright.trilinos_sparsity_pattern()),
                               true);
          mod_B->FillComplete(inputright.domain_partitioner(),
                              inputright.range_partitioner());
          Assert (inputright.local_range() == V.local_range(),
                  ExcMessage ("Parallel distribution of matrix B and vector V "
                              "does not match."));

          const int local_N = inputright.local_size();
          for (int i=0; i<local_N; ++i)
            {
              int N_entries = -1;
              double *new_data, *B_data;
              mod_B->ExtractMyRowView (i, N_entries, new_data);
              inputright.trilinos_matrix().ExtractMyRowView (i, N_entries, B_data);
              double value = V.trilinos_vector()[0][i];
              for (TrilinosWrappers::types::int_type j=0; j<N_entries; ++j)
                new_data[j] = value * B_data[j];
            }
        }

      // use ML built-in method for performing
      // the matrix-matrix product.
      // create ML operators on top of the
      // Epetra matrices. if we use a
      // transposed matrix, let ML know it
      ML_Comm *comm;
      ML_Comm_Create(&comm);
#ifdef ML_MPI
      const Epetra_MpiComm *epcomm = dynamic_cast<const Epetra_MpiComm *>(&(inputleft.trilinos_matrix().Comm()));
      // Get the MPI communicator, as it may not be MPI_COMM_W0RLD, and update the ML comm object
      if (epcomm) ML_Comm_Set_UsrComm(comm,epcomm->Comm());
#endif
      ML_Operator *A_ = ML_Operator_Create(comm);
      ML_Operator *B_ = ML_Operator_Create(comm);
      ML_Operator *C_ = ML_Operator_Create(comm);
      SparseMatrix transposed_mat;

      if (transpose_left == false)
        ML_Operator_WrapEpetraCrsMatrix
        (const_cast<Epetra_CrsMatrix *>(&inputleft.trilinos_matrix()),A_,
         false);
      else
        {
          // create transposed matrix
          SparsityPattern sparsity_transposed (inputleft.domain_partitioner(),
                                               inputleft.range_partitioner());
          Assert (inputleft.domain_partitioner().LinearMap() == true,
                  ExcMessage("Matrix must be partitioned contiguously between procs."));
          for (unsigned int i=0; i<inputleft.local_size(); ++i)
            {
              int num_entries, * indices;
              inputleft.trilinos_sparsity_pattern().ExtractMyRowView(i, num_entries,
                                                                     indices);
              Assert (num_entries >= 0, ExcInternalError());
#ifndef DEAL_II_WITH_64BIT_INDICES
              const size_type GID = inputleft.row_partitioner().GID(i);
              for (TrilinosWrappers::types::int_type j=0; j<num_entries; ++j)
                sparsity_transposed.add (inputleft.col_partitioner().GID(indices[j]),
                                         GID);
#else
              const size_type GID = inputleft.row_partitioner().GID64(i);
              for (TrilinosWrappers::types::int_type j=0; j<num_entries; ++j)
                sparsity_transposed.add (inputleft.col_partitioner().GID64(indices[j]),
                                         GID);
#endif
            }

          sparsity_transposed.compress();
          transposed_mat.reinit (sparsity_transposed);
          for (unsigned int i=0; i<inputleft.local_size(); ++i)
            {
              int num_entries, * indices;
              double *values;
              inputleft.trilinos_matrix().ExtractMyRowView(i, num_entries,
                                                           values, indices);
              Assert (num_entries >= 0, ExcInternalError());
#ifndef DEAL_II_WITH_64BIT_INDICES
              const size_type GID = inputleft.row_partitioner().GID(i);
              for (TrilinosWrappers::types::int_type j=0; j<num_entries; ++j)
                transposed_mat.set (inputleft.col_partitioner().GID(indices[j]),
                                    GID, values[j]);
#else
              const size_type GID = inputleft.row_partitioner().GID64(i);
              for (TrilinosWrappers::types::int_type j=0; j<num_entries; ++j)
                transposed_mat.set (inputleft.col_partitioner().GID64(indices[j]),
                                    GID, values[j]);
#endif
            }
          transposed_mat.compress();
          ML_Operator_WrapEpetraCrsMatrix
          (const_cast<Epetra_CrsMatrix *>(&transposed_mat.trilinos_matrix()),
           A_,false);
        }
      ML_Operator_WrapEpetraCrsMatrix(mod_B.get(),B_,false);

      // We implement the multiplication by
      // hand in a similar way as is done in
      // ml/src/Operator/ml_rap.c for a triple
      // matrix product. This means that the
      // code is very similar to the one found
      // in ml/src/Operator/ml_rap.c

      // import data if necessary
      ML_Operator *Btmp, *Ctmp, *Ctmp2, *tptr;
      ML_CommInfoOP *getrow_comm;
      int max_per_proc;
      TrilinosWrappers::types::int_type N_input_vector = B_->invec_leng;
      getrow_comm = B_->getrow->pre_comm;
      if ( getrow_comm != NULL)
        for (TrilinosWrappers::types::int_type i = 0; i < getrow_comm->N_neighbors; i++)
          for (TrilinosWrappers::types::int_type j = 0; j < getrow_comm->neighbors[i].N_send; j++)
            AssertThrow (getrow_comm->neighbors[i].send_list[j] < N_input_vector,
                         ExcInternalError());

      ML_create_unique_col_id(N_input_vector, &(B_->getrow->loc_glob_map),
                              getrow_comm, &max_per_proc, B_->comm);
      B_->getrow->use_loc_glob_map = ML_YES;
      if (A_->getrow->pre_comm != NULL)
        ML_exchange_rows( B_, &Btmp, A_->getrow->pre_comm);
      else Btmp = B_;

      // perform matrix-matrix product
      ML_matmat_mult(A_, Btmp , &Ctmp);

      // release temporary structures we needed
      // for multiplication
      ML_free(B_->getrow->loc_glob_map);
      B_->getrow->loc_glob_map = NULL;
      B_->getrow->use_loc_glob_map = ML_NO;
      if (A_->getrow->pre_comm != NULL)
        {
          tptr = Btmp;
          while ( (tptr!= NULL) && (tptr->sub_matrix != B_))
            tptr = tptr->sub_matrix;
          if (tptr != NULL) tptr->sub_matrix = NULL;
          ML_RECUR_CSR_MSRdata_Destroy(Btmp);
          ML_Operator_Destroy(&Btmp);
        }

      // make correct data structures
      if (A_->getrow->post_comm != NULL)
        ML_exchange_rows(Ctmp, &Ctmp2, A_->getrow->post_comm);
      else
        Ctmp2 = Ctmp;

      ML_back_to_csrlocal(Ctmp2, C_, max_per_proc);

      ML_RECUR_CSR_MSRdata_Destroy (Ctmp);
      ML_Operator_Destroy (&Ctmp);

      if (A_->getrow->post_comm != NULL)
        {
          ML_RECUR_CSR_MSRdata_Destroy(Ctmp2);
          ML_Operator_Destroy (&Ctmp2);
        }

      // create an Epetra matrix from the ML
      // matrix that we got as a result.
      Epetra_CrsMatrix *C_mat;
      ML_Operator2EpetraCrsMatrix(C_, C_mat);
      C_mat->FillComplete();
      C_mat->OptimizeStorage();
      result.reinit (*C_mat);

      // destroy allocated memory
      delete C_mat;
      ML_Operator_Destroy (&A_);
      ML_Operator_Destroy (&B_);
      ML_Operator_Destroy (&C_);
      ML_Comm_Destroy (&comm);
    }
  }


  void
  SparseMatrix::mmult (SparseMatrix       &C,
                       const SparseMatrix &B,
                       const VectorBase   &V) const
  {
#ifdef DEAL_II_WITH_64BIT_INDICES
    Assert(false,ExcNotImplemented())
#endif
    internals::perform_mmult (*this, B, C, V, false);
  }



  void
  SparseMatrix::Tmmult (SparseMatrix       &C,
                        const SparseMatrix &B,
                        const VectorBase   &V) const
  {
#ifdef DEAL_II_WITH_64BIT_INDICES
    Assert(false,ExcNotImplemented())
#endif
    internals::perform_mmult (*this, B, C, V, true);
  }



  void
  SparseMatrix::add (const TrilinosScalar  factor,
                     const SparseMatrix   &rhs)
  {
    Assert (rhs.m() == m(), ExcDimensionMismatch (rhs.m(), m()));
    Assert (rhs.n() == n(), ExcDimensionMismatch (rhs.n(), n()));

    const std::pair<size_type, size_type>
    local_range = rhs.local_range();

    int ierr;

    // If both matrices have been transformed to local index space (in
    // Trilinos speak: they are filled) and the matrices are based on the same
    // sparsity pattern, we can extract views of the column data on both
    // matrices and simply manipulate the values that are addressed by the
    // pointers.
    if (matrix->Filled() == true &&
        rhs.matrix->Filled() == true &&
        &matrix->Graph() == &rhs.matrix->Graph())
      for (size_type row=local_range.first;
           row < local_range.second; ++row)
        {
          Assert (matrix->NumGlobalEntries(row) ==
                  rhs.matrix->NumGlobalEntries(row),
                  ExcDimensionMismatch(matrix->NumGlobalEntries(row),
                                       rhs.matrix->NumGlobalEntries(row)));

          const TrilinosWrappers::types::int_type row_local =
            matrix->RowMap().LID(static_cast<TrilinosWrappers::types::int_type>(row));
          int n_entries, rhs_n_entries;
          TrilinosScalar *value_ptr, *rhs_value_ptr;

          // In debug mode, we want to check whether the indices really are
          // the same in the calling matrix and the input matrix. The reason
          // for doing this only in debug mode is that both extracting indices
          // and comparing indices is relatively slow compared to just working
          // with the values.
#ifdef DEBUG
          int *index_ptr, *rhs_index_ptr;
          ierr = rhs.matrix->ExtractMyRowView (row_local, rhs_n_entries,
                                               rhs_value_ptr, rhs_index_ptr);
          Assert (ierr == 0, ExcTrilinosError(ierr));

          ierr = matrix->ExtractMyRowView (row_local, n_entries, value_ptr,
                                           index_ptr);
          Assert (ierr == 0, ExcTrilinosError(ierr));
#else
          rhs.matrix->ExtractMyRowView (row_local, rhs_n_entries,rhs_value_ptr);
          matrix->ExtractMyRowView (row_local, n_entries, value_ptr);
#endif

          AssertDimension (n_entries, rhs_n_entries);

          for (TrilinosWrappers::types::int_type i=0; i<n_entries; ++i)
            {
              *value_ptr++ += *rhs_value_ptr++ * factor;
#ifdef DEBUG
              Assert (*index_ptr++ == *rhs_index_ptr++,
                      ExcInternalError());
#endif
            }
        }
    // If we have different sparsity patterns, we have to be more careful (in
    // particular when we use multiple processors) and extract a copy of the
    // row data, multiply it by the factor and then add it to the matrix using
    // the respective add() function.
    else
      {
        int max_row_length = 0;
        for (size_type row=local_range.first;
             row < local_range.second; ++row)
          max_row_length
            = std::max (max_row_length,rhs.matrix->NumGlobalEntries(row));

        std::vector<TrilinosScalar> values (max_row_length);
        std::vector<TrilinosWrappers::types::int_type> column_indices (max_row_length);
        for (size_type row=local_range.first;
             row < local_range.second; ++row)
          {
            int n_entries;
            ierr = rhs.matrix->Epetra_CrsMatrix::ExtractGlobalRowCopy
                   (static_cast<TrilinosWrappers::types::int_type>(row),
                    max_row_length, n_entries, &values[0], &column_indices[0]);
            Assert (ierr == 0, ExcTrilinosError(ierr));

            // Filter away zero elements
            unsigned int n_actual_entries = 0.;
            for (int i=0; i<n_entries; ++i)
              if (std::abs(values[i]) != 0.)
                {
                  column_indices[n_actual_entries] = column_indices[i];
                  values[n_actual_entries++] = values[i] * factor;
                }

            ierr = matrix->Epetra_CrsMatrix::SumIntoGlobalValues
                   (static_cast<TrilinosWrappers::types::int_type>(row),
                    n_actual_entries, &values[0], &column_indices[0]);
            Assert (ierr == 0,
                    ExcMessage("Adding the entries from the other matrix "
                               "failed, possibly because the sparsity pattern "
                               "of that matrix includes more elements than the "
                               "calling matrix, which is not allowed."));
          }
        compress (VectorOperation::add);

      }
  }



  void
  SparseMatrix::transpose ()
  {
    // This only flips a flag that tells
    // Trilinos that any vmult operation
    // should be done with the
    // transpose. However, the matrix
    // structure is not reset.
    int ierr;

    if (!matrix->UseTranspose())
      {
        ierr = matrix->SetUseTranspose (true);
        AssertThrow (ierr == 0, ExcTrilinosError(ierr));
      }
    else
      {
        ierr = matrix->SetUseTranspose (false);
        AssertThrow (ierr == 0, ExcTrilinosError(ierr));
      }
  }



  void
  SparseMatrix::write_ascii ()
  {
    Assert (false, ExcNotImplemented());
  }



  // As of now, no particularly neat
  // ouput is generated in case of
  // multiple processors.
  void
  SparseMatrix::print (std::ostream &out,
                       const bool    print_detailed_trilinos_information) const
  {
    if (print_detailed_trilinos_information == true)
      out << *matrix;
    else
      {
        double *values;
        int *indices;
        int num_entries;

        for (int i=0; i<matrix->NumMyRows(); ++i)
          {
            matrix->ExtractMyRowView (i, num_entries, values, indices);
            for (TrilinosWrappers::types::int_type j=0; j<num_entries; ++j)
              out << "(" << global_row_index(*matrix,i) << ","
                  << global_column_index(*matrix,indices[j]) << ") "
                  << values[j] << std::endl;
          }
      }

    AssertThrow (out, ExcIO());
  }



  SparseMatrix::size_type
  SparseMatrix::memory_consumption () const
  {
    size_type static_memory = sizeof(this) + sizeof (*matrix)
                              + sizeof(*matrix->Graph().DataPtr());
    return ((sizeof(TrilinosScalar)+sizeof(TrilinosWrappers::types::int_type))*
            matrix->NumMyNonzeros() + sizeof(int)*local_size() + static_memory);
  }
}



// explicit instantiations
#include "trilinos_sparse_matrix.inst"


// TODO: put these instantiations into generic file
namespace TrilinosWrappers
{
  template void
  SparseMatrix::reinit (const dealii::SparsityPattern &);
  template void
  SparseMatrix::reinit (const CompressedSparsityPattern &);
  template void
  SparseMatrix::reinit (const CompressedSetSparsityPattern &);
  template void
  SparseMatrix::reinit (const CompressedSimpleSparsityPattern &);

  template void
  SparseMatrix::reinit (const Epetra_Map &,
                        const dealii::SparsityPattern &,
                        const bool);
  template void
  SparseMatrix::reinit (const Epetra_Map &,
                        const CompressedSparsityPattern &,
                        const bool);
  template void
  SparseMatrix::reinit (const Epetra_Map &,
                        const CompressedSetSparsityPattern &,
                        const bool);
  template void
  SparseMatrix::reinit (const Epetra_Map &,
                        const CompressedSimpleSparsityPattern &,
                        const bool);


  template void
  SparseMatrix::reinit (const Epetra_Map &,
                        const Epetra_Map &,
                        const dealii::SparsityPattern &,
                        const bool);
  template void
  SparseMatrix::reinit (const Epetra_Map &,
                        const Epetra_Map &,
                        const CompressedSparsityPattern &,
                        const bool);
  template void
  SparseMatrix::reinit (const Epetra_Map &,
                        const Epetra_Map &,
                        const CompressedSetSparsityPattern &,
                        const bool);

}

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_TRILINOS
