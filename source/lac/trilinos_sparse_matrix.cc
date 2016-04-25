// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2015 by the deal.II authors
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
#  include <deal.II/lac/dynamic_sparsity_pattern.h>
#  include <deal.II/lac/sparsity_tools.h>
#  include <deal.II/lac/parallel_vector.h>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#  include <Epetra_Export.h>
#  include <ml_epetra_utils.h>
#  include <ml_struct.h>
#  include <Teuchos_RCP.hpp>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

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
      // if we are asked to visit the past-the-end line, then simply
      // release all our caches and go on with life.
      //
      // do the same if the row we're supposed to visit is not locally
      // owned. this is simply going to make non-locally owned rows
      // look like they're empty
      if ((this->a_row == matrix->m())
          ||
          (matrix->in_local_range (this->a_row) == false))
        {
          colnum_cache.reset ();
          value_cache.reset ();

          return;
        }

      // get a representation of the present row
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
    compress(VectorOperation::insert);
  }



  SparseMatrix::~SparseMatrix ()
  {}



  void
  SparseMatrix::copy_from (const SparseMatrix &rhs)
  {
    if (this == &rhs)
      return;

    nonlocal_matrix.reset();
    nonlocal_matrix_exporter.reset();

    // check whether we need to update the whole matrix layout (we have
    // different maps or if we detect a row where the columns of the two
    // matrices do not match)
    bool needs_deep_copy =
      !matrix->RowMap().SameAs(rhs.matrix->RowMap()) ||
      !matrix->ColMap().SameAs(rhs.matrix->ColMap()) ||
      !matrix->DomainMap().SameAs(rhs.matrix->DomainMap()) ||
      n_nonzero_elements() != rhs.n_nonzero_elements();
    if (!needs_deep_copy)
      {
        const std::pair<size_type, size_type>
        local_range = rhs.local_range();

        int ierr;
        // Try to copy all the rows of the matrix one by one. In case of error
        // (i.e., the column indices are different), we need to abort and blow
        // away the matrix.
        for (size_type row=local_range.first; row < local_range.second; ++row)
          {
            const int row_local =
              matrix->RowMap().LID(static_cast<TrilinosWrappers::types::int_type>(row));

            int n_entries, rhs_n_entries;
            TrilinosScalar *value_ptr, *rhs_value_ptr;
            int *index_ptr, *rhs_index_ptr;
            ierr = rhs.matrix->ExtractMyRowView (row_local, rhs_n_entries,
                                                 rhs_value_ptr, rhs_index_ptr);
            (void)ierr;
            Assert (ierr == 0, ExcTrilinosError(ierr));

            ierr = matrix->ExtractMyRowView (row_local, n_entries, value_ptr,
                                             index_ptr);
            Assert (ierr == 0, ExcTrilinosError(ierr));

            if (n_entries != rhs_n_entries ||
                std::memcmp(static_cast<void *>(index_ptr),
                            static_cast<void *>(rhs_index_ptr),
                            sizeof(int)*n_entries) != 0)
              {
                needs_deep_copy = true;
                break;
              }

            for (int i=0; i<n_entries; ++i)
              value_ptr[i] = rhs_value_ptr[i];
          }
      }

    if (needs_deep_copy)
      {
        column_space_map.reset (new Epetra_Map (rhs.domain_partitioner()));

        // release memory before reallocation
        matrix.reset ();
        matrix.reset (new Epetra_FECrsMatrix(*rhs.matrix));

        matrix->FillComplete(*column_space_map, matrix->RowMap());
      }

    if (rhs.nonlocal_matrix.get() != 0)
      nonlocal_matrix.reset(new Epetra_CrsMatrix(Copy, rhs.nonlocal_matrix->Graph()));
  }



  namespace
  {
    typedef SparseMatrix::size_type size_type;

    template <typename SparsityPatternType>
    void
    reinit_matrix (const Epetra_Map                          &input_row_map,
                   const Epetra_Map                          &input_col_map,
                   const SparsityPatternType                 &sparsity_pattern,
                   const bool                                 exchange_data,
                   std_cxx11::shared_ptr<Epetra_Map>         &column_space_map,
                   std_cxx11::shared_ptr<Epetra_FECrsMatrix> &matrix,
                   std_cxx11::shared_ptr<Epetra_CrsMatrix>   &nonlocal_matrix,
                   std_cxx11::shared_ptr<Epetra_Export>      &nonlocal_matrix_exporter)
    {
      // release memory before reallocation
      matrix.reset();
      nonlocal_matrix.reset();
      nonlocal_matrix_exporter.reset();

      if (input_row_map.Comm().MyPID() == 0)
        {
          AssertDimension (sparsity_pattern.n_rows(),
                           static_cast<size_type>(n_global_elements(input_row_map)));
          AssertDimension (sparsity_pattern.n_cols(),
                           static_cast<size_type>(n_global_elements(input_col_map)));
        }

      column_space_map.reset (new Epetra_Map (input_col_map));

      // if we want to exchange data, build a usual Trilinos sparsity pattern
      // and let that handle the exchange. otherwise, manually create a
      // CrsGraph, which consumes considerably less memory because it can set
      // correct number of indices right from the start
      if (exchange_data)
        {
          SparsityPattern trilinos_sparsity;
          trilinos_sparsity.reinit (input_row_map, input_col_map,
                                    sparsity_pattern, exchange_data);
          matrix.reset (new Epetra_FECrsMatrix
                        (Copy, trilinos_sparsity.trilinos_sparsity_pattern(), false));

          return;
        }

      const size_type first_row = min_my_gid(input_row_map),
                      last_row = max_my_gid(input_row_map)+1;
      std::vector<int> n_entries_per_row(last_row-first_row);

      for (size_type row=first_row; row<last_row; ++row)
        n_entries_per_row[row-first_row] = sparsity_pattern.row_length(row);

      // The deal.II notation of a Sparsity pattern corresponds to the Epetra
      // concept of a Graph. Hence, we generate a graph by copying the
      // sparsity pattern into it, and then build up the matrix from the
      // graph. This is considerable faster than directly filling elements
      // into the matrix. Moreover, it consumes less memory, since the
      // internal reordering is done on ints only, and we can leave the
      // doubles aside.

      // for more than one processor, need to specify only row map first and
      // let the matrix entries decide about the column map (which says which
      // columns are present in the matrix, not to be confused with the
      // col_map that tells how the domain dofs of the matrix will be
      // distributed). for only one processor, we can directly assign the
      // columns as well. Compare this with bug # 4123 in the Sandia Bugzilla.
      std_cxx11::shared_ptr<Epetra_CrsGraph> graph;
      if (input_row_map.Comm().NumProc() > 1)
        graph.reset (new Epetra_CrsGraph (Copy, input_row_map,
                                          &n_entries_per_row[0], true));
      else
        graph.reset (new Epetra_CrsGraph (Copy, input_row_map, input_col_map,
                                          &n_entries_per_row[0], true));

      // This functions assumes that the sparsity pattern sits on all
      // processors (completely). The parallel version uses an Epetra graph
      // that is already distributed.

      // now insert the indices
      std::vector<TrilinosWrappers::types::int_type>   row_indices;

      for (size_type row=first_row; row<last_row; ++row)
        {
          const int row_length = sparsity_pattern.row_length(row);
          if (row_length == 0)
            continue;

          row_indices.resize (row_length, -1);
          {
            typename SparsityPatternType::iterator p = sparsity_pattern.begin(row);
            for (size_type col=0; p != sparsity_pattern.end(row); ++p, ++col)
              row_indices[col] = p->column();
          }
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
      AssertDimension (sparsity_pattern.n_cols(),
                       static_cast<size_type>(n_global_cols(*graph)));
      (void)n_global_cols;

      // And now finally generate the matrix.
      matrix.reset (new Epetra_FECrsMatrix(Copy, *graph, false));
    }



    // for the non-local graph, we need to circumvent the problem that some
    // processors will not add into the non-local graph at all: We do not want
    // to insert dummy elements on >5000 processors because that gets very
    // slow. Thus, we set a flag in Epetra_CrsGraph that sets the correct
    // flag. Since it is protected, we need to expose this information by
    // deriving a class from Epetra_CrsGraph for the purpose of creating the
    // data structure
    class Epetra_CrsGraphMod : public Epetra_CrsGraph
    {
    public:
      Epetra_CrsGraphMod (const Epetra_Map &row_map,
                          const int        *n_entries_per_row)
        :
        Epetra_CrsGraph(Copy, row_map, n_entries_per_row, true)
      {};

      void SetIndicesAreGlobal()
      {
        this->Epetra_CrsGraph::SetIndicesAreGlobal(true);
      }
    };


    // specialization for DynamicSparsityPattern which can provide us with
    // more information about the non-locally owned rows
    template <>
    void
    reinit_matrix (const Epetra_Map                          &input_row_map,
                   const Epetra_Map                          &input_col_map,
                   const DynamicSparsityPattern              &sparsity_pattern,
                   const bool                                 exchange_data,
                   std_cxx11::shared_ptr<Epetra_Map>         &column_space_map,
                   std_cxx11::shared_ptr<Epetra_FECrsMatrix> &matrix,
                   std_cxx11::shared_ptr<Epetra_CrsMatrix>   &nonlocal_matrix,
                   std_cxx11::shared_ptr<Epetra_Export>      &nonlocal_matrix_exporter)
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
      // the owned rows. In that case, do not create the nonlocal graph and
      // fill the columns by demand
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

      Epetra_Map off_processor_map(-1, ghost_rows.size(),
                                   (ghost_rows.size()>0)?(&ghost_rows[0]):NULL,
                                   0, input_row_map.Comm());

      std_cxx11::shared_ptr<Epetra_CrsGraph> graph;
      std_cxx11::shared_ptr<Epetra_CrsGraphMod> nonlocal_graph;
      if (input_row_map.Comm().NumProc() > 1)
        {
          graph.reset (new Epetra_CrsGraph (Copy, input_row_map,
                                            (n_entries_per_row.size()>0)?(&n_entries_per_row[0]):NULL,
                                            exchange_data ? false : true));
          if (have_ghost_rows == true)
            nonlocal_graph.reset (new Epetra_CrsGraphMod (off_processor_map,
                                                          &n_entries_per_ghost_row[0]));
        }
      else
        graph.reset (new Epetra_CrsGraph (Copy, input_row_map, input_col_map,
                                          (n_entries_per_row.size()>0)?(&n_entries_per_row[0]):NULL,
                                          true));

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
          for (int col=0; col < row_length; ++col)
            row_indices[col] = sparsity_pattern.column_number(global_row, col);

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
          // must make sure the IndicesAreGlobal flag is set on all processors
          // because some processors might not call InsertGlobalIndices (and
          // we do not want to insert dummy indices on all processors for
          // large-scale simulations due to the bad impact on performance)
          nonlocal_graph->SetIndicesAreGlobal();
          Assert(nonlocal_graph->IndicesAreGlobal() == true,
                 ExcInternalError());
          nonlocal_graph->FillComplete(input_col_map, input_row_map);
          nonlocal_graph->OptimizeStorage();

          // insert data from nonlocal graph into the final sparsity pattern
          if (exchange_data)
            {
              Epetra_Export exporter(nonlocal_graph->RowMap(), input_row_map);
              int ierr = graph->Export(*nonlocal_graph, exporter, Add);
              (void)ierr;
              Assert (ierr==0, ExcTrilinosError(ierr));
            }

          nonlocal_matrix.reset (new Epetra_CrsMatrix(Copy, *nonlocal_graph));
        }

      graph->FillComplete(input_col_map, input_row_map);
      graph->OptimizeStorage();

      AssertDimension (sparsity_pattern.n_cols(),static_cast<size_type>(
                         n_global_cols(*graph)));

      matrix.reset (new Epetra_FECrsMatrix(Copy, *graph, false));
    }
  }



  template <typename SparsityPatternType>
  void
  SparseMatrix::reinit (const SparsityPatternType &sparsity_pattern)
  {
    const Epetra_Map rows (static_cast<TrilinosWrappers::types::int_type>(sparsity_pattern.n_rows()),
                           0,
                           Utilities::Trilinos::comm_self());
    const Epetra_Map columns (static_cast<TrilinosWrappers::types::int_type>(sparsity_pattern.n_cols()),
                              0,
                              Utilities::Trilinos::comm_self());

    reinit_matrix (rows, columns, sparsity_pattern, false,
                   column_space_map, matrix, nonlocal_matrix,
                   nonlocal_matrix_exporter);
  }



  template <typename SparsityPatternType>
  void
  SparseMatrix::reinit (const Epetra_Map          &input_map,
                        const SparsityPatternType &sparsity_pattern,
                        const bool                 exchange_data)
  {
    reinit_matrix (input_map, input_map, sparsity_pattern, exchange_data,
                   column_space_map, matrix, nonlocal_matrix,
                   nonlocal_matrix_exporter);
  }






  template <typename SparsityPatternType>
  inline
  void SparseMatrix::reinit (const IndexSet            &row_parallel_partitioning,
                             const IndexSet            &col_parallel_partitioning,
                             const SparsityPatternType &sparsity_pattern,
                             const MPI_Comm            &communicator,
                             const bool                 exchange_data)
  {
    Epetra_Map row_map =
      row_parallel_partitioning.make_trilinos_map (communicator, false);
    Epetra_Map col_map =
      col_parallel_partitioning.make_trilinos_map (communicator, false);
    reinit_matrix (row_map, col_map, sparsity_pattern, exchange_data,
                   column_space_map, matrix, nonlocal_matrix,
                   nonlocal_matrix_exporter);

    // In the end, the matrix needs to be compressed in order to be really
    // ready.
    last_action = Zero;
    compress(VectorOperation::insert);
  }



  template <typename SparsityPatternType>
  inline
  void SparseMatrix::reinit (const Epetra_Map          &row_map,
                             const Epetra_Map          &col_map,
                             const SparsityPatternType &sparsity_pattern,
                             const bool                 exchange_data)
  {
    reinit_matrix (row_map, col_map, sparsity_pattern, exchange_data,
                   column_space_map, matrix, nonlocal_matrix,
                   nonlocal_matrix_exporter);

    // In the end, the matrix needs to be compressed in order to be really
    // ready.
    last_action = Zero;
    compress(VectorOperation::insert);
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

    last_action = Zero;
    compress(VectorOperation::insert);
  }



  void
  SparseMatrix::reinit (const SparseMatrix &sparse_matrix)
  {
    if (this == &sparse_matrix)
      return;

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

    last_action = Zero;
    compress(VectorOperation::insert);
  }



  template <typename number>
  inline
  void SparseMatrix::reinit (const IndexSet      &row_parallel_partitioning,
                             const IndexSet      &col_parallel_partitioning,
                             const ::dealii::SparseMatrix<number> &dealii_sparse_matrix,
                             const MPI_Comm      &communicator,
                             const double         drop_tolerance,
                             const bool           copy_values,
                             const ::dealii::SparsityPattern *use_this_sparsity)
  {
    if (copy_values == false)
      {
        // in case we do not copy values, just
        // call the other function.
        if (use_this_sparsity == 0)
          reinit (row_parallel_partitioning, col_parallel_partitioning,
                  dealii_sparse_matrix.get_sparsity_pattern(),
                  communicator, false);
        else
          reinit (row_parallel_partitioning, col_parallel_partitioning,
                  *use_this_sparsity, communicator, false);
        return;
      }

    const size_type n_rows = dealii_sparse_matrix.m();

    AssertDimension (row_parallel_partitioning.size(), n_rows);
    AssertDimension (col_parallel_partitioning.size(), dealii_sparse_matrix.n());

    const ::dealii::SparsityPattern &sparsity_pattern =
      (use_this_sparsity!=0)? *use_this_sparsity :
      dealii_sparse_matrix.get_sparsity_pattern();

    if (matrix.get() == 0 ||
        m() != n_rows ||
        n_nonzero_elements() != sparsity_pattern.n_nonzero_elements())
      {
        reinit (row_parallel_partitioning, col_parallel_partitioning,
                sparsity_pattern, communicator, false);
      }

    // fill the values. the same as above: go through all rows of the
    // matrix, and then all columns. since the sparsity patterns of the
    // input matrix and the specified sparsity pattern might be different,
    // need to go through the row for both these sparsity structures
    // simultaneously in order to really set the correct values.
    size_type maximum_row_length = matrix->MaxNumEntries();
    std::vector<size_type> row_indices (maximum_row_length);
    std::vector<TrilinosScalar> values (maximum_row_length);

    for (size_type row=0; row<n_rows; ++row)
      // see if the row is locally stored on this processor
      if (row_parallel_partitioning.is_element(row) == true)
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
    compress(VectorOperation::insert);
  }



  template <typename number>
  void
  SparseMatrix::reinit (const ::dealii::SparseMatrix<number> &dealii_sparse_matrix,
                        const double                          drop_tolerance,
                        const bool                            copy_values,
                        const ::dealii::SparsityPattern      *use_this_sparsity)
  {
    reinit (complete_index_set(dealii_sparse_matrix.m()),
            complete_index_set(dealii_sparse_matrix.n()),
            dealii_sparse_matrix, MPI_COMM_SELF, drop_tolerance,
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
    reinit (IndexSet(input_map), IndexSet(input_map), dealii_sparse_matrix,
            MPI_COMM_SELF, drop_tolerance, copy_values, use_this_sparsity);
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
    reinit (IndexSet(input_row_map), IndexSet(input_col_map),
            dealii_sparse_matrix, MPI_COMM_SELF,
            drop_tolerance, copy_values, use_this_sparsity);
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

    last_action = Zero;
    compress(VectorOperation::insert);
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
          ExcMessage("Operation and argument to compress() do not match"));
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
        (void)ierr;

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
    for (size_type row=0; row<rows.size(); ++row)
      clear_row(rows[row], new_diag_value);
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
        (void)ierr;
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
        (void)ierr;
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



  void
  SparseMatrix::set (const std::vector<size_type>     &row_indices,
                     const std::vector<size_type>     &col_indices,
                     const FullMatrix<TrilinosScalar> &values,
                     const bool                        elide_zero_values)
  {
    Assert (row_indices.size() == values.m(),
            ExcDimensionMismatch(row_indices.size(), values.m()));
    Assert (col_indices.size() == values.n(),
            ExcDimensionMismatch(col_indices.size(), values.n()));

    for (size_type i=0; i<row_indices.size(); ++i)
      set (row_indices[i], col_indices.size(), &col_indices[0], &values(i,0),
           elide_zero_values);
  }



  void
  SparseMatrix::set (const size_type                    row,
                     const std::vector<size_type>      &col_indices,
                     const std::vector<TrilinosScalar> &values,
                     const bool                         elide_zero_values)
  {
    Assert (col_indices.size() == values.size(),
            ExcDimensionMismatch(col_indices.size(), values.size()));

    set (row, col_indices.size(), &col_indices[0], &values[0],
         elide_zero_values);
  }



  void
  SparseMatrix::set (const size_type       row,
                     const size_type       n_cols,
                     const size_type      *col_indices,
                     const TrilinosScalar *values,
                     const bool            elide_zero_values)
  {
    AssertIndexRange(row, this->m());

    int ierr;
    if (last_action == Add)
      {
        ierr = matrix->GlobalAssemble (*column_space_map, matrix->RowMap(),
                                       true);

        Assert (ierr == 0, ExcTrilinosError(ierr));
      }

    last_action = Insert;

    TrilinosWrappers::types::int_type *col_index_ptr;
    TrilinosScalar *col_value_ptr;
    TrilinosWrappers::types::int_type n_columns;

    TrilinosScalar short_val_array[100];
    TrilinosWrappers::types::int_type short_index_array[100];
    std::vector<TrilinosScalar> long_val_array;
    std::vector<TrilinosWrappers::types::int_type> long_index_array;


    // If we don't elide zeros, the pointers are already available... need to
    // cast to non-const pointers as that is the format taken by Trilinos (but
    // we will not modify const data)
    if (elide_zero_values == false)
      {
        col_index_ptr = (TrilinosWrappers::types::int_type *)col_indices;
        col_value_ptr = const_cast<TrilinosScalar *>(values);
        n_columns = n_cols;
      }
    else
      {
        // Otherwise, extract nonzero values in each row and get the
        // respective indices.
        if (n_cols > 100)
          {
            long_val_array.resize(n_cols);
            long_index_array.resize(n_cols);
            col_index_ptr = &long_index_array[0];
            col_value_ptr = &long_val_array[0];
          }
        else
          {
            col_index_ptr = &short_index_array[0];
            col_value_ptr = &short_val_array[0];
          }

        n_columns = 0;
        for (size_type j=0; j<n_cols; ++j)
          {
            const double value = values[j];
            AssertIsFinite(value);
            if (value != 0)
              {
                col_index_ptr[n_columns] = col_indices[j];
                col_value_ptr[n_columns] = value;
                n_columns++;
              }
          }

        Assert(n_columns <= (TrilinosWrappers::types::int_type)n_cols, ExcInternalError());
      }


    // If the calling matrix owns the row to which we want to insert values,
    // we can directly call the Epetra_CrsMatrix input function, which is much
    // faster than the Epetra_FECrsMatrix function. We distinguish between two
    // cases: the first one is when the matrix is not filled (i.e., it is
    // possible to add new elements to the sparsity pattern), and the second
    // one is when the pattern is already fixed. In the former case, we add
    // the possibility to insert new values, and in the second we just replace
    // data.
    if (matrix->RowMap().MyGID(static_cast<TrilinosWrappers::types::int_type>(row)) == true)
      {
        if (matrix->Filled() == false)
          {
            ierr = matrix->Epetra_CrsMatrix::InsertGlobalValues(
                     static_cast<TrilinosWrappers::types::int_type>(row),
                     static_cast<int>(n_columns),const_cast<double *>(col_value_ptr),
                     col_index_ptr);

            // When inserting elements, we do not want to create exceptions in
            // the case when inserting non-local data (since that's what we
            // want to do right now).
            if (ierr > 0)
              ierr = 0;
          }
        else
          ierr = matrix->Epetra_CrsMatrix::ReplaceGlobalValues(row, n_columns,
                                                               col_value_ptr,
                                                               col_index_ptr);
      }
    else
      {
        // When we're at off-processor data, we have to stick with the
        // standard Insert/ReplaceGlobalValues function. Nevertheless, the way
        // we call it is the fastest one (any other will lead to repeated
        // allocation and deallocation of memory in order to call the function
        // we already use, which is very inefficient if writing one element at
        // a time).
        compressed = false;

        if (matrix->Filled() == false)
          {
            ierr = matrix->InsertGlobalValues (1,
                                               (TrilinosWrappers::types::int_type *)&row,
                                               n_columns, col_index_ptr,
                                               &col_value_ptr,
                                               Epetra_FECrsMatrix::ROW_MAJOR);
            if (ierr > 0)
              ierr = 0;
          }
        else
          ierr = matrix->ReplaceGlobalValues (1,
                                              (TrilinosWrappers::types::int_type *)&row,
                                              n_columns, col_index_ptr,
                                              &col_value_ptr,
                                              Epetra_FECrsMatrix::ROW_MAJOR);
        // use the FECrsMatrix facilities for set even in the case when we
        // have explicitly set the off-processor rows because that only works
        // properly when adding elements, not when setting them (since we want
        // to only touch elements that have been set explicitly, and there is
        // no way on the receiving processor to identify them otherwise)
      }

    Assert (ierr <= 0, ExcAccessToNonPresentElement(row, col_index_ptr[0]));
    AssertThrow (ierr >= 0, ExcTrilinosError(ierr));
  }



  void
  SparseMatrix::add (const std::vector<size_type>     &indices,
                     const FullMatrix<TrilinosScalar> &values,
                     const bool                        elide_zero_values)
  {
    Assert (indices.size() == values.m(),
            ExcDimensionMismatch(indices.size(), values.m()));
    Assert (values.m() == values.n(), ExcNotQuadratic());

    for (size_type i=0; i<indices.size(); ++i)
      add (indices[i], indices.size(), &indices[0], &values(i,0),
           elide_zero_values);
  }



  void
  SparseMatrix::add (const std::vector<size_type>  &row_indices,
                     const std::vector<size_type>  &col_indices,
                     const FullMatrix<TrilinosScalar> &values,
                     const bool                        elide_zero_values)
  {
    Assert (row_indices.size() == values.m(),
            ExcDimensionMismatch(row_indices.size(), values.m()));
    Assert (col_indices.size() == values.n(),
            ExcDimensionMismatch(col_indices.size(), values.n()));

    for (size_type i=0; i<row_indices.size(); ++i)
      add (row_indices[i], col_indices.size(), &col_indices[0],
           &values(i,0), elide_zero_values);
  }



  void
  SparseMatrix::add (const size_type                    row,
                     const std::vector<size_type>      &col_indices,
                     const std::vector<TrilinosScalar> &values,
                     const bool                         elide_zero_values)
  {
    Assert (col_indices.size() == values.size(),
            ExcDimensionMismatch(col_indices.size(), values.size()));

    add (row, col_indices.size(), &col_indices[0], &values[0],
         elide_zero_values);
  }



  void
  SparseMatrix::add (const size_type       row,
                     const size_type       n_cols,
                     const size_type      *col_indices,
                     const TrilinosScalar *values,
                     const bool            elide_zero_values,
                     const bool            /*col_indices_are_sorted*/)
  {
    AssertIndexRange(row, this->m());
    int ierr;
    if (last_action == Insert)
      {
        // TODO: this could lead to a dead lock when only one processor
        // calls GlobalAssemble.
        ierr = matrix->GlobalAssemble(*column_space_map,
                                      matrix->RowMap(), false);

        AssertThrow (ierr == 0, ExcTrilinosError(ierr));
      }

    last_action = Add;

    TrilinosWrappers::types::int_type *col_index_ptr;
    TrilinosScalar *col_value_ptr;
    TrilinosWrappers::types::int_type n_columns;

    double short_val_array[100];
    TrilinosWrappers::types::int_type short_index_array[100];
    std::vector<TrilinosScalar> long_val_array;
    std::vector<TrilinosWrappers::types::int_type> long_index_array;

    // If we don't elide zeros, the pointers are already available... need to
    // cast to non-const pointers as that is the format taken by Trilinos (but
    // we will not modify const data)
    if (elide_zero_values == false)
      {
        col_index_ptr = (TrilinosWrappers::types::int_type *)col_indices;
        col_value_ptr = const_cast<TrilinosScalar *>(values);
        n_columns = n_cols;
#ifdef DEBUG
        for (size_type j=0; j<n_cols; ++j)
          AssertIsFinite(values[j]);
#endif
      }
    else
      {
        // Otherwise, extract nonzero values in each row and the corresponding
        // index.
        if (n_cols > 100)
          {
            long_val_array.resize(n_cols);
            long_index_array.resize(n_cols);
            col_index_ptr = &long_index_array[0];
            col_value_ptr = &long_val_array[0];
          }
        else
          {
            col_index_ptr = &short_index_array[0];
            col_value_ptr = &short_val_array[0];
          }

        n_columns = 0;
        for (size_type j=0; j<n_cols; ++j)
          {
            const double value = values[j];

            AssertIsFinite(value);
            if (value != 0)
              {
                col_index_ptr[n_columns] = col_indices[j];
                col_value_ptr[n_columns] = value;
                n_columns++;
              }
          }

        Assert(n_columns <= (TrilinosWrappers::types::int_type)n_cols, ExcInternalError());

      }
    // Exit early if there is nothing to do
    if (n_columns == 0)
      {
        return;
      }

    // If the calling processor owns the row to which we want to add values, we
    // can directly call the Epetra_CrsMatrix input function, which is much
    // faster than the Epetra_FECrsMatrix function.
    if (matrix->RowMap().MyGID(static_cast<TrilinosWrappers::types::int_type>(row)) == true)
      {
        ierr = matrix->Epetra_CrsMatrix::SumIntoGlobalValues(row, n_columns,
                                                             col_value_ptr,
                                                             col_index_ptr);
      }
    else if (nonlocal_matrix.get() != 0)
      {
        compressed = false;
        // this is the case when we have explicitly set the off-processor rows
        // and want to create a separate matrix object for them (to retain
        // thread-safety)
        Assert (nonlocal_matrix->RowMap().LID(static_cast<TrilinosWrappers::types::int_type>(row)) != -1,
                ExcMessage("Attempted to write into off-processor matrix row "
                           "that has not be specified as being writable upon "
                           "initialization"));
        ierr = nonlocal_matrix->SumIntoGlobalValues(row, n_columns,
                                                    col_value_ptr,
                                                    col_index_ptr);
      }
    else
      {
        // When we're at off-processor data, we have to stick with the
        // standard SumIntoGlobalValues function. Nevertheless, the way we
        // call it is the fastest one (any other will lead to repeated
        // allocation and deallocation of memory in order to call the function
        // we already use, which is very inefficient if writing one element at
        // a time).
        compressed = false;

        ierr = matrix->SumIntoGlobalValues (1,
                                            (TrilinosWrappers::types::int_type *)&row, n_columns,
                                            col_index_ptr,
                                            &col_value_ptr,
                                            Epetra_FECrsMatrix::ROW_MAJOR);
      }

#ifdef DEBUG
    if (ierr > 0)
      {
        std::cout << "------------------------------------------"
                  << std::endl;
        std::cout << "Got error " << ierr << " in row " << row
                  << " of proc " << matrix->RowMap().Comm().MyPID()
                  << " when trying to add the columns:" << std::endl;
        for (TrilinosWrappers::types::int_type i=0; i<n_columns; ++i)
          std::cout << col_index_ptr[i] << " ";
        std::cout << std::endl << std::endl;
        std::cout << "Matrix row "
                  << (matrix->RowMap().MyGID(static_cast<TrilinosWrappers::types::int_type>(row)) == false ? "(nonlocal part)" : "")
                  << " has the following indices:" << std::endl;
        std::vector<TrilinosWrappers::types::int_type> indices;
        const Epetra_CrsGraph *graph =
          (nonlocal_matrix.get() != 0 &&
           matrix->RowMap().MyGID(static_cast<TrilinosWrappers::types::int_type>(row)) == false) ?
          &nonlocal_matrix->Graph() : &matrix->Graph();

        indices.resize(graph->NumGlobalIndices(static_cast<TrilinosWrappers::types::int_type>(row)));
        int n_indices = 0;
        graph->ExtractGlobalRowCopy(static_cast<TrilinosWrappers::types::int_type>(row),
                                    indices.size(), n_indices, &indices[0]);
        AssertDimension(static_cast<unsigned int>(n_indices), indices.size());

        for (TrilinosWrappers::types::int_type i=0; i<n_indices; ++i)
          std::cout << indices[i] << " ";
        std::cout << std::endl << std::endl;
        Assert (ierr <= 0,
                ExcAccessToNonPresentElement(row, col_index_ptr[0]));
      }
#endif
    Assert (ierr >= 0, ExcTrilinosError(ierr));
  }



  SparseMatrix &
  SparseMatrix::operator = (const double d)
  {
    Assert (d==0, ExcScalarAssignmentOnlyForZeroValue());
    compress (::dealii::VectorOperation::unknown); // TODO: why do we do this? Should we not check for is_compressed?

    const int ierr = matrix->PutScalar(d);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
    if (nonlocal_matrix.get() != 0)
      nonlocal_matrix->PutScalar(d);

    return *this;
  }



  void
  SparseMatrix::add (const TrilinosScalar  factor,
                     const SparseMatrix   &rhs)
  {
    AssertDimension (rhs.m(), m());
    AssertDimension (rhs.n(), n());
    AssertDimension (rhs.local_range().first, local_range().first);
    AssertDimension (rhs.local_range().second, local_range().second);
    Assert(matrix->RowMap().SameAs(rhs.matrix->RowMap()),
           ExcMessage("Can only add matrices with same distribution of rows"));
    Assert(matrix->Filled() && rhs.matrix->Filled(),
           ExcMessage("Addition of matrices only allowed if matrices are "
                      "filled, i.e., compress() has been called"));

    const std::pair<size_type, size_type>
    local_range = rhs.local_range();
    const bool same_col_map = matrix->ColMap().SameAs(rhs.matrix->ColMap());

    int ierr;
    for (size_type row=local_range.first; row < local_range.second; ++row)
      {
        const int row_local =
          matrix->RowMap().LID(static_cast<TrilinosWrappers::types::int_type>(row));

        // First get a view to the matrix columns of both matrices. Note that
        // the data is in local index spaces so we need to be careful not only
        // to compare column indices in case they are derived from the same
        // map.
        int n_entries, rhs_n_entries;
        TrilinosScalar *value_ptr, *rhs_value_ptr;
        int *index_ptr, *rhs_index_ptr;
        ierr = rhs.matrix->ExtractMyRowView (row_local, rhs_n_entries,
                                             rhs_value_ptr, rhs_index_ptr);
        (void)ierr;
        Assert (ierr == 0, ExcTrilinosError(ierr));

        ierr = matrix->ExtractMyRowView (row_local, n_entries, value_ptr,
                                         index_ptr);
        Assert (ierr == 0, ExcTrilinosError(ierr));
        bool expensive_checks = (n_entries != rhs_n_entries || !same_col_map);
        if (!expensive_checks)
          {
            // check if the column indices are the same. If yes, can simply
            // copy over the data.
            expensive_checks = std::memcmp(static_cast<void *>(index_ptr),
                                           static_cast<void *>(rhs_index_ptr),
                                           sizeof(int)*n_entries) != 0;
            if (!expensive_checks)
              for (int i=0; i<n_entries; ++i)
                value_ptr[i] += rhs_value_ptr[i] * factor;
          }
        // Now to the expensive case where we need to check all column indices
        // against each other (transformed into global index space) and where
        // we need to make sure that all entries we are about to add into the
        // lhs matrix actually exist
        if (expensive_checks)
          {
            for (int i=0; i<rhs_n_entries; ++i)
              {
                if (rhs_value_ptr[i] == 0.)
                  continue;
                const TrilinosWrappers::types::int_type rhs_global_col =
                  global_column_index(*rhs.matrix, rhs_index_ptr[i]);
                int local_col = matrix->ColMap().LID(rhs_global_col);
                int *local_index = Utilities::lower_bound(index_ptr,
                                                          index_ptr+n_entries,
                                                          local_col);
                Assert(local_index != index_ptr + n_entries &&
                       *local_index == local_col,
                       ExcMessage("Adding the entries from the other matrix "
                                  "failed, because the sparsity pattern "
                                  "of that matrix includes more elements than the "
                                  "calling matrix, which is not allowed."));
                value_ptr[local_index-index_ptr] += factor * rhs_value_ptr[i];
              }
          }
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



  SparseMatrix &
  SparseMatrix::operator *= (const TrilinosScalar a)
  {
    const int ierr = matrix->Scale (a);
    Assert (ierr == 0, ExcTrilinosError(ierr));
    (void)ierr; // removes -Wunused-variable in optimized mode

    return *this;
  }



  SparseMatrix &
  SparseMatrix::operator /= (const TrilinosScalar a)
  {
    Assert (a !=0, ExcDivideByZero());

    const TrilinosScalar factor = 1./a;

    const int ierr = matrix->Scale (factor);
    Assert (ierr == 0, ExcTrilinosError(ierr));
    (void)ierr; // removes -Wunused-variable in optimized mode

    return *this;
  }



  TrilinosScalar
  SparseMatrix::l1_norm () const
  {
    Assert (matrix->Filled(), ExcMatrixNotCompressed());
    return matrix->NormOne();
  }



  TrilinosScalar
  SparseMatrix::linfty_norm () const
  {
    Assert (matrix->Filled(), ExcMatrixNotCompressed());
    return matrix->NormInf();
  }



  TrilinosScalar
  SparseMatrix::frobenius_norm () const
  {
    Assert (matrix->Filled(), ExcMatrixNotCompressed());
    return matrix->NormFrobenius();
  }



  namespace internal
  {
    namespace SparseMatrix
    {
      template <typename VectorType>
      inline
      void check_vector_map_equality(const Epetra_CrsMatrix &,
                                     const VectorType &,
                                     const VectorType &)
      {
      }

      inline
      void check_vector_map_equality(const Epetra_CrsMatrix              &m,
                                     const TrilinosWrappers::MPI::Vector &in,
                                     const TrilinosWrappers::MPI::Vector &out)
      {
        Assert (in.vector_partitioner().SameAs(m.DomainMap()) == true,
                ExcMessage ("Column map of matrix does not fit with vector map!"));
        Assert (out.vector_partitioner().SameAs(m.RangeMap()) == true,
                ExcMessage ("Row map of matrix does not fit with vector map!"));
        (void)m;
        (void)in;
        (void)out;
      }
    }
  }


  template <typename VectorType>
  void
  SparseMatrix::vmult (VectorType       &dst,
                       const VectorType &src) const
  {
    Assert (&src != &dst, ExcSourceEqualsDestination());
    Assert (matrix->Filled(), ExcMatrixNotCompressed());
    (void)src;
    (void)dst;

    internal::SparseMatrix::check_vector_map_equality(*matrix, src, dst);
    const size_type dst_local_size = dst.end() - dst.begin();
    AssertDimension (dst_local_size, static_cast<size_type>(matrix->RangeMap().NumMyPoints()));
    const size_type src_local_size = src.end() - src.begin();
    AssertDimension (src_local_size, static_cast<size_type>(matrix->DomainMap().NumMyPoints()));

    Epetra_MultiVector tril_dst (View, matrix->RangeMap(), dst.begin(),
                                 dst_local_size, 1);
    Epetra_MultiVector tril_src (View, matrix->DomainMap(),
                                 const_cast<TrilinosScalar *>(src.begin()),
                                 src_local_size, 1);

    const int ierr = matrix->Multiply (false, tril_src, tril_dst);
    Assert (ierr == 0, ExcTrilinosError(ierr));
    (void)ierr; // removes -Wunused-variable in optimized mode
  }



  template <typename VectorType>
  void
  SparseMatrix::Tvmult (VectorType       &dst,
                        const VectorType &src) const
  {
    Assert (&src != &dst, ExcSourceEqualsDestination());
    Assert (matrix->Filled(), ExcMatrixNotCompressed());

    internal::SparseMatrix::check_vector_map_equality(*matrix, dst, src);
    const size_type dst_local_size = dst.end() - dst.begin();
    AssertDimension (dst_local_size, static_cast<size_type>(matrix->DomainMap().NumMyPoints()));
    const size_type src_local_size = src.end() - src.begin();
    AssertDimension (src_local_size, static_cast<size_type>(matrix->RangeMap().NumMyPoints()));

    Epetra_MultiVector tril_dst (View, matrix->DomainMap(), dst.begin(),
                                 dst_local_size, 1);
    Epetra_MultiVector tril_src (View, matrix->RangeMap(),
                                 const_cast<double *>(src.begin()),
                                 src_local_size, 1);

    const int ierr = matrix->Multiply (true, tril_src, tril_dst);
    Assert (ierr == 0, ExcTrilinosError(ierr));
    (void)ierr; // removes -Wunused-variable in optimized mode
  }



  template <typename VectorType>
  void
  SparseMatrix::vmult_add (VectorType       &dst,
                           const VectorType &src) const
  {
    Assert (&src != &dst, ExcSourceEqualsDestination());

    // Reinit a temporary vector with fast argument set, which does not
    // overwrite the content (to save time). However, the
    // TrilinosWrappers::Vector classes do not support this, so create a
    // deal.II local vector that has this fast setting. It will be accepted in
    // vmult because it only checks the local size.
    dealii::Vector<TrilinosScalar> temp_vector;
    temp_vector.reinit(dst.end()-dst.begin(), true);
    dealii::VectorView<TrilinosScalar> src_view(src.end()-src.begin(), src.begin());
    dealii::VectorView<TrilinosScalar> dst_view(dst.end()-dst.begin(), dst.begin());
    vmult (temp_vector, static_cast<const dealii::Vector<TrilinosScalar>&>(src_view));
    if (dst_view.size() > 0)
      dst_view += temp_vector;
  }



  template <typename VectorType>
  void
  SparseMatrix::Tvmult_add (VectorType       &dst,
                            const VectorType &src) const
  {
    Assert (&src != &dst, ExcSourceEqualsDestination());

    // Reinit a temporary vector with fast argument set, which does not
    // overwrite the content (to save time). However, the
    // TrilinosWrappers::Vector classes do not support this, so create a
    // deal.II local vector that has this fast setting. It will be accepted in
    // vmult because it only checks the local size.
    dealii::Vector<TrilinosScalar> temp_vector;
    temp_vector.reinit(dst.end()-dst.begin(), true);
    dealii::VectorView<TrilinosScalar> src_view(src.end()-src.begin(), src.begin());
    dealii::VectorView<TrilinosScalar> dst_view(dst.end()-dst.begin(), dst.begin());
    Tvmult (temp_vector, static_cast<const dealii::Vector<TrilinosScalar>&>(src_view));
    if (dst_view.size() > 0)
      dst_view += temp_vector;
  }



  TrilinosScalar
  SparseMatrix::matrix_norm_square (const VectorBase &v) const
  {
    Assert (matrix->RowMap().SameAs(matrix->DomainMap()),
            ExcNotQuadratic());

    VectorBase temp_vector;
    temp_vector.reinit(v, true);

    vmult (temp_vector, v);
    return temp_vector*v;
  }



  TrilinosScalar
  SparseMatrix::matrix_scalar_product (const VectorBase &u,
                                       const VectorBase &v) const
  {
    Assert (matrix->RowMap().SameAs(matrix->DomainMap()),
            ExcNotQuadratic());

    VectorBase temp_vector;
    temp_vector.reinit(v, true);

    vmult (temp_vector, v);
    return u*temp_vector;
  }



  TrilinosScalar
  SparseMatrix::residual (VectorBase       &dst,
                          const VectorBase &x,
                          const VectorBase &b) const
  {
    vmult (dst, x);
    dst -= b;
    dst *= -1.;

    return dst.l2_norm();
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
              const size_type GID = inputleft.trilinos_matrix().RowMap().GID(i);
              for (TrilinosWrappers::types::int_type j=0; j<num_entries; ++j)
                sparsity_transposed.add (inputleft.trilinos_matrix().ColMap().GID(indices[j]),
                                         GID);
#else
              const size_type GID = inputleft.trilinos_matrix().RowMap().GID64(i);
              for (TrilinosWrappers::types::int_type j=0; j<num_entries; ++j)
                sparsity_transposed.add (inputleft.trilinos_matrix().ColMap().GID64(indices[j]),
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
              const size_type GID = inputleft.trilinos_matrix().RowMap().GID(i);
              for (TrilinosWrappers::types::int_type j=0; j<num_entries; ++j)
                transposed_mat.set (inputleft.trilinos_matrix().ColMap().GID(indices[j]),
                                    GID, values[j]);
#else
              const size_type GID = inputleft.trilinos_matrix().RowMap().GID64(i);
              for (TrilinosWrappers::types::int_type j=0; j<num_entries; ++j)
                transposed_mat.set (inputleft.trilinos_matrix().ColMap().GID64(indices[j]),
                                    GID, values[j]);
#endif
            }
          transposed_mat.compress(VectorOperation::insert);
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
      C_mat->FillComplete(mod_B->DomainMap(),
                          transpose_left ?
                          inputleft.trilinos_matrix().DomainMap() :
                          inputleft.trilinos_matrix().RangeMap());
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



  const Epetra_Map &
  SparseMatrix::domain_partitioner () const
  {
    return matrix->DomainMap();
  }



  const Epetra_Map &
  SparseMatrix::range_partitioner () const
  {
    return matrix->RangeMap();
  }



  const Epetra_Map &
  SparseMatrix::row_partitioner () const
  {
    return matrix->RowMap();
  }



  const Epetra_Map &
  SparseMatrix::col_partitioner () const
  {
    return matrix->ColMap();
  }



  MPI_Comm SparseMatrix::get_mpi_communicator () const
  {

#ifdef DEAL_II_WITH_MPI

    const Epetra_MpiComm *mpi_comm
      = dynamic_cast<const Epetra_MpiComm *>(&matrix->RangeMap().Comm());
    return mpi_comm->Comm();
#else

    return MPI_COMM_SELF;

#endif

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
  SparseMatrix::reinit (const DynamicSparsityPattern &);

  template void
  SparseMatrix::reinit (const Epetra_Map &,
                        const dealii::SparsityPattern &,
                        const bool);
  template void
  SparseMatrix::reinit (const Epetra_Map &,
                        const DynamicSparsityPattern &,
                        const bool);
  template void
  SparseMatrix::reinit (const Epetra_Map &,
                        const Epetra_Map &,
                        const dealii::SparsityPattern &,
                        const bool);
  template void
  SparseMatrix::reinit (const Epetra_Map &,
                        const Epetra_Map &,
                        const DynamicSparsityPattern &,
                        const bool);
  template void
  SparseMatrix::reinit (const IndexSet &,
                        const IndexSet &,
                        const dealii::SparsityPattern &,
                        const MPI_Comm &,
                        const bool);
  template void
  SparseMatrix::reinit (const IndexSet &,
                        const IndexSet &,
                        const DynamicSparsityPattern &,
                        const MPI_Comm &,
                        const bool);

  template void
  SparseMatrix::vmult (VectorBase &,
                       const VectorBase &) const;
  template void
  SparseMatrix::vmult (Vector &,
                       const Vector &) const;
  template void
  SparseMatrix::vmult (MPI::Vector &,
                       const MPI::Vector &) const;
  template void
  SparseMatrix::vmult (dealii::Vector<double> &,
                       const dealii::Vector<double> &) const;
  template void
  SparseMatrix::vmult (dealii::parallel::distributed::Vector<double> &,
                       const dealii::parallel::distributed::Vector<double> &) const;
  template void
  SparseMatrix::Tvmult (VectorBase &,
                        const VectorBase &) const;
  template void
  SparseMatrix::Tvmult (Vector &,
                        const Vector &) const;
  template void
  SparseMatrix::Tvmult (MPI::Vector &,
                        const MPI::Vector &) const;
  template void
  SparseMatrix::Tvmult (dealii::Vector<double> &,
                        const dealii::Vector<double> &) const;
  template void
  SparseMatrix::Tvmult (dealii::parallel::distributed::Vector<double> &,
                        const dealii::parallel::distributed::Vector<double> &) const;
  template void
  SparseMatrix::vmult_add (VectorBase &,
                           const VectorBase &) const;
  template void
  SparseMatrix::vmult_add (Vector &,
                           const Vector &) const;
  template void
  SparseMatrix::vmult_add (MPI::Vector &,
                           const MPI::Vector &) const;
  template void
  SparseMatrix::vmult_add (dealii::Vector<double> &,
                           const dealii::Vector<double> &) const;
  template void
  SparseMatrix::vmult_add (dealii::parallel::distributed::Vector<double> &,
                           const dealii::parallel::distributed::Vector<double> &) const;
  template void
  SparseMatrix::Tvmult_add (VectorBase &,
                            const VectorBase &) const;
  template void
  SparseMatrix::Tvmult_add (Vector &,
                            const Vector &) const;
  template void
  SparseMatrix::Tvmult_add (MPI::Vector &,
                            const MPI::Vector &) const;
  template void
  SparseMatrix::Tvmult_add (dealii::Vector<double> &,
                            const dealii::Vector<double> &) const;
  template void
  SparseMatrix::Tvmult_add (dealii::parallel::distributed::Vector<double> &,
                            const dealii::parallel::distributed::Vector<double> &) const;
}

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_TRILINOS
