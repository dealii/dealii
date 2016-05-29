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

#include <deal.II/lac/trilinos_sparsity_pattern.h>

#ifdef DEAL_II_WITH_TRILINOS

#  include <deal.II/base/utilities.h>
#  include <deal.II/base/mpi.h>
#  include <deal.II/lac/sparsity_pattern.h>
#  include <deal.II/lac/dynamic_sparsity_pattern.h>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#  include <Epetra_Export.h>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

DEAL_II_NAMESPACE_OPEN

namespace TrilinosWrappers
{
  namespace
  {
    // define a helper function that queries the size of an Epetra_Map object
    // by calling either the 32- or 64-bit function necessary, and returns the
    // result in the correct data type so that we can use it in calling other
    // Epetra member functions that are overloaded by index type
#ifndef DEAL_II_WITH_64BIT_INDICES
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

    int n_global_rows(const Epetra_CrsGraph &graph)
    {
      return graph.NumGlobalRows();
    }

    int n_global_cols(const Epetra_CrsGraph &graph)
    {
      return graph.NumGlobalCols();
    }

    int n_global_entries(const Epetra_CrsGraph &graph)
    {
      return graph.NumGlobalEntries();
    }

    int global_row_index(const Epetra_CrsGraph &graph, int i)
    {
      return graph.GRID(i);
    }
#else
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

    long long int n_global_rows(const Epetra_CrsGraph &graph)
    {
      return graph.NumGlobalRows64();
    }

    long long int n_global_cols(const Epetra_CrsGraph &graph)
    {
      return graph.NumGlobalCols64();
    }

    long long int n_global_entries(const Epetra_CrsGraph &graph)
    {
      return graph.NumGlobalEntries64();
    }

    long long int global_row_index(const Epetra_CrsGraph &graph, int i)
    {
      return graph.GRID64(i);
    }
#endif
  }

  namespace SparsityPatternIterators
  {
    void
    Accessor::visit_present_row ()
    {
      // if we are asked to visit the
      // past-the-end line, then simply
      // release all our caches and go on
      // with life
      if (this->a_row == sparsity_pattern->n_rows())
        {
          colnum_cache.reset ();

          return;
        }
//TODO: Is this thread safe?

      // otherwise first flush Trilinos caches
      sparsity_pattern->compress ();

      // get a representation of the present
      // row
      int ncols;
      // TODO: casting a size_type to an int, could be a problem
      int colnums = sparsity_pattern->n_cols();

      int ierr;
      ierr = sparsity_pattern->graph->ExtractGlobalRowCopy((TrilinosWrappers::types::int_type)this->a_row,
                                                           colnums,
                                                           ncols,
                                                           (TrilinosWrappers::types::int_type *)&(*colnum_cache)[0]);
      AssertThrow (ierr == 0, ExcTrilinosError(ierr));

      // copy it into our caches if the
      // line isn't empty. if it is, then
      // we've done something wrong, since
      // we shouldn't have initialized an
      // iterator for an empty line (what
      // would it point to?)
      Assert (ncols != 0, ExcInternalError());
      colnum_cache.reset (new std::vector<size_type> (colnums,
                                                      colnums+ncols));
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
  SparsityPattern::SparsityPattern ()
  {
    column_space_map.reset(new Epetra_Map (TrilinosWrappers::types::int_type(0),
                                           TrilinosWrappers::types::int_type(0),
                                           Utilities::Trilinos::comm_self()));
    graph.reset (new Epetra_FECrsGraph(View,
                                       *column_space_map,
                                       *column_space_map,
                                       0));
    graph->FillComplete();
  }


  SparsityPattern::SparsityPattern (const Epetra_Map  &input_map,
                                    const size_type n_entries_per_row)
  {
    reinit (input_map, input_map, n_entries_per_row);
  }



  SparsityPattern::SparsityPattern (const Epetra_Map             &input_map,
                                    const std::vector<size_type> &n_entries_per_row)
  {
    reinit (input_map, input_map, n_entries_per_row);
  }



  SparsityPattern::SparsityPattern (const Epetra_Map  &input_row_map,
                                    const Epetra_Map  &input_col_map,
                                    const size_type n_entries_per_row)
  {
    reinit (input_row_map, input_col_map, n_entries_per_row);
  }



  SparsityPattern::SparsityPattern (const Epetra_Map                &input_row_map,
                                    const Epetra_Map                &input_col_map,
                                    const std::vector<size_type> &n_entries_per_row)
  {
    reinit (input_row_map, input_col_map, n_entries_per_row);
  }



  SparsityPattern::SparsityPattern (const size_type m,
                                    const size_type n,
                                    const size_type n_entries_per_row)
  {
    reinit (m, n, n_entries_per_row);
  }



  SparsityPattern::SparsityPattern (const size_type               m,
                                    const size_type               n,
                                    const std::vector<size_type> &n_entries_per_row)
  {
    reinit (m, n, n_entries_per_row);
  }


  // Copy function only works if the
  // sparsity pattern is empty.
  SparsityPattern::SparsityPattern (const SparsityPattern &input_sparsity)
    :
    Subscriptor(),
    column_space_map (new Epetra_Map(TrilinosWrappers::types::int_type(0),
                                     TrilinosWrappers::types::int_type(0),
                                     Utilities::Trilinos::comm_self())),
    graph (new Epetra_FECrsGraph(View,
                                 *column_space_map,
                                 *column_space_map,
                                 0))
  {
    (void)input_sparsity;
    Assert (input_sparsity.n_rows() == 0,
            ExcMessage ("Copy constructor only works for empty sparsity patterns."));
  }



  SparsityPattern::SparsityPattern  (const IndexSet  &parallel_partitioning,
                                     const MPI_Comm  &communicator,
                                     const size_type  n_entries_per_row)
  {
    reinit (parallel_partitioning, parallel_partitioning, communicator,
            n_entries_per_row);
  }



  SparsityPattern::SparsityPattern  (const IndexSet     &parallel_partitioning,
                                     const MPI_Comm     &communicator,
                                     const std::vector<size_type> &n_entries_per_row)
  {
    reinit (parallel_partitioning, parallel_partitioning, communicator,
            n_entries_per_row);
  }



  SparsityPattern::SparsityPattern  (const IndexSet  &row_parallel_partitioning,
                                     const IndexSet  &col_parallel_partitioning,
                                     const MPI_Comm  &communicator,
                                     const size_type  n_entries_per_row)
  {
    reinit (row_parallel_partitioning, col_parallel_partitioning,
            communicator, n_entries_per_row);
  }



  SparsityPattern::
  SparsityPattern  (const IndexSet     &row_parallel_partitioning,
                    const IndexSet     &col_parallel_partitioning,
                    const MPI_Comm     &communicator,
                    const std::vector<size_type> &n_entries_per_row)
  {
    reinit (row_parallel_partitioning, col_parallel_partitioning,
            communicator, n_entries_per_row);
  }



  SparsityPattern::
  SparsityPattern  (const IndexSet     &row_parallel_partitioning,
                    const IndexSet     &col_parallel_partitioning,
                    const IndexSet     &writable_rows,
                    const MPI_Comm     &communicator,
                    const size_type     n_max_entries_per_row)
  {
    reinit (row_parallel_partitioning, col_parallel_partitioning,
            writable_rows, communicator, n_max_entries_per_row);
  }



  SparsityPattern::~SparsityPattern ()
  {}



  void
  SparsityPattern::reinit (const size_type  m,
                           const size_type  n,
                           const size_type  n_entries_per_row)
  {
    reinit (complete_index_set(m), complete_index_set(n), MPI_COMM_SELF,
            n_entries_per_row);
  }



  void
  SparsityPattern::reinit (const size_type  m,
                           const size_type  n,
                           const std::vector<size_type> &n_entries_per_row)
  {
    reinit (complete_index_set(m), complete_index_set(n), MPI_COMM_SELF,
            n_entries_per_row);
  }



  namespace
  {
    typedef SparsityPattern::size_type size_type;

    void
    reinit_sp (const Epetra_Map                         &row_map,
               const Epetra_Map                         &col_map,
               const size_type                           n_entries_per_row,
               std_cxx11::shared_ptr<Epetra_Map>        &column_space_map,
               std_cxx11::shared_ptr<Epetra_FECrsGraph> &graph,
               std_cxx11::shared_ptr<Epetra_CrsGraph>   &nonlocal_graph)
    {
      Assert(row_map.IsOneToOne(),
             ExcMessage("Row map must be 1-to-1, i.e., no overlap between "
                        "the maps of different processors."));
      Assert(col_map.IsOneToOne(),
             ExcMessage("Column map must be 1-to-1, i.e., no overlap between "
                        "the maps of different processors."));

      nonlocal_graph.reset();
      graph.reset ();
      column_space_map.reset (new Epetra_Map (col_map));

      // for more than one processor, need to specify only row map first and
      // let the matrix entries decide about the column map (which says which
      // columns are present in the matrix, not to be confused with the
      // col_map that tells how the domain dofs of the matrix will be
      // distributed). for only one processor, we can directly assign the
      // columns as well. If we use a recent Trilinos version, we can also
      // require building a non-local graph which gives us thread-safe
      // initialization.
      if (row_map.Comm().NumProc() > 1)
        graph.reset (new Epetra_FECrsGraph(Copy, row_map,
                                           n_entries_per_row, false
                                           // TODO: Check which new Trilinos
                                           // version supports this... Remember
                                           // to change tests/trilinos/assemble_matrix_parallel_07
                                           // too.
                                           //#if DEAL_II_TRILINOS_VERSION_GTE(11,14,0)
                                           //, true
                                           //#endif
                                          ));
      else
        graph.reset (new Epetra_FECrsGraph(Copy, row_map, col_map,
                                           n_entries_per_row, false));
    }



    void
    reinit_sp (const Epetra_Map                         &row_map,
               const Epetra_Map                         &col_map,
               const std::vector<size_type>             &n_entries_per_row,
               std_cxx11::shared_ptr<Epetra_Map>        &column_space_map,
               std_cxx11::shared_ptr<Epetra_FECrsGraph> &graph,
               std_cxx11::shared_ptr<Epetra_CrsGraph>   &nonlocal_graph)
    {
      Assert(row_map.IsOneToOne(),
             ExcMessage("Row map must be 1-to-1, i.e., no overlap between "
                        "the maps of different processors."));
      Assert(col_map.IsOneToOne(),
             ExcMessage("Column map must be 1-to-1, i.e., no overlap between "
                        "the maps of different processors."));

      // release memory before reallocation
      nonlocal_graph.reset();
      graph.reset ();
      AssertDimension (n_entries_per_row.size(),
                       static_cast<size_type>(n_global_elements(row_map)));

      column_space_map.reset (new Epetra_Map (col_map));
      std::vector<int> local_entries_per_row(max_my_gid(row_map)-
                                             min_my_gid(row_map));
      for (unsigned int i=0; i<local_entries_per_row.size(); ++i)
        local_entries_per_row[i] = n_entries_per_row[min_my_gid(row_map)+i];

      if (row_map.Comm().NumProc() > 1)
        graph.reset(new Epetra_FECrsGraph(Copy, row_map,
                                          &local_entries_per_row[0],
                                          false
                                          // TODO: Check which new Trilinos
                                          // version supports this... Remember
                                          // to change tests/trilinos/assemble_matrix_parallel_07
                                          // too.
                                          //#if DEAL_II_TRILINOS_VERSION_GTE(11,14,0)
                                          //, true
                                          //#endif
                                         ));
      else
        graph.reset(new Epetra_FECrsGraph(Copy, row_map, col_map,
                                          &local_entries_per_row[0],
                                          false));
    }



    template <typename SparsityPatternType>
    void
    reinit_sp (const Epetra_Map                         &row_map,
               const Epetra_Map                         &col_map,
               const SparsityPatternType                &sp,
               const bool                                exchange_data,
               std_cxx11::shared_ptr<Epetra_Map>        &column_space_map,
               std_cxx11::shared_ptr<Epetra_FECrsGraph> &graph,
               std_cxx11::shared_ptr<Epetra_CrsGraph>   &nonlocal_graph)
    {
      nonlocal_graph.reset ();
      graph.reset ();

      AssertDimension (sp.n_rows(),
                       static_cast<size_type>(n_global_elements(row_map)));
      AssertDimension (sp.n_cols(),
                       static_cast<size_type>(n_global_elements(col_map)));

      column_space_map.reset (new Epetra_Map (col_map));

      Assert (row_map.LinearMap() == true,
              ExcMessage ("This function only works if the row map is contiguous."));

      const size_type first_row = min_my_gid(row_map),
                      last_row = max_my_gid(row_map)+1;
      std::vector<int> n_entries_per_row(last_row - first_row);

      // Trilinos wants the row length as an int this is hopefully never going
      // to be a problem.
      for (size_type row=first_row; row<last_row; ++row)
        n_entries_per_row[row-first_row] = static_cast<int>(sp.row_length(row));

      if (row_map.Comm().NumProc() > 1)
        graph.reset(new Epetra_FECrsGraph(Copy, row_map,
                                          &n_entries_per_row[0],
                                          false));
      else
        graph.reset (new Epetra_FECrsGraph(Copy, row_map, col_map,
                                           &n_entries_per_row[0],
                                           false));

      AssertDimension (sp.n_rows(),
                       static_cast<size_type>(n_global_rows(*graph)));

      std::vector<TrilinosWrappers::types::int_type> row_indices;

      // Include possibility to exchange data since DynamicSparsityPattern is
      // able to do so
      if (exchange_data==false)
        for (size_type row=first_row; row<last_row; ++row)
          {
            const TrilinosWrappers::types::int_type row_length = sp.row_length(row);
            if (row_length == 0)
              continue;

            row_indices.resize (row_length, -1);
            {
              typename SparsityPatternType::iterator p = sp.begin(row);
              // avoid incrementing p over the end of the current row because
              // it is slow for DynamicSparsityPattern in parallel
              for (int col=0; col<row_length; )
                {
                  row_indices[col++] = p->column();
                  if (col < row_length)
                    ++p;
                }
            }
            graph->Epetra_CrsGraph::InsertGlobalIndices (row, row_length,
                                                         &row_indices[0]);
          }
      else
        for (size_type row=0; row<sp.n_rows(); ++row)
          {
            const TrilinosWrappers::types::int_type row_length = sp.row_length(row);
            if (row_length == 0)
              continue;

            row_indices.resize (row_length, -1);
            {
              typename SparsityPatternType::iterator p = sp.begin(row);
              // avoid incrementing p over the end of the current row because
              // it is slow for DynamicSparsityPattern in parallel
              for (int col=0; col<row_length; )
                {
                  row_indices[col++] = p->column();
                  if (col < row_length)
                    ++p;
                }
            }
            graph->InsertGlobalIndices (1,
                                        reinterpret_cast<TrilinosWrappers::types::int_type *>(&row),
                                        row_length, &row_indices[0]);
          }

      int ierr =
        graph->GlobalAssemble (*column_space_map,
                               static_cast<const Epetra_Map &>(graph->RangeMap()),
                               true);
      AssertThrow (ierr == 0, ExcTrilinosError(ierr));

      ierr = graph->OptimizeStorage ();
      AssertThrow (ierr == 0, ExcTrilinosError(ierr));
    }
  }


  void
  SparsityPattern::reinit (const Epetra_Map  &input_map,
                           const size_type    n_entries_per_row)
  {
    reinit_sp (input_map, input_map, n_entries_per_row,
               column_space_map, graph, nonlocal_graph);
  }



  void
  SparsityPattern::reinit (const Epetra_Map  &input_row_map,
                           const Epetra_Map  &input_col_map,
                           const size_type    n_entries_per_row)
  {
    reinit_sp (input_row_map, input_col_map, n_entries_per_row,
               column_space_map, graph, nonlocal_graph);
  }



  void
  SparsityPattern::reinit (const Epetra_Map   &input_map,
                           const std::vector<size_type> &n_entries_per_row)
  {
    reinit_sp (input_map, input_map, n_entries_per_row,
               column_space_map, graph, nonlocal_graph);
  }



  void
  SparsityPattern::reinit (const Epetra_Map   &input_row_map,
                           const Epetra_Map   &input_col_map,
                           const std::vector<size_type> &n_entries_per_row)
  {
    reinit_sp (input_row_map, input_col_map, n_entries_per_row,
               column_space_map, graph, nonlocal_graph);
  }



  void
  SparsityPattern::reinit (const IndexSet  &parallel_partitioning,
                           const MPI_Comm  &communicator,
                           const size_type  n_entries_per_row)
  {
    Epetra_Map map = parallel_partitioning.make_trilinos_map (communicator,
                                                              false);
    reinit_sp (map, map, n_entries_per_row,
               column_space_map, graph, nonlocal_graph);
  }



  void SparsityPattern::reinit (const IndexSet     &parallel_partitioning,
                                const MPI_Comm     &communicator,
                                const std::vector<size_type> &n_entries_per_row)
  {
    Epetra_Map map = parallel_partitioning.make_trilinos_map (communicator,
                                                              false);
    reinit_sp (map, map, n_entries_per_row,
               column_space_map, graph, nonlocal_graph);
  }



  void SparsityPattern::reinit (const IndexSet &row_parallel_partitioning,
                                const IndexSet &col_parallel_partitioning,
                                const MPI_Comm &communicator,
                                const size_type  n_entries_per_row)
  {
    Epetra_Map row_map =
      row_parallel_partitioning.make_trilinos_map (communicator, false);
    Epetra_Map col_map =
      col_parallel_partitioning.make_trilinos_map (communicator, false);
    reinit_sp (row_map, col_map, n_entries_per_row,
               column_space_map, graph, nonlocal_graph);
  }



  void
  SparsityPattern::reinit (const IndexSet     &row_parallel_partitioning,
                           const IndexSet     &col_parallel_partitioning,
                           const MPI_Comm     &communicator,
                           const std::vector<size_type> &n_entries_per_row)
  {
    Epetra_Map row_map =
      row_parallel_partitioning.make_trilinos_map (communicator, false);
    Epetra_Map col_map =
      col_parallel_partitioning.make_trilinos_map (communicator, false);
    reinit_sp (row_map, col_map, n_entries_per_row,
               column_space_map, graph, nonlocal_graph);
  }



  void
  SparsityPattern::reinit (const IndexSet  &row_parallel_partitioning,
                           const IndexSet  &col_parallel_partitioning,
                           const IndexSet  &writable_rows,
                           const MPI_Comm  &communicator,
                           const size_type  n_entries_per_row)
  {
    Epetra_Map row_map =
      row_parallel_partitioning.make_trilinos_map (communicator, false);
    Epetra_Map col_map =
      col_parallel_partitioning.make_trilinos_map (communicator, false);
    reinit_sp (row_map, col_map, n_entries_per_row,
               column_space_map, graph, nonlocal_graph);

    IndexSet nonlocal_partitioner = writable_rows;
    AssertDimension(nonlocal_partitioner.size(), row_parallel_partitioning.size());
#ifdef DEBUG
    {
      IndexSet tmp = writable_rows & row_parallel_partitioning;
      Assert (tmp == row_parallel_partitioning,
              ExcMessage("The set of writable rows passed to this method does not "
                         "contain the locally owned rows, which is not allowed."));
    }
#endif
    nonlocal_partitioner.subtract_set(row_parallel_partitioning);
    if (Utilities::MPI::n_mpi_processes(communicator) > 1)
      {
        Epetra_Map nonlocal_map =
          nonlocal_partitioner.make_trilinos_map(communicator, true);
        nonlocal_graph.reset(new Epetra_CrsGraph(Copy, nonlocal_map, 0));
      }
    else
      Assert(nonlocal_partitioner.n_elements() == 0, ExcInternalError());
  }



  template<typename SparsityPatternType>
  void
  SparsityPattern::reinit (const IndexSet            &row_parallel_partitioning,
                           const IndexSet            &col_parallel_partitioning,
                           const SparsityPatternType &nontrilinos_sparsity_pattern,
                           const MPI_Comm            &communicator,
                           const bool                 exchange_data)
  {
    Epetra_Map row_map =
      row_parallel_partitioning.make_trilinos_map (communicator, false);
    Epetra_Map col_map =
      col_parallel_partitioning.make_trilinos_map (communicator, false);
    reinit_sp (row_map, col_map, nontrilinos_sparsity_pattern, exchange_data,
               column_space_map, graph, nonlocal_graph);
  }



  template<typename SparsityPatternType>
  void
  SparsityPattern::reinit (const IndexSet            &parallel_partitioning,
                           const SparsityPatternType &nontrilinos_sparsity_pattern,
                           const MPI_Comm            &communicator,
                           const bool                 exchange_data)
  {
    Epetra_Map map = parallel_partitioning.make_trilinos_map (communicator,
                                                              false);
    reinit_sp (map, map, nontrilinos_sparsity_pattern, exchange_data,
               column_space_map, graph, nonlocal_graph);
  }



  template <typename SparsityPatternType>
  void
  SparsityPattern::reinit (const Epetra_Map          &input_map,
                           const SparsityPatternType &sp,
                           const bool                 exchange_data)
  {
    reinit_sp (input_map, input_map, sp, exchange_data,
               column_space_map, graph, nonlocal_graph);
  }



  template <typename SparsityPatternType>
  void
  SparsityPattern::reinit (const Epetra_Map          &input_row_map,
                           const Epetra_Map          &input_col_map,
                           const SparsityPatternType &sp,
                           const bool                 exchange_data)
  {
    reinit_sp (input_row_map, input_col_map, sp, exchange_data,
               column_space_map, graph, nonlocal_graph);

    compress();
  }



  SparsityPattern &
  SparsityPattern::operator = (const SparsityPattern &)
  {
    Assert (false, ExcNotImplemented());
    return *this;
  }



  template <typename SparsityPatternType>
  void
  SparsityPattern::copy_from (const SparsityPatternType &sp)
  {
    const Epetra_Map rows (TrilinosWrappers::types::int_type(sp.n_rows()), 0,
                           Utilities::Trilinos::comm_self());
    const Epetra_Map columns (TrilinosWrappers::types::int_type(sp.n_cols()), 0,
                              Utilities::Trilinos::comm_self());

    reinit_sp (rows, columns, sp, false,
               column_space_map, graph, nonlocal_graph);
  }



  void
  SparsityPattern::clear ()
  {
    // When we clear the matrix, reset
    // the pointer and generate an
    // empty sparsity pattern.
    column_space_map.reset (new Epetra_Map (TrilinosWrappers::types::int_type(0),
                                            TrilinosWrappers::types::int_type(0),
                                            Utilities::Trilinos::comm_self()));
    graph.reset (new Epetra_FECrsGraph(View, *column_space_map,
                                       *column_space_map, 0));
    graph->FillComplete();

    nonlocal_graph.reset();
  }



  void
  SparsityPattern::compress ()
  {
    int ierr;
    Assert (column_space_map.get() != 0, ExcInternalError());
    if (nonlocal_graph.get() != 0)
      {
        if (nonlocal_graph->IndicesAreGlobal() == false &&
            nonlocal_graph->RowMap().NumMyElements() > 0)
          {
            // insert dummy element
            TrilinosWrappers::types::int_type row = nonlocal_graph->RowMap().MyGID(
                                                      static_cast<TrilinosWrappers::types::int_type> (0));
            nonlocal_graph->InsertGlobalIndices(row, 1, &row);
          }
        Assert(nonlocal_graph->RowMap().NumMyElements() == 0 ||
               nonlocal_graph->IndicesAreGlobal() == true,
               ExcInternalError());
        nonlocal_graph->FillComplete(*column_space_map,
                                     static_cast<const Epetra_Map &>(graph->RangeMap()));
        nonlocal_graph->OptimizeStorage();
        Epetra_Export exporter(nonlocal_graph->RowMap(), graph->RowMap());
        ierr = graph->Export(*nonlocal_graph, exporter, Add);
        AssertThrow (ierr == 0, ExcTrilinosError(ierr));
        ierr =
          graph->FillComplete(*column_space_map,
                              static_cast<const Epetra_Map &>(graph->RangeMap()));
      }
    else
      ierr = graph->GlobalAssemble (*column_space_map,
                                    static_cast<const Epetra_Map &>(graph->RangeMap()),
                                    true);

    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    ierr = graph->OptimizeStorage ();
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



  bool
  SparsityPattern::exists (const size_type i,
                           const size_type j) const
  {
    // Extract local indices in
    // the matrix.
    int trilinos_i = graph->LRID(static_cast<TrilinosWrappers::types::int_type>(i)),
        trilinos_j = graph->LCID(static_cast<TrilinosWrappers::types::int_type>(j));

    // If the data is not on the
    // present processor, we throw
    // an exception. This is on of
    // the two tiny differences to
    // the el(i,j) call, which does
    // not throw any assertions.
    if (trilinos_i == -1)
      {
        return false;
      }
    else
      {
        // Check whether the matrix
        // already is transformed to
        // local indices.
        if (graph->Filled() == false)
          {
            int nnz_present = graph->NumGlobalIndices(i);
            int nnz_extracted;
            TrilinosWrappers::types::int_type *col_indices;

            // Generate the view and make
            // sure that we have not generated
            // an error.
            // TODO: trilinos_i is the local row index -> it is an int but
            // ExtractGlobalRowView requires trilinos_i to be the global row
            // index and thus it should be a long long int
            int ierr = graph->ExtractGlobalRowView(
                         static_cast<TrilinosWrappers::types::int_type>(trilinos_i),
                         nnz_extracted, col_indices);
            (void)ierr;
            Assert (ierr==0, ExcTrilinosError(ierr));
            Assert (nnz_present == nnz_extracted,
                    ExcDimensionMismatch(nnz_present, nnz_extracted));

            // Search the index
            TrilinosWrappers::types::int_type *el_find =
              std::find(col_indices, col_indices + nnz_present, trilinos_j);

            TrilinosWrappers::types::int_type local_col_index =
              (TrilinosWrappers::types::int_type)(el_find - col_indices);

            if (local_col_index == nnz_present)
              return false;
          }
        else
          {
            // Prepare pointers for extraction
            // of a view of the row.
            int nnz_present = graph->NumGlobalIndices(
                                static_cast<TrilinosWrappers::types::int_type>(i));
            int nnz_extracted;
            int *col_indices;

            // Generate the view and make
            // sure that we have not generated
            // an error.
            int ierr = graph->ExtractMyRowView(trilinos_i,
                                               nnz_extracted, col_indices);
            (void)ierr;
            Assert (ierr==0, ExcTrilinosError(ierr));

            Assert (nnz_present == nnz_extracted,
                    ExcDimensionMismatch(nnz_present, nnz_extracted));

            // Search the index
            int *el_find = std::find(col_indices, col_indices + nnz_present,
                                     static_cast<int>(trilinos_j));

            int local_col_index = (int)(el_find - col_indices);

            if (local_col_index == nnz_present)
              return false;
          }
      }

    return true;
  }



  SparsityPattern::size_type
  SparsityPattern::bandwidth () const
  {
    size_type local_b=0;
    TrilinosWrappers::types::int_type global_b=0;
    for (int i=0; i<(int)local_size(); ++i)
      {
        int *indices;
        int num_entries;
        graph->ExtractMyRowView(i, num_entries, indices);
        for (unsigned int j=0; j<(unsigned int)num_entries; ++j)
          {
            if (static_cast<size_type>(std::abs(static_cast<TrilinosWrappers::types::int_type>(i-indices[j]))) > local_b)
              local_b = std::abs(static_cast<TrilinosWrappers::types::int_type>(i-indices[j]));
          }
      }
    graph->Comm().MaxAll((TrilinosWrappers::types::int_type *)&local_b, &global_b, 1);
    return static_cast<size_type>(global_b);
  }



  SparsityPattern::size_type
  SparsityPattern::n_rows () const
  {
    const TrilinosWrappers::types::int_type n_rows = n_global_rows(*graph);
    return n_rows;
  }



  SparsityPattern::size_type
  SparsityPattern::n_cols () const
  {
    TrilinosWrappers::types::int_type n_cols;
    if (graph->Filled() == true)
      n_cols = n_global_cols(*graph);
    else
      n_cols = n_global_elements(*column_space_map);

    return n_cols;
  }



  unsigned int
  SparsityPattern::local_size () const
  {
    int n_rows = graph -> NumMyRows();

    return n_rows;
  }



  std::pair<SparsityPattern::size_type, SparsityPattern::size_type>
  SparsityPattern::local_range () const
  {
    size_type begin, end;
    begin =  min_my_gid(graph->RowMap());
    end = max_my_gid(graph->RowMap())+1;

    return std::make_pair (begin, end);
  }



  SparsityPattern::size_type
  SparsityPattern::n_nonzero_elements () const
  {
    TrilinosWrappers::types::int_type nnz = n_global_entries(*graph);

    return static_cast<size_type>(nnz);
  }



  unsigned int
  SparsityPattern::max_entries_per_row () const
  {
    int nnz = graph->MaxNumIndices();

    return static_cast<unsigned int>(nnz);
  }



  SparsityPattern::size_type
  SparsityPattern::row_length (const size_type row) const
  {
    Assert (row < n_rows(), ExcInternalError());

    // get a representation of the
    // present row
    TrilinosWrappers::types::int_type ncols = -1;
    TrilinosWrappers::types::int_type local_row =
      graph->LRID(static_cast<TrilinosWrappers::types::int_type>(row));

    // on the processor who owns this
    // row, we'll have a non-negative
    // value.
    if (local_row >= 0)
      ncols = graph->NumMyIndices (local_row);

    return static_cast<size_type>(ncols);
  }



  const Epetra_Map &
  SparsityPattern::domain_partitioner () const
  {
    return static_cast<const Epetra_Map &>(graph->DomainMap());
  }



  const Epetra_Map &
  SparsityPattern::range_partitioner () const
  {
    return static_cast<const Epetra_Map &>(graph->RangeMap());
  }



  const Epetra_Map &
  SparsityPattern::row_partitioner () const
  {
    return static_cast<const Epetra_Map &>(graph->RowMap());
  }



  const Epetra_Map &
  SparsityPattern::col_partitioner () const
  {
    return static_cast<const Epetra_Map &>(graph->ColMap());
  }



  const Epetra_Comm &
  SparsityPattern::trilinos_communicator () const
  {
    return graph->RangeMap().Comm();
  }



  MPI_Comm
  SparsityPattern::get_mpi_communicator () const
  {

#ifdef DEAL_II_WITH_MPI

    const Epetra_MpiComm *mpi_comm
      = dynamic_cast<const Epetra_MpiComm *>(&graph->RangeMap().Comm());
    return mpi_comm->Comm();
#else

    return MPI_COMM_SELF;

#endif

  }



  void
  SparsityPattern::write_ascii ()
  {
    Assert (false, ExcNotImplemented());
  }



  // As of now, no particularly neat
  // ouput is generated in case of
  // multiple processors.
  void
  SparsityPattern::print (std::ostream &out,
                          const bool    write_extended_trilinos_info) const
  {
    if (write_extended_trilinos_info)
      out << *graph;
    else
      {
        int *indices;
        int num_entries;

        for (int i=0; i<graph->NumMyRows(); ++i)
          {
            graph->ExtractMyRowView (i, num_entries, indices);
            for (int j=0; j<num_entries; ++j)
              out << "(" << i << "," << indices[global_row_index(*graph,j)] << ") "
                  << std::endl;
          }
      }

    AssertThrow (out, ExcIO());
  }



  void
  SparsityPattern::print_gnuplot (std::ostream &out) const
  {
    Assert (graph->Filled() == true, ExcInternalError());
    for (unsigned int row=0; row<local_size(); ++row)
      {
        int *indices;
        int num_entries;
        graph->ExtractMyRowView (row, num_entries, indices);

        for (unsigned int j=0; j<(unsigned int)num_entries; ++j)
          // while matrix entries are usually
          // written (i,j), with i vertical and
          // j horizontal, gnuplot output is
          // x-y, that is we have to exchange
          // the order of output
          out << indices[global_row_index(*graph,static_cast<int>(j))]
              << " " << -static_cast<signed int>(row) << std::endl;
      }

    AssertThrow (out, ExcIO());
  }

//TODO: Implement!
  std::size_t
  SparsityPattern::memory_consumption() const
  {
    Assert(false, ExcNotImplemented());
    return 0;
  }


  // explicit instantiations
  //
  template void
  SparsityPattern::copy_from (const dealii::SparsityPattern &);
  template void
  SparsityPattern::copy_from (const dealii::DynamicSparsityPattern &);


  template void
  SparsityPattern::reinit (const Epetra_Map &,
                           const dealii::SparsityPattern &,
                           bool);
  template void
  SparsityPattern::reinit (const Epetra_Map &,
                           const dealii::DynamicSparsityPattern &,
                           bool);

  template void
  SparsityPattern::reinit (const Epetra_Map &,
                           const Epetra_Map &,
                           const dealii::SparsityPattern &,
                           bool);
  template void
  SparsityPattern::reinit (const Epetra_Map &,
                           const Epetra_Map &,
                           const dealii::DynamicSparsityPattern &,
                           bool);


  template void
  SparsityPattern::reinit (const IndexSet &,
                           const dealii::SparsityPattern &,
                           const MPI_Comm &,
                           bool);
  template void
  SparsityPattern::reinit (const IndexSet &,
                           const dealii::DynamicSparsityPattern &,
                           const MPI_Comm &,
                           bool);


  template void
  SparsityPattern::reinit (const IndexSet &,
                           const IndexSet &,
                           const dealii::SparsityPattern &,
                           const MPI_Comm &,
                           bool);
  template void
  SparsityPattern::reinit (const IndexSet &,
                           const IndexSet &,
                           const dealii::DynamicSparsityPattern &,
                           const MPI_Comm &,
                           bool);

}

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_TRILINOS
