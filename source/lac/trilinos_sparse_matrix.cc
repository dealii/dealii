// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2020 by the deal.II authors
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
#include <deal.II/lac/trilinos_sparse_matrix.h>

#ifdef DEAL_II_WITH_TRILINOS

#  include <deal.II/base/utilities.h>

#  include <deal.II/lac/dynamic_sparsity_pattern.h>
#  include <deal.II/lac/la_parallel_vector.h>
#  include <deal.II/lac/sparse_matrix.h>
#  include <deal.II/lac/sparsity_pattern.h>
#  include <deal.II/lac/sparsity_tools.h>
#  include <deal.II/lac/trilinos_precondition.h>
#  include <deal.II/lac/trilinos_sparsity_pattern.h>

#  include <boost/container/small_vector.hpp>

#  ifdef DEAL_II_TRILINOS_WITH_EPETRAEXT
#    include <EpetraExt_MatrixMatrix.h>
#  endif
#  include <Epetra_Export.h>
#  include <Teuchos_RCP.hpp>
#  include <ml_epetra_utils.h>
#  include <ml_struct.h>

#  include <memory>

DEAL_II_NAMESPACE_OPEN

namespace TrilinosWrappers
{
  namespace internal
  {
    template <typename VectorType>
    typename VectorType::value_type *
    begin(VectorType &V)
    {
      return V.begin();
    }

    template <typename VectorType>
    const typename VectorType::value_type *
    begin(const VectorType &V)
    {
      return V.begin();
    }

    template <typename VectorType>
    typename VectorType::value_type *
    end(VectorType &V)
    {
      return V.end();
    }

    template <typename VectorType>
    const typename VectorType::value_type *
    end(const VectorType &V)
    {
      return V.end();
    }

#  ifdef DEAL_II_WITH_MPI
    template <>
    double *
    begin(LinearAlgebra::EpetraWrappers::Vector &V)
    {
      return V.trilinos_vector()[0];
    }

    template <>
    const double *
    begin(const LinearAlgebra::EpetraWrappers::Vector &V)
    {
      return V.trilinos_vector()[0];
    }

    template <>
    double *
    end(LinearAlgebra::EpetraWrappers::Vector &V)
    {
      return V.trilinos_vector()[0] + V.trilinos_vector().MyLength();
    }

    template <>
    const double *
    end(const LinearAlgebra::EpetraWrappers::Vector &V)
    {
      return V.trilinos_vector()[0] + V.trilinos_vector().MyLength();
    }

#    ifdef DEAL_II_TRILINOS_WITH_TPETRA
    template <typename Number>
    Number *
    begin(LinearAlgebra::TpetraWrappers::Vector<Number> &V)
    {
      return V.trilinos_vector().getDataNonConst().get();
    }

    template <typename Number>
    const Number *
    begin(const LinearAlgebra::TpetraWrappers::Vector<Number> &V)
    {
      return V.trilinos_vector().getData().get();
    }

    template <typename Number>
    Number *
    end(LinearAlgebra::TpetraWrappers::Vector<Number> &V)
    {
      return V.trilinos_vector().getDataNonConst().get() +
             V.trilinos_vector().getLocalLength();
    }

    template <typename Number>
    const Number *
    end(const LinearAlgebra::TpetraWrappers::Vector<Number> &V)
    {
      return V.trilinos_vector().getData().get() +
             V.trilinos_vector().getLocalLength();
    }
#    endif
#  endif
  } // namespace internal


  namespace SparseMatrixIterators
  {
    void
    AccessorBase::visit_present_row()
    {
      // if we are asked to visit the past-the-end line, then simply
      // release all our caches and go on with life.
      //
      // do the same if the row we're supposed to visit is not locally
      // owned. this is simply going to make non-locally owned rows
      // look like they're empty
      if ((this->a_row == matrix->m()) ||
          (matrix->in_local_range(this->a_row) == false))
        {
          colnum_cache.reset();
          value_cache.reset();

          return;
        }

      // get a representation of the present row
      int                               ncols;
      TrilinosWrappers::types::int_type colnums = matrix->n();
      if (value_cache.get() == nullptr)
        {
          value_cache =
            std::make_shared<std::vector<TrilinosScalar>>(matrix->n());
          colnum_cache = std::make_shared<std::vector<size_type>>(matrix->n());
        }
      else
        {
          value_cache->resize(matrix->n());
          colnum_cache->resize(matrix->n());
        }

      int ierr = matrix->trilinos_matrix().ExtractGlobalRowCopy(
        this->a_row,
        colnums,
        ncols,
        value_cache->data(),
        reinterpret_cast<TrilinosWrappers::types::int_type *>(
          colnum_cache->data()));
      value_cache->resize(ncols);
      colnum_cache->resize(ncols);
      AssertThrow(ierr == 0, ExcTrilinosError(ierr));

      // copy it into our caches if the
      // line isn't empty. if it is, then
      // we've done something wrong, since
      // we shouldn't have initialized an
      // iterator for an empty line (what
      // would it point to?)
    }
  } // namespace SparseMatrixIterators


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
  SparseMatrix::SparseMatrix()
    : column_space_map(new Epetra_Map(0, 0, Utilities::Trilinos::comm_self()))
    , matrix(
        new Epetra_FECrsMatrix(View, *column_space_map, *column_space_map, 0))
    , last_action(Zero)
    , compressed(true)
  {
    matrix->FillComplete();
  }



  SparseMatrix::SparseMatrix(const size_type    m,
                             const size_type    n,
                             const unsigned int n_max_entries_per_row)
    : column_space_map(
        new Epetra_Map(static_cast<TrilinosWrappers::types::int_type>(n),
                       0,
                       Utilities::Trilinos::comm_self()))
    ,

    // on one processor only, we know how the
    // columns of the matrix will be
    // distributed (everything on one
    // processor), so we can hand in this
    // information to the constructor. we
    // can't do so in parallel, where the
    // information from columns is only
    // available when entries have been added
    matrix(new Epetra_FECrsMatrix(
      Copy,
      Epetra_Map(static_cast<TrilinosWrappers::types::int_type>(m),
                 0,
                 Utilities::Trilinos::comm_self()),
      *column_space_map,
      n_max_entries_per_row,
      false))
    , last_action(Zero)
    , compressed(false)
  {}



  SparseMatrix::SparseMatrix(const size_type                  m,
                             const size_type                  n,
                             const std::vector<unsigned int> &n_entries_per_row)
    : column_space_map(
        new Epetra_Map(static_cast<TrilinosWrappers::types::int_type>(n),
                       0,
                       Utilities::Trilinos::comm_self()))
    , matrix(new Epetra_FECrsMatrix(
        Copy,
        Epetra_Map(static_cast<TrilinosWrappers::types::int_type>(m),
                   0,
                   Utilities::Trilinos::comm_self()),
        *column_space_map,
        reinterpret_cast<int *>(
          const_cast<unsigned int *>(n_entries_per_row.data())),
        false))
    , last_action(Zero)
    , compressed(false)
  {}



  SparseMatrix::SparseMatrix(const IndexSet &   parallel_partitioning,
                             const MPI_Comm &   communicator,
                             const unsigned int n_max_entries_per_row)
    : column_space_map(new Epetra_Map(
        parallel_partitioning.make_trilinos_map(communicator, false)))
    , matrix(new Epetra_FECrsMatrix(Copy,
                                    *column_space_map,
                                    n_max_entries_per_row,
                                    false))
    , last_action(Zero)
    , compressed(false)
  {}



  SparseMatrix::SparseMatrix(const IndexSet &parallel_partitioning,
                             const MPI_Comm &communicator,
                             const std::vector<unsigned int> &n_entries_per_row)
    : column_space_map(new Epetra_Map(
        parallel_partitioning.make_trilinos_map(communicator, false)))
    , matrix(new Epetra_FECrsMatrix(Copy,
                                    *column_space_map,
                                    reinterpret_cast<int *>(
                                      const_cast<unsigned int *>(
                                        n_entries_per_row.data())),
                                    false))
    , last_action(Zero)
    , compressed(false)
  {}



  SparseMatrix::SparseMatrix(const IndexSet &row_parallel_partitioning,
                             const IndexSet &col_parallel_partitioning,
                             const MPI_Comm &communicator,
                             const size_type n_max_entries_per_row)
    : column_space_map(new Epetra_Map(
        col_parallel_partitioning.make_trilinos_map(communicator, false)))
    , matrix(new Epetra_FECrsMatrix(
        Copy,
        row_parallel_partitioning.make_trilinos_map(communicator, false),
        n_max_entries_per_row,
        false))
    , last_action(Zero)
    , compressed(false)
  {}



  SparseMatrix::SparseMatrix(const IndexSet &row_parallel_partitioning,
                             const IndexSet &col_parallel_partitioning,
                             const MPI_Comm &communicator,
                             const std::vector<unsigned int> &n_entries_per_row)
    : column_space_map(new Epetra_Map(
        col_parallel_partitioning.make_trilinos_map(communicator, false)))
    , matrix(new Epetra_FECrsMatrix(
        Copy,
        row_parallel_partitioning.make_trilinos_map(communicator, false),
        reinterpret_cast<int *>(
          const_cast<unsigned int *>(n_entries_per_row.data())),
        false))
    , last_action(Zero)
    , compressed(false)
  {}



  SparseMatrix::SparseMatrix(const SparsityPattern &sparsity_pattern)
    : column_space_map(new Epetra_Map(sparsity_pattern.domain_partitioner()))
    , matrix(
        new Epetra_FECrsMatrix(Copy,
                               sparsity_pattern.trilinos_sparsity_pattern(),
                               false))
    , last_action(Zero)
    , compressed(true)
  {
    Assert(sparsity_pattern.trilinos_sparsity_pattern().Filled() == true,
           ExcMessage(
             "The Trilinos sparsity pattern has not been compressed."));
    compress(VectorOperation::insert);
  }



  SparseMatrix::SparseMatrix(SparseMatrix &&other) noexcept
    : column_space_map(std::move(other.column_space_map))
    , matrix(std::move(other.matrix))
    , nonlocal_matrix(std::move(other.nonlocal_matrix))
    , nonlocal_matrix_exporter(std::move(other.nonlocal_matrix_exporter))
    , last_action(other.last_action)
    , compressed(other.compressed)
  {
    other.last_action = Zero;
    other.compressed  = false;
  }



  void
  SparseMatrix::copy_from(const SparseMatrix &rhs)
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
        // Try to copy all the rows of the matrix one by one. In case of error
        // (i.e., the column indices are different), we need to abort and blow
        // away the matrix.
        for (const auto row : locally_owned_range_indices())
          {
            const int row_local = matrix->RowMap().LID(
              static_cast<TrilinosWrappers::types::int_type>(row));
            Assert((row_local >= 0), ExcAccessToNonlocalRow(row));

            int             n_entries, rhs_n_entries;
            TrilinosScalar *value_ptr, *rhs_value_ptr;
            int *           index_ptr, *rhs_index_ptr;
            int             ierr = rhs.matrix->ExtractMyRowView(row_local,
                                                    rhs_n_entries,
                                                    rhs_value_ptr,
                                                    rhs_index_ptr);
            (void)ierr;
            Assert(ierr == 0, ExcTrilinosError(ierr));

            ierr = matrix->ExtractMyRowView(row_local,
                                            n_entries,
                                            value_ptr,
                                            index_ptr);
            Assert(ierr == 0, ExcTrilinosError(ierr));

            if (n_entries != rhs_n_entries ||
                std::memcmp(static_cast<void *>(index_ptr),
                            static_cast<void *>(rhs_index_ptr),
                            sizeof(int) * n_entries) != 0)
              {
                needs_deep_copy = true;
                break;
              }

            for (int i = 0; i < n_entries; ++i)
              value_ptr[i] = rhs_value_ptr[i];
          }
      }

    if (needs_deep_copy)
      {
        column_space_map =
          std::make_unique<Epetra_Map>(rhs.trilinos_matrix().DomainMap());

        // release memory before reallocation
        matrix = std::make_unique<Epetra_FECrsMatrix>(*rhs.matrix);

        matrix->FillComplete(*column_space_map, matrix->RowMap());
      }

    if (rhs.nonlocal_matrix.get() != nullptr)
      nonlocal_matrix =
        std::make_unique<Epetra_CrsMatrix>(Copy, rhs.nonlocal_matrix->Graph());
  }



  namespace
  {
    using size_type = SparseMatrix::size_type;

    template <typename SparsityPatternType>
    void
    reinit_matrix(const IndexSet &             row_parallel_partitioning,
                  const IndexSet &             column_parallel_partitioning,
                  const SparsityPatternType &  sparsity_pattern,
                  const bool                   exchange_data,
                  const MPI_Comm &             communicator,
                  std::unique_ptr<Epetra_Map> &column_space_map,
                  std::unique_ptr<Epetra_FECrsMatrix> &matrix,
                  std::unique_ptr<Epetra_CrsMatrix> &  nonlocal_matrix,
                  std::unique_ptr<Epetra_Export> &     nonlocal_matrix_exporter)
    {
      // release memory before reallocation
      matrix.reset();
      nonlocal_matrix.reset();
      nonlocal_matrix_exporter.reset();

      column_space_map = std::make_unique<Epetra_Map>(
        column_parallel_partitioning.make_trilinos_map(communicator, false));

      if (column_space_map->Comm().MyPID() == 0)
        {
          AssertDimension(sparsity_pattern.n_rows(),
                          row_parallel_partitioning.size());
          AssertDimension(sparsity_pattern.n_cols(),
                          column_parallel_partitioning.size());
        }

      Epetra_Map row_space_map =
        row_parallel_partitioning.make_trilinos_map(communicator, false);

      // if we want to exchange data, build a usual Trilinos sparsity pattern
      // and let that handle the exchange. otherwise, manually create a
      // CrsGraph, which consumes considerably less memory because it can set
      // correct number of indices right from the start
      if (exchange_data)
        {
          SparsityPattern trilinos_sparsity;
          trilinos_sparsity.reinit(row_parallel_partitioning,
                                   column_parallel_partitioning,
                                   sparsity_pattern,
                                   communicator,
                                   exchange_data);
          matrix = std::make_unique<Epetra_FECrsMatrix>(
            Copy, trilinos_sparsity.trilinos_sparsity_pattern(), false);

          return;
        }

      const size_type first_row = TrilinosWrappers::min_my_gid(row_space_map),
                      last_row =
                        TrilinosWrappers::max_my_gid(row_space_map) + 1;
      std::vector<int> n_entries_per_row(last_row - first_row);

      for (size_type row = first_row; row < last_row; ++row)
        n_entries_per_row[row - first_row] = sparsity_pattern.row_length(row);

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
      std::unique_ptr<Epetra_CrsGraph> graph;
      if (row_space_map.Comm().NumProc() > 1)
        graph = std::make_unique<Epetra_CrsGraph>(Copy,
                                                  row_space_map,
                                                  n_entries_per_row.data(),
                                                  true);
      else
        graph = std::make_unique<Epetra_CrsGraph>(Copy,
                                                  row_space_map,
                                                  *column_space_map,
                                                  n_entries_per_row.data(),
                                                  true);

      // This functions assumes that the sparsity pattern sits on all
      // processors (completely). The parallel version uses an Epetra graph
      // that is already distributed.

      // now insert the indices
      std::vector<TrilinosWrappers::types::int_type> row_indices;

      for (size_type row = first_row; row < last_row; ++row)
        {
          const int row_length = sparsity_pattern.row_length(row);
          if (row_length == 0)
            continue;

          row_indices.resize(row_length, -1);
          {
            typename SparsityPatternType::iterator p =
              sparsity_pattern.begin(row);
            for (size_type col = 0; p != sparsity_pattern.end(row); ++p, ++col)
              row_indices[col] = p->column();
          }
          graph->Epetra_CrsGraph::InsertGlobalIndices(row,
                                                      row_length,
                                                      row_indices.data());
        }

      // Eventually, optimize the graph structure (sort indices, make memory
      // contiguous, etc). note that the documentation of the function indeed
      // states that we first need to provide the column (domain) map and then
      // the row (range) map
      graph->FillComplete(*column_space_map, row_space_map);
      graph->OptimizeStorage();

      // check whether we got the number of columns right.
      AssertDimension(sparsity_pattern.n_cols(),
                      TrilinosWrappers::n_global_cols(*graph));
      (void)n_global_cols;

      // And now finally generate the matrix.
      matrix = std::make_unique<Epetra_FECrsMatrix>(Copy, *graph, false);
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
      Epetra_CrsGraphMod(const Epetra_Map &row_map,
                         const int *       n_entries_per_row)
        : Epetra_CrsGraph(Copy, row_map, n_entries_per_row, true)
      {}

      void
      SetIndicesAreGlobal()
      {
        this->Epetra_CrsGraph::SetIndicesAreGlobal(true);
      }
    };



    // specialization for DynamicSparsityPattern which can provide us with
    // more information about the non-locally owned rows
    template <>
    void
    reinit_matrix(const IndexSet &              row_parallel_partitioning,
                  const IndexSet &              column_parallel_partitioning,
                  const DynamicSparsityPattern &sparsity_pattern,
                  const bool                    exchange_data,
                  const MPI_Comm &              communicator,
                  std::unique_ptr<Epetra_Map> & column_space_map,
                  std::unique_ptr<Epetra_FECrsMatrix> &matrix,
                  std::unique_ptr<Epetra_CrsMatrix> &  nonlocal_matrix,
                  std::unique_ptr<Epetra_Export> &     nonlocal_matrix_exporter)
    {
      matrix.reset();
      nonlocal_matrix.reset();
      nonlocal_matrix_exporter.reset();

      column_space_map = std::make_unique<Epetra_Map>(
        column_parallel_partitioning.make_trilinos_map(communicator, false));

      AssertDimension(sparsity_pattern.n_rows(),
                      row_parallel_partitioning.size());
      AssertDimension(sparsity_pattern.n_cols(),
                      column_parallel_partitioning.size());

      Epetra_Map row_space_map =
        row_parallel_partitioning.make_trilinos_map(communicator, false);

      IndexSet relevant_rows(sparsity_pattern.row_index_set());
      // serial case
      if (relevant_rows.size() == 0)
        {
          relevant_rows.set_size(
            TrilinosWrappers::n_global_elements(row_space_map));
          relevant_rows.add_range(
            0, TrilinosWrappers::n_global_elements(row_space_map));
        }
      relevant_rows.compress();
      Assert(relevant_rows.n_elements() >=
               static_cast<unsigned int>(row_space_map.NumMyElements()),
             ExcMessage(
               "Locally relevant rows of sparsity pattern must contain "
               "all locally owned rows"));

      // check whether the relevant rows correspond to exactly the same map as
      // the owned rows. In that case, do not create the nonlocal graph and
      // fill the columns by demand
      const bool have_ghost_rows = [&]() {
        std::vector<dealii::types::global_dof_index> indices;
        relevant_rows.fill_index_vector(indices);
        Epetra_Map relevant_map(
          TrilinosWrappers::types::int_type(-1),
          TrilinosWrappers::types::int_type(relevant_rows.n_elements()),
          (indices.empty() ?
             nullptr :
             reinterpret_cast<TrilinosWrappers::types::int_type *>(
               indices.data())),
          0,
          row_space_map.Comm());
        return !relevant_map.SameAs(row_space_map);
      }();

      const unsigned int n_rows = relevant_rows.n_elements();
      std::vector<TrilinosWrappers::types::int_type> ghost_rows;
      std::vector<int> n_entries_per_row(row_space_map.NumMyElements());
      std::vector<int> n_entries_per_ghost_row;
      for (unsigned int i = 0, own = 0; i < n_rows; ++i)
        {
          const TrilinosWrappers::types::int_type global_row =
            relevant_rows.nth_index_in_set(i);
          if (row_space_map.MyGID(global_row))
            n_entries_per_row[own++] = sparsity_pattern.row_length(global_row);
          else if (sparsity_pattern.row_length(global_row) > 0)
            {
              ghost_rows.push_back(global_row);
              n_entries_per_ghost_row.push_back(
                sparsity_pattern.row_length(global_row));
            }
        }

      Epetra_Map off_processor_map(-1,
                                   ghost_rows.size(),
                                   (ghost_rows.size() > 0) ?
                                     (ghost_rows.data()) :
                                     nullptr,
                                   0,
                                   row_space_map.Comm());

      std::unique_ptr<Epetra_CrsGraph>    graph;
      std::unique_ptr<Epetra_CrsGraphMod> nonlocal_graph;
      if (row_space_map.Comm().NumProc() > 1)
        {
          graph =
            std::make_unique<Epetra_CrsGraph>(Copy,
                                              row_space_map,
                                              (n_entries_per_row.size() > 0) ?
                                                (n_entries_per_row.data()) :
                                                nullptr,
                                              exchange_data ? false : true);
          if (have_ghost_rows == true)
            nonlocal_graph = std::make_unique<Epetra_CrsGraphMod>(
              off_processor_map, n_entries_per_ghost_row.data());
        }
      else
        graph =
          std::make_unique<Epetra_CrsGraph>(Copy,
                                            row_space_map,
                                            *column_space_map,
                                            (n_entries_per_row.size() > 0) ?
                                              (n_entries_per_row.data()) :
                                              nullptr,
                                            true);

      // now insert the indices, select between the right matrix
      std::vector<TrilinosWrappers::types::int_type> row_indices;

      for (unsigned int i = 0; i < n_rows; ++i)
        {
          const TrilinosWrappers::types::int_type global_row =
            relevant_rows.nth_index_in_set(i);
          const int row_length = sparsity_pattern.row_length(global_row);
          if (row_length == 0)
            continue;

          row_indices.resize(row_length, -1);
          for (int col = 0; col < row_length; ++col)
            row_indices[col] = sparsity_pattern.column_number(global_row, col);

          if (row_space_map.MyGID(global_row))
            graph->InsertGlobalIndices(global_row,
                                       row_length,
                                       row_indices.data());
          else
            {
              Assert(nonlocal_graph.get() != nullptr, ExcInternalError());
              nonlocal_graph->InsertGlobalIndices(global_row,
                                                  row_length,
                                                  row_indices.data());
            }
        }

      // finalize nonlocal graph and create nonlocal matrix
      if (nonlocal_graph.get() != nullptr)
        {
          // must make sure the IndicesAreGlobal flag is set on all processors
          // because some processors might not call InsertGlobalIndices (and
          // we do not want to insert dummy indices on all processors for
          // large-scale simulations due to the bad impact on performance)
          nonlocal_graph->SetIndicesAreGlobal();
          Assert(nonlocal_graph->IndicesAreGlobal() == true,
                 ExcInternalError());
          nonlocal_graph->FillComplete(*column_space_map, row_space_map);
          nonlocal_graph->OptimizeStorage();

          // insert data from nonlocal graph into the final sparsity pattern
          if (exchange_data)
            {
              Epetra_Export exporter(nonlocal_graph->RowMap(), row_space_map);
              int ierr = graph->Export(*nonlocal_graph, exporter, Add);
              (void)ierr;
              Assert(ierr == 0, ExcTrilinosError(ierr));
            }

          nonlocal_matrix =
            std::make_unique<Epetra_CrsMatrix>(Copy, *nonlocal_graph);
        }

      graph->FillComplete(*column_space_map, row_space_map);
      graph->OptimizeStorage();

      AssertDimension(sparsity_pattern.n_cols(),
                      TrilinosWrappers::n_global_cols(*graph));

      matrix = std::make_unique<Epetra_FECrsMatrix>(Copy, *graph, false);
    }
  } // namespace



  template <typename SparsityPatternType>
  void
  SparseMatrix::reinit(const SparsityPatternType &sparsity_pattern)
  {
    reinit_matrix(complete_index_set(sparsity_pattern.n_rows()),
                  complete_index_set(sparsity_pattern.n_cols()),
                  sparsity_pattern,
                  false,
                  MPI_COMM_SELF,
                  column_space_map,
                  matrix,
                  nonlocal_matrix,
                  nonlocal_matrix_exporter);
  }



  template <typename SparsityPatternType>
  inline typename std::enable_if<
    !std::is_same<SparsityPatternType,
                  dealii::SparseMatrix<double>>::value>::type
  SparseMatrix::reinit(const IndexSet &           row_parallel_partitioning,
                       const IndexSet &           col_parallel_partitioning,
                       const SparsityPatternType &sparsity_pattern,
                       const MPI_Comm &           communicator,
                       const bool                 exchange_data)
  {
    reinit_matrix(row_parallel_partitioning,
                  col_parallel_partitioning,
                  sparsity_pattern,
                  exchange_data,
                  communicator,
                  column_space_map,
                  matrix,
                  nonlocal_matrix,
                  nonlocal_matrix_exporter);

    // In the end, the matrix needs to be compressed in order to be really
    // ready.
    last_action = Zero;
    compress(VectorOperation::insert);
  }



  void
  SparseMatrix::reinit(const SparsityPattern &sparsity_pattern)
  {
    matrix.reset();
    nonlocal_matrix_exporter.reset();

    // reinit with a (parallel) Trilinos sparsity pattern.
    column_space_map =
      std::make_unique<Epetra_Map>(sparsity_pattern.domain_partitioner());
    matrix = std::make_unique<Epetra_FECrsMatrix>(
      Copy, sparsity_pattern.trilinos_sparsity_pattern(), false);

    if (sparsity_pattern.nonlocal_graph.get() != nullptr)
      nonlocal_matrix =
        std::make_unique<Epetra_CrsMatrix>(Copy,
                                           *sparsity_pattern.nonlocal_graph);
    else
      nonlocal_matrix.reset();

    last_action = Zero;
    compress(VectorOperation::insert);
  }



  void
  SparseMatrix::reinit(const SparseMatrix &sparse_matrix)
  {
    if (this == &sparse_matrix)
      return;

    column_space_map =
      std::make_unique<Epetra_Map>(sparse_matrix.trilinos_matrix().DomainMap());
    matrix.reset();
    nonlocal_matrix_exporter.reset();
    matrix = std::make_unique<Epetra_FECrsMatrix>(
      Copy, sparse_matrix.trilinos_sparsity_pattern(), false);

    if (sparse_matrix.nonlocal_matrix != nullptr)
      nonlocal_matrix = std::make_unique<Epetra_CrsMatrix>(
        Copy, sparse_matrix.nonlocal_matrix->Graph());
    else
      nonlocal_matrix.reset();

    last_action = Zero;
    compress(VectorOperation::insert);
  }



  template <typename number>
  inline void
  SparseMatrix::reinit(
    const IndexSet &                      row_parallel_partitioning,
    const IndexSet &                      col_parallel_partitioning,
    const ::dealii::SparseMatrix<number> &dealii_sparse_matrix,
    const MPI_Comm &                      communicator,
    const double                          drop_tolerance,
    const bool                            copy_values,
    const ::dealii::SparsityPattern *     use_this_sparsity)
  {
    if (copy_values == false)
      {
        // in case we do not copy values, just
        // call the other function.
        if (use_this_sparsity == nullptr)
          reinit(row_parallel_partitioning,
                 col_parallel_partitioning,
                 dealii_sparse_matrix.get_sparsity_pattern(),
                 communicator,
                 false);
        else
          reinit(row_parallel_partitioning,
                 col_parallel_partitioning,
                 *use_this_sparsity,
                 communicator,
                 false);
        return;
      }

    const size_type n_rows = dealii_sparse_matrix.m();

    AssertDimension(row_parallel_partitioning.size(), n_rows);
    AssertDimension(col_parallel_partitioning.size(), dealii_sparse_matrix.n());

    const ::dealii::SparsityPattern &sparsity_pattern =
      (use_this_sparsity != nullptr) ?
        *use_this_sparsity :
        dealii_sparse_matrix.get_sparsity_pattern();

    if (matrix.get() == nullptr || m() != n_rows ||
        n_nonzero_elements() != sparsity_pattern.n_nonzero_elements())
      {
        reinit(row_parallel_partitioning,
               col_parallel_partitioning,
               sparsity_pattern,
               communicator,
               false);
      }

    // fill the values. the same as above: go through all rows of the
    // matrix, and then all columns. since the sparsity patterns of the
    // input matrix and the specified sparsity pattern might be different,
    // need to go through the row for both these sparsity structures
    // simultaneously in order to really set the correct values.
    size_type                   maximum_row_length = matrix->MaxNumEntries();
    std::vector<size_type>      row_indices(maximum_row_length);
    std::vector<TrilinosScalar> values(maximum_row_length);

    for (size_type row = 0; row < n_rows; ++row)
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



  template <typename number>
  void
  SparseMatrix::reinit(
    const ::dealii::SparseMatrix<number> &dealii_sparse_matrix,
    const double                          drop_tolerance,
    const bool                            copy_values,
    const ::dealii::SparsityPattern *     use_this_sparsity)
  {
    reinit(complete_index_set(dealii_sparse_matrix.m()),
           complete_index_set(dealii_sparse_matrix.n()),
           dealii_sparse_matrix,
           MPI_COMM_SELF,
           drop_tolerance,
           copy_values,
           use_this_sparsity);
  }



  void
  SparseMatrix::reinit(const Epetra_CrsMatrix &input_matrix,
                       const bool              copy_values)
  {
    Assert(input_matrix.Filled() == true,
           ExcMessage("Input CrsMatrix has not called FillComplete()!"));

    column_space_map = std::make_unique<Epetra_Map>(input_matrix.DomainMap());

    const Epetra_CrsGraph *graph = &input_matrix.Graph();

    nonlocal_matrix.reset();
    nonlocal_matrix_exporter.reset();
    matrix.reset();
    matrix = std::make_unique<Epetra_FECrsMatrix>(Copy, *graph, false);

    matrix->FillComplete(*column_space_map, input_matrix.RangeMap(), true);

    if (copy_values == true)
      {
        // point to the first data entry in the two
        // matrices and copy the content
        const TrilinosScalar *in_values   = input_matrix[0];
        TrilinosScalar *      values      = (*matrix)[0];
        const size_type       my_nonzeros = input_matrix.NumMyNonzeros();
        std::memcpy(values, in_values, my_nonzeros * sizeof(TrilinosScalar));
      }

    last_action = Zero;
    compress(VectorOperation::insert);
  }



  void
  SparseMatrix::compress(::dealii::VectorOperation::values operation)
  {
    Epetra_CombineMode mode = last_action;
    if (last_action == Zero)
      {
        if ((operation == ::dealii::VectorOperation::add) ||
            (operation == ::dealii::VectorOperation::unknown))
          mode = Add;
        else if (operation == ::dealii::VectorOperation::insert)
          mode = Insert;
        else
          Assert(
            false,
            ExcMessage(
              "compress() can only be called with VectorOperation add, insert, or unknown"));
      }
    else
      {
        Assert(((last_action == Add) &&
                (operation != ::dealii::VectorOperation::insert)) ||
                 ((last_action == Insert) &&
                  (operation != ::dealii::VectorOperation::add)),
               ExcMessage("Operation and argument to compress() do not match"));
      }

    // flush buffers
    int ierr;
    if (nonlocal_matrix.get() != nullptr && mode == Add)
      {
        // do only export in case of an add() operation, otherwise the owning
        // processor must have set the correct entry
        nonlocal_matrix->FillComplete(*column_space_map, matrix->RowMap());
        if (nonlocal_matrix_exporter.get() == nullptr)
          nonlocal_matrix_exporter =
            std::make_unique<Epetra_Export>(nonlocal_matrix->RowMap(),
                                            matrix->RowMap());
        ierr =
          matrix->Export(*nonlocal_matrix, *nonlocal_matrix_exporter, mode);
        AssertThrow(ierr == 0, ExcTrilinosError(ierr));
        ierr = matrix->FillComplete(*column_space_map, matrix->RowMap());
        nonlocal_matrix->PutScalar(0);
      }
    else
      ierr =
        matrix->GlobalAssemble(*column_space_map, matrix->RowMap(), true, mode);

    AssertThrow(ierr == 0, ExcTrilinosError(ierr));

    ierr = matrix->OptimizeStorage();
    AssertThrow(ierr == 0, ExcTrilinosError(ierr));

    last_action = Zero;

    compressed = true;
  }



  void
  SparseMatrix::clear()
  {
    // When we clear the matrix, reset
    // the pointer and generate an
    // empty matrix.
    column_space_map =
      std::make_unique<Epetra_Map>(0, 0, Utilities::Trilinos::comm_self());
    matrix = std::make_unique<Epetra_FECrsMatrix>(View, *column_space_map, 0);
    nonlocal_matrix.reset();
    nonlocal_matrix_exporter.reset();

    matrix->FillComplete();

    compressed = true;
  }



  void
  SparseMatrix::clear_row(const size_type      row,
                          const TrilinosScalar new_diag_value)
  {
    Assert(matrix->Filled() == true, ExcMatrixNotCompressed());

    // Only do this on the rows owned
    // locally on this processor.
    int local_row =
      matrix->LRID(static_cast<TrilinosWrappers::types::int_type>(row));
    if (local_row >= 0)
      {
        TrilinosScalar *values;
        int *           col_indices;
        int             num_entries;
        const int       ierr =
          matrix->ExtractMyRowView(local_row, num_entries, values, col_indices);
        (void)ierr;

        Assert(ierr == 0, ExcTrilinosError(ierr));

        const std::ptrdiff_t diag_index =
          std::find(col_indices, col_indices + num_entries, local_row) -
          col_indices;

        for (TrilinosWrappers::types::int_type j = 0; j < num_entries; ++j)
          if (diag_index != j || new_diag_value == 0)
            values[j] = 0.;

        if (diag_index != num_entries && std::fabs(values[diag_index]) == 0.0 &&
            new_diag_value != 0.0)
          values[diag_index] = new_diag_value;
      }
  }



  void
  SparseMatrix::clear_rows(const std::vector<size_type> &rows,
                           const TrilinosScalar          new_diag_value)
  {
    for (const auto row : rows)
      clear_row(row, new_diag_value);
  }



  TrilinosScalar
  SparseMatrix::operator()(const size_type i, const size_type j) const
  {
    // Extract local indices in
    // the matrix.
    int trilinos_i =
          matrix->LRID(static_cast<TrilinosWrappers::types::int_type>(i)),
        trilinos_j =
          matrix->LCID(static_cast<TrilinosWrappers::types::int_type>(j));
    TrilinosScalar value = 0.;

    // If the data is not on the
    // present processor, we throw
    // an exception. This is one of
    // the two tiny differences to
    // the el(i,j) call, which does
    // not throw any assertions.
    if (trilinos_i == -1)
      {
        Assert(false,
               ExcAccessToNonLocalElement(
                 i, j, local_range().first, local_range().second));
      }
    else
      {
        // Check whether the matrix has
        // already been transformed to local
        // indices.
        Assert(matrix->Filled(), ExcMatrixNotCompressed());

        // Prepare pointers for extraction
        // of a view of the row.
        int             nnz_present = matrix->NumMyEntries(trilinos_i);
        int             nnz_extracted;
        int *           col_indices;
        TrilinosScalar *values;

        // Generate the view and make
        // sure that we have not generated
        // an error.
        // TODO Check that col_indices are int and not long long
        int ierr = matrix->ExtractMyRowView(trilinos_i,
                                            nnz_extracted,
                                            values,
                                            col_indices);
        (void)ierr;
        Assert(ierr == 0, ExcTrilinosError(ierr));

        Assert(nnz_present == nnz_extracted,
               ExcDimensionMismatch(nnz_present, nnz_extracted));

        // Search the index where we
        // look for the value, and then
        // finally get it.
        const std::ptrdiff_t local_col_index =
          std::find(col_indices, col_indices + nnz_present, trilinos_j) -
          col_indices;

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
            Assert(false, ExcInvalidIndex(i, j));
          }
        else
          value = values[local_col_index];
      }

    return value;
  }



  TrilinosScalar
  SparseMatrix::el(const size_type i, const size_type j) const
  {
    // Extract local indices in
    // the matrix.
    int trilinos_i =
          matrix->LRID(static_cast<TrilinosWrappers::types::int_type>(i)),
        trilinos_j =
          matrix->LCID(static_cast<TrilinosWrappers::types::int_type>(j));
    TrilinosScalar value = 0.;

    // If the data is not on the
    // present processor, we can't
    // continue. Just print out zero
    // as discussed in the
    // documentation of this
    // function. if you want error
    // checking, use operator().
    if ((trilinos_i == -1) || (trilinos_j == -1))
      return 0.;
    else
      {
        // Check whether the matrix
        // already is transformed to
        // local indices.
        Assert(matrix->Filled(), ExcMatrixNotCompressed());

        // Prepare pointers for extraction
        // of a view of the row.
        int             nnz_present = matrix->NumMyEntries(trilinos_i);
        int             nnz_extracted;
        int *           col_indices;
        TrilinosScalar *values;

        // Generate the view and make
        // sure that we have not generated
        // an error.
        int ierr = matrix->ExtractMyRowView(trilinos_i,
                                            nnz_extracted,
                                            values,
                                            col_indices);
        (void)ierr;
        Assert(ierr == 0, ExcTrilinosError(ierr));

        Assert(nnz_present == nnz_extracted,
               ExcDimensionMismatch(nnz_present, nnz_extracted));

        // Search the index where we
        // look for the value, and then
        // finally get it.
        const std::ptrdiff_t local_col_index =
          std::find(col_indices, col_indices + nnz_present, trilinos_j) -
          col_indices;

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
  SparseMatrix::diag_element(const size_type i) const
  {
    Assert(m() == n(), ExcNotQuadratic());

#  ifdef DEBUG
    // use operator() in debug mode because
    // it checks if this is a valid element
    // (in parallel)
    return operator()(i, i);
#  else
    // Trilinos doesn't seem to have a
    // more efficient way to access the
    // diagonal than by just using the
    // standard el(i,j) function.
    return el(i, i);
#  endif
  }



  unsigned int
  SparseMatrix::row_length(const size_type row) const
  {
    Assert(row < m(), ExcInternalError());

    // get a representation of the
    // present row
    int ncols = -1;
    int local_row =
      matrix->LRID(static_cast<TrilinosWrappers::types::int_type>(row));
    Assert((local_row >= 0), ExcAccessToNonlocalRow(row));

    // on the processor who owns this
    // row, we'll have a non-negative
    // value.
    if (local_row >= 0)
      {
        int ierr = matrix->NumMyRowEntries(local_row, ncols);
        AssertThrow(ierr == 0, ExcTrilinosError(ierr));
      }

    return static_cast<unsigned int>(ncols);
  }



  void
  SparseMatrix::set(const std::vector<size_type> &    row_indices,
                    const std::vector<size_type> &    col_indices,
                    const FullMatrix<TrilinosScalar> &values,
                    const bool                        elide_zero_values)
  {
    Assert(row_indices.size() == values.m(),
           ExcDimensionMismatch(row_indices.size(), values.m()));
    Assert(col_indices.size() == values.n(),
           ExcDimensionMismatch(col_indices.size(), values.n()));

    for (size_type i = 0; i < row_indices.size(); ++i)
      set(row_indices[i],
          col_indices.size(),
          col_indices.data(),
          &values(i, 0),
          elide_zero_values);
  }



  void
  SparseMatrix::set(const size_type                    row,
                    const std::vector<size_type> &     col_indices,
                    const std::vector<TrilinosScalar> &values,
                    const bool                         elide_zero_values)
  {
    Assert(col_indices.size() == values.size(),
           ExcDimensionMismatch(col_indices.size(), values.size()));

    set(row,
        col_indices.size(),
        col_indices.data(),
        values.data(),
        elide_zero_values);
  }



  template <>
  void
  SparseMatrix::set<TrilinosScalar>(const size_type       row,
                                    const size_type       n_cols,
                                    const size_type *     col_indices,
                                    const TrilinosScalar *values,
                                    const bool            elide_zero_values)
  {
    AssertIndexRange(row, this->m());

    int ierr;
    if (last_action == Add)
      {
        ierr =
          matrix->GlobalAssemble(*column_space_map, matrix->RowMap(), true);

        Assert(ierr == 0, ExcTrilinosError(ierr));
      }

    last_action = Insert;

    const TrilinosWrappers::types::int_type *col_index_ptr;
    const TrilinosScalar *                   col_value_ptr;
    const TrilinosWrappers::types::int_type  trilinos_row = row;
    TrilinosWrappers::types::int_type        n_columns;

    boost::container::small_vector<TrilinosScalar, 200> local_value_array(
      n_cols);
    boost::container::small_vector<TrilinosWrappers::types::int_type, 200>
      local_index_array(n_cols);

    // If we don't elide zeros, the pointers are already available... need to
    // cast to non-const pointers as that is the format taken by Trilinos (but
    // we will not modify const data)
    if (elide_zero_values == false)
      {
        col_index_ptr =
          reinterpret_cast<const TrilinosWrappers::types::int_type *>(
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
                n_columns++;
              }
          }

        AssertIndexRange(n_columns, n_cols + 1);
      }


    // If the calling matrix owns the row to which we want to insert values,
    // we can directly call the Epetra_CrsMatrix input function, which is much
    // faster than the Epetra_FECrsMatrix function. We distinguish between two
    // cases: the first one is when the matrix is not filled (i.e., it is
    // possible to add new elements to the sparsity pattern), and the second
    // one is when the pattern is already fixed. In the former case, we add
    // the possibility to insert new values, and in the second we just replace
    // data.
    if (matrix->RowMap().MyGID(
          static_cast<TrilinosWrappers::types::int_type>(row)) == true)
      {
        if (matrix->Filled() == false)
          {
            ierr = matrix->Epetra_CrsMatrix::InsertGlobalValues(
              row, static_cast<int>(n_columns), col_value_ptr, col_index_ptr);

            // When inserting elements, we do not want to create exceptions in
            // the case when inserting non-local data (since that's what we
            // want to do right now).
            if (ierr > 0)
              ierr = 0;
          }
        else
          ierr = matrix->Epetra_CrsMatrix::ReplaceGlobalValues(row,
                                                               n_columns,
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
            ierr = matrix->InsertGlobalValues(1,
                                              &trilinos_row,
                                              n_columns,
                                              col_index_ptr,
                                              &col_value_ptr,
                                              Epetra_FECrsMatrix::ROW_MAJOR);
            if (ierr > 0)
              ierr = 0;
          }
        else
          ierr = matrix->ReplaceGlobalValues(1,
                                             &trilinos_row,
                                             n_columns,
                                             col_index_ptr,
                                             &col_value_ptr,
                                             Epetra_FECrsMatrix::ROW_MAJOR);
        // use the FECrsMatrix facilities for set even in the case when we
        // have explicitly set the off-processor rows because that only works
        // properly when adding elements, not when setting them (since we want
        // to only touch elements that have been set explicitly, and there is
        // no way on the receiving processor to identify them otherwise)
      }

    Assert(ierr <= 0, ExcAccessToNonPresentElement(row, col_index_ptr[0]));
    AssertThrow(ierr >= 0, ExcTrilinosError(ierr));
  }



  void
  SparseMatrix::add(const std::vector<size_type> &    indices,
                    const FullMatrix<TrilinosScalar> &values,
                    const bool                        elide_zero_values)
  {
    Assert(indices.size() == values.m(),
           ExcDimensionMismatch(indices.size(), values.m()));
    Assert(values.m() == values.n(), ExcNotQuadratic());

    for (size_type i = 0; i < indices.size(); ++i)
      add(indices[i],
          indices.size(),
          indices.data(),
          &values(i, 0),
          elide_zero_values);
  }



  void
  SparseMatrix::add(const std::vector<size_type> &    row_indices,
                    const std::vector<size_type> &    col_indices,
                    const FullMatrix<TrilinosScalar> &values,
                    const bool                        elide_zero_values)
  {
    Assert(row_indices.size() == values.m(),
           ExcDimensionMismatch(row_indices.size(), values.m()));
    Assert(col_indices.size() == values.n(),
           ExcDimensionMismatch(col_indices.size(), values.n()));

    for (size_type i = 0; i < row_indices.size(); ++i)
      add(row_indices[i],
          col_indices.size(),
          col_indices.data(),
          &values(i, 0),
          elide_zero_values);
  }



  void
  SparseMatrix::add(const size_type                    row,
                    const std::vector<size_type> &     col_indices,
                    const std::vector<TrilinosScalar> &values,
                    const bool                         elide_zero_values)
  {
    Assert(col_indices.size() == values.size(),
           ExcDimensionMismatch(col_indices.size(), values.size()));

    add(row,
        col_indices.size(),
        col_indices.data(),
        values.data(),
        elide_zero_values);
  }



  void
  SparseMatrix::add(const size_type       row,
                    const size_type       n_cols,
                    const size_type *     col_indices,
                    const TrilinosScalar *values,
                    const bool            elide_zero_values,
                    const bool /*col_indices_are_sorted*/)
  {
    AssertIndexRange(row, this->m());
    int ierr;
    if (last_action == Insert)
      {
        // TODO: this could lead to a dead lock when only one processor
        // calls GlobalAssemble.
        ierr =
          matrix->GlobalAssemble(*column_space_map, matrix->RowMap(), false);

        AssertThrow(ierr == 0, ExcTrilinosError(ierr));
      }

    last_action = Add;

    const TrilinosWrappers::types::int_type *col_index_ptr;
    const TrilinosScalar *                   col_value_ptr;
    const TrilinosWrappers::types::int_type  trilinos_row = row;
    TrilinosWrappers::types::int_type        n_columns;

    boost::container::small_vector<TrilinosScalar, 100> local_value_array(
      n_cols);
    boost::container::small_vector<TrilinosWrappers::types::int_type, 100>
      local_index_array(n_cols);

    // If we don't elide zeros, the pointers are already available... need to
    // cast to non-const pointers as that is the format taken by Trilinos (but
    // we will not modify const data)
    if (elide_zero_values == false)
      {
        col_index_ptr =
          reinterpret_cast<const TrilinosWrappers::types::int_type *>(
            col_indices);
        col_value_ptr = values;
        n_columns     = n_cols;
#  ifdef DEBUG
        for (size_type j = 0; j < n_cols; ++j)
          AssertIsFinite(values[j]);
#  endif
      }
    else
      {
        // Otherwise, extract nonzero values in each row and the corresponding
        // index.
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
    // Exit early if there is nothing to do
    if (n_columns == 0)
      {
        return;
      }

    // If the calling processor owns the row to which we want to add values, we
    // can directly call the Epetra_CrsMatrix input function, which is much
    // faster than the Epetra_FECrsMatrix function.
    if (matrix->RowMap().MyGID(
          static_cast<TrilinosWrappers::types::int_type>(row)) == true)
      {
        ierr = matrix->Epetra_CrsMatrix::SumIntoGlobalValues(row,
                                                             n_columns,
                                                             col_value_ptr,
                                                             col_index_ptr);
        AssertThrow(ierr == 0, ExcTrilinosError(ierr));
      }
    else if (nonlocal_matrix.get() != nullptr)
      {
        compressed = false;
        // this is the case when we have explicitly set the off-processor rows
        // and want to create a separate matrix object for them (to retain
        // thread-safety)
        Assert(nonlocal_matrix->RowMap().LID(
                 static_cast<TrilinosWrappers::types::int_type>(row)) != -1,
               ExcMessage("Attempted to write into off-processor matrix row "
                          "that has not be specified as being writable upon "
                          "initialization"));
        ierr = nonlocal_matrix->SumIntoGlobalValues(row,
                                                    n_columns,
                                                    col_value_ptr,
                                                    col_index_ptr);
        AssertThrow(ierr == 0, ExcTrilinosError(ierr));
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

        ierr = matrix->SumIntoGlobalValues(1,
                                           &trilinos_row,
                                           n_columns,
                                           col_index_ptr,
                                           &col_value_ptr,
                                           Epetra_FECrsMatrix::ROW_MAJOR);
        AssertThrow(ierr == 0, ExcTrilinosError(ierr));
      }

#  ifdef DEBUG
    if (ierr > 0)
      {
        std::cout << "------------------------------------------" << std::endl;
        std::cout << "Got error " << ierr << " in row " << row << " of proc "
                  << matrix->RowMap().Comm().MyPID()
                  << " when trying to add the columns:" << std::endl;
        for (TrilinosWrappers::types::int_type i = 0; i < n_columns; ++i)
          std::cout << col_index_ptr[i] << " ";
        std::cout << std::endl << std::endl;
        std::cout << "Matrix row "
                  << (matrix->RowMap().MyGID(
                        static_cast<TrilinosWrappers::types::int_type>(row)) ==
                          false ?
                        "(nonlocal part)" :
                        "")
                  << " has the following indices:" << std::endl;
        std::vector<TrilinosWrappers::types::int_type> indices;
        const Epetra_CrsGraph *                        graph =
          (nonlocal_matrix.get() != nullptr &&
           matrix->RowMap().MyGID(
             static_cast<TrilinosWrappers::types::int_type>(row)) == false) ?
            &nonlocal_matrix->Graph() :
            &matrix->Graph();

        indices.resize(graph->NumGlobalIndices(row));
        int n_indices = 0;
        graph->ExtractGlobalRowCopy(row,
                                    indices.size(),
                                    n_indices,
                                    indices.data());
        AssertDimension(n_indices, indices.size());

        for (TrilinosWrappers::types::int_type i = 0; i < n_indices; ++i)
          std::cout << indices[i] << " ";
        std::cout << std::endl << std::endl;
        Assert(ierr <= 0, ExcAccessToNonPresentElement(row, col_index_ptr[0]));
      }
#  endif
    Assert(ierr >= 0, ExcTrilinosError(ierr));
  }



  SparseMatrix &
  SparseMatrix::operator=(const double d)
  {
    Assert(d == 0, ExcScalarAssignmentOnlyForZeroValue());
    compress(
      ::dealii::VectorOperation::unknown); // TODO: why do we do this? Should we
                                           // not check for is_compressed?

    const int ierr = matrix->PutScalar(d);
    AssertThrow(ierr == 0, ExcTrilinosError(ierr));
    if (nonlocal_matrix.get() != nullptr)
      nonlocal_matrix->PutScalar(d);

    return *this;
  }



  void
  SparseMatrix::add(const TrilinosScalar factor, const SparseMatrix &rhs)
  {
    AssertDimension(rhs.m(), m());
    AssertDimension(rhs.n(), n());
    AssertDimension(rhs.local_range().first, local_range().first);
    AssertDimension(rhs.local_range().second, local_range().second);
    Assert(matrix->RowMap().SameAs(rhs.matrix->RowMap()),
           ExcMessage("Can only add matrices with same distribution of rows"));
    Assert(matrix->Filled() && rhs.matrix->Filled(),
           ExcMessage("Addition of matrices only allowed if matrices are "
                      "filled, i.e., compress() has been called"));

    const bool same_col_map = matrix->ColMap().SameAs(rhs.matrix->ColMap());

    for (const auto row : locally_owned_range_indices())
      {
        const int row_local = matrix->RowMap().LID(
          static_cast<TrilinosWrappers::types::int_type>(row));
        Assert((row_local >= 0), ExcAccessToNonlocalRow(row));

        // First get a view to the matrix columns of both matrices. Note that
        // the data is in local index spaces so we need to be careful not only
        // to compare column indices in case they are derived from the same
        // map.
        int             n_entries, rhs_n_entries;
        TrilinosScalar *value_ptr, *rhs_value_ptr;
        int *           index_ptr, *rhs_index_ptr;
        int             ierr = rhs.matrix->ExtractMyRowView(row_local,
                                                rhs_n_entries,
                                                rhs_value_ptr,
                                                rhs_index_ptr);
        (void)ierr;
        Assert(ierr == 0, ExcTrilinosError(ierr));

        ierr =
          matrix->ExtractMyRowView(row_local, n_entries, value_ptr, index_ptr);
        Assert(ierr == 0, ExcTrilinosError(ierr));
        bool expensive_checks = (n_entries != rhs_n_entries || !same_col_map);
        if (!expensive_checks)
          {
            // check if the column indices are the same. If yes, can simply
            // copy over the data.
            expensive_checks = std::memcmp(static_cast<void *>(index_ptr),
                                           static_cast<void *>(rhs_index_ptr),
                                           sizeof(int) * n_entries) != 0;
            if (!expensive_checks)
              for (int i = 0; i < n_entries; ++i)
                value_ptr[i] += rhs_value_ptr[i] * factor;
          }
        // Now to the expensive case where we need to check all column indices
        // against each other (transformed into global index space) and where
        // we need to make sure that all entries we are about to add into the
        // lhs matrix actually exist
        if (expensive_checks)
          {
            for (int i = 0; i < rhs_n_entries; ++i)
              {
                if (rhs_value_ptr[i] == 0.)
                  continue;
                const TrilinosWrappers::types::int_type rhs_global_col =
                  global_column_index(*rhs.matrix, rhs_index_ptr[i]);
                int  local_col   = matrix->ColMap().LID(rhs_global_col);
                int *local_index = Utilities::lower_bound(index_ptr,
                                                          index_ptr + n_entries,
                                                          local_col);
                Assert(local_index != index_ptr + n_entries &&
                         *local_index == local_col,
                       ExcMessage(
                         "Adding the entries from the other matrix "
                         "failed, because the sparsity pattern "
                         "of that matrix includes more elements than the "
                         "calling matrix, which is not allowed."));
                value_ptr[local_index - index_ptr] += factor * rhs_value_ptr[i];
              }
          }
      }
  }



  void
  SparseMatrix::transpose()
  {
    // This only flips a flag that tells
    // Trilinos that any vmult operation
    // should be done with the
    // transpose. However, the matrix
    // structure is not reset.
    int ierr;

    if (!matrix->UseTranspose())
      {
        ierr = matrix->SetUseTranspose(true);
        AssertThrow(ierr == 0, ExcTrilinosError(ierr));
      }
    else
      {
        ierr = matrix->SetUseTranspose(false);
        AssertThrow(ierr == 0, ExcTrilinosError(ierr));
      }
  }



  SparseMatrix &
  SparseMatrix::operator*=(const TrilinosScalar a)
  {
    const int ierr = matrix->Scale(a);
    Assert(ierr == 0, ExcTrilinosError(ierr));
    (void)ierr; // removes -Wunused-variable in optimized mode

    return *this;
  }



  SparseMatrix &
  SparseMatrix::operator/=(const TrilinosScalar a)
  {
    Assert(a != 0, ExcDivideByZero());

    const TrilinosScalar factor = 1. / a;

    const int ierr = matrix->Scale(factor);
    Assert(ierr == 0, ExcTrilinosError(ierr));
    (void)ierr; // removes -Wunused-variable in optimized mode

    return *this;
  }



  TrilinosScalar
  SparseMatrix::l1_norm() const
  {
    Assert(matrix->Filled(), ExcMatrixNotCompressed());
    return matrix->NormOne();
  }



  TrilinosScalar
  SparseMatrix::linfty_norm() const
  {
    Assert(matrix->Filled(), ExcMatrixNotCompressed());
    return matrix->NormInf();
  }



  TrilinosScalar
  SparseMatrix::frobenius_norm() const
  {
    Assert(matrix->Filled(), ExcMatrixNotCompressed());
    return matrix->NormFrobenius();
  }



  namespace internal
  {
    namespace SparseMatrixImplementation
    {
      template <typename VectorType>
      inline void
      check_vector_map_equality(const Epetra_CrsMatrix &,
                                const VectorType &,
                                const VectorType &)
      {}

      inline void
      check_vector_map_equality(const Epetra_CrsMatrix &             m,
                                const TrilinosWrappers::MPI::Vector &in,
                                const TrilinosWrappers::MPI::Vector &out)
      {
        Assert(in.trilinos_partitioner().SameAs(m.DomainMap()) == true,
               ExcMessage(
                 "Column map of matrix does not fit with vector map!"));
        Assert(out.trilinos_partitioner().SameAs(m.RangeMap()) == true,
               ExcMessage("Row map of matrix does not fit with vector map!"));
        (void)m;
        (void)in;
        (void)out;
      }
    } // namespace SparseMatrixImplementation
  }   // namespace internal


  template <typename VectorType>
  typename std::enable_if<
    std::is_same<typename VectorType::value_type, TrilinosScalar>::value>::type
  SparseMatrix::vmult(VectorType &dst, const VectorType &src) const
  {
    Assert(&src != &dst, ExcSourceEqualsDestination());
    Assert(matrix->Filled(), ExcMatrixNotCompressed());
    (void)src;
    (void)dst;

    internal::SparseMatrixImplementation::check_vector_map_equality(*matrix,
                                                                    src,
                                                                    dst);
    const size_type dst_local_size = internal::end(dst) - internal::begin(dst);
    AssertDimension(dst_local_size, matrix->RangeMap().NumMyPoints());
    const size_type src_local_size = internal::end(src) - internal::begin(src);
    AssertDimension(src_local_size, matrix->DomainMap().NumMyPoints());

    Epetra_MultiVector tril_dst(
      View, matrix->RangeMap(), internal::begin(dst), dst_local_size, 1);
    Epetra_MultiVector tril_src(View,
                                matrix->DomainMap(),
                                const_cast<TrilinosScalar *>(
                                  internal::begin(src)),
                                src_local_size,
                                1);

    const int ierr = matrix->Multiply(false, tril_src, tril_dst);
    Assert(ierr == 0, ExcTrilinosError(ierr));
    (void)ierr; // removes -Wunused-variable in optimized mode
  }



  template <typename VectorType>
  typename std::enable_if<
    !std::is_same<typename VectorType::value_type, TrilinosScalar>::value>::type
  SparseMatrix::vmult(VectorType & /*dst*/, const VectorType & /*src*/) const
  {
    AssertThrow(false, ExcNotImplemented());
  }



  template <typename VectorType>
  typename std::enable_if<
    std::is_same<typename VectorType::value_type, TrilinosScalar>::value>::type
  SparseMatrix::Tvmult(VectorType &dst, const VectorType &src) const
  {
    Assert(&src != &dst, ExcSourceEqualsDestination());
    Assert(matrix->Filled(), ExcMatrixNotCompressed());

    internal::SparseMatrixImplementation::check_vector_map_equality(*matrix,
                                                                    dst,
                                                                    src);
    const size_type dst_local_size = internal::end(dst) - internal::begin(dst);
    AssertDimension(dst_local_size, matrix->DomainMap().NumMyPoints());
    const size_type src_local_size = internal::end(src) - internal::begin(src);
    AssertDimension(src_local_size, matrix->RangeMap().NumMyPoints());

    Epetra_MultiVector tril_dst(
      View, matrix->DomainMap(), internal::begin(dst), dst_local_size, 1);
    Epetra_MultiVector tril_src(View,
                                matrix->RangeMap(),
                                const_cast<double *>(internal::begin(src)),
                                src_local_size,
                                1);

    const int ierr = matrix->Multiply(true, tril_src, tril_dst);
    Assert(ierr == 0, ExcTrilinosError(ierr));
    (void)ierr; // removes -Wunused-variable in optimized mode
  }



  template <typename VectorType>
  typename std::enable_if<
    !std::is_same<typename VectorType::value_type, TrilinosScalar>::value>::type
  SparseMatrix::Tvmult(VectorType & /*dst*/, const VectorType & /*src*/) const
  {
    AssertThrow(false, ExcNotImplemented());
  }



  template <typename VectorType>
  void
  SparseMatrix::vmult_add(VectorType &dst, const VectorType &src) const
  {
    Assert(&src != &dst, ExcSourceEqualsDestination());

    // Reinit a temporary vector with fast argument set, which does not
    // overwrite the content (to save time).
    VectorType tmp_vector;
    tmp_vector.reinit(dst, true);
    vmult(tmp_vector, src);
    dst += tmp_vector;
  }



  template <typename VectorType>
  void
  SparseMatrix::Tvmult_add(VectorType &dst, const VectorType &src) const
  {
    Assert(&src != &dst, ExcSourceEqualsDestination());

    // Reinit a temporary vector with fast argument set, which does not
    // overwrite the content (to save time).
    VectorType tmp_vector;
    tmp_vector.reinit(dst, true);
    Tvmult(tmp_vector, src);
    dst += tmp_vector;
  }



  TrilinosScalar
  SparseMatrix::matrix_norm_square(const MPI::Vector &v) const
  {
    Assert(matrix->RowMap().SameAs(matrix->DomainMap()), ExcNotQuadratic());

    MPI::Vector temp_vector;
    temp_vector.reinit(v, true);

    vmult(temp_vector, v);
    return temp_vector * v;
  }



  TrilinosScalar
  SparseMatrix::matrix_scalar_product(const MPI::Vector &u,
                                      const MPI::Vector &v) const
  {
    Assert(matrix->RowMap().SameAs(matrix->DomainMap()), ExcNotQuadratic());

    MPI::Vector temp_vector;
    temp_vector.reinit(v, true);

    vmult(temp_vector, v);
    return u * temp_vector;
  }



  TrilinosScalar
  SparseMatrix::residual(MPI::Vector &      dst,
                         const MPI::Vector &x,
                         const MPI::Vector &b) const
  {
    vmult(dst, x);
    dst -= b;
    dst *= -1.;

    return dst.l2_norm();
  }



  namespace internals
  {
    using size_type = dealii::types::global_dof_index;

    void
    perform_mmult(const SparseMatrix &inputleft,
                  const SparseMatrix &inputright,
                  SparseMatrix &      result,
                  const MPI::Vector & V,
                  const bool          transpose_left)
    {
      const bool use_vector = (V.size() == inputright.m() ? true : false);
      if (transpose_left == false)
        {
          Assert(inputleft.n() == inputright.m(),
                 ExcDimensionMismatch(inputleft.n(), inputright.m()));
          Assert(inputleft.trilinos_matrix().DomainMap().SameAs(
                   inputright.trilinos_matrix().RangeMap()),
                 ExcMessage("Parallel partitioning of A and B does not fit."));
        }
      else
        {
          Assert(inputleft.m() == inputright.m(),
                 ExcDimensionMismatch(inputleft.m(), inputright.m()));
          Assert(inputleft.trilinos_matrix().RangeMap().SameAs(
                   inputright.trilinos_matrix().RangeMap()),
                 ExcMessage("Parallel partitioning of A and B does not fit."));
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
          mod_B = Teuchos::rcp(const_cast<Epetra_CrsMatrix *>(
                                 &inputright.trilinos_matrix()),
                               false);
        }
      else
        {
          mod_B = Teuchos::rcp(
            new Epetra_CrsMatrix(Copy, inputright.trilinos_sparsity_pattern()),
            true);
          mod_B->FillComplete(inputright.trilinos_matrix().DomainMap(),
                              inputright.trilinos_matrix().RangeMap());
          Assert(inputright.local_range() == V.local_range(),
                 ExcMessage("Parallel distribution of matrix B and vector V "
                            "does not match."));

          const int local_N = inputright.local_size();
          for (int i = 0; i < local_N; ++i)
            {
              int     N_entries = -1;
              double *new_data, *B_data;
              mod_B->ExtractMyRowView(i, N_entries, new_data);
              inputright.trilinos_matrix().ExtractMyRowView(i,
                                                            N_entries,
                                                            B_data);
              double value = V.trilinos_vector()[0][i];
              for (TrilinosWrappers::types::int_type j = 0; j < N_entries; ++j)
                new_data[j] = value * B_data[j];
            }
        }


      SparseMatrix tmp_result(transpose_left ?
                                inputleft.locally_owned_domain_indices() :
                                inputleft.locally_owned_range_indices(),
                              inputright.locally_owned_domain_indices(),
                              inputleft.get_mpi_communicator());

#  ifdef DEAL_II_TRILINOS_WITH_EPETRAEXT
      EpetraExt::MatrixMatrix::Multiply(inputleft.trilinos_matrix(),
                                        transpose_left,
                                        *mod_B,
                                        false,
                                        const_cast<Epetra_CrsMatrix &>(
                                          tmp_result.trilinos_matrix()));
#  else
      Assert("false", ExcMessage("This function requires EpetraExt."));
#  endif
      result.reinit(tmp_result.trilinos_matrix());
    }
  } // namespace internals


  void
  SparseMatrix::mmult(SparseMatrix &      C,
                      const SparseMatrix &B,
                      const MPI::Vector & V) const
  {
    internals::perform_mmult(*this, B, C, V, false);
  }



  void
  SparseMatrix::Tmmult(SparseMatrix &      C,
                       const SparseMatrix &B,
                       const MPI::Vector & V) const
  {
    internals::perform_mmult(*this, B, C, V, true);
  }



  void
  SparseMatrix::write_ascii()
  {
    Assert(false, ExcNotImplemented());
  }



  // As of now, no particularly neat
  // output is generated in case of
  // multiple processors.
  void
  SparseMatrix::print(std::ostream &out,
                      const bool    print_detailed_trilinos_information) const
  {
    if (print_detailed_trilinos_information == true)
      out << *matrix;
    else
      {
        double *values;
        int *   indices;
        int     num_entries;

        for (int i = 0; i < matrix->NumMyRows(); ++i)
          {
            const int ierr =
              matrix->ExtractMyRowView(i, num_entries, values, indices);
            (void)ierr;
            Assert(ierr == 0, ExcTrilinosError(ierr));

            for (TrilinosWrappers::types::int_type j = 0; j < num_entries; ++j)
              out << "(" << TrilinosWrappers::global_row_index(*matrix, i)
                  << ","
                  << TrilinosWrappers::global_column_index(*matrix, indices[j])
                  << ") " << values[j] << std::endl;
          }
      }

    AssertThrow(out, ExcIO());
  }



  SparseMatrix::size_type
  SparseMatrix::memory_consumption() const
  {
    size_type static_memory =
      sizeof(*this) + sizeof(*matrix) + sizeof(*matrix->Graph().DataPtr());
    return (
      (sizeof(TrilinosScalar) + sizeof(TrilinosWrappers::types::int_type)) *
        matrix->NumMyNonzeros() +
      sizeof(int) * local_size() + static_memory);
  }



  MPI_Comm
  SparseMatrix::get_mpi_communicator() const
  {
#  ifdef DEAL_II_WITH_MPI

    const Epetra_MpiComm *mpi_comm =
      dynamic_cast<const Epetra_MpiComm *>(&matrix->RangeMap().Comm());
    Assert(mpi_comm != nullptr, ExcInternalError());
    return mpi_comm->Comm();
#  else

    return MPI_COMM_SELF;

#  endif
  }
} // namespace TrilinosWrappers


namespace TrilinosWrappers
{
  namespace internal
  {
    namespace
    {
#  ifndef DEAL_II_WITH_MPI
      Epetra_Map
      make_serial_Epetra_map(const IndexSet &serial_partitioning)
      {
        // See IndexSet::make_trilinos_map
        return Epetra_Map(
          TrilinosWrappers::types::int_type(serial_partitioning.size()),
          TrilinosWrappers::types::int_type(serial_partitioning.n_elements()),
          0,
          Epetra_SerialComm());
      }
#  endif
    } // namespace

    namespace LinearOperatorImplementation
    {
      TrilinosPayload::TrilinosPayload()
        : use_transpose(false)
        ,
#  ifdef DEAL_II_WITH_MPI
        communicator(MPI_COMM_SELF)
        , domain_map(IndexSet().make_trilinos_map(communicator.Comm()))
        , range_map(IndexSet().make_trilinos_map(communicator.Comm()))
#  else
        domain_map(internal::make_serial_Epetra_map(IndexSet()))
        , range_map(internal::make_serial_Epetra_map(IndexSet()))
#  endif
      {
        vmult = [](Range &, const Domain &) {
          Assert(false,
                 ExcMessage("Uninitialized TrilinosPayload::vmult called "
                            "(Default constructor)"));
        };

        Tvmult = [](Domain &, const Range &) {
          Assert(false,
                 ExcMessage("Uninitialized TrilinosPayload::Tvmult called "
                            "(Default constructor)"));
        };

        inv_vmult = [](Domain &, const Range &) {
          Assert(false,
                 ExcMessage("Uninitialized TrilinosPayload::inv_vmult called "
                            "(Default constructor)"));
        };

        inv_Tvmult = [](Range &, const Domain &) {
          Assert(false,
                 ExcMessage("Uninitialized TrilinosPayload::inv_Tvmult called "
                            "(Default constructor)"));
        };
      }



      TrilinosPayload::TrilinosPayload(
        const TrilinosWrappers::SparseMatrix &matrix_exemplar,
        const TrilinosWrappers::SparseMatrix &matrix)
        : use_transpose(matrix_exemplar.trilinos_matrix().UseTranspose())
        ,
#  ifdef DEAL_II_WITH_MPI
        communicator(matrix_exemplar.get_mpi_communicator())
        , domain_map(
            matrix_exemplar.locally_owned_domain_indices().make_trilinos_map(
              communicator.Comm()))
        , range_map(
            matrix_exemplar.locally_owned_range_indices().make_trilinos_map(
              communicator.Comm()))
#  else
        domain_map(internal::make_serial_Epetra_map(
          matrix_exemplar.locally_owned_domain_indices()))
        , range_map(internal::make_serial_Epetra_map(
            matrix_exemplar.locally_owned_range_indices()))
#  endif
      {
        vmult = [&matrix_exemplar, &matrix](Range &       tril_dst,
                                            const Domain &tril_src) {
          // Duplicated from TrilinosWrappers::SparseMatrix::vmult
          Assert(&tril_src != &tril_dst,
                 TrilinosWrappers::SparseMatrix::ExcSourceEqualsDestination());
          Assert(matrix.trilinos_matrix().Filled(),
                 TrilinosWrappers::SparseMatrix::ExcMatrixNotCompressed());
          internal::check_vector_map_equality(
            matrix_exemplar.trilinos_matrix(),
            tril_src,
            tril_dst,
            matrix_exemplar.trilinos_matrix().UseTranspose());
          internal::check_vector_map_equality(
            matrix.trilinos_matrix(),
            tril_src,
            tril_dst,
            matrix.trilinos_matrix().UseTranspose());

          const int ierr = matrix.trilinos_matrix().Apply(tril_src, tril_dst);
          AssertThrow(ierr == 0, ExcTrilinosError(ierr));
        };

        Tvmult = [&matrix_exemplar, &matrix](Domain &     tril_dst,
                                             const Range &tril_src) {
          // Duplicated from TrilinosWrappers::SparseMatrix::Tvmult
          Assert(&tril_src != &tril_dst,
                 TrilinosWrappers::SparseMatrix::ExcSourceEqualsDestination());
          Assert(matrix.trilinos_matrix().Filled(),
                 TrilinosWrappers::SparseMatrix::ExcMatrixNotCompressed());
          internal::check_vector_map_equality(
            matrix_exemplar.trilinos_matrix(),
            tril_src,
            tril_dst,
            !matrix_exemplar.trilinos_matrix().UseTranspose());
          internal::check_vector_map_equality(
            matrix.trilinos_matrix(),
            tril_src,
            tril_dst,
            !matrix.trilinos_matrix().UseTranspose());

          Epetra_CrsMatrix &tril_mtrx_non_const =
            const_cast<Epetra_CrsMatrix &>(matrix.trilinos_matrix());
          tril_mtrx_non_const.SetUseTranspose(
            !matrix.trilinos_matrix().UseTranspose());
          const int ierr = matrix.trilinos_matrix().Apply(tril_src, tril_dst);
          AssertThrow(ierr == 0, ExcTrilinosError(ierr));
          tril_mtrx_non_const.SetUseTranspose(
            !matrix.trilinos_matrix().UseTranspose());
        };

        inv_vmult = [](Domain &, const Range &) {
          Assert(false,
                 ExcMessage("Uninitialized TrilinosPayload::inv_vmult called "
                            "(Matrix constructor with matrix exemplar)"));
        };

        inv_Tvmult = [](Range &, const Domain &) {
          Assert(false,
                 ExcMessage("Uninitialized TrilinosPayload::inv_Tvmult called "
                            "(Matrix constructor with matrix exemplar)"));
        };
      }



      TrilinosPayload::TrilinosPayload(
        const TrilinosWrappers::SparseMatrix &    matrix_exemplar,
        const TrilinosWrappers::PreconditionBase &preconditioner)
        : use_transpose(matrix_exemplar.trilinos_matrix().UseTranspose())
        ,
#  ifdef DEAL_II_WITH_MPI
        communicator(matrix_exemplar.get_mpi_communicator())
        , domain_map(
            matrix_exemplar.locally_owned_domain_indices().make_trilinos_map(
              communicator.Comm()))
        , range_map(
            matrix_exemplar.locally_owned_range_indices().make_trilinos_map(
              communicator.Comm()))
#  else
        domain_map(internal::make_serial_Epetra_map(
          matrix_exemplar.locally_owned_domain_indices()))
        , range_map(internal::make_serial_Epetra_map(
            matrix_exemplar.locally_owned_range_indices()))
#  endif
      {
        vmult = [&matrix_exemplar, &preconditioner](Range &       tril_dst,
                                                    const Domain &tril_src) {
          // Duplicated from TrilinosWrappers::PreconditionBase::vmult
          // as well as from TrilinosWrappers::SparseMatrix::Tvmult
          Assert(&tril_src != &tril_dst,
                 TrilinosWrappers::SparseMatrix::ExcSourceEqualsDestination());
          internal::check_vector_map_equality(
            matrix_exemplar.trilinos_matrix(),
            tril_src,
            tril_dst,
            matrix_exemplar.trilinos_matrix().UseTranspose());
          internal::check_vector_map_equality(
            preconditioner.trilinos_operator(),
            tril_src,
            tril_dst,
            preconditioner.trilinos_operator().UseTranspose());

          const int ierr =
            preconditioner.trilinos_operator().ApplyInverse(tril_src, tril_dst);
          AssertThrow(ierr == 0, ExcTrilinosError(ierr));
        };

        Tvmult = [&matrix_exemplar, &preconditioner](Domain &     tril_dst,
                                                     const Range &tril_src) {
          // Duplicated from TrilinosWrappers::PreconditionBase::vmult
          // as well as from TrilinosWrappers::SparseMatrix::Tvmult
          Assert(&tril_src != &tril_dst,
                 TrilinosWrappers::SparseMatrix::ExcSourceEqualsDestination());
          internal::check_vector_map_equality(
            matrix_exemplar.trilinos_matrix(),
            tril_src,
            tril_dst,
            !matrix_exemplar.trilinos_matrix().UseTranspose());
          internal::check_vector_map_equality(
            preconditioner.trilinos_operator(),
            tril_src,
            tril_dst,
            !preconditioner.trilinos_operator().UseTranspose());

          preconditioner.trilinos_operator().SetUseTranspose(
            !preconditioner.trilinos_operator().UseTranspose());
          const int ierr =
            preconditioner.trilinos_operator().ApplyInverse(tril_src, tril_dst);
          AssertThrow(ierr == 0, ExcTrilinosError(ierr));
          preconditioner.trilinos_operator().SetUseTranspose(
            !preconditioner.trilinos_operator().UseTranspose());
        };

        inv_vmult = [&matrix_exemplar, &preconditioner](Domain &     tril_dst,
                                                        const Range &tril_src) {
          // Duplicated from TrilinosWrappers::PreconditionBase::vmult
          // as well as from TrilinosWrappers::SparseMatrix::Tvmult
          Assert(&tril_src != &tril_dst,
                 TrilinosWrappers::SparseMatrix::ExcSourceEqualsDestination());
          internal::check_vector_map_equality(
            matrix_exemplar.trilinos_matrix(),
            tril_src,
            tril_dst,
            !matrix_exemplar.trilinos_matrix().UseTranspose());
          internal::check_vector_map_equality(
            preconditioner.trilinos_operator(),
            tril_src,
            tril_dst,
            !preconditioner.trilinos_operator().UseTranspose());

          const int ierr =
            preconditioner.trilinos_operator().ApplyInverse(tril_src, tril_dst);
          AssertThrow(ierr == 0, ExcTrilinosError(ierr));
        };

        inv_Tvmult = [&matrix_exemplar,
                      &preconditioner](Range &       tril_dst,
                                       const Domain &tril_src) {
          // Duplicated from TrilinosWrappers::PreconditionBase::vmult
          // as well as from TrilinosWrappers::SparseMatrix::Tvmult
          Assert(&tril_src != &tril_dst,
                 TrilinosWrappers::SparseMatrix::ExcSourceEqualsDestination());
          internal::check_vector_map_equality(
            matrix_exemplar.trilinos_matrix(),
            tril_src,
            tril_dst,
            matrix_exemplar.trilinos_matrix().UseTranspose());
          internal::check_vector_map_equality(
            preconditioner.trilinos_operator(),
            tril_src,
            tril_dst,
            preconditioner.trilinos_operator().UseTranspose());

          preconditioner.trilinos_operator().SetUseTranspose(
            !preconditioner.trilinos_operator().UseTranspose());
          const int ierr =
            preconditioner.trilinos_operator().ApplyInverse(tril_src, tril_dst);
          AssertThrow(ierr == 0, ExcTrilinosError(ierr));
          preconditioner.trilinos_operator().SetUseTranspose(
            !preconditioner.trilinos_operator().UseTranspose());
        };
      }



      TrilinosPayload::TrilinosPayload(
        const TrilinosWrappers::PreconditionBase &preconditioner_exemplar,
        const TrilinosWrappers::PreconditionBase &preconditioner)
        : use_transpose(
            preconditioner_exemplar.trilinos_operator().UseTranspose())
        ,
#  ifdef DEAL_II_WITH_MPI
        communicator(preconditioner_exemplar.get_mpi_communicator())
        , domain_map(preconditioner_exemplar.locally_owned_domain_indices()
                       .make_trilinos_map(communicator.Comm()))
        , range_map(preconditioner_exemplar.locally_owned_range_indices()
                      .make_trilinos_map(communicator.Comm()))
#  else
        domain_map(internal::make_serial_Epetra_map(
          preconditioner_exemplar.locally_owned_domain_indices()))
        , range_map(internal::make_serial_Epetra_map(
            preconditioner_exemplar.locally_owned_range_indices()))
#  endif
      {
        vmult = [&preconditioner_exemplar,
                 &preconditioner](Range &tril_dst, const Domain &tril_src) {
          // Duplicated from TrilinosWrappers::PreconditionBase::vmult
          // as well as from TrilinosWrappers::SparseMatrix::Tvmult
          Assert(&tril_src != &tril_dst,
                 TrilinosWrappers::SparseMatrix::ExcSourceEqualsDestination());
          internal::check_vector_map_equality(
            preconditioner_exemplar.trilinos_operator(),
            tril_src,
            tril_dst,
            preconditioner_exemplar.trilinos_operator().UseTranspose());
          internal::check_vector_map_equality(
            preconditioner.trilinos_operator(),
            tril_src,
            tril_dst,
            preconditioner.trilinos_operator().UseTranspose());

          const int ierr =
            preconditioner.trilinos_operator().Apply(tril_src, tril_dst);
          AssertThrow(ierr == 0, ExcTrilinosError(ierr));
        };

        Tvmult = [&preconditioner_exemplar,
                  &preconditioner](Domain &tril_dst, const Range &tril_src) {
          // Duplicated from TrilinosWrappers::PreconditionBase::vmult
          // as well as from TrilinosWrappers::SparseMatrix::Tvmult
          Assert(&tril_src != &tril_dst,
                 TrilinosWrappers::SparseMatrix::ExcSourceEqualsDestination());
          internal::check_vector_map_equality(
            preconditioner_exemplar.trilinos_operator(),
            tril_src,
            tril_dst,
            !preconditioner_exemplar.trilinos_operator().UseTranspose());
          internal::check_vector_map_equality(
            preconditioner.trilinos_operator(),
            tril_src,
            tril_dst,
            !preconditioner.trilinos_operator().UseTranspose());

          preconditioner.trilinos_operator().SetUseTranspose(
            !preconditioner.trilinos_operator().UseTranspose());
          const int ierr =
            preconditioner.trilinos_operator().Apply(tril_src, tril_dst);
          AssertThrow(ierr == 0, ExcTrilinosError(ierr));
          preconditioner.trilinos_operator().SetUseTranspose(
            !preconditioner.trilinos_operator().UseTranspose());
        };

        inv_vmult = [&preconditioner_exemplar,
                     &preconditioner](Domain &tril_dst, const Range &tril_src) {
          // Duplicated from TrilinosWrappers::PreconditionBase::vmult
          // as well as from TrilinosWrappers::SparseMatrix::Tvmult
          Assert(&tril_src != &tril_dst,
                 TrilinosWrappers::SparseMatrix::ExcSourceEqualsDestination());
          internal::check_vector_map_equality(
            preconditioner_exemplar.trilinos_operator(),
            tril_src,
            tril_dst,
            !preconditioner_exemplar.trilinos_operator().UseTranspose());
          internal::check_vector_map_equality(
            preconditioner.trilinos_operator(),
            tril_src,
            tril_dst,
            !preconditioner.trilinos_operator().UseTranspose());

          const int ierr =
            preconditioner.trilinos_operator().ApplyInverse(tril_src, tril_dst);
          AssertThrow(ierr == 0, ExcTrilinosError(ierr));
        };

        inv_Tvmult = [&preconditioner_exemplar,
                      &preconditioner](Range &       tril_dst,
                                       const Domain &tril_src) {
          // Duplicated from TrilinosWrappers::PreconditionBase::vmult
          // as well as from TrilinosWrappers::SparseMatrix::Tvmult
          Assert(&tril_src != &tril_dst,
                 TrilinosWrappers::SparseMatrix::ExcSourceEqualsDestination());
          internal::check_vector_map_equality(
            preconditioner_exemplar.trilinos_operator(),
            tril_src,
            tril_dst,
            preconditioner_exemplar.trilinos_operator().UseTranspose());
          internal::check_vector_map_equality(
            preconditioner.trilinos_operator(),
            tril_src,
            tril_dst,
            preconditioner.trilinos_operator().UseTranspose());

          preconditioner.trilinos_operator().SetUseTranspose(
            !preconditioner.trilinos_operator().UseTranspose());
          const int ierr =
            preconditioner.trilinos_operator().ApplyInverse(tril_src, tril_dst);
          AssertThrow(ierr == 0, ExcTrilinosError(ierr));
          preconditioner.trilinos_operator().SetUseTranspose(
            !preconditioner.trilinos_operator().UseTranspose());
        };
      }



      TrilinosPayload::TrilinosPayload(const TrilinosPayload &payload)
        : vmult(payload.vmult)
        , Tvmult(payload.Tvmult)
        , inv_vmult(payload.inv_vmult)
        , inv_Tvmult(payload.inv_Tvmult)
        , use_transpose(payload.use_transpose)
        , communicator(payload.communicator)
        , domain_map(payload.domain_map)
        , range_map(payload.range_map)
      {}



      // Composite copy constructor
      // This is required for PackagedOperations
      TrilinosPayload::TrilinosPayload(const TrilinosPayload &first_op,
                                       const TrilinosPayload &second_op)
        : use_transpose(false)
        , // The combination of operators provides the exact
          // definition of the operation
        communicator(first_op.communicator)
        , domain_map(second_op.domain_map)
        , range_map(first_op.range_map)
      {}



      TrilinosPayload
      TrilinosPayload::identity_payload() const
      {
        TrilinosPayload return_op(*this);

        return_op.vmult = [](Range &tril_dst, const Range &tril_src) {
          tril_dst = tril_src;
        };

        return_op.Tvmult = [](Range &tril_dst, const Range &tril_src) {
          tril_dst = tril_src;
        };

        return_op.inv_vmult = [](Range &tril_dst, const Range &tril_src) {
          tril_dst = tril_src;
        };

        return_op.inv_Tvmult = [](Range &tril_dst, const Range &tril_src) {
          tril_dst = tril_src;
        };

        return return_op;
      }



      TrilinosPayload
      TrilinosPayload::null_payload() const
      {
        TrilinosPayload return_op(*this);

        return_op.vmult = [](Range &tril_dst, const Domain &) {
          const int ierr = tril_dst.PutScalar(0.0);

          AssertThrow(ierr == 0, ExcTrilinosError(ierr));
        };

        return_op.Tvmult = [](Domain &tril_dst, const Range &) {
          const int ierr = tril_dst.PutScalar(0.0);

          AssertThrow(ierr == 0, ExcTrilinosError(ierr));
        };

        return_op.inv_vmult = [](Domain &tril_dst, const Range &) {
          AssertThrow(false,
                      ExcMessage("Cannot compute inverse of null operator"));

          const int ierr = tril_dst.PutScalar(0.0);

          AssertThrow(ierr == 0, ExcTrilinosError(ierr));
        };

        return_op.inv_Tvmult = [](Range &tril_dst, const Domain &) {
          AssertThrow(false,
                      ExcMessage("Cannot compute inverse of null operator"));

          const int ierr = tril_dst.PutScalar(0.0);

          AssertThrow(ierr == 0, ExcTrilinosError(ierr));
        };

        return return_op;
      }



      TrilinosPayload
      TrilinosPayload::transpose_payload() const
      {
        TrilinosPayload return_op(*this);
        return_op.transpose();
        return return_op;
      }



      IndexSet
      TrilinosPayload::locally_owned_domain_indices() const
      {
        return IndexSet(domain_map);
      }



      IndexSet
      TrilinosPayload::locally_owned_range_indices() const
      {
        return IndexSet(range_map);
      }



      MPI_Comm
      TrilinosPayload::get_mpi_communicator() const
      {
#  ifdef DEAL_II_WITH_MPI
        return communicator.Comm();
#  else
        return MPI_COMM_SELF;
#  endif
      }



      void
      TrilinosPayload::transpose()
      {
        SetUseTranspose(!use_transpose);
      }



      bool
      TrilinosPayload::UseTranspose() const
      {
        return use_transpose;
      }



      int
      TrilinosPayload::SetUseTranspose(bool UseTranspose)
      {
        if (use_transpose != UseTranspose)
          {
            use_transpose = UseTranspose;
            std::swap(domain_map, range_map);
            std::swap(vmult, Tvmult);
            std::swap(inv_vmult, inv_Tvmult);
          }
        return 0;
      }



      int
      TrilinosPayload::Apply(const VectorType &X, VectorType &Y) const
      {
        // The transposedness of the operations is taken care of
        // when we hit the transpose flag.
        vmult(Y, X);
        return 0;
      }



      int
      TrilinosPayload::ApplyInverse(const VectorType &Y, VectorType &X) const
      {
        // The transposedness of the operations is taken care of
        // when we hit the transpose flag.
        inv_vmult(X, Y);
        return 0;
      }



      const char *
      TrilinosPayload::Label() const
      {
        return "TrilinosPayload";
      }



      const Epetra_Comm &
      TrilinosPayload::Comm() const
      {
        return communicator;
      }



      const Epetra_Map &
      TrilinosPayload::OperatorDomainMap() const
      {
        return domain_map;
      }



      const Epetra_Map &
      TrilinosPayload::OperatorRangeMap() const
      {
        return range_map;
      }



      bool
      TrilinosPayload::HasNormInf() const
      {
        return false;
      }



      double
      TrilinosPayload::NormInf() const
      {
        AssertThrow(false, ExcNotImplemented());
        return 0.0;
      }



      TrilinosPayload
      operator+(const TrilinosPayload &first_op,
                const TrilinosPayload &second_op)
      {
        using Domain        = typename TrilinosPayload::Domain;
        using Range         = typename TrilinosPayload::Range;
        using Intermediate  = typename TrilinosPayload::VectorType;
        using GVMVectorType = TrilinosWrappers::MPI::Vector;

        Assert(first_op.locally_owned_domain_indices() ==
                 second_op.locally_owned_domain_indices(),
               ExcMessage(
                 "Operators are set to work on incompatible IndexSets."));
        Assert(first_op.locally_owned_range_indices() ==
                 second_op.locally_owned_range_indices(),
               ExcMessage(
                 "Operators are set to work on incompatible IndexSets."));

        TrilinosPayload return_op(first_op, second_op);

        // Capture by copy so the payloads are always valid
        return_op.vmult = [first_op, second_op](Range &       tril_dst,
                                                const Domain &tril_src) {
          // Duplicated from LinearOperator::operator*
          // TODO: Template the constructor on GrowingVectorMemory vector type?
          GrowingVectorMemory<GVMVectorType>   vector_memory;
          VectorMemory<GVMVectorType>::Pointer i(vector_memory);

          // Initialize intermediate vector
          const Epetra_Map &first_op_init_map = first_op.OperatorDomainMap();
          i->reinit(IndexSet(first_op_init_map),
                    first_op.get_mpi_communicator(),
                    /*bool omit_zeroing_entries =*/true);

          // Duplicated from TrilinosWrappers::SparseMatrix::vmult
          const size_type i_local_size = i->end() - i->begin();
          AssertDimension(i_local_size, first_op_init_map.NumMyPoints());
          const Epetra_Map &second_op_init_map = second_op.OperatorDomainMap();
          AssertDimension(i_local_size, second_op_init_map.NumMyPoints());
          (void)second_op_init_map;
          Intermediate tril_int(View,
                                first_op_init_map,
                                const_cast<TrilinosScalar *>(i->begin()),
                                i_local_size,
                                1);

          // These operators may themselves be transposed or not, so we let them
          // decide what the intended outcome is
          second_op.Apply(tril_src, tril_int);
          first_op.Apply(tril_src, tril_dst);
          const int ierr = tril_dst.Update(1.0, tril_int, 1.0);
          AssertThrow(ierr == 0, ExcTrilinosError(ierr));
        };

        return_op.Tvmult = [first_op, second_op](Domain &     tril_dst,
                                                 const Range &tril_src) {
          // Duplicated from LinearOperator::operator*
          // TODO: Template the constructor on GrowingVectorMemory vector type?
          GrowingVectorMemory<GVMVectorType>   vector_memory;
          VectorMemory<GVMVectorType>::Pointer i(vector_memory);

          // These operators may themselves be transposed or not, so we let them
          // decide what the intended outcome is
          // We must first transpose the operators to get the right IndexSets
          // for the input, intermediate and result vectors
          const_cast<TrilinosPayload &>(first_op).transpose();
          const_cast<TrilinosPayload &>(second_op).transpose();

          // Initialize intermediate vector
          const Epetra_Map &first_op_init_map = first_op.OperatorRangeMap();
          i->reinit(IndexSet(first_op_init_map),
                    first_op.get_mpi_communicator(),
                    /*bool omit_zeroing_entries =*/true);

          // Duplicated from TrilinosWrappers::SparseMatrix::vmult
          const size_type i_local_size = i->end() - i->begin();
          AssertDimension(i_local_size, first_op_init_map.NumMyPoints());
          const Epetra_Map &second_op_init_map = second_op.OperatorRangeMap();
          AssertDimension(i_local_size, second_op_init_map.NumMyPoints());
          (void)second_op_init_map;
          Intermediate tril_int(View,
                                first_op_init_map,
                                const_cast<TrilinosScalar *>(i->begin()),
                                i_local_size,
                                1);

          // These operators may themselves be transposed or not, so we let them
          // decide what the intended outcome is
          second_op.Apply(tril_src, tril_int);
          first_op.Apply(tril_src, tril_dst);
          const int ierr = tril_dst.Update(1.0, tril_int, 1.0);
          AssertThrow(ierr == 0, ExcTrilinosError(ierr));

          // Reset transpose flag
          const_cast<TrilinosPayload &>(first_op).transpose();
          const_cast<TrilinosPayload &>(second_op).transpose();
        };

        return_op.inv_vmult = [first_op, second_op](Domain &     tril_dst,
                                                    const Range &tril_src) {
          // Duplicated from LinearOperator::operator*
          // TODO: Template the constructor on GrowingVectorMemory vector type?
          GrowingVectorMemory<GVMVectorType>   vector_memory;
          VectorMemory<GVMVectorType>::Pointer i(vector_memory);

          // Initialize intermediate vector
          const Epetra_Map &first_op_init_map = first_op.OperatorRangeMap();
          i->reinit(IndexSet(first_op_init_map),
                    first_op.get_mpi_communicator(),
                    /*bool omit_zeroing_entries =*/true);

          // Duplicated from TrilinosWrappers::SparseMatrix::vmult
          const size_type i_local_size = i->end() - i->begin();
          AssertDimension(i_local_size, first_op_init_map.NumMyPoints());
          const Epetra_Map &second_op_init_map = second_op.OperatorRangeMap();
          AssertDimension(i_local_size, second_op_init_map.NumMyPoints());
          (void)second_op_init_map;
          Intermediate tril_int(View,
                                first_op_init_map,
                                const_cast<TrilinosScalar *>(i->begin()),
                                i_local_size,
                                1);

          // These operators may themselves be transposed or not, so we let them
          // decide what the intended outcome is
          second_op.ApplyInverse(tril_src, tril_int);
          first_op.ApplyInverse(tril_src, tril_dst);
          const int ierr = tril_dst.Update(1.0, tril_int, 1.0);
          AssertThrow(ierr == 0, ExcTrilinosError(ierr));
        };

        return_op.inv_Tvmult = [first_op, second_op](Range &       tril_dst,
                                                     const Domain &tril_src) {
          // Duplicated from LinearOperator::operator*
          // TODO: Template the constructor on GrowingVectorMemory vector type?
          GrowingVectorMemory<GVMVectorType>   vector_memory;
          VectorMemory<GVMVectorType>::Pointer i(vector_memory);

          // These operators may themselves be transposed or not, so we let them
          // decide what the intended outcome is
          // We must first transpose the operators to get the right IndexSets
          // for the input, intermediate and result vectors
          const_cast<TrilinosPayload &>(first_op).transpose();
          const_cast<TrilinosPayload &>(second_op).transpose();

          // Initialize intermediate vector
          const Epetra_Map &first_op_init_map = first_op.OperatorDomainMap();
          i->reinit(IndexSet(first_op_init_map),
                    first_op.get_mpi_communicator(),
                    /*bool omit_zeroing_entries =*/true);

          // Duplicated from TrilinosWrappers::SparseMatrix::vmult
          const size_type i_local_size = i->end() - i->begin();
          AssertDimension(i_local_size, first_op_init_map.NumMyPoints());
          const Epetra_Map &second_op_init_map = second_op.OperatorDomainMap();
          AssertDimension(i_local_size, second_op_init_map.NumMyPoints());
          (void)second_op_init_map;
          Intermediate tril_int(View,
                                first_op_init_map,
                                const_cast<TrilinosScalar *>(i->begin()),
                                i_local_size,
                                1);

          // These operators may themselves be transposed or not, so we let them
          // decide what the intended outcome is
          second_op.ApplyInverse(tril_src, tril_int);
          first_op.ApplyInverse(tril_src, tril_dst);
          const int ierr = tril_dst.Update(1.0, tril_int, 1.0);
          AssertThrow(ierr == 0, ExcTrilinosError(ierr));

          // Reset transpose flag
          const_cast<TrilinosPayload &>(first_op).transpose();
          const_cast<TrilinosPayload &>(second_op).transpose();
        };

        return return_op;
      }



      TrilinosPayload operator*(const TrilinosPayload &first_op,
                                const TrilinosPayload &second_op)
      {
        using Domain        = typename TrilinosPayload::Domain;
        using Range         = typename TrilinosPayload::Range;
        using Intermediate  = typename TrilinosPayload::VectorType;
        using GVMVectorType = TrilinosWrappers::MPI::Vector;

        AssertThrow(first_op.locally_owned_domain_indices() ==
                      second_op.locally_owned_range_indices(),
                    ExcMessage(
                      "Operators are set to work on incompatible IndexSets."));

        TrilinosPayload return_op(first_op, second_op);

        // Capture by copy so the payloads are always valid
        return_op.vmult = [first_op, second_op](Range &       tril_dst,
                                                const Domain &tril_src) {
          // Duplicated from LinearOperator::operator*
          // TODO: Template the constructor on GrowingVectorMemory vector type?
          GrowingVectorMemory<GVMVectorType>   vector_memory;
          VectorMemory<GVMVectorType>::Pointer i(vector_memory);

          // Initialize intermediate vector
          const Epetra_Map &first_op_init_map = first_op.OperatorDomainMap();
          i->reinit(IndexSet(first_op_init_map),
                    first_op.get_mpi_communicator(),
                    /*bool omit_zeroing_entries =*/true);

          // Duplicated from TrilinosWrappers::SparseMatrix::vmult
          const size_type i_local_size = i->end() - i->begin();
          AssertDimension(i_local_size, first_op_init_map.NumMyPoints());
          const Epetra_Map &second_op_init_map = second_op.OperatorRangeMap();
          AssertDimension(i_local_size, second_op_init_map.NumMyPoints());
          (void)second_op_init_map;
          Intermediate tril_int(View,
                                first_op_init_map,
                                const_cast<TrilinosScalar *>(i->begin()),
                                i_local_size,
                                1);

          // These operators may themselves be transposed or not, so we let them
          // decide what the intended outcome is
          second_op.Apply(tril_src, tril_int);
          first_op.Apply(tril_int, tril_dst);
        };

        return_op.Tvmult = [first_op, second_op](Domain &     tril_dst,
                                                 const Range &tril_src) {
          // Duplicated from LinearOperator::operator*
          // TODO: Template the constructor on GrowingVectorMemory vector type?
          GrowingVectorMemory<GVMVectorType>   vector_memory;
          VectorMemory<GVMVectorType>::Pointer i(vector_memory);

          // These operators may themselves be transposed or not, so we let them
          // decide what the intended outcome is
          // We must first transpose the operators to get the right IndexSets
          // for the input, intermediate and result vectors
          const_cast<TrilinosPayload &>(first_op).transpose();
          const_cast<TrilinosPayload &>(second_op).transpose();

          // Initialize intermediate vector
          const Epetra_Map &first_op_init_map = first_op.OperatorRangeMap();
          i->reinit(IndexSet(first_op_init_map),
                    first_op.get_mpi_communicator(),
                    /*bool omit_zeroing_entries =*/true);

          // Duplicated from TrilinosWrappers::SparseMatrix::vmult
          const size_type i_local_size = i->end() - i->begin();
          AssertDimension(i_local_size, first_op_init_map.NumMyPoints());
          const Epetra_Map &second_op_init_map = second_op.OperatorDomainMap();
          AssertDimension(i_local_size, second_op_init_map.NumMyPoints());
          (void)second_op_init_map;
          Intermediate tril_int(View,
                                first_op_init_map,
                                const_cast<TrilinosScalar *>(i->begin()),
                                i_local_size,
                                1);

          // Apply the operators in the reverse order to vmult
          first_op.Apply(tril_src, tril_int);
          second_op.Apply(tril_int, tril_dst);

          // Reset transpose flag
          const_cast<TrilinosPayload &>(first_op).transpose();
          const_cast<TrilinosPayload &>(second_op).transpose();
        };

        return_op.inv_vmult = [first_op, second_op](Domain &     tril_dst,
                                                    const Range &tril_src) {
          // Duplicated from LinearOperator::operator*
          // TODO: Template the constructor on GrowingVectorMemory vector type?
          GrowingVectorMemory<GVMVectorType>   vector_memory;
          VectorMemory<GVMVectorType>::Pointer i(vector_memory);

          // Initialize intermediate vector
          const Epetra_Map &first_op_init_map = first_op.OperatorRangeMap();
          i->reinit(IndexSet(first_op_init_map),
                    first_op.get_mpi_communicator(),
                    /*bool omit_zeroing_entries =*/true);

          // Duplicated from TrilinosWrappers::SparseMatrix::vmult
          const size_type i_local_size = i->end() - i->begin();
          AssertDimension(i_local_size, first_op_init_map.NumMyPoints());
          const Epetra_Map &second_op_init_map = second_op.OperatorDomainMap();
          AssertDimension(i_local_size, second_op_init_map.NumMyPoints());
          (void)second_op_init_map;
          Intermediate tril_int(View,
                                first_op_init_map,
                                const_cast<TrilinosScalar *>(i->begin()),
                                i_local_size,
                                1);

          // Apply the operators in the reverse order to vmult
          // and the same order as Tvmult
          first_op.ApplyInverse(tril_src, tril_int);
          second_op.ApplyInverse(tril_int, tril_dst);
        };

        return_op.inv_Tvmult = [first_op, second_op](Range &       tril_dst,
                                                     const Domain &tril_src) {
          // Duplicated from LinearOperator::operator*
          // TODO: Template the constructor on GrowingVectorMemory vector type?
          GrowingVectorMemory<GVMVectorType>   vector_memory;
          VectorMemory<GVMVectorType>::Pointer i(vector_memory);

          // These operators may themselves be transposed or not, so we let them
          // decide what the intended outcome is
          // We must first transpose the operators to get the right IndexSets
          // for the input, intermediate and result vectors
          const_cast<TrilinosPayload &>(first_op).transpose();
          const_cast<TrilinosPayload &>(second_op).transpose();

          // Initialize intermediate vector
          const Epetra_Map &first_op_init_map = first_op.OperatorDomainMap();
          i->reinit(IndexSet(first_op_init_map),
                    first_op.get_mpi_communicator(),
                    /*bool omit_zeroing_entries =*/true);

          // Duplicated from TrilinosWrappers::SparseMatrix::vmult
          const size_type i_local_size = i->end() - i->begin();
          AssertDimension(i_local_size, first_op_init_map.NumMyPoints());
          const Epetra_Map &second_op_init_map = second_op.OperatorRangeMap();
          AssertDimension(i_local_size, second_op_init_map.NumMyPoints());
          (void)second_op_init_map;
          Intermediate tril_int(View,
                                first_op_init_map,
                                const_cast<TrilinosScalar *>(i->begin()),
                                i_local_size,
                                1);

          // These operators may themselves be transposed or not, so we let them
          // decide what the intended outcome is
          // Apply the operators in the reverse order to Tvmult
          // and the same order as vmult
          second_op.ApplyInverse(tril_src, tril_int);
          first_op.ApplyInverse(tril_int, tril_dst);

          // Reset transpose flag
          const_cast<TrilinosPayload &>(first_op).transpose();
          const_cast<TrilinosPayload &>(second_op).transpose();
        };

        return return_op;
      }

    } // namespace LinearOperatorImplementation
  }   /* namespace internal */
} /* namespace TrilinosWrappers */



// explicit instantiations
#  include "trilinos_sparse_matrix.inst"

#  ifndef DOXYGEN
// TODO: put these instantiations into generic file
namespace TrilinosWrappers
{
  template void
  SparseMatrix::reinit(const dealii::SparsityPattern &);

  template void
  SparseMatrix::reinit(const DynamicSparsityPattern &);

  template void
  SparseMatrix::reinit(const IndexSet &,
                       const IndexSet &,
                       const dealii::SparsityPattern &,
                       const MPI_Comm &,
                       const bool);

  template void
  SparseMatrix::reinit(const IndexSet &,
                       const IndexSet &,
                       const DynamicSparsityPattern &,
                       const MPI_Comm &,
                       const bool);

  template void
  SparseMatrix::vmult(MPI::Vector &, const MPI::Vector &) const;

  template void
  SparseMatrix::vmult(dealii::Vector<double> &,
                      const dealii::Vector<double> &) const;

  template void
  SparseMatrix::vmult(
    dealii::LinearAlgebra::distributed::Vector<double> &,
    const dealii::LinearAlgebra::distributed::Vector<double> &) const;

#    ifdef DEAL_II_WITH_MPI
#      ifdef DEAL_II_TRILINOS_WITH_TPETRA
  template void
  SparseMatrix::vmult(
    dealii::LinearAlgebra::TpetraWrappers::Vector<double> &,
    const dealii::LinearAlgebra::TpetraWrappers::Vector<double> &) const;

  template void
  SparseMatrix::vmult(
    dealii::LinearAlgebra::TpetraWrappers::Vector<float> &,
    const dealii::LinearAlgebra::TpetraWrappers::Vector<float> &) const;
#      endif

  template void
  SparseMatrix::vmult(
    dealii::LinearAlgebra::EpetraWrappers::Vector &,
    const dealii::LinearAlgebra::EpetraWrappers::Vector &) const;
#    endif

  template void
  SparseMatrix::Tvmult(MPI::Vector &, const MPI::Vector &) const;

  template void
  SparseMatrix::Tvmult(dealii::Vector<double> &,
                       const dealii::Vector<double> &) const;

  template void
  SparseMatrix::Tvmult(
    dealii::LinearAlgebra::distributed::Vector<double> &,
    const dealii::LinearAlgebra::distributed::Vector<double> &) const;

#    ifdef DEAL_II_WITH_MPI
#      ifdef DEAL_II_TRILINOS_WITH_TPETRA
  template void
  SparseMatrix::Tvmult(
    dealii::LinearAlgebra::TpetraWrappers::Vector<double> &,
    const dealii::LinearAlgebra::TpetraWrappers::Vector<double> &) const;

  template void
  SparseMatrix::Tvmult(
    dealii::LinearAlgebra::TpetraWrappers::Vector<float> &,
    const dealii::LinearAlgebra::TpetraWrappers::Vector<float> &) const;
#      endif

  template void
  SparseMatrix::Tvmult(
    dealii::LinearAlgebra::EpetraWrappers::Vector &,
    const dealii::LinearAlgebra::EpetraWrappers::Vector &) const;
#    endif

  template void
  SparseMatrix::vmult_add(MPI::Vector &, const MPI::Vector &) const;

  template void
  SparseMatrix::vmult_add(dealii::Vector<double> &,
                          const dealii::Vector<double> &) const;

  template void
  SparseMatrix::vmult_add(
    dealii::LinearAlgebra::distributed::Vector<double> &,
    const dealii::LinearAlgebra::distributed::Vector<double> &) const;

#    ifdef DEAL_II_WITH_MPI
#      ifdef DEAL_II_TRILINOS_WITH_TPETRA
  template void
  SparseMatrix::vmult_add(
    dealii::LinearAlgebra::TpetraWrappers::Vector<double> &,
    const dealii::LinearAlgebra::TpetraWrappers::Vector<double> &) const;

  template void
  SparseMatrix::vmult_add(
    dealii::LinearAlgebra::TpetraWrappers::Vector<float> &,
    const dealii::LinearAlgebra::TpetraWrappers::Vector<float> &) const;
#      endif

  template void
  SparseMatrix::vmult_add(
    dealii::LinearAlgebra::EpetraWrappers::Vector &,
    const dealii::LinearAlgebra::EpetraWrappers::Vector &) const;
#    endif

  template void
  SparseMatrix::Tvmult_add(MPI::Vector &, const MPI::Vector &) const;

  template void
  SparseMatrix::Tvmult_add(dealii::Vector<double> &,
                           const dealii::Vector<double> &) const;

  template void
  SparseMatrix::Tvmult_add(
    dealii::LinearAlgebra::distributed::Vector<double> &,
    const dealii::LinearAlgebra::distributed::Vector<double> &) const;

#    ifdef DEAL_II_WITH_MPI
#      ifdef DEAL_II_TRILINOS_WITH_TPETRA
  template void
  SparseMatrix::Tvmult_add(
    dealii::LinearAlgebra::TpetraWrappers::Vector<double> &,
    const dealii::LinearAlgebra::TpetraWrappers::Vector<double> &) const;

  template void
  SparseMatrix::Tvmult_add(
    dealii::LinearAlgebra::TpetraWrappers::Vector<float> &,
    const dealii::LinearAlgebra::TpetraWrappers::Vector<float> &) const;
#      endif

  template void
  SparseMatrix::Tvmult_add(
    dealii::LinearAlgebra::EpetraWrappers::Vector &,
    const dealii::LinearAlgebra::EpetraWrappers::Vector &) const;
#    endif
} // namespace TrilinosWrappers
#  endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_TRILINOS
