// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2023 by the deal.II authors
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

#ifndef dealii_trilinos_tpetra_sparse_matrix_templates_h
#define dealii_trilinos_tpetra_sparse_matrix_templates_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_TRILINOS_WITH_TPETRA

#  include <deal.II/lac/dynamic_sparsity_pattern.h>
#  include <deal.II/lac/trilinos_tpetra_sparse_matrix.h>
#  include <deal.II/lac/trilinos_tpetra_sparsity_pattern.h>

DEAL_II_NAMESPACE_OPEN

namespace LinearAlgebra
{

  namespace TpetraWrappers
  {
    // reinit_matrix():
    namespace
    {
      using size_type = dealii::types::signed_global_dof_index;

      template <typename NodeType>
      using MapType =
        Tpetra::Map<int, dealii::types::signed_global_dof_index, NodeType>;

      template <typename Number, typename NodeType>
      using MatrixType =
        Tpetra::CrsMatrix<Number,
                          int,
                          dealii::types::signed_global_dof_index,
                          NodeType>;

      template <typename NodeType>
      using GraphType =
        Tpetra::CrsGraph<int, dealii::types::signed_global_dof_index, NodeType>;

      template <typename Number,
                typename NodeType,
                typename SparsityPatternType>
      void
      reinit_matrix(const IndexSet            &row_parallel_partitioning,
                    const IndexSet            &column_parallel_partitioning,
                    const SparsityPatternType &sparsity_pattern,
                    const bool                 exchange_data,
                    const MPI_Comm             communicator,
                    Teuchos::RCP<MapType<NodeType>> &column_space_map,
                    Teuchos::RCP<MatrixType<Number, NodeType>> &matrix)
      {
        // release memory before reallocation
        matrix.reset();

        // Get the Tpetra::Maps
        Teuchos::RCP<MapType<NodeType>> row_space_map =
          row_parallel_partitioning.make_tpetra_map_rcp(communicator, false);

        column_space_map =
          column_parallel_partitioning.make_tpetra_map_rcp(communicator, false);

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
            SparsityPattern trilinos_sparsity;
            trilinos_sparsity.reinit(row_parallel_partitioning,
                                     column_parallel_partitioning,
                                     sparsity_pattern,
                                     communicator,
                                     exchange_data);
            matrix = Utilities::Trilinos::internal::make_rcp<
              MatrixType<Number, NodeType>>(
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
        Teuchos::RCP<GraphType<NodeType>> graph;

#  if DEAL_II_TRILINOS_VERSION_GTE(12, 16, 0)
        graph = Utilities::Trilinos::internal::make_rcp<GraphType<NodeType>>(
          row_space_map, n_entries_per_row);
#  else
        graph = Utilities::Trilinos::internal::make_rcp<GraphType<NodeType>>(
          row_space_map, Teuchos::arcpFromArray(n_entries_per_row));
#  endif

        // This functions assumes that the sparsity pattern sits on all
        // processors (completely). The parallel version uses a Tpetra graph
        // that is already distributed.

        // now insert the indices
        std::vector<TrilinosWrappers::types::int_type> row_indices;

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
            graph->insertGlobalIndices(global_row,
                                       row_length,
                                       row_indices.data());
          }

        // Eventually, optimize the graph structure (sort indices, make memory
        // contiguous, etc.). note that the documentation of the function indeed
        // states that we first need to provide the column (domain) map and then
        // the row (range) map
        graph->fillComplete(column_space_map, row_space_map);

        // check whether we got the number of columns right.
        AssertDimension(sparsity_pattern.n_cols(), graph->getGlobalNumCols());

        // And now finally generate the matrix.
        matrix =
          Utilities::Trilinos::internal::make_rcp<MatrixType<Number, NodeType>>(
            graph);
      }
    } // namespace



    // Constructors and initialization:

    // The constructor is actually the only point where we have to check
    // whether we build a serial or a parallel Trilinos matrix.
    // Actually, it does not even matter how many threads there are, but
    // only if we use an MPI compiler or a standard compiler. So, even one
    // thread on a configuration with MPI will still get a parallel interface.
    template <typename Number, typename MemorySpace>
    SparseMatrix<Number, MemorySpace>::SparseMatrix()
      : column_space_map(Utilities::Trilinos::internal::make_rcp<MapType>(
          0,
          0,
          Utilities::Trilinos::tpetra_comm_self()))
    {
      // Prepare the graph
      Teuchos::RCP<GraphType> graph =
        Utilities::Trilinos::internal::make_rcp<GraphType>(column_space_map,
                                                           column_space_map,
                                                           0);
      graph->fillComplete();

      // Create the matrix from the graph
      matrix = Utilities::Trilinos::internal::make_rcp<MatrixType>(graph);

      compressed = false;
    }



    template <typename Number, typename MemorySpace>
    SparseMatrix<Number, MemorySpace>::SparseMatrix(
      const SparsityPattern<MemorySpace> &sparsity_pattern)
      : matrix(Utilities::Trilinos::internal::make_rcp<MatrixType>(
          sparsity_pattern.trilinos_sparsity_pattern()))
    {
      column_space_map =
        Teuchos::rcp_const_cast<MapType>(sparsity_pattern.domain_partitioner());
      compressed = false;
      compress(VectorOperation::add);
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
      reinit_matrix<Number, NodeType, SparsityPatternType>(
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
        Teuchos::rcp_const_cast<MapType>(sparsity_pattern.domain_partitioner());
      matrix = Utilities::Trilinos::internal::make_rcp<MatrixType>(
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
          parallel_partitioning.make_tpetra_map_rcp(communicator, false))
      , matrix(Utilities::Trilinos::internal::make_rcp<MatrixType>(
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
          parallel_partitioning.make_tpetra_map_rcp(communicator, false))
      , compressed(false)
    {
      Teuchos::Array<size_t> n_entries_per_row_array(n_entries_per_row.begin(),
                                                     n_entries_per_row.end());
#  if DEAL_II_TRILINOS_VERSION_GTE(12, 16, 0)
      matrix = Utilities::Trilinos::internal::make_rcp<MatrixType>(
        column_space_map, n_entries_per_row_array);
#  else
      matrix = Utilities::Trilinos::internal::make_rcp<MatrixType>(
        column_space_map, Teuchos::arcpFromArray(n_entries_per_row_array));
#  endif
    }



    template <typename Number, typename MemorySpace>
    SparseMatrix<Number, MemorySpace>::SparseMatrix(
      const IndexSet &row_parallel_partitioning,
      const IndexSet &col_parallel_partitioning,
      const MPI_Comm  communicator,
      const size_type n_max_entries_per_row)
      : column_space_map(
          col_parallel_partitioning.make_tpetra_map_rcp(communicator, false))
      , matrix(Utilities::Trilinos::internal::make_rcp<MatrixType>(
          row_parallel_partitioning.make_tpetra_map_rcp(communicator, false),
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
          col_parallel_partitioning.make_tpetra_map_rcp(communicator, false))
      , compressed(false)
    {
      Teuchos::Array<size_t> n_entries_per_row_array(n_entries_per_row.begin(),
                                                     n_entries_per_row.end());
#  if DEAL_II_TRILINOS_VERSION_GTE(12, 16, 0)
      matrix = Utilities::Trilinos::internal::make_rcp<MatrixType>(
        row_parallel_partitioning.make_tpetra_map_rcp(communicator, false),
        n_entries_per_row_array);
#  else
      matrix = Utilities::Trilinos::internal::make_rcp<MatrixType>(
        row_parallel_partitioning.make_tpetra_map_rcp(communicator, false),
        Teuchos::arcpFromArray(n_entries_per_row_array));
#  endif
    }



    template <typename Number, typename MemorySpace>
    template <typename SparsityPatternType>
    inline void
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
    void
    SparseMatrix<Number, MemorySpace>::reinit(
      const IndexSet &row_parallel_partitioning,

      const IndexSet            &col_parallel_partitioning,
      const SparsityPatternType &sparsity_pattern,
      const MPI_Comm             communicator,
      const bool                 exchange_data)
    {
      reinit_matrix<Number, NodeType, SparsityPatternType>(
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
      (void)d;
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
      matrix->scale(a);
      return *this;
    }



    template <typename Number, typename MemorySpace>
    SparseMatrix<Number, MemorySpace> &
    SparseMatrix<Number, MemorySpace>::operator/=(const Number a)
    {
      Assert(a != 0, ExcDivideByZero());

      const Number factor = 1.0 / a;
      matrix->scale(factor);
      return *this;
    }



    template <typename Number, typename MemorySpace>
    void
    SparseMatrix<Number, MemorySpace>::add(
      const size_type       row,
      const size_type       n_cols,
      const size_type      *col_indices,
      const TrilinosScalar *values,
      const bool            elide_zero_values,
      const bool /*col_indices_are_sorted*/)
    {
      AssertIndexRange(row, this->m());

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
                  AssertIndexRange(col_indices[i], n());
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
            AssertIndexRange(col_indices[i], n());
            col_indices_array[i] = col_indices[i];
            values_array[i]      = values[i];
          }

      // Sum the values into the global matrix.
      matrix->sumIntoGlobalValues(row, col_indices_array, values_array);
    }



    // Multiplications

    template <typename Number, typename MemorySpace>
    void
    SparseMatrix<Number, MemorySpace>::vmult(
      Vector<Number, MemorySpace>       &dst,
      const Vector<Number, MemorySpace> &src) const
    {
      Assert(&src != &dst, ExcSourceEqualsDestination());
      Assert(matrix->isFillComplete(), ExcMatrixNotCompressed());
      Assert(src.trilinos_vector().getMap()->isSameAs(*matrix->getDomainMap()),
             ExcColMapMissmatch());
      Assert(dst.trilinos_vector().getMap()->isSameAs(*matrix->getRangeMap()),
             ExcDomainMapMissmatch());
      matrix->apply(src.trilinos_vector(), dst.trilinos_vector());
    }



    template <typename Number, typename MemorySpace>
    void
    SparseMatrix<Number, MemorySpace>::Tvmult(
      Vector<Number, MemorySpace>       &dst,
      const Vector<Number, MemorySpace> &src) const
    {
      Assert(&src != &dst, ExcSourceEqualsDestination());
      Assert(matrix->isFillComplete(), ExcMatrixNotCompressed());
      Assert(dst.trilinos_vector().getMap()->isSameAs(*matrix->getDomainMap()),
             ExcColMapMissmatch());
      Assert(src.trilinos_vector().getMap()->isSameAs(*matrix->getRangeMap()),
             ExcDomainMapMissmatch());
      matrix->apply(src.trilinos_vector(),
                    dst.trilinos_vector(),
                    Teuchos::TRANS);
    }



    template <typename Number, typename MemorySpace>
    void
    SparseMatrix<Number, MemorySpace>::vmult_add(
      Vector<Number, MemorySpace>       &dst,
      const Vector<Number, MemorySpace> &src) const
    {
      Assert(&src != &dst, ExcSourceEqualsDestination());
      Assert(matrix->isFillComplete(), ExcMatrixNotCompressed());
      Assert(src.trilinos_vector().getMap()->isSameAs(*matrix->getDomainMap()),
             ExcColMapMissmatch());
      Assert(dst.trilinos_vector().getMap()->isSameAs(*matrix->getRangeMap()),
             ExcDomainMapMissmatch());
      matrix->apply(src.trilinos_vector(),
                    dst.trilinos_vector(),
                    Teuchos::NO_TRANS,
                    Teuchos::ScalarTraits<Number>::one(),
                    Teuchos::ScalarTraits<Number>::one());
    }



    template <typename Number, typename MemorySpace>
    void
    SparseMatrix<Number, MemorySpace>::Tvmult_add(
      Vector<Number, MemorySpace>       &dst,
      const Vector<Number, MemorySpace> &src) const
    {
      Assert(&src != &dst, ExcSourceEqualsDestination());
      Assert(matrix->isFillComplete(), ExcMatrixNotCompressed());
      Assert(dst.trilinos_vector().getMap()->isSameAs(*matrix->getDomainMap()),
             ExcColMapMissmatch());
      Assert(src.trilinos_vector().getMap()->isSameAs(*matrix->getRangeMap()),
             ExcDomainMapMissmatch());
      matrix->apply(src.trilinos_vector(),
                    dst.trilinos_vector(),
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
#  if DEAL_II_TRILINOS_VERSION_GTE(13, 2, 0)
          typename MatrixType::values_host_view_type     values;
          typename MatrixType::local_inds_host_view_type indices;
#  else
          Teuchos::ArrayView<const Number> values;
          Teuchos::ArrayView<const int>    indices;
#  endif

#  if DEAL_II_TRILINOS_VERSION_GTE(14, 0, 0)
          for (size_t i = 0; i < matrix->getLocalNumRows(); ++i)
#  else
          for (size_t i = 0; i < matrix->getNodeNumRows(); ++i)
#  endif
            {
              matrix->getLocalRowView(i, indices, values);

              for (size_type j = 0; j < indices.size(); ++j)
                out << "(" << matrix->getRowMap()->getGlobalElement(i) << ","
                    << matrix->getColMap()->getGlobalElement(indices[j]) << ") "
                    << values[j] << std::endl;
            }
        }

      AssertThrow(out.fail() == false, ExcIO());
    }



    template <typename Number, typename MemorySpace>
    void
    SparseMatrix<Number, MemorySpace>::compress(
      [[maybe_unused]] VectorOperation::values operation)
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

  } // namespace TpetraWrappers

} // namespace LinearAlgebra

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_TRILINOS_WITH_TPETRA

#endif // dealii_trilinos_tpetra_sparse_matrix_templates_h
