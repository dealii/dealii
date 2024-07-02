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


#include <deal.II/base/logstream.h>

#include <deal.II/lac/ginkgo_sparse_matrix.h>

#ifdef DEAL_II_WITH_GINKGO

#  include <deal.II/lac/exceptions.h>

#  include <ginkgo/core/base/mtx_io.hpp>
#  include <ginkgo/core/base/types.hpp>

#  include <cmath>
#  include <memory>


DEAL_II_NAMESPACE_OPEN

namespace GinkgoWrappers
{
  namespace detail
  {
    // helper type for the visitor #4
    template <class... Ts>
    struct overloaded : Ts...
    {
      using Ts::operator()...;
    };
    // explicit deduction guide (not needed as of C++20)
    template <class... Ts>
    overloaded(Ts...) -> overloaded<Ts...>;

  } // namespace detail

  template <typename Number, typename IndexType>
  AbstractMatrix<Number, IndexType>::AbstractMatrix(
    std::shared_ptr<const gko::Executor> exec,
    const AbstractMatrix::size_type      m,
    const AbstractMatrix::size_type      n)
    : exec_(std::move(exec))
    , data_{AbstractMatrix::assembly{gko::dim<2>{m, n}}}
  {}

  template <typename Number, typename IndexType>
  AbstractMatrix<Number, IndexType>::AbstractMatrix(
    std::unique_ptr<ginkgo_type> M)
    : exec_(M->get_executor())
    , data_(std::move(M))
  {}

  template <typename Number, typename IndexType>
  AbstractMatrix<Number, IndexType> &
  AbstractMatrix<Number, IndexType>::operator=(AbstractMatrix &&m) noexcept(
    false)
  {
    if (this != &m)
      {
        std::visit(detail::overloaded{
                     [this](assembly &data) {
                       this->data_ =
                         std::exchange(data, assembly{data.get_size()});
                     },
                     [this](matrix_ptr data) {
                       this->data_ = data->create_default();
                       std_cxx17::get<matrix_ptr>(this->data_)->move_from(data);
                     }},
                   m.data_);
      }
    return *this;
  }

  template <typename Number, typename IndexType>
  AbstractMatrix<Number, IndexType>::AbstractMatrix(AbstractMatrix &&m) noexcept
    : exec_(m.exec_)
    , data_(matrix_ptr{})
  {
    *this = std::move(m);
  }

  template <typename Number, typename IndexType>
  const typename AbstractMatrix<Number, IndexType>::ginkgo_type *
  AbstractMatrix<Number, IndexType>::get_gko_object() const
  {
    return std_cxx17::get<matrix_ptr>(this->data_).get();
  }

  template <typename Number, typename IndexType>
  std::shared_ptr<const typename AbstractMatrix<Number, IndexType>::ginkgo_type>
  AbstractMatrix<Number, IndexType>::get_shared_gko_object() const
  {
    return std_cxx17::get<matrix_ptr>(this->data_);
  }

  template <typename Number, typename IndexType>
  void
  AbstractMatrix<Number, IndexType>::set(const AbstractMatrix::size_type i,
                                         const AbstractMatrix::size_type j,
                                         const Number                    value)
  {
    std_cxx17::get<assembly>(data_).set_value(static_cast<IndexType>(i),
                                              static_cast<IndexType>(j),
                                              value);
  }

  template <typename Number, typename IndexType>
  void
  AbstractMatrix<Number, IndexType>::set(
    const AbstractMatrix::size_type  row,
    const AbstractMatrix::size_type  n_cols,
    const AbstractMatrix::size_type *col_indices,
    const Number *                   values,
    const bool                       elide_zero_values)
  {
    for (size_type k = 0; k < n_cols; ++k)
      {
        if (values[k] != Number{0} || !elide_zero_values)
          {
            set(row, col_indices[k], values[k]);
          }
      }
  }

  template <typename Number, typename IndexType>
  void
  AbstractMatrix<Number, IndexType>::set(const std::vector<size_type> &indices,
                                         const FullMatrix<Number> &full_matrix,
                                         const bool elide_zero_values)
  {
    set(indices, indices, full_matrix, elide_zero_values);
  }

  template <typename Number, typename IndexType>
  void
  AbstractMatrix<Number, IndexType>::set(
    const std::vector<size_type> &row_indices,
    const std::vector<size_type> &col_indices,
    const FullMatrix<Number> &    full_matrix,
    const bool                    elide_zero_values)
  {
    for (size_type row = 0; row < row_indices.size(); ++row)
      {
        set(row_indices[row],
            col_indices.size(),
            col_indices.data(),
            full_matrix[row].begin(),
            elide_zero_values);
      }
  }

  template <typename Number, typename IndexType>
  void
  AbstractMatrix<Number, IndexType>::set(
    const AbstractMatrix::size_type row,
    const std::vector<size_type> &  col_indices,
    const std::vector<Number> &     values,
    const bool                      elide_zero_values)
  {
    set(row,
        col_indices.size(),
        col_indices.data(),
        values.data(),
        elide_zero_values);
  }


  template <typename Number, typename IndexType>
  void
  AbstractMatrix<Number, IndexType>::add(const AbstractMatrix::size_type i,
                                         const AbstractMatrix::size_type j,
                                         const Number                    value)
  {
    std_cxx17::get<assembly>(data_).add_value(static_cast<IndexType>(i),
                                              static_cast<IndexType>(j),
                                              value);
  }

  template <typename Number, typename IndexType>
  void
  AbstractMatrix<Number, IndexType>::add(
    const AbstractMatrix::size_type  row,
    const AbstractMatrix::size_type  n_cols,
    const AbstractMatrix::size_type *col_indices,
    const Number *                   values,
    const bool                       elide_zero_values,
    const bool                       col_indices_are_sorted [[maybe_unused]])
  {
    for (size_type k = 0; k < n_cols; ++k)
      {
        if (values[k] != Number{0} || !elide_zero_values)
          {
            add(row, col_indices[k], values[k]);
          }
      }
  }

  template <typename Number, typename IndexType>
  void
  AbstractMatrix<Number, IndexType>::add(const std::vector<size_type> &indices,
                                         const FullMatrix<Number> &full_matrix,
                                         const bool elide_zero_values)
  {
    add(indices, indices, full_matrix, elide_zero_values);
  }

  template <typename Number, typename IndexType>
  void
  AbstractMatrix<Number, IndexType>::add(
    const std::vector<size_type> &row_indices,
    const std::vector<size_type> &col_indices,
    const FullMatrix<Number> &    full_matrix,
    const bool                    elide_zero_values)
  {
    for (size_type i = 0; i < row_indices.size(); ++i)
      {
        add(row_indices[i],
            col_indices.size(),
            col_indices.data(),
            full_matrix[i].begin(),
            elide_zero_values);
      }
  }

  template <typename Number, typename IndexType>
  void
  AbstractMatrix<Number, IndexType>::add(
    const AbstractMatrix::size_type row,
    const std::vector<size_type> &  col_indices,
    const std::vector<Number> &     values,
    const bool                      elide_zero_values)
  {
    add(row,
        col_indices.size(),
        col_indices.data(),
        values.data(),
        elide_zero_values);
  }


  template <typename Number,
            typename IndexType,
            template <class, class>
            class GinkgoType>
  Matrix<Number, IndexType, GinkgoType>::Matrix(
    std::shared_ptr<const gko::Executor> exec,
    const Matrix::size_type              m,
    const Matrix::size_type              n)
    : AbstractMatrix<Number, IndexType>(exec, m, n)
  {}

  template <typename Number,
            typename IndexType,
            template <class, class>
            class GinkgoType>
  Matrix<Number, IndexType, GinkgoType>::Matrix(std::unique_ptr<ginkgo_type> M)
    : AbstractMatrix<Number, IndexType>(std::move(M))
  {}

  template <typename Number,
            typename IndexType,
            template <class, class>
            class GinkgoType>
  void
  Matrix<Number, IndexType, GinkgoType>::compress(
    const VectorOperation::values /* unused */)
  {
    auto assembly_data =
      std::move(std_cxx17::get<typename Base::assembly>(this->data_));
    auto assembled_matrix = ginkgo_type::create(this->exec_);
    assembled_matrix->read(assembly_data.get_ordered_data());
    this->data_ = std::move(assembled_matrix);
  }


  template <typename Number, typename IndexType>
  Csr<Number, IndexType>::Csr(std::shared_ptr<const gko::Executor> exec,
                              const SparseMatrix<Number> &         other)
    : Csr::Base{std::move(exec), 0, 0}
  {
    const size_type M = other.m();
    const size_type N = other.n();

    auto       data_host    = ginkgo_type::create(this->exec_->get_master(),
                                         gko::dim<2>(M, N),
                                         other.n_nonzero_elements());
    Number *   mat_values   = data_host->get_values();
    IndexType *mat_row_ptrs = data_host->get_row_ptrs();
    IndexType *mat_col_idxs = data_host->get_col_idxs();

    // Copy over the data from the matrix to the data structures Ginkgo needs.
    //
    // Final note: if the matrix has entries in the sparsity pattern that are
    // actually occupied by entries that have a zero numerical value, then we
    // keep them anyway. people are supposed to provide accurate sparsity
    // patterns.

    // first fill row lengths array
    mat_row_ptrs[0] = 0;
    for (size_type row = 1; row <= M; ++row)
      mat_row_ptrs[row] = mat_row_ptrs[row - 1] + other.get_row_length(row - 1);

    // Copy over matrix elements. note that for sparse matrices,
    // iterators are sorted so that they traverse each row from start to end
    // before moving on to the next row. however, this isn't true for block
    // matrices, so we have to do a bit of bookkeeping
    {
      // Have an array that for each row points to the first entry not yet
      // written to
      std::vector<IndexType> row_pointers(M + 1);
      std::copy(data_host->get_row_ptrs(),
                data_host->get_row_ptrs() + M + 1,
                row_pointers.begin());

      // Loop over the elements of the matrix row by row, as suggested in the
      // documentation of the sparse matrix iterator class
      for (size_type row = 0; row < M; ++row)
        {
          for (typename ::dealii::SparseMatrix<Number>::const_iterator p =
                 other.begin(row);
               p != other.end(row);
               ++p)
            {
              // Write entry into the first free one for this row
              mat_col_idxs[row_pointers[row]] = p->column();
              mat_values[row_pointers[row]]   = p->value();

              // Then move pointer ahead
              ++row_pointers[row];
            }
        }

      // At the end, we should have written all rows completely
      for (size_type i = 0; i < M - 1; ++i)
        Assert(row_pointers[i] == mat_row_ptrs[i + 1], ExcInternalError());
    }

    auto data = ginkgo_type::create(this->exec_);
    data->move_from(data_host.get());
    data->sort_by_column_index();

    *this = Csr{std::move(data)};
  }

  using gko::int32;
  using gko::int64;

#  define DECLARE_ABSTRACT_MATRIX(_value_type, _index_type) \
    class AbstractMatrix<_value_type, _index_type>

  GKO_INSTANTIATE_FOR_EACH_VALUE_AND_INDEX_TYPE(DECLARE_ABSTRACT_MATRIX);


#  define DECLARE_CSR(_value_type, _index_type) \
    class Csr<_value_type, _index_type>

  GKO_INSTANTIATE_FOR_EACH_VALUE_AND_INDEX_TYPE(DECLARE_CSR);


#  define INSTANTIATE_CONCRETE_MATRIX(_type)                   \
    template class Matrix<double, int32, _type>;               \
    template class Matrix<double, int64, _type>;               \
    template class Matrix<float, int32, _type>;                \
    template class Matrix<float, int64, _type>;                \
    template class Matrix<std::complex<double>, int32, _type>; \
    template class Matrix<std::complex<double>, int64, _type>; \
    template class Matrix<std::complex<float>, int32, _type>;  \
    template class Matrix<std::complex<float>, int64, _type>

  INSTANTIATE_CONCRETE_MATRIX(gko::matrix::Csr);
  INSTANTIATE_CONCRETE_MATRIX(gko::matrix::Coo);
  INSTANTIATE_CONCRETE_MATRIX(gko::matrix::Ell);
  //    INSTANTIATE_CONCRETE_MATRIX(gko::matrix::Fbcsr); // need to fix ginkgo
  //    constructor
  INSTANTIATE_CONCRETE_MATRIX(gko::matrix::Hybrid);
  INSTANTIATE_CONCRETE_MATRIX(gko::matrix::Sellp);
} // namespace GinkgoWrappers

DEAL_II_NAMESPACE_CLOSE

#endif
