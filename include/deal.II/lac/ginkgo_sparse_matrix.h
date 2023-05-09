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

#ifndef dealii_ginkgo_sparse_matrix_h
#define dealii_ginkgo_sparse_matrix_h


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_GINKGO

#  include <deal.II/base/index_set.h>
#  include <deal.II/base/std_cxx17/variant.h>
#  include <deal.II/base/subscriptor.h>

#  include <deal.II/lac/ginkgo_vector.h>
#  include <deal.II/lac/sparse_matrix.h>

#  include <ginkgo/core/matrix/coo.hpp>
#  include <ginkgo/core/matrix/csr.hpp>
#  include <ginkgo/core/matrix/ell.hpp>
#  include <ginkgo/core/matrix/fbcsr.hpp>
#  include <ginkgo/core/matrix/hybrid.hpp>
#  include <ginkgo/core/matrix/sellp.hpp>

#  include <iomanip>
#  include <ios>
#  include <memory>

#  include "full_matrix.h"

DEAL_II_NAMESPACE_OPEN

namespace GinkgoWrappers
{
  template <typename Number, typename IndexType = gko::int32>
  class AbstractMatrix : public Subscriptor
  {
  public:
    using value_type  = Number;
    using real_type   = typename numbers::NumberTraits<Number>::real_type;
    using size_type   = types::global_dof_index;
    using ginkgo_type = gko::LinOp;

    AbstractMatrix() = delete;

    AbstractMatrix(const AbstractMatrix &) = delete;
    AbstractMatrix(AbstractMatrix &&m) noexcept;

    AbstractMatrix &
    operator=(const AbstractMatrix &) = delete;
    AbstractMatrix &
    operator=(AbstractMatrix &&m) noexcept(false);

    size_type
    m() const;
    size_type
    n() const;

    template <typename OtherNumber>
    void
    vmult(Vector<OtherNumber> &u, const Vector<OtherNumber> &v) const;

    template <typename OtherNumber>
    void
    vmult_add(Vector<OtherNumber> &u, const Vector<OtherNumber> &v) const;

    /**
     * @note Since there is no native Ginkgo support for this, this creates a temporary transpose.
     */
    template <typename OtherNumber>
    void
    Tvmult(Vector<OtherNumber> &u, const Vector<OtherNumber> &v) const;

    /**
     * @note Since there is no native Ginkgo support for this, this creates a temporary transpose.
     */
    template <typename OtherNumber>
    void
    Tvmult_add(Vector<OtherNumber> &u, const Vector<OtherNumber> &v) const;

    void
    set(const size_type i, const size_type j, const Number value);

    /**
     * @note Since there is no native Ginkgo support for this, this results in repeated calls to set(const size_type, const_size_type, const Number).
     */
    void
    set(const std::vector<size_type> &indices,
        const FullMatrix<Number> &    full_matrix,
        const bool                    elide_zero_values = false);

    /**
     * @note Since there is no native Ginkgo support for this, this results in repeated calls to set(const size_type, const_size_type, const Number).
     */
    void
    set(const std::vector<size_type> &row_indices,
        const std::vector<size_type> &col_indices,
        const FullMatrix<Number> &    full_matrix,
        const bool                    elide_zero_values = false);

    /**
     * @note Since there is no native Ginkgo support for this, this results in repeated calls to set(const size_type, const_size_type, const Number).
     */
    void
    set(const size_type               row,
        const std::vector<size_type> &col_indices,
        const std::vector<Number> &   values,
        const bool                    elide_zero_values = false);

    /**
     * @note Since there is no native Ginkgo support for this, this results in repeated calls to set(const size_type, const_size_type, const Number).
     */
    void
    set(const size_type  row,
        const size_type  n_cols,
        const size_type *col_indices,
        const Number *   values,
        const bool       elide_zero_values = false);

    void
    add(const size_type i, const size_type j, const Number value);

    /**
     * @note Since there is no native Ginkgo support for this, this results in repeated calls to add(const size_type, const_size_type, const Number).
     */
    void
    add(const std::vector<size_type> &indices,
        const FullMatrix<Number> &    full_matrix,
        const bool                    elide_zero_values = true);

    /**
     * @note Since there is no native Ginkgo support for this, this results in repeated calls to add(const size_type, const_size_type, const Number).
     */
    void
    add(const std::vector<size_type> &row_indices,
        const std::vector<size_type> &col_indices,
        const FullMatrix<Number> &    full_matrix,
        const bool                    elide_zero_values = true);

    /**
     * @note Since there is no native Ginkgo support for this, this results in repeated calls to add(const size_type, const_size_type, const Number).
     */
    void
    add(const size_type               row,
        const std::vector<size_type> &col_indices,
        const std::vector<Number> &   values,
        const bool                    elide_zero_values = true);

    /**
     * @note Since there is no native Ginkgo support for this, this results in repeated calls to add(const size_type, const_size_type, const Number).
     */
    void
    add(const size_type  row,
        const size_type  n_cols,
        const size_type *col_indices,
        const Number *   values,
        const bool       elide_zero_values      = true,
        const bool       col_indices_are_sorted = false);

    virtual void
    compress(
      const VectorOperation::values operation = VectorOperation::insert) = 0;

    const ginkgo_type *
    get_gko_object() const;

    std::shared_ptr<const ginkgo_type>
    get_shared_gko_object() const;

  protected:
    AbstractMatrix(std::shared_ptr<const gko::Executor> exec,
                   const size_type                      m,
                   const size_type                      n);

    AbstractMatrix(std::unique_ptr<ginkgo_type> M);

    std::shared_ptr<const gko::Executor> exec_;

    using assembly   = gko::matrix_assembly_data<Number, IndexType>;
    using matrix_ptr = std::shared_ptr<ginkgo_type>;
    std_cxx17::variant<assembly, matrix_ptr> data_;
  };


  template <typename Number,
            typename IndexType,
            template <typename, typename>
            class GinkgoType>
  class Matrix : public AbstractMatrix<Number, IndexType>
  {
    using Base = AbstractMatrix<Number, IndexType>;

  public:
    using value_type  = typename Base::value_type;
    using real_type   = typename Base::real_type;
    using size_type   = typename Base::size_type;
    using ginkgo_type = GinkgoType<Number, IndexType>;

    Matrix(std::shared_ptr<const gko::Executor> exec,
           const size_type                      m,
           const size_type                      n);
    Matrix(std::unique_ptr<ginkgo_type> M);

    const ginkgo_type *
    get_gko_object() const
    {
      return static_cast<const ginkgo_type *>(Base::get_gko_object());
    }

    std::shared_ptr<const ginkgo_type>
    get_shared_gko_object() const
    {
      return std::dynamic_pointer_cast<const ginkgo_type>(
        Base::get_shared_gko_object());
    }

    void
    compress(const VectorOperation::values operation =
               VectorOperation::insert) override;
  };

  template <typename Number, typename IndexType = gko::int32>
  class Csr : public Matrix<Number, IndexType, gko::matrix::Csr>
  {
    using Base = Matrix<Number, IndexType, gko::matrix::Csr>;

  public:
    using value_type  = typename Base::value_type;
    using real_type   = typename Base::real_type;
    using size_type   = typename Base::size_type;
    using ginkgo_type = typename Base::ginkgo_type;

    using Base::Base;

    Csr(std::shared_ptr<const gko::Executor> exec,
        const SparseMatrix<Number> &         other);
  };

  template <typename Number, typename IndexType = gko::int32>
  using Coo = Matrix<Number, IndexType, gko::matrix::Coo>;

  template <typename Number, typename IndexType = gko::int32>
  using Ell = Matrix<Number, IndexType, gko::matrix::Ell>;

  //  template<typename Number, typename IndexType = gko::int32>
  //  using Fbcsr = Matrix<Number, IndexType, gko::matrix::Fbcsr>;

  template <typename Number, typename IndexType = gko::int32>
  using Hybrid = Matrix<Number, IndexType, gko::matrix::Hybrid>;

  template <typename Number, typename IndexType = gko::int32>
  using Sellp = Matrix<Number, IndexType, gko::matrix::Sellp>;


  template <typename Number, typename IndexType>
  typename AbstractMatrix<Number, IndexType>::size_type
  AbstractMatrix<Number, IndexType>::m() const
  {
    if (std_cxx17::holds_alternative<assembly>(data_))
      {
        return std_cxx17::get<assembly>(data_).get_size()[0];
      }
    else
      {
        return std_cxx17::get<matrix_ptr>(data_)->get_size()[0];
      }
  }

  template <typename Number, typename IndexType>
  typename AbstractMatrix<Number, IndexType>::size_type
  AbstractMatrix<Number, IndexType>::n() const
  {
    if (std_cxx17::holds_alternative<assembly>(data_))
      {
        return std_cxx17::get<assembly>(data_).get_size()[1];
      }
    else
      {
        return std_cxx17::get<matrix_ptr>(data_)->get_size()[1];
      }
  }

  template <typename Number, typename IndexType>
  template <typename OtherNumber>
  void
  AbstractMatrix<Number, IndexType>::vmult_add(
    Vector<OtherNumber> &      u,
    const Vector<OtherNumber> &v) const
  {
    auto one = gko::initialize<gko::matrix::Dense<Number>>({1.0}, exec_);
    std_cxx17::get<matrix_ptr>(data_)->apply(one.get(),
                                             v.get_gko_object(),
                                             one.get(),
                                             u.get_gko_object());
  }

  template <typename Number, typename IndexType>
  template <typename OtherNumber>
  void
  AbstractMatrix<Number, IndexType>::vmult(Vector<OtherNumber> &      u,
                                           const Vector<OtherNumber> &v) const
  {
    std_cxx17::get<matrix_ptr>(data_)->apply(v.get_gko_object(),
                                             u.get_gko_object());
  }

  template <typename Number, typename IndexType>
  template <typename OtherNumber>
  void
  AbstractMatrix<Number, IndexType>::Tvmult_add(
    Vector<OtherNumber> &      u,
    const Vector<OtherNumber> &v) const
  {
    auto one = gko::initialize<gko::matrix::Dense<Number>>({1.0}, exec_);
    gko::as<gko::Transposable>(std_cxx17::get<matrix_ptr>(data_))
      ->transpose()
      ->apply(one.get(), v.get_gko_object(), one.get(), u.get_gko_object());
  }

  template <typename Number, typename IndexType>
  template <typename OtherNumber>
  void
  AbstractMatrix<Number, IndexType>::Tvmult(Vector<OtherNumber> &      u,
                                            const Vector<OtherNumber> &v) const
  {
    gko::as<gko::Transposable>(std_cxx17::get<matrix_ptr>(data_))
      ->transpose()
      ->apply(v.get_gko_object(), u.get_gko_object());
  }

} // namespace GinkgoWrappers

DEAL_II_NAMESPACE_CLOSE

#endif

#endif // dealii_ginkgo_sparse_matrix_h
