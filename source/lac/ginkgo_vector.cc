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

#include <deal.II/lac/ginkgo_vector.h>

#ifdef DEAL_II_WITH_GINKGO

#  include <deal.II/lac/exceptions.h>
#  include <deal.II/lac/vector.h>

#  include <ginkgo/core/base/mtx_io.hpp>
#  include <ginkgo/core/base/types.hpp>

#  include <cmath>
#  include <memory>


DEAL_II_NAMESPACE_OPEN

namespace GinkgoWrappers
{
  template <typename Number>
  Vector<Number>::Vector(std::shared_ptr<const gko::Executor> exec)
    : data_(ginkgo_type::create(std::move(exec)))
  {}

  template <typename Number>
  Vector<Number>::Vector(std::unique_ptr<ginkgo_type> v)
    : data_(std::move(v))
  {}

  template <typename Number>
  Vector<Number>::Vector(std::shared_ptr<const gko::Executor> exec,
                         const size_type                      size)
    : data_(ginkgo_type::create(std::move(exec), gko::dim<2>{size, 1}))
  {}

  template <typename Number>
  Vector<Number>::Vector(std::shared_ptr<const gko::Executor> exec,
                         const std::initializer_list<Number> &list)
    : data_(gko::initialize<ginkgo_type>(list, std::move(exec)))
  {}

  template <typename Number>
  Vector<Number>::Vector(const Vector<Number> &other)
    : Subscriptor()
    , data_(gko::clone(other.data_))
  {}

  template <typename Number>
  Vector<Number>::Vector(Vector<Number> &&other) noexcept
    : Vector(other.data_->get_executor())
  {
    data_->move_from(other.data_.get());
  }

  template <typename Number>
  Vector<Number>
  Vector<Number>::create_with_same_size(const Vector<Number> &V,
                                        const bool omit_zeroing_entries)
  {
    Vector result{V.get_gko_object()->get_executor(), V.size()};
    if (!omit_zeroing_entries)
      {
        result = Number(0.0);
      }
    return result;
  }

  template <typename Number>
  Vector<Number> &
  Vector<Number>::operator=(const Vector<Number> &other)
  {
    if (this != &other)
      {
        data_->copy_from(other.data_.get());
      }
    return *this;
  }

  template <typename Number>
  Vector<Number> &
  Vector<Number>::operator=(Vector<Number> &&other) noexcept(false)
  {
    if (this != &other)
      {
        data_->move_from(other.data_.get());
      }
    return *this;
  }

  template <typename Number>
  std::unique_ptr<const typename Vector<Number>::ginkgo_type>
  Vector<Number>::get_gko_object() const noexcept
  {
    return gko::make_const_dense_view(data_);
  }

  template <typename Number>
  std::unique_ptr<typename Vector<Number>::ginkgo_type>
  Vector<Number>::get_gko_object() noexcept
  {
    return gko::make_dense_view(data_);
  }

  template <typename Number>
  std::size_t
  Vector<Number>::memory_consumption() const
  {
    return sizeof(Vector<Number>) + sizeof(ginkgo_type) +
           sizeof(Number) * size();
  }

  template <typename Number>
  void
  Vector<Number>::print(std::ostream &     out,
                        const unsigned int precision,
                        const bool         scientific,
                        const bool         accross [[maybe_unused]])
  {
    // TODO: figure out the meaning of accross
    const auto default_precision = out.precision();

    out << std::setprecision(precision);
    if (scientific)
      out << std::scientific;

    gko::write(out, data_.get());

    out << std::setprecision(default_precision);
    if (scientific)
      out << std::defaultfloat;
  }

  template <typename Number>
  IndexSet
  Vector<Number>::locally_owned_elements() const
  {
    return IndexSet(size());
  }

  template <typename Number>
  Number
  Vector<Number>::add_and_dot(const Number a, const Vector &V, const Vector &W)
  {
    this->add(a, V);
    return *this * W;
  }

  template <typename Number>
  typename Vector<Number>::real_type
  Vector<Number>::linfty_norm() const
  {
    throw ExcNotImplemented();
  }

  template <typename Number>
  typename Vector<Number>::real_type
  Vector<Number>::l2_norm() const
  {
    auto result = ginkgo_norm_type::create(data_->get_executor()->get_master(),
                                           gko::dim<2>{1, 1});
    data_->compute_norm2(result.get());
    return result->at(0);
  }

  template <typename Number>
  typename Vector<Number>::real_type
  Vector<Number>::l1_norm() const
  {
    auto result = ginkgo_norm_type ::create(data_->get_executor()->get_master(),
                                            gko::dim<2>{1, 1});
    data_->compute_norm1(result.get());
    return result->at(0);
  }

  template <typename Number>
  Number
  Vector<Number>::mean_value() const
  {
    throw ExcNotImplemented();
  }

  template <typename Number>
  bool
  Vector<Number>::all_zero() const
  {
    auto norm = l1_norm();
    return norm <=
           1e2 * std::numeric_limits<gko::remove_complex<Number>>::min();
  }

  template <typename Number>
  Vector<Number> &
  Vector<Number>::operator=(const Number s)
  {
    data_->fill(s);
    return *this;
  }

  template <typename Number>
  void
  Vector<Number>::equ(const Number a, const Vector &V)
  {
    AssertDimension(V.size(), size());
    *this = V;
    *this *= a;
  }

  template <typename Number>
  void
  Vector<Number>::scale(const Vector &scaling_factors)
  {
    AssertDimension(scaling_factors.size(), size());
    auto exec = data_->get_executor();
    auto col_view =
      ginkgo_type::create(exec,
                          gko::dim<2>{1, size()},
                          gko::make_array_view(exec, size(), begin()),
                          size());
    col_view->scale(scaling_factors.data_.get());
  }

  template <typename Number>
  void
  Vector<Number>::sadd(const Number s, const Number a, const Vector &V)
  {
    AssertDimension(V.size(), size());
    *this *= s;
    this->add(a, V);
  }

  template <typename Number>
  void
  Vector<Number>::add(const Number  a,
                      const Vector &V,
                      const Number  b,
                      const Vector &W)
  {
    AssertDimension(V.size(), size());
    AssertDimension(W.size(), size());
    Vector tmp(V);
    tmp.add(b / a, W);
    this->add(a, tmp);
  }

  template <typename Number>
  void
  Vector<Number>::add(const Number a, const Vector &V)
  {
    Assert(V.size() == size(), ExcDimensionMismatch(V.size(), size()));
    auto a_dense = gko::initialize<ginkgo_type>({a}, data_->get_executor());
    data_->add_scaled(a_dense.get(), V.data_.get());
  }

  template <typename Number>
  void
  Vector<Number>::add(const Number a)
  {
    auto a_vec = Vector(data_->get_executor(), size());
    a_vec.data_->fill(a);
    *this += a_vec;
  }

  template <typename Number>
  Vector<Number> &
  Vector<Number>::operator*=(const Number factor)
  {
    auto factor_dense =
      gko::initialize<ginkgo_type>({factor}, data_->get_executor());
    data_->scale(factor_dense.get());
    return *this;
  }

  template <typename Number>
  Vector<Number> &
  Vector<Number>::operator/=(const Number factor)
  {
    auto factor_dense =
      gko::initialize<ginkgo_type>({factor}, data_->get_executor());
    data_->inv_scale(factor_dense.get());
    return *this;
  }

  template <typename Number>
  Vector<Number> &
  Vector<Number>::operator+=(const Vector &V)
  {
    AssertDimension(V.size(), size());
    auto one = gko::initialize<ginkgo_type>({1.0}, data_->get_executor());
    data_->add_scaled(one.get(), V.data_.get());
    return *this;
  }

  template <typename Number>
  Vector<Number> &
  Vector<Number>::operator-=(const Vector &V)
  {
    AssertDimension(V.size(), size());
    auto neg_one = gko::initialize<ginkgo_type>({-1.0}, data_->get_executor());
    data_->add_scaled(neg_one.get(), V.data_.get());
    return *this;
  }

  template <typename Number>
  Number
  Vector<Number>::operator*(const Vector &V) const
  {
    AssertDimension(V.size(), size());
    auto result = ginkgo_type ::create(data_->get_executor()->get_master(),
                                       gko::dim<2>{1, 1});
    data_->compute_conj_dot(V.data_.get(), result.get());
    return result->at(0);
  }

  template <typename GinkgoType, typename ValueType>
  std::unique_ptr<GinkgoType>
  create_view_impl(std::shared_ptr<const gko::Executor> exec,
                   types::global_dof_index              size,
                   ValueType *                          data)
  {
    auto host_exec = exec->get_master();
    return GinkgoType::create(host_exec,
                              gko::dim<2>{size, 1},
                              gko::make_array_view(host_exec, size, data),
                              1);
  }

  template <typename Number>
  std::unique_ptr<Vector<Number>>
  Vector<Number>::create_view(
    std::shared_ptr<const gko::Executor>              exec,
    ::dealii::LinearAlgebra::ReadWriteVector<Number> &other)
  {
    return std::make_unique<Vector<Number>>(create_view_impl<ginkgo_type>(
      std::move(exec), other.size(), other.begin()));
  }

  template <typename Number>
  std::unique_ptr<Vector<Number>>
  Vector<Number>::create_view(std::shared_ptr<const gko::Executor> exec,
                              ::dealii::Vector<Number> &           other)
  {
    if (!exec->memory_accessible(exec->get_master()))
      {
        throw ExcInternalError(
          "Can't create a view on CPU data, since the CPU memory is not accessible from the passed in executor.");
      }
    return std::make_unique<Vector<Number>>(create_view_impl<ginkgo_type>(
      std::move(exec), other.size(), other.begin()));
  }

  template <typename Number>
  std::unique_ptr<const Vector<Number>>
  Vector<Number>::create_view(
    std::shared_ptr<const gko::Executor>                    exec,
    const ::dealii::LinearAlgebra::ReadWriteVector<Number> &other)
  {
    if (!exec->memory_accessible(exec->get_master()))
      {
        throw ExcInternalError(
          "Can't create a view on CPU data, since the CPU memory is not accessible from the passed in executor.");
      }
    return std::make_unique<Vector<Number>>(create_view_impl<ginkgo_type>(
      std::move(exec),
      other.size(),
      const_cast<
        typename ::dealii::LinearAlgebra::ReadWriteVector<Number>::iterator>(
        other.begin())));
  }

  template <typename Number>
  std::unique_ptr<const Vector<Number>>
  Vector<Number>::create_view(std::shared_ptr<const gko::Executor> exec,
                              const ::dealii::Vector<Number> &     other)
  {
    return std::make_unique<Vector<Number>>(create_view_impl<ginkgo_type>(
      std::move(exec),
      other.size(),
      const_cast<
        typename ::dealii::LinearAlgebra::ReadWriteVector<Number>::iterator>(
        other.begin())));
  }

  template <typename Number>
  void
  Vector<Number>::add(std::vector<unsigned int> &indices,
                      std::vector<Number> &      values)
  {
    Assert(indices.size() == values.size(),
           ExcDimensionMismatch(indices.size(), values.size()));
    add(indices.size(), indices.data(), values.data());
  }

  template <typename Number>
  void
  Vector<Number>::add(const Vector::size_type  n_elements,
                      const Vector::size_type *indices,
                      const Number *           values)
  {
    auto host_data =
      gko::make_temporary_clone(data_->get_executor()->get_master(), data_);
    for (size_type i = 0; i < n_elements; ++i)
      {
        host_data->get_values()[indices[i]] = values[i];
      }
  }

  template <typename Number>
  bool
  Vector<Number>::in_local_range(const Vector::size_type index) const
  {
    return index < size();
  }
  template <typename Number>
  void
  Vector<Number>::add(const std::vector<size_type> &indices,
                      const dealii::Vector<Number> &values)
  {
    Assert(indices.size() == values.size(),
           ExcDimensionMismatch(indices.size(), values.size()));
    add(indices.size(), indices.data(), values.begin());
  }


#  define DECLARE_GINKGO_VECTOR(_type) class Vector<_type>

  GKO_INSTANTIATE_FOR_EACH_VALUE_TYPE(DECLARE_GINKGO_VECTOR);

#  undef DECLARE_GINKGO_VECTOR


} // namespace GinkgoWrappers


DEAL_II_NAMESPACE_CLOSE

#endif
