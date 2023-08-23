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

#ifndef dealii_ginkgo_vector_h
#define dealii_ginkgo_vector_h


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_GINKGO

#  include <deal.II/base/index_set.h>
#  include <deal.II/base/subscriptor.h>

#  include <deal.II/lac/read_vector.h>
#  include <deal.II/lac/read_write_vector.h>
#  include <deal.II/lac/vector.h>

#  include <ginkgo/core/base/version.hpp>
#  include <ginkgo/core/matrix/dense.hpp>

#  include <iomanip>
#  include <ios>
#  include <memory>


DEAL_II_NAMESPACE_OPEN


/**
 * A namespace in which wrapper classes for Ginkgo objects reside.
 *
 * @ingroup GinkgoWrappers
 */
namespace GinkgoWrappers
{
  /**
   * Class wrapping a gko::matrix::Dense (multi) vector. This class implements
   * most of the methods common to other vector types, so it can be used as a
   * drop-in replacement.
   * @tparam Number The numerical data type. Both real and complex types are supported
   * @ingroup Vectors
   */
  template <typename Number>
  class Vector : public Subscriptor, public ReadVector<Number>
  {
  public:
    using value_type       = Number;
    using real_type        = typename numbers::NumberTraits<Number>::real_type;
    using iterator         = value_type *;
    using const_iterator   = const value_type *;
    using size_type        = typename ReadVector<Number>::size_type;
    using ginkgo_type      = gko::matrix::Dense<Number>;
    using ginkgo_norm_type = typename ginkgo_type::absolute_type;

    // Disable default constructor, because we always need to know the executor.
    Vector() = delete;

    /**
     * Creates an empty vector on the given executor @p exec.
     */
    Vector(std::shared_ptr<const gko::Executor> exec);

    /**
     * Directly wraps an Ginkgo vector.
     */
    Vector(std::unique_ptr<ginkgo_type> v);

    /**
     * Creates a vector of size @p size on the executor @p exec.
     */
    explicit Vector(std::shared_ptr<const gko::Executor> exec,
                    const size_type                      size);


    /**
     * Creates a vector on the executor @p exec and initializes
     * it with the values of @p list
     */
    explicit Vector(std::shared_ptr<const gko::Executor> exec,
                    const std::initializer_list<Number> &list);

    /**
     * Copy constructor.
     */
    Vector(const Vector &other);

    /**
     * Move constructor.
     */
    Vector(Vector &&other) noexcept;

    /**
     * Creates a vector with the same size and executor as @p V. If
     * @p omit_zeroing_entries is false, then all elements of the vector
     * are set to zero.
     */
    static Vector
    create_with_same_size(const Vector &V,
                          const bool    omit_zeroing_entries = false);

    /**
     * Copy assignment.
     * @note The executor of this is not changed. If other has a
     *       different executor, then the data copied between the
     *       executors accordingly.
     */
    Vector &
    operator=(const Vector &other);

    /**
     * Move assignment. Copies in case of different executors.
     * @note The executor of this is not changed. If other has a
     *       different executor, then the data copied between the
     *       executors accordingly.
     */
    Vector &
    operator=(Vector &&other) noexcept(false);

    /**
     * Creates a vector that is a view on a
     * dealii::LinearAlgebra::ReadWriteVector.
     */
    static std::unique_ptr<Vector>
    create_view(std::shared_ptr<const gko::Executor>              exec,
                ::dealii::LinearAlgebra::ReadWriteVector<Number> &other);


    /**
     * Creates a vector that is a view on a
     * dealii::LinearAlgebra::ReadWriteVector.
     */
    static std::unique_ptr<Vector>
    create_view(std::shared_ptr<const gko::Executor> exec,
                ::dealii::Vector<Number> &           other);


    /**
     * Creates a constant vector that is a view on a constant
     * dealii::LinearAlgebra::ReadWriteVector.
     */
    static std::unique_ptr<const Vector>
    create_view(std::shared_ptr<const gko::Executor>                    exec,
                const ::dealii::LinearAlgebra::ReadWriteVector<Number> &other);

    /**
     * Creates a constant vector that is a view on a constant
     * dealii::LinearAlgebra::ReadWriteVector.
     */
    static std::unique_ptr<const Vector>
    create_view(std::shared_ptr<const gko::Executor> exec,
                const ::dealii::Vector<Number> &     other);

    /**
     * Returns iterator to the begin of the stored data.
     * @warning Depending on the executor, this might be a pointer to device memory.
     */
    iterator
    begin() noexcept;

    /**
     * Returns constant iterator to the begin of the stored data.
     * @warning Depending on the executor, this might be a pointer to device memory.
     */
    const_iterator
    begin() const noexcept;

    /**
     * Returns iterator to the end of the stored data.
     * @warning Depending on the executor, this might be a pointer to device memory.
     */
    iterator
    end() noexcept;

    /**
     * Returns constant iterator to the end of the stored data.
     * @warning Depending on the executor, this might be a pointer to device memory.
     */
    const_iterator
    end() const noexcept;


    /**
     * Access the element @p i of the vector.
     * @warning Depending on the executor, this might access device memory.
     */
    Number
    operator()(const size_type i) const noexcept;

    /**
     * Mutable access the element @p i of the vector.
     * @warning Depending on the executor, this might access device memory.
     */
    Number &
    operator()(const size_type i) noexcept;

    // @copydoc operator()(const size_type) const
    Number
    operator[](const size_type i) const noexcept;

    // @copydoc operator()(const size_type)
    Number &
    operator[](const size_type i) noexcept;

    /**
     * Returns the size of the vector.
     */
    size_type
    size() const noexcept override;

    /**
     * Returns the index set [0, size()). Necessary for compatability
     * with other vector types.
     */
    IndexSet
    locally_owned_elements() const;

    std::unique_ptr<const ginkgo_type>
    get_gko_object() const noexcept;

    std::unique_ptr<ginkgo_type>
    get_gko_object() noexcept;

    Vector &
    operator=(const Number s);

    /**
     * Multiplies all elements with the value @p factor.
     */
    Vector &
    operator*=(const Number factor);

    /**
     * Divides all elements with the value @p factor.
     */
    Vector &
    operator/=(const Number factor);

    /**
     * Adds the vector @p V to this.
     */
    Vector &
    operator+=(const Vector &V);

    /**
     * Subtracts the vector @p V from this.
     */
    Vector &
    operator-=(const Vector &V);

    /**
     * Returns the scalar product between this and @p V. For complex types, @p V is conjugated.
     */
    Number
    operator*(const Vector &V) const;

    /**
     * Adds @p a to all elements of this.
     * @note ince there is no native Ginkgo support for this, this uses a fill and an add operation
     */
    void
    add(const Number a);

    /**
     * Computes the axpy operation `*this = *this + a * V`.
     */
    void
    add(const Number a, const Vector &V);

    /**
     * Computes the operation `*this = *this + a * V + b * W`.
     * @note Since there is no native Ginkgo support for this, this uses a copy and two add operations
     */
    void
    add(const Number a, const Vector &V, const Number b, const Vector &W);

    /**
     * @copydoc add(coconst size_type , const size_type *, const Number *)
     */
    void
    add(const std::vector<size_type> &indices,
        const dealii::Vector<Number> &values);

    /**
     * @copydoc add(coconst size_type , const size_type *, const Number *)
     */
    void
    add(std::vector<unsigned int> &indices, std::vector<Number> &values);

    /**
     * Adds the values from @p values to the elements at indices @p indices, i. e.
     * `*this[indices[i]] += values[i]`.
     * @note Since there is no native Ginkgo support for this, this might involve copies
     *       between host and device memory.
     */
    void
    add(const size_type  n_elements,
        const size_type *indices,
        const Number *   values);


    /**
     * Computes the operation `*this = s * (*this) + a * V`.
     * @note Since there is no native Ginkgo support for this, this uses a scale and add operation.
     */
    void
    sadd(const Number s, const Number a, const Vector &V);

    /**
     * Scales the elements componentwise with the values in @p scaling_factors.
     */
    void
    scale(const Vector &scaling_factors);

    /**
     * Computes the operation `*this = a * V`
     * @note Since there is no native Ginkgo support for this, this uses a copy and a scale operation.
     */
    void
    equ(const Number a, const Vector &V);

    /**
     * Checks if all values of the vector are zero.
     * @note Since there is no native Ginkgo support for this, this computes the L1 norm and compares against 100 * minimal_number
     */
    bool
    all_zero() const;

    /**
     * @warning Not implemented
     */
    Number
    mean_value() const;

    /**
     * Returns the l1 norm of this vector.
     */
    real_type
    l1_norm() const;

    /**
     * Returns the l2 norm of this vector.
     * @return
     */
    real_type
    l2_norm() const;

    /**
     * @warning Not implemented
     */
    real_type
    linfty_norm() const;

    /**
     * Returns the result of scalar product `<*this + a * V, W>`
     * @note Since there is no native Ginkgo support for this, this uses an add and a scalar product operation
     */
    Number
    add_and_dot(const Number a, const Vector &V, const Vector &W);

    /**
     * Prints the vector to the stream @p out. The parameter @p accross is ignored, it
     * will always print the vector using the matrix market format.
     */
    void
    print(std::ostream & out,
          const unsigned precision  = 3,
          const bool     scientific = true,
          const bool     accross    = true);

    /**
     * The memory consumption in bytes.
     */
    std::size_t
    memory_consumption() const;

    /**
     * No-op, only necessary for compatibility.
     */
    void
    compress(const VectorOperation::values = VectorOperation::insert) const
    {}

    /**
     * No-op, only necessary for compatibility.
     */
    bool
    has_ghost_elements() const
    {
      return false;
    }

    void
    extract_subvector_to(const ArrayView<const size_type> &indices,
                         ArrayView<Number> &values) const override;

    /**
     * Checks if the index @p index is within the local range ([0, size())) of this vector.
     */
    bool
    in_local_range(const size_type index) const;

  private:
    std::unique_ptr<ginkgo_type> data_;
  };

  template <typename Number>
  typename Vector<Number>::iterator
  Vector<Number>::begin() noexcept
  {
    return data_->get_values();
  }

  template <typename Number>
  typename Vector<Number>::const_iterator
  Vector<Number>::begin() const noexcept
  {
    return data_->get_const_values();
  }


  template <typename Number>
  typename Vector<Number>::iterator
  Vector<Number>::end() noexcept
  {
    return data_->get_values() + size();
  }

  template <typename Number>
  typename Vector<Number>::const_iterator
  Vector<Number>::end() const noexcept
  {
    return data_->get_const_values() + size();
  }

  template <typename Number>
  Number
  Vector<Number>::operator()(const size_type i) const noexcept
  {
    return data_->at(i);
  }

  template <typename Number>
  Number &
  Vector<Number>::operator()(const size_type i) noexcept
  {
    return data_->at(i);
  }

  template <typename Number>
  Number
  Vector<Number>::operator[](const size_type i) const noexcept
  {
    return data_->at(i);
  }

  template <typename Number>
  Number &
  Vector<Number>::operator[](const size_type i) noexcept
  {
    return data_->at(i);
  }

  template <typename Number>
  typename Vector<Number>::size_type
  Vector<Number>::size() const noexcept
  {
    return data_->get_size()[0];
  }
} // namespace GinkgoWrappers

/**
 * Declare dealii::GinkgoWrappers::Vector< Number > as serial vector.
 *
 * @relatesalso Vector
 */
template <typename Number>
struct is_serial_vector<GinkgoWrappers::Vector<Number>> : std::true_type
{};

DEAL_II_NAMESPACE_CLOSE

#endif

#endif // dealii_ginkgo_vector_h
