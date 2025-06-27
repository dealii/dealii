// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_linear_index_iterator_h
#define dealii_linear_index_iterator_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

#include <iterator>


DEAL_II_NAMESPACE_OPEN
/**
 * Many classes in deal.II, such as FullMatrix, TransposeTable, and
 * SparseMatrix, store their data in contiguous buffers (though the
 * <em>interpretation</em> of what the elements of these buffers represent
 * can, of course, be complex). For example, FullMatrix and TransposeTable
 * store their data in row major and column major order respectively, whereas
 * for SparseMatrix the mapping from buffer location to matrix entry
 * $\mathbf{A}(i, j)$ is more complicated. In any case, however, the
 * contiguous arrangements of elements enables random access iteration.
 *
 * LinearIndexIterator provides most of the functionality needed to write
 * iterators for these classes. LinearIndexIterator is essentially a
 * simplified version of <code>boost::iterator_facade</code> that assumes
 * <code>AccessorType</code> provides certain members (documented below) that
 * completely describe the state of the iterator. The intended use of this
 * class is for containers to define their own accessor classes and then use
 * the curiously recurring template pattern (CRTP) technique to define their
 * iterators. For example, here is a container that uses LinearIndexIterator
 * to define its own iterator classes:
 *
 * @code
 * template <typename T>
 * class Container
 * {
 * protected:
 *   // forward declaration for friendship
 *   template <bool Constness>
 *   class Iterator;
 *
 *   template <bool Constness>
 *   class Accessor
 *   {
 *   public:
 *     // const iterators store a const pointer
 *     using container_pointer_type
 *       = std::conditional_t<Constness,
 *                            const Container<T>*,
 *                            Container<T>*>;
 *
 *     // This alias is assumed to exist.
 *     using size_type = std::size_t;
 *
 *     // constructor.
 *     Accessor(const container_pointer_type container,
 *              const std::ptrdiff_t index);
 *
 *     // constructor.
 *     Accessor();
 *
 *     // get a constant reference to the current value.
 *     const T& value() const;
 *
 *   protected:
 *     container_pointer_type container;
 *     std::ptrdiff_t linear_index;
 *
 *     // LinearIndexIterator needs access to linear_index and container.
 *     friend class LinearIndexIterator<Iterator<Constness>,
 *                                      Accessor<Constness>>;
 *   };
 *
 *   template <bool Constness>
 *   class Iterator : public LinearIndexIterator<Iterator<Constness>,
 *                                               Accessor<Constness>>
 *   {
 *     // Constructor.
 *     Iterator(Container<T> * const container, const std::ptrdiff_t index);
 *
 *     // implement additional constructors here, but all state should be
 *     // contained in the Accessor, which is a member of the base class.
 *   };
 *
 * public:
 *   using size_type = std::size_t;
 *   using const_iterator = Iterator<true>;
 *   using iterator = Iterator<false>;
 *
 *   iterator begin ();
 *   iterator end ();
 *
 *   const_iterator begin () const;
 *   const_iterator end () const;
 * };
 * @endcode
 *
 * @tparam DerivedIterator As shown in the example above, concrete iterator
 * classes should use this class with the CRTP technique: this provides the
 * boiler-plate comparison and arithmetic operators for iterators. This is
 * necessary for, e.g., LinearIndexIterator::operator++() to return the
 * correct type.
 *
 * @tparam AccessorType LinearIndexIterator assumes that the
 * <code>AccessorType</code> template parameter has the following members
 * which completely describe the current state of the iterator:
 * <ol>
 *   <li>A pointer named <code>container</code> to the original container (e.g.,
 *   the relevant SparseMatrix). This should be a <code>const</code> pointer
 *   for <code>const</code> iterators.</li>
 *   <li>An array index named <code>linear_index</code> that stores the current
 *   position in the container's storage buffer. <code>linear_index</code> does
 *   not need to be an integer: it could be a class type (convertible to the
 *   correct index type of the container) that implements
 *   <code>operator+=</code>, <code>operator&lt;</code>, and
 *   <code>operator==</code>. For example, one could implement a strided
 *   iterator by implementing <code>operator+=</code> and
 *   <code>operator-</code> with multiplicative factors.</li>
 * </ol>
 * In addition, <code>AccessorType</code> should declare the relevant
 * LinearIndexIterator instantiation to be a <code>friend</code> and define a
 * <code>size_type</code> type.
 *
 * @note TransposeTable uses this template to implement its iterators.
 */
template <typename DerivedIterator, typename AccessorType>
class LinearIndexIterator
{
public:
  /**
   * Iterator category.
   */
#ifdef DEAL_II_HAVE_CXX20
  using iterator_category = std::contiguous_iterator_tag;
#else
  using iterator_category = std::random_access_iterator_tag;
#endif

  /**
   * An alias for the type you get when you dereference an iterator of the
   * current kind.
   */
  using value_type = AccessorType;

  /**
   * Difference type.
   */
  using difference_type = std::ptrdiff_t;

  /**
   * Reference type.
   */
  using reference = const value_type &;

  /**
   * Pointer type.
   */
  using pointer = const value_type *;

  /**
   * Size type used by the underlying container.
   */
  using size_type = typename value_type::size_type;

  /**
   * Copy operator.
   */
  DerivedIterator &
  operator=(const DerivedIterator &it);

  /**
   * Prefix increment.
   */
  DerivedIterator &
  operator++();

  /**
   * Postfix increment.
   */
  DerivedIterator
  operator++(int);

  /**
   * Prefix decrement.
   */
  DerivedIterator &
  operator--();

  /**
   * Postfix decrement.
   */
  DerivedIterator
  operator--(int);

  /**
   * Return an iterator that is @p n entries ahead of the current one.
   */
  DerivedIterator
  operator+(const difference_type n) const;

  /**
   * Return an iterator that is @p n entries behind the current one.
   */
  DerivedIterator
  operator-(const difference_type n) const;

  /**
   * Increment the iterator position by @p n.
   */
  DerivedIterator &
  operator+=(const difference_type n);

  /**
   * Decrement the iterator position by @p n.
   */
  DerivedIterator &
  operator-=(const difference_type n);

  /**
   * Return the distance between the current iterator and the argument. The
   * distance is given by how many times one has to apply operator++() to the
   * current iterator to get the argument (for a positive return value), or
   * operator--() (for a negative return value).
   */
  difference_type
  operator-(const DerivedIterator &p) const;

  /**
   * Dereferencing operator.
   */
  reference
  operator*() const;

  /**
   * Dereferencing operator.
   */
  pointer
  operator->() const;

  /**
   * Comparison operator. Returns <code>true</code> if both iterators point to
   * the same entry in the same container.
   */
  template <typename OtherIterator>
  std::enable_if_t<std::is_convertible_v<OtherIterator, DerivedIterator>, bool>
  operator==(const OtherIterator &right) const
  {
    const auto &right_2 = static_cast<const DerivedIterator &>(right);
    return this->accessor == right_2.accessor;
  }

  /**
   * Opposite of operator==().
   */
  template <typename OtherIterator>
  std::enable_if_t<std::is_convertible_v<OtherIterator, DerivedIterator>, bool>
  operator!=(const OtherIterator &right) const
  {
    return !(*this == right);
  }

  /**
   * Comparison operator: uses the same ordering as operator<(), but also
   * checks for equality.
   *
   * This function is only valid if both iterators point into the same
   * container.
   */
  bool
  operator<=(const DerivedIterator &) const;

  /**
   * Comparison operator: uses the same ordering as operator>(), but also
   * checks for equality.
   *
   * This function is only valid if both iterators point into the same
   * container.
   */
  bool
  operator>=(const DerivedIterator &) const;

  /**
   * Comparison operator. Result is true if either the first row number is
   * smaller or if the row numbers are equal and the first index is smaller.
   *
   * This function is only valid if both iterators point into the same
   * container.
   */
  bool
  operator<(const DerivedIterator &) const;

  /**
   * Comparison operator. Works in the same way as operator<(), just the other
   * way round.
   */
  bool
  operator>(const DerivedIterator &) const;

protected:
  /*
   * The inheriting class should have a default constructor.
   */
  LinearIndexIterator() = default; // NOLINT

  /**
   * Constructor that copies an accessor.
   */
  LinearIndexIterator(const AccessorType accessor);

protected:
  /**
   * Store an object of the accessor class.
   */
  AccessorType accessor;
};



template <typename DerivedIterator, typename AccessorType>
inline DerivedIterator &
LinearIndexIterator<DerivedIterator, AccessorType>::operator=(
  const DerivedIterator &it)
{
  accessor.container    = it.container;
  accessor.linear_index = it.linear_index;
  return static_cast<DerivedIterator &>(*this);
}



template <typename DerivedIterator, typename AccessorType>
inline DerivedIterator &
LinearIndexIterator<DerivedIterator, AccessorType>::operator++()
{
  return operator+=(1);
}



template <typename DerivedIterator, typename AccessorType>
inline DerivedIterator
LinearIndexIterator<DerivedIterator, AccessorType>::operator++(int)
{
  const DerivedIterator copy(this->accessor);
  operator+=(1);
  return copy;
}



template <typename DerivedIterator, typename AccessorType>
inline DerivedIterator &
LinearIndexIterator<DerivedIterator, AccessorType>::operator--()
{
  return operator+=(-1);
}



template <typename DerivedIterator, typename AccessorType>
inline DerivedIterator
LinearIndexIterator<DerivedIterator, AccessorType>::operator--(int)
{
  const DerivedIterator copy(this->accessor);
  operator+=(-1);
  return copy;
}



template <typename DerivedIterator, typename AccessorType>
inline DerivedIterator
LinearIndexIterator<DerivedIterator, AccessorType>::operator+(
  const difference_type n) const
{
  DerivedIterator copy(this->accessor);
  copy += n;
  return copy;
}



template <typename DerivedIterator, typename AccessorType>
inline DerivedIterator
LinearIndexIterator<DerivedIterator, AccessorType>::operator-(
  const difference_type n) const
{
  DerivedIterator copy(this->accessor);
  copy += -n;
  return copy;
}



template <typename DerivedIterator, typename AccessorType>
inline DerivedIterator &
LinearIndexIterator<DerivedIterator, AccessorType>::operator+=(
  const difference_type n)
{
  accessor.linear_index += n;
  return static_cast<DerivedIterator &>(*this);
}



template <typename DerivedIterator, typename AccessorType>
inline DerivedIterator &
LinearIndexIterator<DerivedIterator, AccessorType>::operator-=(
  const difference_type n)
{
  return operator+=(-n);
}



template <typename DerivedIterator, typename AccessorType>
inline
  typename LinearIndexIterator<DerivedIterator, AccessorType>::difference_type
  LinearIndexIterator<DerivedIterator, AccessorType>::operator-(
    const DerivedIterator &other) const
{
  Assert(this->accessor.container == other.accessor.container,
         ExcMessage(
           "Only iterators pointing to the same container can be compared."));
  return this->accessor.linear_index - other.accessor.linear_index;
}



template <typename DerivedIterator, typename AccessorType>
inline typename LinearIndexIterator<DerivedIterator, AccessorType>::reference
LinearIndexIterator<DerivedIterator, AccessorType>::operator*() const
{
  return accessor;
}



template <typename DerivedIterator, typename AccessorType>
inline typename LinearIndexIterator<DerivedIterator, AccessorType>::pointer
LinearIndexIterator<DerivedIterator, AccessorType>::operator->() const
{
  return &accessor;
}



template <typename DerivedIterator, typename AccessorType>
inline bool
LinearIndexIterator<DerivedIterator, AccessorType>::operator<=(
  const DerivedIterator &other) const
{
  return (*this == other) || (*this < other);
}



template <typename DerivedIterator, typename AccessorType>
inline bool
LinearIndexIterator<DerivedIterator, AccessorType>::operator>=(
  const DerivedIterator &other) const
{
  return !(*this < other);
}



template <typename DerivedIterator, typename AccessorType>
inline bool
LinearIndexIterator<DerivedIterator, AccessorType>::operator<(
  const DerivedIterator &other) const
{
  Assert(this->accessor.container == other.accessor.container,
         ExcMessage(
           "Only iterators pointing to the same container can be compared."));
  return this->accessor.linear_index < other.accessor.linear_index;
}



template <typename DerivedIterator, typename AccessorType>
inline bool
LinearIndexIterator<DerivedIterator, AccessorType>::operator>(
  const DerivedIterator &other) const
{
  return other < static_cast<const DerivedIterator &>(*this);
}



template <typename DerivedIterator, typename AccessorType>
inline LinearIndexIterator<DerivedIterator, AccessorType>::LinearIndexIterator(
  const AccessorType accessor)
  : accessor(accessor)
{}


DEAL_II_NAMESPACE_CLOSE

#endif
