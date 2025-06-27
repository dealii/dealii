// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2009 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_synchronous_iterator_h
#define dealii_synchronous_iterator_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

#include <iterator>
#include <tuple>

DEAL_II_NAMESPACE_OPEN

/**
 * A class that represents a set of iterators each of which are incremented by
 * one at the same time. This is typically used in calls like
 * <code>std::transform(a.begin(), a.end(), b.begin(), functor);</code> where
 * we have synchronous iterators marching through the containers
 * <code>a,b</code>. If an object of this type represents the end of a range,
 * only the first element is considered (we only have <code>a.end()</code>,
 * not <code>b.end()</code>). An example of how this class is used is given
 * in step-35.
 *
 * The template argument of the current class shall be of type
 * <code>std::tuple</code> with arguments equal to the iterator types.
 *
 * The individual iterators can be accessed using
 * <code>std::get<X>(synchronous_iterator.iterators)</code> where X is
 * the number corresponding to the desired iterator.
 *
 * This type, and the helper functions associated with it, are used as the
 * Value concept for the blocked_range type of the Threading Building Blocks.
 */
template <typename Iterators>
struct SynchronousIterators
{
  /**
   * Constructor.
   */
  SynchronousIterators(const Iterators &i);

  /**
   * Dereference const operator. Returns a const reference to the iterators
   * represented by the current class.
   */
  const Iterators &
  operator*() const;

  /**
   * Dereference operator. Returns a reference to the iterators
   * represented by the current class.
   */
  Iterators &
  operator*();

  /**
   * The following traits are necessary for std::distance() to
   * work with SynchronousIterators.
   */
  using difference_type   = std::size_t;
  using iterator_category = std::bidirectional_iterator_tag;
  using value_type        = Iterators;
  using pointer           = Iterators *;
  using reference         = Iterators &;

private:
  /**
   * Storage for the iterators represented by the current class.
   */
  Iterators iterators;
};



template <typename Iterators>
inline SynchronousIterators<Iterators>::SynchronousIterators(const Iterators &i)
  : iterators(i)
{}



template <typename Iterators>
inline const Iterators &
SynchronousIterators<Iterators>::operator*() const
{
  return iterators;
}



template <typename Iterators>
inline Iterators &
SynchronousIterators<Iterators>::operator*()
{
  return iterators;
}



/**
 * Return whether the first element of the first argument is less than the
 * first element of the second argument. Since the objects compared march
 * forward all elements at the same time, comparing the first element is
 * sufficient.
 *
 * @relatesalso SynchronousIterators
 */
template <typename Iterators>
inline bool
operator<(const SynchronousIterators<Iterators> &a,
          const SynchronousIterators<Iterators> &b)
{
  return std::get<0>(*a) < std::get<0>(*b);
}



/**
 * Return the distance between the first and the second argument. Since the
 * objects compared march forward all elements at the same time, differencing
 * the first element is sufficient.
 *
 * @relatesalso SynchronousIterators
 */
template <typename Iterators>
inline std::size_t
operator-(const SynchronousIterators<Iterators> &a,
          const SynchronousIterators<Iterators> &b)
{
  Assert(std::distance(std::get<0>(*b), std::get<0>(*a)) >= 0,
         ExcInternalError());
  return std::distance(std::get<0>(*b), std::get<0>(*a));
}


/**
 * Advance a tuple of iterators by $n$.
 *
 * @relatesalso SynchronousIterators
 */
template <typename I1, typename I2>
inline void
advance(std::tuple<I1, I2> &t, const unsigned int n)
{
  std::advance(std::get<0>(t), n);
  std::advance(std::get<1>(t), n);
}

/**
 * Advance a tuple of iterators by $n$.
 *
 * @relatesalso SynchronousIterators
 */
template <typename I1, typename I2, typename I3>
inline void
advance(std::tuple<I1, I2, I3> &t, const unsigned int n)
{
  std::advance(std::get<0>(t), n);
  std::advance(std::get<1>(t), n);
  std::advance(std::get<2>(t), n);
}

/**
 * Advance a tuple of iterators by $n$.
 *
 * @relatesalso SynchronousIterators
 */
template <typename I1, typename I2, typename I3, typename I4>
inline void
advance(std::tuple<I1, I2, I3, I4> &t, const unsigned int n)
{
  std::advance(std::get<0>(t), n);
  std::advance(std::get<1>(t), n);
  std::advance(std::get<2>(t), n);
  std::advance(std::get<3>(t), n);
}

/**
 * Reverse a tuple of iterators by 1.
 *
 * @relatesalso SynchronousIterators
 */
template <typename I1, typename I2>
inline void
prev(std::tuple<I1, I2> &t)
{
  --std::get<0>(t);
  --std::get<1>(t);
}

/**
 * Reverse a tuple of iterators by 1.
 *
 * @relatesalso SynchronousIterators
 */
template <typename I1, typename I2, typename I3>
inline void
prev(std::tuple<I1, I2, I3> &t)
{
  --std::get<0>(t);
  --std::get<1>(t);
  --std::get<2>(t);
}

/**
 * Reverse a tuple of iterators by 1.
 *
 * @relatesalso SynchronousIterators
 */
template <typename I1, typename I2, typename I3, typename I4>
inline void
prev(std::tuple<I1, I2, I3, I4> &t)
{
  --std::get<0>(t);
  --std::get<1>(t);
  --std::get<2>(t);
  --std::get<3>(t);
}

/**
 * Advance a tuple of iterators by 1.
 *
 * @relatesalso SynchronousIterators
 */
template <typename I1, typename I2>
inline void
advance_by_one(std::tuple<I1, I2> &t)
{
  ++std::get<0>(t);
  ++std::get<1>(t);
}

/**
 * Advance a tuple of iterators by 1.
 *
 * @relatesalso SynchronousIterators
 */
template <typename I1, typename I2, typename I3>
inline void
advance_by_one(std::tuple<I1, I2, I3> &t)
{
  ++std::get<0>(t);
  ++std::get<1>(t);
  ++std::get<2>(t);
}

/**
 * Advance a tuple of iterators by 1.
 *
 * @relatesalso SynchronousIterators
 */
template <typename I1, typename I2, typename I3, typename I4>
inline void
advance_by_one(std::tuple<I1, I2, I3, I4> &t)
{
  ++std::get<0>(t);
  ++std::get<1>(t);
  ++std::get<2>(t);
  ++std::get<3>(t);
}

/**
 * Advance the elements of this iterator by $n$.
 *
 * @relatesalso SynchronousIterators
 */
template <typename Iterators>
inline SynchronousIterators<Iterators>
operator+(const SynchronousIterators<Iterators> &a, const std::size_t n)
{
  SynchronousIterators<Iterators> x(a);
  dealii::advance(*x, n);
  return x;
}

/**
 * Advance the elements of this iterator by 1 (pre-increment).
 *
 * @relatesalso SynchronousIterators
 */
template <typename Iterators>
inline SynchronousIterators<Iterators> &
operator++(SynchronousIterators<Iterators> &a)
{
  dealii::advance_by_one(*a);
  return a;
}

/**
 * Reverse the elements of this iterator by 1 (pre-decrement).
 *
 * @relatesalso SynchronousIterators
 */
template <typename Iterators>
inline SynchronousIterators<Iterators> &
operator--(SynchronousIterators<Iterators> &a)
{
  dealii::prev(*a);
  return a;
}

/**
 * Advance the elements of this iterator by 1 (post-increment).
 *
 * @relatesalso SynchronousIterators
 */
template <typename Iterators>
inline SynchronousIterators<Iterators>
operator++(SynchronousIterators<Iterators> &a, int)
{
  SynchronousIterators<Iterators> b = a;
  dealii::advance_by_one(*a);
  return b;
}

/**
 * Compare synch iterators for inequality. Since they march in synch,
 * comparing only the first element is sufficient.
 *
 * @relatesalso SynchronousIterators
 */
template <typename Iterators>
inline bool
operator!=(const SynchronousIterators<Iterators> &a,
           const SynchronousIterators<Iterators> &b)
{
  return (std::get<0>(*a) != std::get<0>(*b));
}

DEAL_II_NAMESPACE_CLOSE

#endif
