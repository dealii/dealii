// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2013 by the deal.II authors
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

#ifndef __deal2__synchronous_iterator_h
#define __deal2__synchronous_iterator_h


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>

#include <deal.II/base/std_cxx11/tuple.h>

#include <iterator>

DEAL_II_NAMESPACE_OPEN

/**
 * A class that represents a set of
 * iterators each of which are
 * incremented by one at the same
 * time. This is typically used in calls
 * like <code>std::transform(a.begin(),
 * a.end(), b.begin(), functor);</code>
 * where we have synchronous iterators
 * marching through the containers
 * <code>a,b</code>. If an object of this
 * type represents the end of a range,
 * only the first element is considered
 * (we only have <code>a.end()</code>,
 * not <code>b.end()</code>)
 *
 * The template argument of the current
 * class shall be of type
 * <code>std_cxx11::tuple</code> with
 * arguments equal to the iterator types.
 *
 * This type, and the helper functions
 * associated with it, are used as the
 * Value concept for the blocked_range
 * type of the Threading Building Blocks.
 *
 * @author Wolfgang Bangerth, 2008
 */
template <typename Iterators>
struct SynchronousIterators
{
  /**
   * Constructor.
   */
  SynchronousIterators (const Iterators &i);

  /**
   * Copy constructor.
   */
  SynchronousIterators (const SynchronousIterators &i);

  /**
   * Storage for the iterators
   * represented by the current class.
   */
  Iterators iterators;
};



template <typename Iterators>
inline
SynchronousIterators<Iterators>::
SynchronousIterators (const Iterators &i)
  :
  iterators (i)
{}


template <typename Iterators>
inline
SynchronousIterators<Iterators>::
SynchronousIterators (const SynchronousIterators &i)
  :
  iterators (i.iterators)
{}



/**
 * Return whether the first element of
 * the first argument is less than the
 * first element of the second
 * argument. Since the objects compared
 * march forward all elements at the same
 * time, comparing the first element is
 * sufficient.
 */
template <typename Iterators>
inline
bool
operator< (const SynchronousIterators<Iterators> &a,
           const SynchronousIterators<Iterators> &b)
{
  return std_cxx11::get<0>(a.iterators) < std_cxx11::get<0>(b.iterators);
}



/**
 * Return the distance between the first
 * and the second argument. Since the
 * objects compared march forward all
 * elements at the same time,
 * differencing the first element is
 * sufficient.
 */
template <typename Iterators>
inline
std::size_t
operator- (const SynchronousIterators<Iterators> &a,
           const SynchronousIterators<Iterators> &b)
{
  Assert (std::distance (std_cxx11::get<0>(b.iterators),
                         std_cxx11::get<0>(a.iterators)) >= 0,
          ExcInternalError());
  return std::distance (std_cxx11::get<0>(b.iterators),
                        std_cxx11::get<0>(a.iterators));
}


/**
 * Advance a tuple of iterators by $n$.
 */
template <typename I1, typename I2>
inline
void advance (std_cxx11::tuple<I1,I2> &t,
              const unsigned int       n)
{
  std::advance (std_cxx11::get<0>(t), n);
  std::advance (std_cxx11::get<1>(t), n);
}

/**
 * Advance a tuple of iterators by $n$.
 */
template <typename I1, typename I2, typename I3>
inline
void advance (std_cxx11::tuple<I1,I2,I3> &t,
              const unsigned int          n)
{
  std::advance (std_cxx11::get<0>(t), n);
  std::advance (std_cxx11::get<1>(t), n);
  std::advance (std_cxx11::get<2>(t), n);
}

/**
 * Advance a tuple of iterators by $n$.
 */
template <typename I1, typename I2,
          typename I3, typename I4>
inline
void advance (std_cxx11::tuple<I1,I2,I3, I4> &t,
              const unsigned int              n)
{
  std::advance (std_cxx11::get<0>(t), n);
  std::advance (std_cxx11::get<1>(t), n);
  std::advance (std_cxx11::get<2>(t), n);
  std::advance (std_cxx11::get<3>(t), n);
}



/**
 * Advance a tuple of iterators by 1.
 */
template <typename I1, typename I2>
inline
void advance_by_one (std_cxx11::tuple<I1,I2> &t)
{
  ++std_cxx11::get<0>(t);
  ++std_cxx11::get<1>(t);
}

/**
 * Advance a tuple of iterators by 1.
 */
template <typename I1, typename I2, typename I3>
inline
void advance_by_one (std_cxx11::tuple<I1,I2,I3> &t)
{
  ++std_cxx11::get<0>(t);
  ++std_cxx11::get<1>(t);
  ++std_cxx11::get<2>(t);
}

/**
 * Advance a tuple of iterators by 1.
 */
template <typename I1, typename I2,
          typename I3, typename I4>
inline
void advance_by_one (std_cxx11::tuple<I1,I2,I3,I4> &t)
{
  ++std_cxx11::get<0>(t);
  ++std_cxx11::get<1>(t);
  ++std_cxx11::get<2>(t);
  ++std_cxx11::get<3>(t);
}



/**
 * Advance the elements of this iterator
 * by $n$.
 */
template <typename Iterators>
inline
SynchronousIterators<Iterators>
operator + (const SynchronousIterators<Iterators> &a,
            const std::size_t                      n)
{
  SynchronousIterators<Iterators> x (a);
  dealii::advance (x.iterators, n);
  return x;
}

/**
 * Advance the elements of this iterator
 * by 1.
 */
template <typename Iterators>
inline
SynchronousIterators<Iterators>
operator ++ (SynchronousIterators<Iterators> &a)
{
  dealii::advance_by_one (a.iterators);
  return a;
}


/**
 * Compare synch iterators for
 * inequality. Since they march in synch,
 * comparing only the first element is
 * sufficient.
 */
template <typename Iterators>
inline
bool
operator != (const SynchronousIterators<Iterators> &a,
             const SynchronousIterators<Iterators> &b)
{
  return (std_cxx11::get<0>(a.iterators) !=
          std_cxx11::get<0>(b.iterators));
}

DEAL_II_NAMESPACE_CLOSE

#endif
