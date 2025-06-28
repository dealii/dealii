// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_iterator_range_h
#define dealii_iterator_range_h


#include <deal.II/base/config.h>

#include <deal.II/base/std_cxx20/type_traits.h>
#include <deal.II/base/template_constraints.h>

#include <iterator>

DEAL_II_NAMESPACE_OPEN


// Forward declaration
#ifndef DOXYGEN
template <typename Iterator>
class IteratorOverIterators;
#endif


/**
 * A class that is used to denote a collection of iterators that can be
 * expressed in terms of a range of iterators characterized by a begin and an
 * end iterator. As is common in C++, these ranges are specified as half open
 * intervals defined by a begin iterator and a one-past-the-end iterator.
 *
 * The purpose of this class is so that classes such as Triangulation and
 * DoFHandler can return ranges of cell iterators using an object of the
 * current type from functions such as
 * Triangulation::cell_iterators() and that such an
 * object can then be used in a range-based for loop as supported by C++11,
 * see also
 * @ref CPP11 "C++11 standard".
 *
 * For example, such a loop could look like this if the goal is to set the
 * user flag on every active cell:
 * @code
 *   Triangulation<dim> triangulation;
 *   ...
 *   for (auto &cell : triangulation.active_cell_iterators())
 *     cell->set_user_flag();
 * @endcode
 * In other words, the <code>cell</code> objects are iterators, and the range
 * object returned by Triangulation::active_cell_iterators() and similar
 * functions are conceptually thought of as <i>collections of iterators</i>.
 *
 * Of course, the class may also be used to denote other iterator ranges using
 * different kinds of iterators into other containers.
 *
 *
 * <h3>Class design: Motivation</h3>
 *
 * Informally, the way the C++11 standard describes <a
 * href="http://en.wikipedia.org/wiki/C%2B%2B11#Range-based_for_loop">range-based
 * for loops</a> works as follows: A <i>range-based for loop</i> of the
 * form
 * @code
 *   Container c;
 *   for (auto v : c)
 *     statement;
 * @endcode
 * where <code>c</code> is a container or collection, is equivalent to the
 * following loop:
 * @code
 *   Container c;
 *   for (auto tmp=c.begin(); tmp!=c.end(); ++tmp)
 *     {
 *       auto v = *tmp;
 *       statement;
 *     }
 * @endcode
 * (The precise definition can be found here:
 * https://en.cppreference.com/w/cpp/language/range-for .)
 * In other words, the compiler introduces a temporary variable that
 * <i>iterates</i> over the elements of the container or collection, and the
 * original variable <code>v</code> that appeared in the range-based for loop
 * represents the <i>dereferenced</i> state of these iterators -- namely,
 * the <i>elements</i> of the collection.
 *
 * In the context of loops over cells, we typically want to retain the fact
 * that the loop variable is an iterator, not a value. This is because in
 * deal.II, we never actually use the <i>dereferenced state</i> of a cell
 * iterator: conceptually, it would represent a cell, and technically it is
 * implemented by classes such as CellAccessor and DoFCellAccessor, but these
 * classes are never used explicitly. Consequently, what we would like is that
 * a call such as Triangulation::active_cell_iterators() returns an object
 * that represents a <i>collection of iterators</i> of the kind <code>{begin,
 * begin+1, ..., end-1}</code>. This is conveniently expressed as the half
 * open interval <code>[begin,end)</code>. The loop variable in the
 * range-based for loop would then take on each of these iterators in turn.
 *
 *
 * <h3>Class design: Implementation</h3>
 *
 * To represent the desired semantics as outlined above, this class stores a
 * half-open range of iterators <code>[b,e)</code> of the given template type.
 * Secondly, the class needs to provide begin() and end() functions in such a
 * way that if you <i>dereference</i> the result of IteratorRange::begin(),
 * you get the <code>b</code> iterator. Furthermore, you must be able to
 * increment the object returned by IteratorRange::begin() so that
 * <code>*(std::next(begin())) == b+1</code>. In other words,
 * IteratorRange::begin() must return an iterator that when dereferenced returns
 * an iterator of the template type <code>Iterator</code>: It is an iterator
 * over iterators in the same sense as if you had a pointer into an array of
 * pointers.
 *
 * This is implemented in the form of the IteratorRange::IteratorOverIterators
 * class.
 *
 * @ingroup CPP11
 */
template <typename Iterator>
class IteratorRange
{
public:
  /**
   * Typedef for the iterator type that iterates over other iterators.
   */
  using IteratorOverIterators = dealii::IteratorOverIterators<Iterator>;


  /**
   * Typedef for the iterator type represent by this class.
   */
  using iterator = Iterator;

  /**
   * Default constructor. Create a range represented by two default
   * constructed iterators. This range is likely (depending on the type of the
   * iterators) empty.
   */
  IteratorRange();

  /**
   * Constructor. Constructs a range given the begin and end iterators.
   *
   * @param[in] begin An iterator pointing to the first element of the range
   * @param[in] end   An iterator pointing past the last element represented
   * by this range.
   */
  IteratorRange(const iterator begin, const iterator end);

  /**
   * Return the iterator pointing to the first element of this range.
   */
  IteratorOverIterators
  begin();

  /**
   * Return the iterator pointing to the first element of this range.
   */
  IteratorOverIterators
  begin() const;

  /**
   * Return the iterator pointing to the element past the last element of this
   * range.
   */
  IteratorOverIterators
  end() const;

  /**
   * Return the iterator pointing to the element past the last element of this
   * range.
   */
  IteratorOverIterators
  end();

private:
  /**
   * Iterators characterizing the begin and end of the range.
   */
  const IteratorOverIterators it_begin;
  const IteratorOverIterators it_end;
};



/**
 * A class that implements the semantics of iterators over iterators as
 * discussed in the design sections of the IteratorRange class.
 */
template <typename Iterator>
class IteratorOverIterators
{
public:
  /**
   * Typedef the elements of the collection to give them a name that is more
   * distinct.
   */
  using BaseIterator = Iterator;

  /**
   * Constructor. Initialize this iterator-over-iterator in such a way that
   * it points to the given argument.
   *
   * @param iterator An iterator to which this object is supposed to point.
   */
  explicit IteratorOverIterators(const BaseIterator &iterator);

  /**
   * Dereferencing operator.
   * @return The iterator within the collection currently pointed to.
   */
  const BaseIterator &
  operator*() const;

  /**
   * Dereferencing operator.
   * @return The iterator within the collection currently pointed to.
   */
  const BaseIterator *
  operator->() const;

  /**
   * Prefix increment operator. Move the current iterator to the next
   * element of the collection and return the new value.
   */
  IteratorOverIterators &
  operator++();

  /**
   * Postfix increment operator. Move the current iterator to the next
   * element of the collection, but return the previous value of the
   * iterator.
   */
  IteratorOverIterators
  operator++(int);

  /**
   * Comparison operator
   * @param i_o_i Another iterator over iterators.
   * @return Returns whether the current iterator points to a different
   * object than the iterator represented by the argument.
   */
  bool
  operator!=(const IteratorOverIterators &i_o_i) const;

  /**
   * Implicit conversion operator.
   *
   * @warning When you call this conversion operator (i.e., you convert this
   * iterator-over-iterators to the iterator we are currently pointing to),
   * you obtain a `const` reference to this underlying iterator. The only
   * thing you can really do with this result is dereferencing itself: it
   * presumably points to something useful, but since you don't know where
   * the pointed to object lives, you shouldn't increment or decrement the
   * iterator you get from this operator. As a consequence, the returned
   * iterator is marked as `const`, as this should prevent you from doing
   * anything other than dereference it.
   */
  operator const BaseIterator &() const;

  /**
   * Mark the class as forward iterator and declare some alias which are
   * standard for iterators and are used by algorithms to enquire about the
   * specifics of the iterators they work on.
   */
  using iterator_category = std::forward_iterator_tag;
  using value_type        = Iterator;
  using difference_type   = typename Iterator::difference_type;
  using pointer           = Iterator *;
  using reference         = Iterator &;

private:
  /**
   * The object this iterator currently points to.
   */
  BaseIterator element_of_iterator_collection;
};



/**
 * Create an object of type IteratorRange given the beginning and
 * end iterator.
 */
template <typename BaseIterator>
IteratorRange<BaseIterator>
make_iterator_range(const BaseIterator                             &begin,
                    const std_cxx20::type_identity_t<BaseIterator> &end)
{
  IteratorRange<BaseIterator> ir(begin, end);
  return ir;
}


// ------------------- template member functions


template <typename Iterator>
inline IteratorOverIterators<Iterator>::IteratorOverIterators(
  const BaseIterator &iterator)
  : element_of_iterator_collection(iterator)
{}



template <typename Iterator>
inline const typename IteratorOverIterators<Iterator>::BaseIterator &
IteratorOverIterators<Iterator>::operator*() const
{
  return element_of_iterator_collection;
}



template <typename Iterator>
inline const typename IteratorOverIterators<Iterator>::BaseIterator *
IteratorOverIterators<Iterator>::operator->() const
{
  return &element_of_iterator_collection;
}



template <typename Iterator>
inline IteratorOverIterators<Iterator> &
IteratorOverIterators<Iterator>::operator++()
{
  ++element_of_iterator_collection;
  return *this;
}



template <typename Iterator>
inline IteratorOverIterators<Iterator>
IteratorOverIterators<Iterator>::operator++(int)
{
  const IteratorOverIterators old_value = *this;
  ++element_of_iterator_collection;
  return *old_value;
}



template <typename Iterator>
inline bool
IteratorOverIterators<Iterator>::operator!=(
  const IteratorOverIterators &i_o_i) const
{
  return element_of_iterator_collection != i_o_i.element_of_iterator_collection;
}



template <typename Iterator>
inline IteratorOverIterators<Iterator>::operator const BaseIterator &() const
{
  return element_of_iterator_collection;
}



template <typename Iterator>
inline IteratorRange<Iterator>::IteratorRange()
  : it_begin()
  , it_end()
{}



template <typename Iterator>
inline IteratorRange<Iterator>::IteratorRange(const iterator b,
                                              const iterator e)
  : it_begin(b)
  , it_end(e)
{}


template <typename Iterator>
inline typename IteratorRange<Iterator>::IteratorOverIterators
IteratorRange<Iterator>::begin()
{
  return it_begin;
}


template <typename Iterator>
inline typename IteratorRange<Iterator>::IteratorOverIterators
IteratorRange<Iterator>::begin() const
{
  return it_begin;
}


template <typename Iterator>
inline typename IteratorRange<Iterator>::IteratorOverIterators
IteratorRange<Iterator>::end()
{
  return it_end;
}


template <typename Iterator>
inline typename IteratorRange<Iterator>::IteratorOverIterators
IteratorRange<Iterator>::end() const
{
  return it_end;
}


DEAL_II_NAMESPACE_CLOSE

#endif
