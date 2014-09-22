// ---------------------------------------------------------------------
//
// Copyright (C) 2014 by the deal.II authors
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

#ifndef __deal2__iterator_range_h
#define __deal2__iterator_range_h


#include <deal.II/base/config.h>

#include <iterator>


DEAL_II_NAMESPACE_OPEN


/**
 * A class that is used to denote a collection of iterators that can
 * be expressed in terms of a range of iterators characterized by a begin
 * and an end iterator. As is common in C++, these ranges are specified as
 * half open intervals defined by a begin iterator and a one-past-the-end
 * iterator.
 *
 * The purpose of this class is so that classes such as Triangulation and
 * DoFHandler can return ranges of cell iterators using an object of the
 * current type from functions such as Triangulation::cells() and that such an
 * object can then be used in a range-based for loop as supported by C++11,
 * see also @ref CPP11 "C++11 standard".
 *
 * For example, such a loop could look like this if the goal is to set
 * the user flag on every active cell:
 * @code
 *   Triangulation<dim> triangulation;
 *   ...
 *   for (auto cell : triangulation.active_cell_iterators())
 *     cell->set_user_flag();
 * @endcode
 * In other words, the <code>cell</code> objects are iterators, and the
 * range object returned by Triangulation::active_cell_iterators() and
 * similar functions are conceptually thought of as <i>collections of
 * iterators</i>.
 *
 * Of course, the class may also be used to denote other iterator
 * ranges using different kinds of iterators into other containers.
 *
 *
 * <h3>Class design: Motivation</h3>
 *
 * Informally, the way the C++11 standard describes
 * <a href="http://en.wikipedia.org/wiki/C%2B%2B11#Range-based_for_loop">range-based
 * for loops</a> works as follows: A <i>range-based for loop</i> of the form
 * @code
 *   Container c;
 *   for (auto v : c)
 *     statement;
 * @endcode
 * where <code>c</code> is a container or collection, is equivalent to the following
 * loop:
 * @code
 *   Container c;
 *   for (auto tmp=c.begin(); tmp!=c.end(); ++tmp)
 *     {
 *       auto v = *tmp;
 *       statement;
 *     }
 * @endcode
 * In other words, the compiler introduces a temporary variable that <i>iterates</i>
 * over the elements of the container or collection, and the original variable
 * <code>v</code> that appeared in the range-based for loop represents the
 * <i>dereferenced</i> state of these iterators -- in other words, the
 * <i>elements</i> of the collection.
 *
 * In the context of loops over cells, we typically want to retain the fact that
 * the loop variable is an iterator, not a value. This is because in deal.II,
 * we never actually use the <i>dereferenced state</i> of a cell iterator:
 * conceptually, it would represent a cell, and technically it is implemented
 * by classes such as CellAccessor and DoFCellAccessor, but these classes are
 * never used explicitly. Consequently, what we would like is that a call
 * such as Triangulation::active_cell_iterators() returns an object that
 * represents a <i>collection of iterators</i> of the kind
 * <code>{begin, begin+1, ..., end-1}</code>. This is conveniently expressed
 * as the half open interval <code>[begin,end)</code>. The loop variable in the
 * range-based for loop would then take on each of these iterators in turn.
 *
 *
 * <h3>Class design: Implementation</h3>
 *
 * To represent the desired semantics as outlined above, this class
 * stores a half-open range of iterators <code>[b,e)</code> of
 * the given template type. Secondly, the class needs to provide begin()
 * and end() functions in such a way that if you <i>dereference</i> the
 * result of IteratorRange::begin(), you get the <code>b</code> iterator.
 * Furthermore, you must be able to increment the object returned by
 * IteratorRange::begin() so that <code>*(++begin()) == b+1</code>.
 * In other words, IteratorRange::begin() must return an iterator that
 * when dereferenced returns an iterator of the template type
 * <code>Iterator</code>: It is an iterator over iterators in the same
 * sense as if you had a pointer into an array of pointers.
 *
 * This is implemented in the form of the IteratorRange::IteratorOverIterators
 * class.
 *
 * @ingroup CPP11
 * @author Wolfgang Bangerth, 2014
 */
template <typename Iterator>
class IteratorRange
{
public:
  /**
   * A class that implements the semantics of iterators over iterators
   * as discussed in the design sections of the IteratorRange class.
   */
  class IteratorOverIterators : public std::iterator<std::forward_iterator_tag, Iterator,
    typename Iterator::difference_type>
  {
  public:
    /**
     * Typedef the elements of the collection to give them a name that is
     * more distinct.
     */
    typedef Iterator BaseIterator;

    /**
     * Constructor. Initialize this iterator-over-iterator in
     * such a way that it points to the given argument.
     *
     * @param iterator An iterator to which this object
     *   is supposed to point.
     */
    IteratorOverIterators (const BaseIterator &iterator);

    /**
     * Dereferencing operator.
     * @return The iterator within the collection currently pointed to.
     */
    BaseIterator operator* () const;

    /**
     * Dereferencing operator.
     * @return The iterator within the collection currently pointed to.
     */
    const BaseIterator *operator-> () const;

    /**
     * Prefix increment operator. Move the current iterator to the next
     * element of the collection and return the new value.
     */
    IteratorOverIterators &operator ++ ();

    /**
     * Postfix increment operator. Move the current iterator to the next
     * element of the collection, but return the previous value of the
     * iterator.
     */
    IteratorOverIterators operator ++ (int);

    /**
     * Comparison operator
     * @param i_o_i Another iterator over iterators.
     * @return Returns whether the current iterator points to a
     *   different object than the iterator represented by the
     *   argument.
     */
    bool operator != (const IteratorOverIterators &i_o_i);

  private:
    /**
     * The object this iterator currently points to.
     */
    BaseIterator element_of_iterator_collection;
  };


  /**
   * Typedef for the iterator type represent by this class.
   */
  typedef Iterator iterator;

  /**
   * Default constructor. Create a range represented by two
   * default constructed iterators. This range is likely (depending
   * on the type of the iterators) empty.
   */
  IteratorRange();

  /**
   * Constructor. Constructs a range given the begin and end iterators.
   *
   * @param[in] begin An iterator pointing to the first element of the range
   * @param[in] end   An iterator pointing past the last element represented
   *   by this range.
   */
  IteratorRange (const iterator begin,
                 const iterator end);

  /**
   * Return the iterator pointing to the first element of this range.
   */
  IteratorOverIterators begin();

  /**
   * Return the iterator pointing to the element past the last
   * element of this range.
   */
  IteratorOverIterators end();

private:
  /**
   * Iterators characterizing the begin and end of the range.
   */
  const iterator it_begin;
  const iterator it_end;
};


// ------------------- template member functions


template <typename Iterator>
inline
IteratorRange<Iterator>::IteratorOverIterators::
IteratorOverIterators (const BaseIterator &iterator)
  :
  element_of_iterator_collection (iterator)
{}



template <typename Iterator>
inline
typename IteratorRange<Iterator>::IteratorOverIterators::BaseIterator
IteratorRange<Iterator>::IteratorOverIterators::operator* () const
{
  return element_of_iterator_collection;
}



template <typename Iterator>
inline
const typename IteratorRange<Iterator>::IteratorOverIterators::BaseIterator *
IteratorRange<Iterator>::IteratorOverIterators::operator-> () const
{
  return &element_of_iterator_collection;
}



template <typename Iterator>
inline
typename IteratorRange<Iterator>::IteratorOverIterators &
IteratorRange<Iterator>::IteratorOverIterators::operator ++ ()
{
  ++element_of_iterator_collection;
  return *this;
}



template <typename Iterator>
inline
typename IteratorRange<Iterator>::IteratorOverIterators
IteratorRange<Iterator>::IteratorOverIterators::operator ++ (int)
{
  const IteratorOverIterators old_value = *this;
  ++element_of_iterator_collection;
  return *old_value;
}



template <typename Iterator>
inline
bool
IteratorRange<Iterator>::IteratorOverIterators::operator != (const IteratorOverIterators &i_o_i)
{
  return element_of_iterator_collection != i_o_i.element_of_iterator_collection;
}


template <typename Iterator>
inline
IteratorRange<Iterator>::IteratorRange ()
  :
  it_begin(),
  it_end()
{}



template <typename Iterator>
inline
IteratorRange<Iterator>::IteratorRange (const iterator b,
                                        const iterator e)
  :
  it_begin(b),
  it_end(e)
{}


template <typename Iterator>
inline
typename IteratorRange<Iterator>::IteratorOverIterators
IteratorRange<Iterator>::begin()
{
  return IteratorOverIterators(it_begin);
}


template <typename Iterator>
inline
typename IteratorRange<Iterator>::IteratorOverIterators
IteratorRange<Iterator>::end()
{
  return IteratorOverIterators(it_end);
}


DEAL_II_NAMESPACE_CLOSE

#endif
