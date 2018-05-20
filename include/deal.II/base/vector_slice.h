// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2018 by the deal.II authors
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

#ifndef dealii_vector_slice_h
#define dealii_vector_slice_h

#include <deal.II/base/array_view.h>
#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>

DEAL_II_NAMESPACE_OPEN

/**
 * Filter a range out of any object having a random access <tt>operator[]
 * (unsigned int)</tt> and a function <tt>size() const</tt>.
 *
 * The use of this object is straightforward. It duplicates the random access
 * operator of the <tt>VectorType</tt> and adds an offset to every index.
 *
 * Some precautions have to be taken if it is used for a constant vector: the
 * VectorSlice object has to be constant, too. The appropriate initialization
 * sequence is like this:
 *
 * @code
 *   void f(const std::vector<int>& v)
 *   {
 *     const VectorSlice<const std::vector<int> > slice(v,...);
 *     ...
 *   }
 * @endcode
 *
 * @ingroup data
 * @author Guido Kanschat, 2004
 */
template <typename VectorType>
class VectorSlice
{
public:
  /**
   * Construct a vector slice containing the whole vector. Comes handy, if you
   * did not want to have a slice at all, but the function you call wants it:
   * just put in the vector itself as argument and let this constructor make a
   * slice for you.
   */
  VectorSlice(VectorType& v);
  /**
   * The real constructor for a vector slice, allowing you to specify the
   * start index and the length of the slice.
   */
  VectorSlice(VectorType& v, unsigned int start, unsigned int length);

  /**
   * Conversion operator to an ArrayView object that represents an array of
   * non-const elements pointing to the same location as the current object.
   */
  operator ArrayView<typename VectorType::value_type*>();

  /**
   * Conversion operator to an ArrayView object that represents an array of
   * const elements pointing to the same location as the current object.
   */
  operator ArrayView<const typename VectorType::value_type*>() const;

  /**
   * Return the length of the slice using the same interface as
   * <tt>std::vector</tt>.
   */
  unsigned int
  size() const;

  /**
   * Return a reference to the $i$th element of the range represented by the
   * current object.
   */
  typename VectorType::reference operator[](unsigned int i);

  /**
   * Return a @p const reference to the $i$th element of the range represented
   * by the current object.
   */
  typename VectorType::const_reference operator[](unsigned int i) const;

  /**
   * Standard-conforming iterator function.
   */
  typename VectorType::iterator
  begin();

  /**
   * Standard-conforming iterator function.
   */
  typename VectorType::const_iterator
  begin() const;

  /**
   * Standard-conforming iterator function.
   */
  typename VectorType::iterator
  end();

  /**
   * Standard-conforming iterator function.
   */
  typename VectorType::const_iterator
  end() const;

private:
  /**
   * The vector we extract from.
   */
  VectorType& v;
  /**
   * The start index of the slice.
   */
  const unsigned int start;
  /**
   * The length of the slice.
   */
  const unsigned int length;
};

/**
 * Helper function for creating temporary objects without typing template
 * arguments.
 *
 * @relatesalso VectorSlice
 * @author Guido Kanschat, 2004
 */
template <typename VectorType>
inline const VectorSlice<const VectorType>
make_slice(VectorType& v)
{
  const VectorSlice<const VectorType> r(v);
  return r;
}

/**
 * Helper function for creating temporary objects without typing template
 * arguments.
 *
 * @relatesalso VectorSlice
 * @author Guido Kanschat, 2004
 */
template <typename VectorType>
inline const VectorSlice<const VectorType>
make_slice(VectorType& v, const unsigned int start, const unsigned int length)
{
  const VectorSlice<const VectorType> r(v, start, length);
  return r;
}

//---------------------------------------------------------------------------

template <typename VectorType>
inline VectorSlice<VectorType>::VectorSlice(VectorType& v)
  : v(v), start(0), length(v.size())
{}

template <typename VectorType>
inline VectorSlice<VectorType>::VectorSlice(VectorType&  v,
                                            unsigned int start,
                                            unsigned int length)
  : v(v), start(start), length(length)
{
  Assert((start + length <= v.size()),
         ExcIndexRange(length, 0, v.size() - start + 1));
}

template <typename VectorType>
inline unsigned int
VectorSlice<VectorType>::size() const
{
  return length;
}

template <typename VectorType>
VectorSlice<VectorType>::operator ArrayView<typename VectorType::value_type*>()
{
  return ArrayView<typename VectorType::value_type*>(&v[start], length);
}

template <typename VectorType>
VectorSlice<VectorType>::
operator ArrayView<const typename VectorType::value_type*>() const
{
  return ArrayView<const typename VectorType::value_type*>(&v[start], length);
}

template <typename VectorType>
inline typename VectorType::reference VectorSlice<VectorType>::
                                      operator[](unsigned int i)
{
  Assert((i < length), ExcIndexRange(i, 0, length));

  return v[start + i];
}

template <typename VectorType>
inline typename VectorType::const_reference VectorSlice<VectorType>::
                                            operator[](unsigned int i) const
{
  Assert((i < length), ExcIndexRange(i, 0, length));

  return v[start + i];
}

template <typename VectorType>
inline typename VectorType::const_iterator
VectorSlice<VectorType>::begin() const
{
  return v.begin() + start;
}

template <typename VectorType>
inline typename VectorType::iterator
VectorSlice<VectorType>::begin()
{
  return v.begin() + start;
}

template <typename VectorType>
inline typename VectorType::const_iterator
VectorSlice<VectorType>::end() const
{
  return v.begin() + start + length;
}

template <typename VectorType>
inline typename VectorType::iterator
VectorSlice<VectorType>::end()
{
  return v.begin() + start + length;
}

DEAL_II_NAMESPACE_CLOSE

#endif
