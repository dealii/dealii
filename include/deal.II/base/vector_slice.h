// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2020 by the deal.II authors
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

#ifndef dealii_vector_slice_h
#define dealii_vector_slice_h

#include <deal.II/base/config.h>

#include <deal.II/base/array_view.h>
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
 *
 * @deprecated Use the more general ArrayView class instead.
 */
template <typename VectorType>
class DEAL_II_DEPRECATED VectorSlice
  : public ArrayView<
      typename std::conditional<std::is_const<VectorType>::value,
                                const typename VectorType::value_type,
                                typename VectorType::value_type>::type>
{
public:
  /**
   * Construct a vector slice containing the whole vector. Comes handy, if you
   * did not want to have a slice at all, but the function you call wants it:
   * just put in the vector itself as argument and let this constructor make a
   * slice for you.
   */
  VectorSlice(VectorType &v);
  /**
   * The real constructor for a vector slice, allowing you to specify the
   * start index and the length of the slice.
   */
  VectorSlice(VectorType &v, unsigned int start, unsigned int length);

protected:
  /**
   * Alias for the base class name.
   */
  using ArrayViewType =
    ArrayView<typename std::conditional<std::is_const<VectorType>::value,
                                        const typename VectorType::value_type,
                                        typename VectorType::value_type>::type>;
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
make_slice(VectorType &v)
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
make_slice(VectorType &v, const unsigned int start, const unsigned int length)
{
  const VectorSlice<const VectorType> r(v, start, length);
  return r;
}



//---------------------------------------------------------------------------

template <typename VectorType>
inline VectorSlice<VectorType>::VectorSlice(VectorType &v)
  : ArrayViewType(make_array_view(std::begin(v), std::end(v)))
{}



template <typename VectorType>
inline VectorSlice<VectorType>::VectorSlice(VectorType & v,
                                            unsigned int start,
                                            unsigned int length)
  : ArrayViewType(
      make_array_view(std::begin(v) + start, std::begin(v) + start + length))
{
  AssertIndexRange(length, v.size() - start + 1);
}

DEAL_II_NAMESPACE_CLOSE

#endif
