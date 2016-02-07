// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2016 by the deal.II authors
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

#ifndef dealii__array_view_h
#define dealii__array_view_h

#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/table.h>

#include <boost/type_traits/remove_cv.hpp>

#include <vector>

DEAL_II_NAMESPACE_OPEN


/**
 * A class that represents a window of memory locations of type @p ElementType
 * and presents it as if it was an array that can be accessed via an
 * <code>operator[]</code>. In essence, this class is nothing more than just a
 * pointer to the first location and an integer that represents the length of
 * the array in elements. The memory remains owned by whoever allocated it, as
 * this class does not take over ownership.
 *
 * The advantage of using this class is that you don't have to pass around
 * pairs of pointers and that <code>operator[]</code> checks for the validity
 * of the index with which you subscript this array view.
 *
 * This class can handle views to both non-constant and constant memory
 * locations. If you want to represent a view of a constant array, then the
 * template argument type of this class needs to be @p const as well. The
 * following code snippet gives an example:
 * @code
 *   std::vector<int>       array       = get_data();  // a writable array
 *
 *   ArrayView<int> view (&array[5], 5);               // a view of elements 5..9 (inclusive)
 *   view[2] = 42;                                     // array[7] is set to 42
 *
 *   ArrayView<const int> const_view (&array[5], 5);   // same view, but read-only
 *   int element_7 = const_view[2];                    // returns 42
 *   const_view[2] = 42;                               // error, can't write into this view
 * @endcode
 * In either case, accessing an element of a view does not change the
 * ArrayView object itself, and consequently ArrayView::operator[] is a @p
 * const function. This corresponds to the notion that a view simply
 * represents a, well, "view" of memory that is owned by someone else. Thus,
 * accessing elements of the view changes the memory managed by some other
 * object, but not the view itself, allowing us to make ArrayView::operator[]
 * a @p const member function. This is in contrast to, say, std::vector, which
 * manages the memory it points to and changing an element of the std::vector
 * therefore changes the std::vector object itself -- consequently, the
 * std::vector::operator[] is non-@p const.
 *
 * @ingroup data
 * @author Wolfgang Bangerth, 2015
 */
template <typename ElementType>
class ArrayView
{
public:
  /**
   * A typedef that denotes the "value_type" of this container-like class,
   * i.e., the type of the element it "stores" or points to.
   */
  typedef ElementType value_type;

  /**
   * Constructor.
   *
   * @param[in] starting_element A pointer to the first element of the array
   * this object should represent.
   * @param[in] n_elements The length (in elements) of the chunk of memory
   * this object should represent.
   *
   * @note The object that is constructed from these arguments has no
   * knowledge how large the object into which it points really is. As a
   * consequence, whenever you call ArrayView::operator[], the array view can
   * check that the given index is within the range of the view, but it can't
   * check that the view is indeed a subset of the valid range of elements of
   * the underlying object that allocated that range. In other words, you need
   * to ensure that the range of the view specified by the two arguments to
   * this constructor is in fact a subset of the elements of the array into
   * which it points. The appropriate way to do this is to use the
   * make_array_view() functions.
   */
  ArrayView (value_type       *starting_element,
             const std::size_t n_elements);

  /**
   * Copy constructor from array views that point to non-@p const elements. If
   * the current object will point to non-@p const elements, then this is a
   * straight forward copy constructor. On the other hand, if the current
   * type's @p ElementType template argument is a @p const qualified type,
   * then the current constructor is a conversion constructor that converts a
   * non-@p const view to a @p const view, akin to converting a non-@p const
   * pointer to a @p const pointer.
   */
  ArrayView (const ArrayView<typename boost::remove_cv<value_type>::type> &view);

  /**
   * Return the size (in elements) of the view of memory this object
   * represents.
   */
  std::size_t size() const;

  /**
   * Return a reference to the $i$th element of the range represented by the
   * current object.
   *
   * This function is marked as @p const because it does not change the
   * <em>view object</em>. It may however return a reference to a non-@p const
   * memory location depending on whether the template type of the class is @p
   * const or not.
   */
  value_type &operator[] (const std::size_t i) const;

private:
  /**
   * A pointer to the first element of the range of locations in memory that
   * this object represents.
   */
  value_type  *const starting_element;

  /**
   * The length of the array this object represents.
   */
  const std::size_t  n_elements;
};



//---------------------------------------------------------------------------


template <typename ElementType>
inline
ArrayView<ElementType>::ArrayView(value_type        *starting_element,
                                  const std::size_t  n_elements)
  :
  starting_element (starting_element),
  n_elements(n_elements)
{}



template <typename ElementType>
inline
ArrayView<ElementType>::ArrayView(const ArrayView<typename boost::remove_cv<value_type>::type> &view)
  :
  starting_element (&view[0]),
  n_elements(view.size())
{}



template <typename ElementType>
inline
std::size_t
ArrayView<ElementType>::size() const
{
  return n_elements;
}


template <typename ElementType>
inline
typename ArrayView<ElementType>::value_type &
ArrayView<ElementType>::operator[](const std::size_t i) const
{
  Assert (i<n_elements, ExcIndexRange(i, 0, n_elements));

  return *(starting_element + i);
}


/**
 * Create a view to an entire std::vector object. This is equivalent to
 * initializing an ArrayView object with a pointer to the first element and
 * the size of the given argument.
 *
 * This function is used for non-@p const references to objects of vector
 * type. Such objects contain elements that can be written to. Consequently,
 * the return type of this function is a view to a set of writable objects.
 *
 * @param[in] vector The vector for which we want to have an array view
 * object. The array view corresponds to the <em>entire</em> vector.
 *
 * @relates ArrayView
 */
template <typename ElementType>
inline
ArrayView<ElementType>
make_array_view (std::vector<ElementType> &vector)
{
  return ArrayView<ElementType> (&vector[0], vector.size());
}



/**
 * Create a view to an entire std::vector object. This is equivalent to
 * initializing an ArrayView object with a pointer to the first element and
 * the size of the given argument.
 *
 * This function is used for @p const references to objects of vector type
 * because they contain immutable elements. Consequently, the return type of
 * this function is a view to a set of @p const objects.
 *
 * @param[in] vector The vector for which we want to have an array view
 * object. The array view corresponds to the <em>entire</em> vector.
 *
 * @relates ArrayView
 */
template <typename ElementType>
inline
ArrayView<const ElementType>
make_array_view (const std::vector<ElementType> &vector)
{
  return ArrayView<const ElementType> (&vector[0], vector.size());
}



/**
 * Create a view to a part of a std::vector object. This is equivalent to
 * initializing the ArrayView object with a pointer to the @p starting_index-
 * th element and the @p size_of_view as the length of the view.
 *
 * This function is used for non-@p const references to objects of vector
 * type. Such objects contain elements that can be written to. Consequently,
 * the return type of this function is a view to a set of writable objects.
 *
 * @param[in] vector The vector for which we want to have an array view
 * object.
 * @param[in] starting_index The index of the first element of the vector that
 * will be part of this view.
 * @param[in] size_of_view
 *
 * @pre <code>starting_index + size_of_view <= vector.size()</code>
 *
 * @relates ArrayView
 */
template <typename ElementType>
inline
ArrayView<ElementType>
make_array_view (std::vector<ElementType> &vector,
                 const std::size_t         starting_index,
                 const std::size_t         size_of_view)
{
  Assert (starting_index + size_of_view <= vector.size(),
          ExcMessage ("The starting index and size of the view you want to "
                      "create would lead to a view that extends beyond the end "
                      "of the given vector."));
  return ArrayView<ElementType> (&vector[starting_index], size_of_view);
}



/**
 * Create a view to a part of a std::vector object. This is equivalent to
 * initializing the ArrayView object with a pointer to the @p starting_index-
 * th element and the @p size_of_view as the length of the view.
 *
 * This function is used for @p const references to objects of vector type
 * because they contain immutable elements. Consequently, the return type of
 * this function is a view to a set of @p const objects.
 *
 * @param[in] vector The vector for which we want to have an array view
 * object.
 * @param[in] starting_index The index of the first element of the vector that
 * will be part of this view.
 * @param[in] size_of_view
 *
 * @pre <code>starting_index + size_of_view <= vector.size()</code>
 *
 * @relates ArrayView
 */
template <typename ElementType>
inline
ArrayView<const ElementType>
make_array_view (const std::vector<ElementType> &vector,
                 const std::size_t         starting_index,
                 const std::size_t         size_of_view)
{
  Assert (starting_index + size_of_view <= vector.size(),
          ExcMessage ("The starting index and size of the view you want to "
                      "create would lead to a view that extends beyond the end "
                      "of the given vector."));
  return ArrayView<const ElementType> (&vector[starting_index], size_of_view);
}



/**
 * Create a view to an entire row of a Table<2> object. This is equivalent to
 * initializing an ArrayView object with a pointer to the first element of the
 * given row, and the length of the row as the length of the view.
 *
 * This function is used for non-@p const references to objects of Table type.
 * Such objects contain elements that can be written to. Consequently, the
 * return type of this function is a view to a set of writable objects.
 *
 * @param[in] table The Table for which we want to have an array view object.
 * The array view corresponds to an <em>entire</em> row.
 * @param[in] row The index of the row into the table to which this view
 * should correspond.
 *
 * @relates ArrayView
 */
template <typename ElementType>
inline
ArrayView<ElementType>
make_array_view (Table<2,ElementType>                           &table,
                 const typename Table<2,ElementType>::size_type  row)
{
  AssertIndexRange (row, table.size()[0]);
  return ArrayView<ElementType> (&table[row][0], table.size()[1]);
}



/**
 * Create a view to an entire row of a Table<2> object. This is equivalent to
 * initializing an ArrayView object with a pointer to the first element of the
 * given row, and the length of the row as the length of the view.
 *
 * This function is used for @p const references to objects of Table type
 * because they contain immutable elements. Consequently, the return type of
 * this function is a view to a set of @p const objects.
 *
 * @param[in] table The Table for which we want to have an array view object.
 * The array view corresponds to an <em>entire</em> row.
 * @param[in] row The index of the row into the table to which this view
 * should correspond.
 *
 * @relates ArrayView
 */
template <typename ElementType>
inline
ArrayView<const ElementType>
make_array_view (const Table<2,ElementType>                     &table,
                 const typename Table<2,ElementType>::size_type  row)
{
  AssertIndexRange (row, table.size()[0]);
  return ArrayView<const ElementType> (&table[row][0], table.size()[1]);
}



/**
 * Create a view to (a part of) a row of a Table<2> object.
 *
 * This function is used for non-@p const references to objects of Table type.
 * Such objects contain elements that can be written to. Consequently, the
 * return type of this function is a view to a set of writable objects.
 *
 * @param[in] table The Table for which we want to have an array view object.
 * The array view corresponds to an <em>entire</em> row.
 * @param[in] row The index of the row into the table to which this view
 * should correspond.
 * @param[in] starting_column The index of the column into the given row of
 * the table that corresponds to the first element of this view.
 * @param[in] size_of_view The number of elements this view should have. This
 * corresponds to the number of columns in the current row to which the view
 * should correspond.
 *
 * @relates ArrayView
 */
template <typename ElementType>
inline
ArrayView<ElementType>
make_array_view (Table<2,ElementType>                           &table,
                 const typename Table<2,ElementType>::size_type  row,
                 const typename Table<2,ElementType>::size_type  starting_column,
                 const std::size_t                               size_of_view)
{
  AssertIndexRange (row, table.size()[0]);
  AssertIndexRange (starting_column, table.size()[1]);
  Assert (starting_column + size_of_view <= table.size()[1],
          ExcMessage ("The starting index and size of the view you want to "
                      "create would lead to a view that extends beyond the end "
                      "of a column of the given table."));
  return ArrayView<ElementType> (&table[row][starting_column], size_of_view);
}



/**
 * Create a view to (a part of) a row of a Table<2> object.
 *
 * This function is used for @p const references to objects of Table type
 * because they contain immutable elements. Consequently, the return type of
 * this function is a view to a set of @p const objects.
 *
 * @param[in] table The Table for which we want to have an array view object.
 * The array view corresponds to an <em>entire</em> row.
 * @param[in] row The index of the row into the table to which this view
 * should correspond.
 * @param[in] starting_column The index of the column into the given row of
 * the table that corresponds to the first element of this view.
 * @param[in] size_of_view The number of elements this view should have. This
 * corresponds to the number of columns in the current row to which the view
 * should correspond.
 *
 * @relates ArrayView
 */
template <typename ElementType>
inline
ArrayView<const ElementType>
make_array_view (const Table<2,ElementType>                     &table,
                 const typename Table<2,ElementType>::size_type  row,
                 const typename Table<2,ElementType>::size_type  starting_column,
                 const std::size_t                               size_of_view)
{
  AssertIndexRange (row, table.size()[0]);
  AssertIndexRange (starting_column, table.size()[1]);
  Assert (starting_column + size_of_view <= table.size()[1],
          ExcMessage ("The starting index and size of the view you want to "
                      "create would lead to a view that extends beyond the end "
                      "of a column of the given table."));
  return ArrayView<const ElementType> (&table[row][starting_column], size_of_view);
}



DEAL_II_NAMESPACE_CLOSE

#endif
