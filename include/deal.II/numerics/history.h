// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2020 by the deal.II authors
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

#ifndef dealii_storage_h
#define dealii_storage_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

#include <deque>
#include <memory>
#include <type_traits>

DEAL_II_NAMESPACE_OPEN

/**
 * A helper class to store a finite-size collection of objects of type `T`.
 * If the number of elements exceeds the specified maximum size of the
 * container, the oldest element is removed. Additionally, random access and
 * removal of elements is implemented. Indexing is done relative to the last
 * added element.
 *
 * In order to optimize the container for usage with
 * memory-demanding objects (i.e. linear algebra vectors), the removal of an
 * element does not free the memory. Instead the element
 * is being kept in a separate cache so that subsequent addition does not
 * require re-allocation of memory.
 *
 * The primary usage of this class is in solvers to store a history of vectors.
 * That is, if at the iteration $k$ we store $m$ vectors from previous
 * iterations
 * $\{k-1,k-2,...,k-m\}$, then addition of the new element will make
 * the object contain elements from iterations $\{k,k-1,k-2,...,k-m+1\}$.
 *
 * @author Denis Davydov, 2018
 */
template <typename T>
class FiniteSizeHistory
{
public:
  static_assert(
    std::is_default_constructible<T>::value,
    "This class requires that the elements of type T are default constructible.");

  /**
   * Constructor.
   *
   * @param max_elements maximum number of elements to be stored in the
   * history.
   */
  FiniteSizeHistory(const std::size_t max_elements = 0);

  /**
   * Add new object by copying.
   * If the maximum number of elements is reached, the oldest element is
   * removed.
   */
  void
  add(const T &element);

  /**
   * Remove an element with index @p index,
   * counting from the last added element.
   * `index==0` therefore corresponds to removing
   * the newset element.
   */
  void
  remove(const std::size_t index);

  /**
   * Read/write access to an element with index @p index,
   * counting from the last added element.
   * `index==0` therefore corresponds to the newset element.
   */
  T &operator[](const std::size_t index);

  /**
   * Read access to an element with index @p index,
   * counting from the last added element.
   * `index==0` therefore corresponds to the newset element.
   */
  const T &operator[](const std::size_t index) const;

  /**
   * Return the current size of the history.
   */
  std::size_t
  size() const;

  /**
   * Return the maximum allowed size of the history.
   */
  std::size_t
  max_size() const;

  /**
   * Clear the contents, including the cache.
   */
  void
  clear();

private:
  /**
   * Maximum number of elements to be stored.
   */
  std::size_t max_n_elements;

  /**
   * A deque which stores the data.
   */
  std::deque<std::unique_ptr<T>> data;

  /**
   * A deque to cache data, in particular for the case when
   * removal is followed by addition.
   */
  std::deque<std::unique_ptr<T>> cache;
};



// -------------------  inline and template functions ----------------
#ifndef DOXYGEN



template <typename T>
FiniteSizeHistory<T>::FiniteSizeHistory(const std::size_t max_elements)
  : max_n_elements(max_elements)
{}



template <typename T>
void
FiniteSizeHistory<T>::remove(const std::size_t ind)
{
  AssertIndexRange(ind, data.size());
  auto el = std::move(data[ind]);
  data.erase(data.begin() + ind);

  cache.push_back(std::move(el));

  // whatever we do, we shall not store more than the maximum number of
  // elements
  Assert(data.size() + cache.size() <= max_n_elements, ExcInternalError());
}



template <typename T>
void
FiniteSizeHistory<T>::add(const T &element)
{
  std::unique_ptr<T> new_el;
  if (data.size() < max_n_elements)
    // have not reached the maximum number of elements yet
    {
      if (cache.size() == 0)
        // nothing is cached, just copy a given element
        {
          new_el = std::make_unique<T>(element);
        }
      else
        // something is cached, take one element and copy
        // the user provided one there.
        {
          new_el    = std::move(cache.back());
          (*new_el) = element;

          cache.pop_back(); // removes a pointer that is now a nullptr anyway
        }
    }
  else
    // we reached the maximum number of elements and
    // thus have to re-order/cycle elements currently stored
    {
      new_el    = std::move(data.back());
      (*new_el) = element;

      data.pop_back(); // removes a pointer that is now a nullptr anyway
    }

  // finally insert the new one where appropriate
  data.push_front(std::move(new_el));

  // whatever we do, we shall not store more than the maximum number of
  // elements
  Assert(data.size() + cache.size() <= max_n_elements, ExcInternalError());
}



template <typename T>
T &FiniteSizeHistory<T>::operator[](const std::size_t ind)
{
  AssertIndexRange(ind, data.size());
  return *data[ind];
}



template <typename T>
const T &FiniteSizeHistory<T>::operator[](const std::size_t ind) const
{
  AssertIndexRange(ind, data.size());
  return *data[ind];
}



template <typename T>
std::size_t
FiniteSizeHistory<T>::size() const
{
  return data.size();
}



template <typename T>
std::size_t
FiniteSizeHistory<T>::max_size() const
{
  return max_n_elements;
}



template <typename T>
void
FiniteSizeHistory<T>::clear()
{
  data.clear();
  cache.clear();
}

#endif // Doxygen

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_storage_h
