// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_hp_collection_h
#define dealii_hp_collection_h

#include <deal.II/base/config.h>

#include <deal.II/base/enable_observer_pointer.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/memory_consumption.h>

#include <iterator>
#include <memory>
#include <vector>

DEAL_II_NAMESPACE_OPEN

namespace hp
{
  /**
   * Exception thrown when comparing hp::Collection iterators into different
   * objects.
   *
   * @ingroup Exceptions
   */
  DeclExceptionMsg(ExcDifferentCollection,
                   "You are trying to compare iterators into different "
                   "hp::Collection objects.");

  /**
   * An iterator for hp::Collection.
   */
  template <typename T>
  class CollectionIterator
  {
  public:
    /**
     * Constructor.
     *
     * @param data The actual data of hp::Collection.
     * @param index The current index.
     */
    CollectionIterator(const std::vector<std::shared_ptr<const T>> &data,
                       const std::size_t                            index)
      : data(&data)
      , index(index)
    {}

    /**
     * Copy constructor.
     */
    CollectionIterator(const CollectionIterator<T> &other) = default;

    /**
     * Copy assignment.
     */
    CollectionIterator<T> &
    operator=(const CollectionIterator<T> &other) = default;

    /**
     * Compare for equality.
     */
    bool
    operator==(const CollectionIterator<T> &other) const
    {
      Assert(this->data == other.data, ExcDifferentCollection());
      return this->index == other.index;
    }

    /**
     * Compare for inequality.
     */
    bool
    operator!=(const CollectionIterator<T> &other) const
    {
      Assert(this->data == other.data, ExcDifferentCollection());
      return this->index != other.index;
    }

    /**
     * Compare indices.
     */
    bool
    operator<(const CollectionIterator<T> &other) const
    {
      Assert(this->data == other.data, ExcDifferentCollection());
      return this->index < other.index;
    }

    /**
     * Compare indices.
     */
    bool
    operator<=(const CollectionIterator<T> &other) const
    {
      Assert(this->data == other.data, ExcDifferentCollection());
      return this->index <= other.index;
    }

    /**
     * Compare indices.
     */
    bool
    operator>(const CollectionIterator<T> &other) const
    {
      Assert(this->data == other.data, ExcDifferentCollection());
      return this->index > other.index;
    }

    /**
     * Compare indices.
     */
    bool
    operator>=(const CollectionIterator<T> &other) const
    {
      Assert(this->data == other.data, ExcDifferentCollection());
      return this->index >= other.index;
    }

    /**
     * Dereferencing operator: returns the value of the current index.
     */
    const T &
    operator*() const
    {
      AssertIndexRange(index, data->size());
      return *(*data)[index];
    }

    /**
     * Prefix <tt>++</tt> operator: <tt>++iterator</tt>. This operator advances
     * the iterator to the next index and returns a reference to
     * <tt>*this</tt>.
     */
    CollectionIterator<T> &
    operator++()
    {
      AssertIndexRange(index + 1, data->size() + 1);
      ++index;
      return *this;
    }

    /**
     * This operator advances the iterator by @p offset and returns a
     * reference to <tt>*this</tt>.
     */
    CollectionIterator<T> &
    operator+=(const std::size_t offset)
    {
      AssertIndexRange(index + offset, data->size() + 1);
      index += offset;
      return *this;
    }

    /**
     * Prefix <tt>--</tt> operator: <tt>--iterator</tt>. This operator advances
     * the iterator to the previous index and returns a reference to
     * <tt>*this</tt>.
     */
    CollectionIterator<T> &
    operator--()
    {
      Assert(
        index > 0,
        ExcMessage(
          "You can't decrement an iterator that is already at the beginning of the range."));
      --index;
      return *this;
    }

    /**
     * Create new iterator, which is shifted by @p offset.
     */
    CollectionIterator<T>
    operator+(const std::size_t &offset) const
    {
      AssertIndexRange(index + offset, T::size() + 1);
      return CollectionIterator<T>(*data, index + offset);
    }

    /**
     * Compute distance between this iterator and iterator @p other.
     */
    std::ptrdiff_t
    operator-(const CollectionIterator<T> &other) const
    {
      return static_cast<std::ptrdiff_t>(index) -
             static_cast<std::ptrdiff_t>(other.index);
    }

  private:
    /**
     * Pointer to the actual data of hp::Collection.
     */
    const std::vector<std::shared_ptr<const T>> *data;

    /**
     * Current index.
     */
    std::size_t index;
  };

  /**
   * This class implements a collection of objects.
   *
   * It implements the concepts stated in the
   * @ref hpcollection
   * topic described in the doxygen documentation.
   *
   * @ingroup hp hpcollection
   */
  template <typename T>
  class Collection : public EnableObserverPointer
  {
  public:
    /**
     * Default constructor. Leads to an empty collection that can later be
     * filled using push_back().
     */
    Collection() = default;

    /**
     * Add a new object.
     */
    void
    push_back(const std::shared_ptr<const T> &new_entry);

    /**
     * Return the object which was specified by the user for the
     * active FE index which is provided as a parameter to this method.
     *
     * @pre @p index must be between zero and the number of elements of the
     * collection.
     */
    const T &
    operator[](const unsigned int index) const;

    /**
     * Return the number of objects stored in this container.
     */
    unsigned int
    size() const;

    /**
     * The return value of this function is equivalent to
     * <code>size() == 0</code>
     */
    bool
    empty() const;

    /**
     * Determine an estimate for the memory consumption (in bytes) of this
     * object.
     */
    std::size_t
    memory_consumption() const;

    /**
     * @return An iterator pointing to the beginning of the underlying data (`const`
     * version).
     */
    CollectionIterator<T>
    begin() const;

    /**
     * @return An iterator pointing to the end of the underlying data (`const`
     * version).
     */
    CollectionIterator<T>
    end() const;

  private:
    /**
     * The real container, which stores pointers to the different objects.
     */
    std::vector<std::shared_ptr<const T>> entries;
  };


  /* --------------- inline functions ------------------- */



  template <typename T>
  std::size_t
  Collection<T>::memory_consumption() const
  {
    return (sizeof(*this) + MemoryConsumption::memory_consumption(entries));
  }



  template <typename T>
  void
  Collection<T>::push_back(const std::shared_ptr<const T> &new_entry)
  {
    entries.push_back(new_entry);
  }



  template <typename T>
  inline unsigned int
  Collection<T>::size() const
  {
    return entries.size();
  }



  template <typename T>
  inline bool
  Collection<T>::empty() const
  {
    return this->size() == 0;
  }



  template <typename T>
  inline const T &
  Collection<T>::operator[](const unsigned int index) const
  {
    AssertIndexRange(index, entries.size());
    return *entries[index];
  }



  template <typename T>
  CollectionIterator<T>
  Collection<T>::begin() const
  {
    return CollectionIterator<T>(entries, 0);
  }



  template <typename T>
  CollectionIterator<T>
  Collection<T>::end() const
  {
    return CollectionIterator<T>(entries, entries.size());
  }

} // namespace hp


DEAL_II_NAMESPACE_CLOSE

namespace std
{
  /**
   * Iterator traits for hp::CollectionIterator.
   */
  template <class T>
  struct iterator_traits<dealii::hp::CollectionIterator<T>>
    : public iterator_traits<
        typename std::vector<std::shared_ptr<const T>>::iterator>
  {};
} // namespace std

#endif
