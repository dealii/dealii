// ---------------------------------------------------------------------
//
// Copyright (C) 2021 by the deal.II authors
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

#ifndef dealii_hp_collection_h
#define dealii_hp_collection_h

#include <deal.II/base/config.h>

#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/subscriptor.h>

#include <memory>
#include <vector>

DEAL_II_NAMESPACE_OPEN

namespace hp
{
  /**
   * This class implements a collection of objects.
   *
   * It implements the concepts stated in the @ref hpcollection
   * module described in the doxygen documentation.
   *
   * @ingroup hp hpcollection
   */
  template <typename T, int N>
  class CollectionBase : public Subscriptor
  {
  public:
    /**
     * Default constructor. Leads to an empty collection that can later be
     * filled using push_back().
     */
    CollectionBase() = default;

    /**
     * TODO.
     */
    CollectionBase(const Table<N, std::shared_ptr<const T>> &entries)
      : entries(entries)
    {}

    /**
     * Determine an estimate for the memory consumption (in bytes) of this
     * object.
     */
    std::size_t
    memory_consumption() const;

    const Table<N, std::shared_ptr<const T>> &
    get_entries() const
    {
      return entries;
    }

  protected:
    /**
     * The real container, which stores pointers to the different objects.
     */
    Table<N, std::shared_ptr<const T>> entries;
  };



  template <typename T, int N, typename U = T>
  class Collection;



  template <typename T, typename U>
  class Collection<T, 0, U>
  {
  public:
    /**
     * TODO.
     */
    Collection() = default;
  };



  template <typename T, typename U>
  class Collection<T, 1, U> : public CollectionBase<T, 1>
  {
  public:
    /**
     * Default constructor. Leads to an empty collection that can later be
     * filled using push_back().
     */
    Collection() = default;

    /**
     * TODO.
     */
    Collection(const Table<1, std::shared_ptr<const T>> &entries)
      : CollectionBase<T, 1>(entries)
    {}

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
    const T &operator[](const unsigned int index) const;

    /**
     * Return the number of objects stored in this container.
     */
    unsigned int
    size() const;
  };



  template <typename T, typename U>
  class Collection<T, 2, U> : public CollectionBase<T, 2>
  {
  public:
    /**
     * TODO
     */
    Collection() = default;

    /**
     * TODO.
     */
    Collection(const Table<2, std::shared_ptr<const T>> &entries)
      : CollectionBase<T, 2>(entries)
    {}

    /**
     * TODO
     */
    const U operator[](const unsigned int index) const
    {
      Table<1, std::shared_ptr<const T>> new_enties(this->entries.size()[1]);

      for (unsigned int i = 0; i < this->entries.size()[1]; ++i)
        new_enties[i] = this->entries[index][i];

      return U(new_enties);
    }

    /**
     * TODO
     */
    unsigned int
    size() const
    {
      return this->entries.size()[0];
    }
  };


  /* --------------- inline functions ------------------- */



  template <typename T, int N>
  std::size_t
  CollectionBase<T, N>::memory_consumption() const
  {
    return (sizeof(*this) + MemoryConsumption::memory_consumption(entries));
  }



  template <typename T, typename U>
  void
  Collection<T, 1, U>::push_back(const std::shared_ptr<const T> &new_entry)
  {
    const auto         temp     = this->entries;
    const unsigned int old_size = this->size();

    this->entries = Table<1, std::shared_ptr<const T>>(old_size + 1);

    for (unsigned int i = 0; i < old_size; ++i)
      this->entries[i] = temp[i];

    this->entries[old_size] = new_entry;
  }



  template <typename T, typename U>
  inline unsigned int
  Collection<T, 1, U>::size() const
  {
    return this->entries.size()[0];
  }



  template <typename T, typename U>
  inline const T &Collection<T, 1, U>::
                  operator[](const unsigned int index) const
  {
    AssertIndexRange(index, this->size());
    return *this->entries[index];
  }

} // namespace hp


DEAL_II_NAMESPACE_CLOSE

#endif
