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
  template <typename T, int N = 1>
  class CollectionBase : public Subscriptor
  {
  public:
    /**
     * Default constructor. Leads to an empty collection that can later be
     * filled using push_back().
     */
    CollectionBase() = default;

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

    /**
     * Determine an estimate for the memory consumption (in bytes) of this
     * object.
     */
    std::size_t
    memory_consumption() const;

  protected:
    /**
     * The real container, which stores pointers to the different objects.
     */
    Table<N, std::shared_ptr<const T>> entries;
  };



  template <typename T, int N>
  class Collection;



  template <typename T>
  class Collection<T, 1> : public CollectionBase<T, 1>
  {
  public:
    Collection() = default;
  };


  /* --------------- inline functions ------------------- */



  template <typename T, int N>
  std::size_t
  CollectionBase<T, N>::memory_consumption() const
  {
    return (sizeof(*this) + MemoryConsumption::memory_consumption(entries));
  }



  template <typename T, int N>
  void
  CollectionBase<T, N>::push_back(const std::shared_ptr<const T> &new_entry)
  {
    const auto         temp     = entries;
    const unsigned int old_size = this->size();

    AssertDimension(N, 1);

    entries = Table<N, std::shared_ptr<const T>>(old_size + 1);

    for (unsigned int i = 0; i < old_size; ++i)
      entries[i] = temp[i];

    entries[old_size] = new_entry;
  }



  template <typename T, int N>
  inline unsigned int
  CollectionBase<T, N>::size() const
  {
    return entries.size()[0];
  }



  template <typename T, int N>
  inline const T &CollectionBase<T, N>::
                  operator[](const unsigned int index) const
  {
    AssertIndexRange(index, this->size());
    return *entries[index];
  }

} // namespace hp


DEAL_II_NAMESPACE_CLOSE

#endif
