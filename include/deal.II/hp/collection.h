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
  template <typename T>
  class Collection : public Subscriptor
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
  inline const T &Collection<T>::operator[](const unsigned int index) const
  {
    AssertIndexRange(index, entries.size());
    return *entries[index];
  }

} // namespace hp


DEAL_II_NAMESPACE_CLOSE

#endif
