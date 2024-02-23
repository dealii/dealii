// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_mutex_h
#define dealii_mutex_h


#include <deal.II/base/config.h>

#include <mutex>

DEAL_II_NAMESPACE_OPEN

/**
 * @addtogroup threads
 * @{
 */

namespace Threads
{
  /**
   * A class implementing a <a
   * href="https://en.wikipedia.org/wiki/Lock_(computer_science)">mutex</a>.
   * Mutexes are used to lock data structures to ensure that only a
   * single thread of execution can access them at the same time.
   *
   * This class is a thin wrapper around `std::mutex`. The only difference
   * is that this class is copyable when `std::mutex` is not.  Indeed, when
   * copied, the receiving object does not copy any state from the object
   * being copied, i.e. an entirely new mutex is created. These semantics
   * are consistent with the common use case if a mutex is used as a member
   * variable to lock the other member variables of a class: in that case,
   * the mutex of the copied-to object should only guard the members of the
   * copied-to object, not the members of both the copied-to and
   * copied-from object. Since at the time when the class is copied, the
   * destination's member variable is not used yet, its corresponding mutex
   * should also remain in its original state.
   */
  class Mutex : public std::mutex
  {
  public:
    /**
     * Default constructor.
     */
    Mutex() = default;

    /**
     * Copy constructor. As discussed in this class's documentation, no state
     * is copied from the object given as argument.
     */
    Mutex(const Mutex &)
      : std::mutex()
    {}

    /**
     * Copy operators. As discussed in this class's documentation, no state
     * is copied from the object given as argument.
     */
    Mutex &
    operator=(const Mutex &)
    {
      return *this;
    }
  };
} // namespace Threads

/**
 * @}
 */

DEAL_II_NAMESPACE_CLOSE
#endif
