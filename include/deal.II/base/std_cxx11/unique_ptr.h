// ---------------------------------------------------------------------
//
// Copyright (C) 2015, 2016 by the deal.II authors
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

#ifndef dealii_std_cxx11_unique_ptr_h
#define dealii_std_cxx11_unique_ptr_h


#include <deal.II/base/config.h>


#  include <memory>
DEAL_II_NAMESPACE_OPEN
namespace std_cxx11
{
  using std::unique_ptr;
}
DEAL_II_NAMESPACE_CLOSE

// then allow using the old namespace name instead of the new one
DEAL_II_NAMESPACE_OPEN
namespace std_cxx1x = std_cxx11;
DEAL_II_NAMESPACE_CLOSE


DEAL_II_NAMESPACE_OPEN

/**
 * Convert an object of type `std::unique_ptr<From>` to an object of
 * type `std::unique_ptr<To>`, where it is assumed that we can cast
 * the pointer to `From` to a pointer to `To` using a `dynamic_cast`
 * -- in other words, we assume that `From` and `To` are connected
 * through a class hierarchy, and that the object pointed to is in
 * fact of a type that contains both a `From` and a `To`. An example
 * is if either `To` is derived from `From` or the other way around.
 *
 * The function throws an exception of type `std::bad_cast` if the
 * `dynamic_cast` does not succeed. This is the same exception you
 * would get if a regular `dynamic_cast` between object types (but not
 * pointer types) does not succeed.
 *
 * An example of how this function works is as follows:
 * @code
 *   class B { ... };            // A base class. Assume that it has virtual
 *                               // functions so that dynamic_cast can work.
 *   class D : public B { ... }; // A derived class
 *
 *
 *   std::unique_ptr<B> create_object (...) {...}  // A factory function
 *
 *   void foo (...)
 *   {
 *     std::unique_ptr<B> b = create_object (...);
 *
 *     // Assume that we know for some reason that the object above must
 *     // have created a D object but returned it as a std::unique_ptr<B>.
 *     // In order to access the D functionality, we need to cast the
 *     // pointer. Use the equivalent to dynamic_cast:
 *     std::unique_ptr<D> d = dynamic_unique_cast<D>(std::move(b));
 *
 *     // If the object really was a D, then 'd' now points to it. Note
 *     // also that in accordance with the semantics of std::unique_ptr,
 *     // it was necessary to std::move the 'b' object, and indeed 'b'
 *     // now no longer points to anything -- ownership has been
 *     // transferred to 'd'!
 * @endcode
 *
 * @note This function does not try to convert the `Deleter` objects stored
 *   by `std::unique_ptr` objects. The function therefore only works if the
 *   deleter objects are at their defaults, i.e., if they are of type
 *   `std::default_delete<To>` and `std::default_delete<From>`.
 */
template <typename To, typename From>
std::unique_ptr<To>
dynamic_unique_cast (std::unique_ptr<From> &&p)
{
  // Let's see if we can cast from 'From' to 'To'. If so, do the cast,
  // and then release the pointer from the old
  // owner
  if (To *cast = dynamic_cast<To *>(p.get()))
    {
      std::unique_ptr<To> result(cast);
      p.release();
      return result;
    }
  else
    throw std::bad_cast();
}

DEAL_II_NAMESPACE_CLOSE

#endif
