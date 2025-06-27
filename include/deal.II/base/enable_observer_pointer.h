// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1998 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_enable_observer_pointer_h
#define dealii_enable_observer_pointer_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

#include <atomic>
#include <map>
#include <mutex>
#include <ostream>
#include <string>
#include <typeinfo>
#include <vector>

DEAL_II_NAMESPACE_OPEN


// Forward declaration
template <typename T, typename P>
class ObserverPointer;


/**
 * This class supports the functioning of the ObserverPointer class.
 *
 * Used as a base class, this class allows the ObserverPointer class to register
 * with an object that it is observing that object. As a consequence, the object
 * can report an error when its destructor is called if there are still
 * ObserverPointers pointing to it, since destroying the object would lead to
 * dangling pointers that point to memory that's no longer valid.
 *
 * In actual practice, the mechanism works as follows: Assume you have an
 * `ObserverPointer<X>` object that points objects of type `X` (where `X` is
 * derived from EnableObserverPointer). When you write
 * @code
 *   ObserverPointer<X> ptr;
 *   X                  object;
 *   ptr = &x;
 * @endcode
 * then `ObserverPointer::operator=()` called in the last line of the code above
 * will call the EnableObserverPointer::subscribe() function (of the class
 * documented here) with arguments that allow both
 * sides to track each other. If `object` goes out of scope while `ptr` is still
 * pointing to it, this would create a dangling pointer, and the destructor of
 * `X` (i.e., actually the destructor of the EnableObserverPointer base class of
 * `X`) will abort the program with an error.
 *
 * In practice, when such an error happens, it is often difficult to tell
 * *which* pointer is still pointing to the `object` that is being destroyed. To
 * this end, ObserverPointer can also pass a string to the subscribe() member
 * function that identifies which observer is still alive. This string is passed
 * to the ObserverPointer object upon construction and is user-defined so that
 * code developed can give names to their observer pointers.
 *
 * @ingroup memory
 */
class EnableObserverPointer
{
public:
  /**
   * Constructor setting the counter to zero.
   */
  EnableObserverPointer();

  /**
   * Copy-constructor.
   *
   * The counter of the copy is zero, since references point to the original
   * object.
   */
  EnableObserverPointer(const EnableObserverPointer &);

  /**
   * Move constructor.
   *
   * An object inheriting from EnableObserverPointer can only be
   * moved if no other objects are subscribing to it.
   */
  EnableObserverPointer(EnableObserverPointer &&) noexcept;

  /**
   * Destructor, asserting that the counter is zero.
   */
  virtual ~EnableObserverPointer();

  /**
   * Assignment operator.
   *
   * This has to be handled with care, too, because the counter has to remain
   * the same. It therefore does nothing more than returning <tt>*this</tt>.
   */
  EnableObserverPointer &
  operator=(const EnableObserverPointer &);

  /**
   * Move assignment operator. Only invalidates the object moved from.
   */
  EnableObserverPointer &
  operator=(EnableObserverPointer &&) noexcept;

  /**
   * @name Querying the observer pointers an object has.
   *
   * @{
   */

  /**
   * Return the present number of subscriptions to this object. This allows to
   * use this class for reference counted lifetime determination where the
   * last one to unsubscribe also deletes the object.
   */
  unsigned int
  n_subscriptions() const;

  /**
   * List the subscribers to the input @p stream.
   */
  template <typename StreamType>
  void
  list_subscribers(StreamType &stream) const;

  /**
   * List the subscribers to @p deallog.
   */
  void
  list_subscribers() const;

  /** @} */

  /**
   * @addtogroup Exceptions
   * @{
   */

  /**
   * Exception: Object may not be deleted, since it is used.
   */
  DeclException3(ExcInUse,
                 int,
                 std::string,
                 std::string,
                 << "Object of class " << arg2 << " is still used by " << arg1
                 << " other objects."
                 << "\n\n"
                 << "(Additional information: " << arg3 << ")\n\n"
                 << "See the entry in the Frequently Asked Questions of "
                 << "deal.II (linked to from http://www.dealii.org/) for "
                 << "a lot more information on what this error means and "
                 << "how to fix programs in which it happens.");

  /**
   * A subscriber with the identification string given to
   * EnableObserverPointer::unsubscribe() did not subscribe to the
   * object.
   */
  DeclException2(ExcNoSubscriber,
                 std::string,
                 std::string,
                 << "No subscriber with identifier <" << arg2
                 << "> subscribes to this object of class " << arg1
                 << ". Consequently, it cannot be unsubscribed.");
  /** @} */

  /**
   * Read or write the data of this object to or from a stream for the purpose
   * of serialization using the [BOOST serialization
   * library](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html).
   *
   * This function does not actually serialize any of the member variables of
   * this class. The reason is that what this class stores is only who
   * subscribes to this object, but who does so at the time of storing the
   * contents of this object does not necessarily have anything to do with who
   * subscribes to the object when it is restored. Consequently, we do not
   * want to overwrite the subscribers at the time of restoring, and then
   * there is no reason to write the subscribers out in the first place.
   */
  template <class Archive>
  void
  serialize(Archive &ar, const unsigned int version);

private:
  /**
   * Store the number of objects which subscribed to this object. Initially,
   * this number is zero, and upon destruction it shall be zero again (i.e.
   * all objects which subscribed should have unsubscribed again).
   *
   * The creator (and owner) of an object is counted in the map below if HE
   * manages to supply identification.
   *
   * We use the <tt>mutable</tt> keyword in order to allow subscription to
   * constant objects also.
   *
   * This counter may be read from and written to concurrently in
   * multithreaded code: hence we use the <code>std::atomic</code> class
   * template.
   */
  mutable std::atomic<unsigned int> counter;

  /**
   * In this map, we count subscriptions for each different identification
   * string supplied to subscribe().
   */
  mutable std::map<std::string, unsigned int> counter_map;

  /**
   * The data type used in #counter_map.
   */
  using map_value_type = decltype(counter_map)::value_type;

  /**
   * The iterator type used in #counter_map.
   */
  using map_iterator = decltype(counter_map)::iterator;

  /**
   * In this vector, we store pointers to the validity bool in the
   * ObserverPointer objects that subscribe to this class.
   */
  mutable std::vector<std::atomic<bool> *> validity_pointers;

  /**
   * Pointer to the typeinfo object of this object, from which we can later
   * deduce the class name. Since this information on the derived class is
   * neither available in the destructor, nor in the constructor, we obtain it
   * in between and store it here.
   */
  mutable const std::type_info *object_info;

  /**
   * A mutex used to ensure data consistency when accessing the `mutable`
   * members of this class. This lock is used in the subscribe() and
   * unsubscribe() functions, as well as in `list_subscribers()`.
   */
  static std::mutex mutex;

  /**
   * @name EnableObserverPointer functionality
   *
   * Classes derived from EnableObserverPointer provide a facility
   * to subscribe to this object. This is mostly used by the ObserverPointer
   * class.
   * @{
   */

  /**
   * Subscribes a user of the object by storing the pointer @p validity. The
   * subscriber may be identified by text supplied as @p identifier.
   */
  void
  subscribe(std::atomic<bool> *const validity,
            const std::string       &identifier = "") const;

  /**
   * Unsubscribes a user from the object.
   *
   * @note The @p identifier and the @p validity pointer must be the same as
   * the one supplied to subscribe().
   */
  void
  unsubscribe(std::atomic<bool> *const validity,
              const std::string       &identifier = "") const;

  /**
   * Check that there are no objects subscribing to this object. If this check
   * passes then it is safe to destroy the current object. It this check fails
   * then this function will either abort or print an error message to deallog
   * (by using the AssertNothrow mechanism), but will not throw an exception.
   *
   * @note Since this function is just a consistency check it does nothing in
   * release mode.
   *
   * @note If this function is called when there is an uncaught exception
   * then, rather than aborting, this function prints an error message to the
   * standard error stream and returns.
   */
  void
  check_no_subscribers() const noexcept;

  template <typename, typename>
  friend class ObserverPointer;

  /** @} */
};


/**
 * A type alias for the EnableObserverPointer class that makes sure
 * the previous name of the class, Subscriptor, continues to be available.
 *
 * @deprecated Use the new name of the class, ObserverPointer, instead.
 */
using Subscriptor DEAL_II_DEPRECATED_EARLY_WITH_COMMENT(
  "Use the new name of the class, EnableObserverPointer.") =
  EnableObserverPointer;


//---------------------------------------------------------------------------

inline EnableObserverPointer::EnableObserverPointer()
  : counter(0)
  , object_info(nullptr)
{}



inline EnableObserverPointer::EnableObserverPointer(
  const EnableObserverPointer &)
  : counter(0)
  , object_info(nullptr)
{}



inline EnableObserverPointer &
EnableObserverPointer::operator=(const EnableObserverPointer &s)
{
  object_info = s.object_info;
  return *this;
}



inline unsigned int
EnableObserverPointer::n_subscriptions() const
{
  return counter;
}



template <class Archive>
inline void
EnableObserverPointer::serialize(Archive &, const unsigned int)
{
  // do nothing, as explained in the
  // documentation of this function
}

template <typename StreamType>
inline void
EnableObserverPointer::list_subscribers(StreamType &stream) const
{
  std::lock_guard<std::mutex> lock(mutex);

  for (const auto &it : counter_map)
    stream << it.second << '/' << counter << " subscriptions from \""
           << it.first << '\"' << std::endl;
}

DEAL_II_NAMESPACE_CLOSE

#endif
