// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2022 by the deal.II authors
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

#ifndef dealii_subscriptor_h
#define dealii_subscriptor_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

#include <atomic>
#include <map>
#include <mutex>
#include <string>
#include <typeinfo>
#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * Handling of subscriptions.
 *
 * This class, as a base class, allows to keep track of other objects using a
 * specific object. It is used to avoid that pointers that point to an object of
 * a class derived from Subscriptor are referenced after that object has been
 * invalidated. Here, invalidation is assumed to happen when the object is
 * moved from or destroyed.
 * The mechanism works as follows: The member function subscribe() accepts a
 * pointer to a boolean that is modified on invalidation. The object that owns
 * this pointer (usually an object of class type SmartPointer) is then expected
 * to check the state of the boolean before trying to access this class.
 *
 * The utility of this class is even enhanced by providing identifying strings
 * to the functions subscribe() and unsubscribe(). These strings are represented
 * as <code>const char</code> pointers since the underlying buffer comes from
 * (and is managed by) the run-time type information system: more exactly, these
 * pointers are the result the function call <code>typeid(x).name()</code> where
 * <code>x</code> is some object. Therefore, the pointers provided to
 * subscribe() and to unsubscribe() must be the same. Strings with equal
 * contents will not be recognized to be the same. The handling in
 * SmartPointer will take care of this.
 * The current subscribers to this class can be obtained by calling
 * list_subscribers().
 *
 * @ingroup memory
 */
class Subscriptor
{
public:
  /**
   * Constructor setting the counter to zero.
   */
  Subscriptor();

  /**
   * Copy-constructor.
   *
   * The counter of the copy is zero, since references point to the original
   * object.
   */
  Subscriptor(const Subscriptor &);

  /**
   * Move constructor.
   *
   * An object inheriting from Subscriptor can only be moved if no other
   * objects are subscribing to it.
   */
  Subscriptor(Subscriptor &&) noexcept;

  /**
   * Destructor, asserting that the counter is zero.
   */
  virtual ~Subscriptor();

  /**
   * Assignment operator.
   *
   * This has to be handled with care, too, because the counter has to remain
   * the same. It therefore does nothing more than returning <tt>*this</tt>.
   */
  Subscriptor &
  operator=(const Subscriptor &);

  /**
   * Move assignment operator. Only invalidates the object moved from.
   */
  Subscriptor &
  operator=(Subscriptor &&) noexcept;

  /**
   * @name Subscriptor functionality
   *
   * Classes derived from Subscriptor provide a facility to subscribe to this
   * object. This is mostly used by the SmartPointer class.
   */
  // @{

  /**
   * Subscribes a user of the object by storing the pointer @p validity. The
   * subscriber may be identified by text supplied as @p identifier.
   */
  void
  subscribe(std::atomic<bool> *const validity,
            const std::string &      identifier = "") const;

  /**
   * Unsubscribes a user from the object.
   *
   * @note The @p identifier and the @p validity pointer must be the same as
   * the one supplied to subscribe().
   */
  void
  unsubscribe(std::atomic<bool> *const validity,
              const std::string &      identifier = "") const;

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

  // @}

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
   * Subscriptor::unsubscribe() did not subscribe to the object.
   */
  DeclException2(ExcNoSubscriber,
                 std::string,
                 std::string,
                 << "No subscriber with identifier <" << arg2
                 << "> subscribes to this object of class " << arg1
                 << ". Consequently, it cannot be unsubscribed.");
  //@}

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
   * In this vector, we store pointers to the validity bool in the SmartPointer
   * objects that subscribe to this class.
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

  /**
   * A mutex used to ensure data consistency when printing out the list
   * of subscribers.
   */
  static std::mutex mutex;
};

//---------------------------------------------------------------------------

inline Subscriptor::Subscriptor()
  : counter(0)
  , object_info(nullptr)
{}



inline Subscriptor::Subscriptor(const Subscriptor &)
  : counter(0)
  , object_info(nullptr)
{}



inline Subscriptor &
Subscriptor::operator=(const Subscriptor &s)
{
  object_info = s.object_info;
  return *this;
}



inline unsigned int
Subscriptor::n_subscriptions() const
{
  return counter;
}



template <class Archive>
inline void
Subscriptor::serialize(Archive &, const unsigned int)
{
  // do nothing, as explained in the
  // documentation of this function
}

template <typename StreamType>
inline void
Subscriptor::list_subscribers(StreamType &stream) const
{
  std::lock_guard<std::mutex> lock(mutex);

  for (const auto &it : counter_map)
    stream << it.second << '/' << counter << " subscriptions from \""
           << it.first << '\"' << std::endl;
}

DEAL_II_NAMESPACE_CLOSE

#endif
