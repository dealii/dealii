// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2018 by the deal.II authors
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

#ifndef dealii_subscriptor_h
#define dealii_subscriptor_h


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>

#include <atomic>
#include <map>
#include <string>
#include <typeinfo>

DEAL_II_NAMESPACE_OPEN

/**
 * Handling of subscriptions.
 *
 * This class, as a base class, allows to keep track of other objects using a
 * specific object. It is used, when an object, given to a constructor by
 * reference, is stored. Then, the original object may not be deleted before
 * the dependent object is deleted. You can assert this constraint by letting
 * the object passed be derived from this class and let the user subscribe()
 * to this object. The destructor the used object inherits from the
 * Subscriptor class then will lead to an error when destruction is attempted
 * while there are still subscriptions.
 *
 * The utility of this class is even enhanced by providing identifying strings
 * to the functions subscribe() and unsubscribe(). In case of a hanging
 * subscription during destruction, this string will be listed in the
 * exception's message. These strings are represented as <code>const char
 * *</code> pointers since the underlying buffer comes from (and is managed
 * by) the run-time type information system: more exactly, these pointers are
 * the result the function call <code>typeid(x).name()</code> where
 * <code>x</code> is some object. Therefore, the pointers provided to
 * subscribe() and to unsubscribe() must be the same. Strings with equal
 * contents will not be recognized to be the same. The handling in
 * SmartPointer will take care of this.
 *
 * @note This feature is switched off if multithreading is used (i.e., if
 * <code>DEAL_II_WITH_THREADS</code> is on).
 *
 * @ingroup memory
 * @author Guido Kanschat, 1998 - 2005
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
  virtual
  ~Subscriptor();

  /**
   * Assignment operator.
   *
   * This has to be handled with care, too, because the counter has to remain
   * the same. It therefore does nothing more than returning <tt>*this</tt>.
   */
  Subscriptor &
  operator = (const Subscriptor &);

  /**
   * Move assignment operator.
   *
   * Asserts that the counter for the moved object is zero.
   */
  Subscriptor &
  operator = (Subscriptor &&) noexcept;

  /**
   * Subscribes a user of the object. The subscriber may be identified by text
   * supplied as <tt>identifier</tt>.
   */
  void
  subscribe (const char *identifier = nullptr) const;

  /**
   * Unsubscribes a user from the object.
   *
   * @note The <tt>identifier</tt> must be the <b>same pointer</b> as the one
   * supplied to subscribe(), not just the same text.
   */
  void
  unsubscribe (const char *identifier = nullptr) const;

  /**
   * Return the present number of subscriptions to this object. This allows to
   * use this class for reference counted lifetime determination where the
   * last one to unsubscribe also deletes the object.
   */
  unsigned int
  n_subscriptions () const;

  /**
   * List the subscribers to @p deallog.
   */
  void
  list_subscribers () const;

  /**
   * @addtogroup Exceptions
   * @{
   */

  /**
   * Exception: Object may not be deleted, since it is used.
   */
  DeclException3(ExcInUse,
                 int, std::string, std::string,
                 << "Object of class " << arg2
                 << " is still used by " << arg1 << " other objects."
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
  DeclException2(ExcNoSubscriber, std::string, std::string,
                 << "No subscriber with identifier <" << arg2
                 << "> subscribes to this object of class " << arg1
                 << ". Consequently, it cannot be unsubscribed.");
  //@}

  /**
   * Read or write the data of this object to or from a stream for the purpose
   * of serialization.
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
   * The data type used in #counter_map.
   */
  typedef std::map<const char *, unsigned int>::value_type
  map_value_type;

  /**
   * The iterator type used in #counter_map.
   */
  typedef std::map<const char *, unsigned int>::iterator
  map_iterator;

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
   *
   */
  mutable std::atomic<unsigned int> counter;

  /**
   * In this map, we count subscriptions for each different identification
   * string supplied to subscribe().
   */
  mutable std::map<const char *, unsigned int> counter_map;

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
  check_no_subscribers () const noexcept;
};

//---------------------------------------------------------------------------

template <class Archive>
inline
void
Subscriptor::serialize(Archive &,
                       const unsigned int)
{
  // do nothing, as explained in the
  // documentation of this function
}

DEAL_II_NAMESPACE_CLOSE

#endif
