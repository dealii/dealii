// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2013 by the deal.II authors
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

#ifndef __deal2__smartpointer_h
#define __deal2__smartpointer_h


#include <deal.II/base/config.h>
#include <deal.II/base/subscriptor.h>

DEAL_II_NAMESPACE_OPEN

/**
 * Smart pointers avoid destruction of an object in use. They can be used just
 * like a pointer (i.e. using the <tt>*</tt> and <tt>-></tt> operators and through casting)
 * but make sure that the object pointed to is not deleted in the course of
 * use of the pointer by signalling the pointee its use.
 *
 * Objects pointed to, i.e. the class T, should inherit
 * Subscriptor or must implement the same functionality. Null pointers
 * are an exception from this rule and are allowed, too.
 *
 * The second template argument P only serves a single purpose: if a
 * constructor without a debug string is used, then the name of P is
 * used as the debug string.
 *
 * SmartPointer does NOT implement any memory handling! Especially,
 * deleting a SmartPointer does not delete the object. Writing
 * @code
 * SmartPointer<T,P> dont_do_this = new T;
 * @endcode
 * is a sure way to program a memory leak! The secure version is
 * @code
 * T* p = new T;
 * {
 *   SmartPointer<T,P> t(p);
 *   ...
 * }
 * delete p;
 * @endcode
 *
 * Note that a smart pointer can handle <tt>const</tt>ness of an object, i.e.
 * a <tt>SmartPointer<const ABC></tt> really behaves as if it were a pointer to
 * a constant object (disallowing write access when dereferenced), while
 * <tt>SmartPointer<ABC></tt> is a mutable pointer.
 *
 * @ingroup memory
 * @author Guido Kanschat, Wolfgang Bangerth, 1998 - 2009
 */
template<typename T, typename P = void>
class SmartPointer
{
public:
  /**
   * Standard constructor for null
   * pointer. The id of this
   * pointer is set to the name of
   * the class P.
   */
  SmartPointer ();

  /**
   * Copy constructor for
   * SmartPointer. We do now
   * copy the object subscribed to
   * from <tt>tt</tt>, but subscribe
   * ourselves to it again.
   */
  template <class Q>
  SmartPointer (const SmartPointer<T,Q> &tt);

  /**
   * Copy constructor for
   * SmartPointer. We do now
   * copy the object subscribed to
   * from <tt>tt</tt>, but subscribe
   * ourselves to it again.
   */
  SmartPointer (const SmartPointer<T,P> &tt);

  /**
   * Constructor taking a normal
   * pointer.  If possible, i.e. if
   * the pointer is not a null
   * pointer, the constructor
   * subscribes to the given object
   * to lock it, i.e. to prevent
   * its destruction before the end
   * of its use.
   *
   * The <tt>id</tt> is used in the
   * call to
   * Subscriptor::subscribe(id) and
   * by ~SmartPointer() in the call
   * to Subscriptor::unsubscribe().
   */
  SmartPointer (T *t, const char *id);

  /**
   * Constructor taking a normal
   * pointer.  If possible, i.e. if
   * the pointer is not a null
   * pointer, the constructor
   * subscribes to the given object
   * to lock it, i.e. to prevent
   * its destruction before the end
   * of its use. The id of this
   * pointer is set to the name of
   * the class P.
   */
  SmartPointer (T *t);


  /**
   * Destructor, removing the
   * subscription.
   */
  ~SmartPointer();

  /**
   * Assignment operator for normal
   * pointers. The pointer
   * subscribes to the new object
   * automatically and unsubscribes
   * to an old one if it exists. It
   * will not try to subscribe to a
   * null-pointer, but still
   * delete the old subscription.
   */
  SmartPointer<T,P> &operator= (T *tt);

  /**
   * Assignment operator for
   * SmartPointer.  The pointer
   * subscribes to the new object
   * automatically and unsubscribes
   * to an old one if it exists.
   */
  template <class Q>
  SmartPointer<T,P> &operator= (const SmartPointer<T,Q> &tt);

  /**
   * Assignment operator for
   * SmartPointer.  The pointer
   * subscribes to the new object
   * automatically and unsubscribes
   * to an old one if it exists.
   */
  SmartPointer<T,P> &operator= (const SmartPointer<T,P> &tt);

  /**
   * Delete the object pointed to
   * and set the pointer to zero.
   */
  void clear ();

  /**
   * Conversion to normal pointer.
   */
  operator T *() const;

  /**
   * Dereferencing operator. This
   * operator throws an
   * ExcNotInitialized if the
   * pointer is a null pointer.
   */
  T &operator * () const;

  /**
   * Dereferencing operator. This
   * operator throws an
   * ExcNotInitialized if the
   * pointer is a null pointer.
   */
  T *operator -> () const;

  /**
   * Exchange the pointers of this
   * object and the argument. Since
   * both the objects to which is
   * pointed are subscribed to
   * before and after, we do not
   * have to change their
   * subscription counters.
   *
   * Note that this function (with
   * two arguments) and the
   * respective functions where one
   * of the arguments is a pointer
   * and the other one is a C-style
   * pointer are implemented in
   * global namespace.
   */
  template <class Q>
  void swap (SmartPointer<T,Q> &tt);

  /**
   * Swap pointers between this
   * object and the pointer
   * given. As this releases the
   * object pointed to presently,
   * we reduce its subscription
   * count by one, and increase it
   * at the object which we will
   * point to in the future.
   *
   * Note that we indeed need a
   * reference of a pointer, as we
   * want to change the pointer
   * variable which we are given.
   */
  void swap (T *&tt);

  /**
   * Return an estimate of the
   * amount of memory (in bytes)
   * used by this class. Note in
   * particular, that this only
   * includes the amount of memory
   * used by <b>this</b> object, not
   * by the object pointed to.
   */
  std::size_t memory_consumption () const;

private:
  /**
   * Pointer to the object we want
   * to subscribt to. Since it is
   * often necessary to follow this
   * pointer when debugging, we
   * have deliberately chosen a
   * short name.
   */
  T *t;

  /**
   * The identification for the
   * subscriptor.
   */
  const char *const id;
};


/* --------------------------- inline Template functions ------------------------------*/


template <typename T, typename P>
inline
SmartPointer<T,P>::SmartPointer ()
  :
  t (0), id(typeid(P).name())
{}



template <typename T, typename P>
inline
SmartPointer<T,P>::SmartPointer (T *t)
  :
  t (t), id(typeid(P).name())
{
  if (t != 0)
    t->subscribe(id);
}



template <typename T, typename P>
inline
SmartPointer<T,P>::SmartPointer (T *t, const char *id)
  :
  t (t), id(id)
{
  if (t != 0)
    t->subscribe(id);
}



template <typename T, typename P>
template <class Q>
inline
SmartPointer<T,P>::SmartPointer (const SmartPointer<T,Q> &tt)
  :
  t (tt.t), id(tt.id)
{
  if (t != 0)
    t->subscribe(id);
}



template <typename T, typename P>
inline
SmartPointer<T,P>::SmartPointer (const SmartPointer<T,P> &tt)
  :
  t (tt.t), id(tt.id)
{
  if (t != 0)
    t->subscribe(id);
}



template <typename T, typename P>
inline
SmartPointer<T,P>::~SmartPointer ()
{
  if (t != 0)
    t->unsubscribe(id);
}



template <typename T, typename P>
inline
void
SmartPointer<T,P>::clear ()
{
  if (t != 0)
    {
      t->unsubscribe(id);
      delete t;
      t = 0;
    }
}



template <typename T, typename P>
inline
SmartPointer<T,P> &SmartPointer<T,P>::operator = (T *tt)
{
  // optimize if no real action is
  // requested
  if (t == tt)
    return *this;

  if (t != 0)
    t->unsubscribe(id);
  t = tt;
  if (tt != 0)
    tt->subscribe(id);
  return *this;
}



template <typename T, typename P>
template <class Q>
inline
SmartPointer<T,P> &
SmartPointer<T,P>::operator = (const SmartPointer<T,Q> &tt)
{
  // if objects on the left and right
  // hand side of the operator= are
  // the same, then this is a no-op
  if (&tt == this)
    return *this;

  if (t != 0)
    t->unsubscribe(id);
  t = static_cast<T *>(tt);
  if (tt != 0)
    tt->subscribe(id);
  return *this;
}



template <typename T, typename P>
inline
SmartPointer<T,P> &
SmartPointer<T,P>::operator = (const SmartPointer<T,P> &tt)
{
  // if objects on the left and right
  // hand side of the operator= are
  // the same, then this is a no-op
  if (&tt == this)
    return *this;

  if (t != 0)
    t->unsubscribe(id);
  t = static_cast<T *>(tt);
  if (tt != 0)
    tt->subscribe(id);
  return *this;
}



template <typename T, typename P>
inline
SmartPointer<T,P>::operator T *() const
{
  return t;
}



template <typename T, typename P>
inline
T &SmartPointer<T,P>::operator * () const
{
  Assert(t != 0, ExcNotInitialized());
  return *t;
}



template <typename T, typename P>
inline
T *SmartPointer<T,P>::operator -> () const
{
  Assert(t != 0, ExcNotInitialized());
  return t;
}



template <typename T, typename P>
template <class Q>
inline
void SmartPointer<T,P>::swap (SmartPointer<T,Q> &tt)
{
#ifdef DEBUG
  SmartPointer<T,P> aux(t,id);
  *this = tt;
  tt = aux;
#else
  std::swap (t, tt.t);
#endif
}



template <typename T, typename P>
inline
void SmartPointer<T,P>::swap (T *&tt)
{
  if (t != 0)
    t->unsubscribe (id);

  std::swap (t, tt);

  if (t != 0)
    t->subscribe (id);
}



template <typename T, typename P>
inline
std::size_t
SmartPointer<T,P>::memory_consumption () const
{
  return sizeof(SmartPointer<T,P>);
}



// The following function is not strictly necessary but is an optimization
// for places where you call swap(p1,p2) with SmartPointer objects p1, p2.
// Unfortunately, MS Visual Studio (at least up to the 2013 edition) trips
// over it when calling std::swap(v1,v2) where v1,v2 are std::vectors of
// SmartPointer objects: it can't determine whether it should call std::swap
// or dealii::swap on the individual elements (see bug #184 on our Google Code
// site. Consequently, just take this function out of the competition for this
// compiler.
#ifndef _MSC_VER
/**
 * Global function to swap the contents of two smart pointers. As both
 * objects to which the pointers point retain to be subscribed to, we
 * do not have to change their subscription count.
 */
template <typename T, typename P, class Q>
inline
void swap (SmartPointer<T,P> &t1, SmartPointer<T,Q> &t2)
{
  t1.swap (t2);
}
#endif


/**
 * Global function to swap the contents of a smart pointer and a
 * C-style pointer.
 *
 * Note that we indeed need a reference of a pointer, as we want to
 * change the pointer variable which we are given.
 */
template <typename T, typename P>
inline
void swap (SmartPointer<T,P> &t1, T *&t2)
{
  t1.swap (t2);
}



/**
 * Global function to swap the contents of a C-style pointer and a
 * smart pointer.
 *
 * Note that we indeed need a reference of a pointer, as we want to
 * change the pointer variable which we are given.
 */
template <typename T, typename P>
inline
void swap (T *&t1, SmartPointer<T,P> &t2)
{
  t2.swap (t1);
}

DEAL_II_NAMESPACE_CLOSE

#endif
