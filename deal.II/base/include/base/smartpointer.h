//----------------------------  smartpointer.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  smartpointer.h  ---------------------------
#ifndef __deal2__smartpointer_h
#define __deal2__smartpointer_h


#ifndef __subscriptor_H
#include <base/subscriptor.h>
#endif

/**
 * Smart pointers avoid destruction of an object in use. They can be used just
 * like a pointer (i.e. using the #*# and #-># operators and through casting)
 * but make sure that the object pointed to is not deleted in the course of
 * use of the pointer by signalling the pointee its use.
 *
 * Objects pointed to should inherit #Subscriptor# or must implement
 * the same functionality. Null pointers are an exception from this
 * rule and are allowed, too.
 *
 * #SmartPointer# does NOT implement any memory handling! Especially,
 * deleting a #SmartPointer# does not delete the object. Writing
 * \begin{verbatim}
 * SmartPointer<T> t = new T;
 * \end{verbatim}
 * is a sure way to program a memory leak! The secure version is
 * \begin{verbatim}
 * T* p = new T;
 * {
 *   SmartPointer<T> t = p;
 *   ...
 * }
 * delete p;
 * \end{verbatim}
 *
 * Note that a smart pointer can handle #const#ness of an object, i.e.
 * a #SmartPointer<const ABC># really behaves as if it were a pointer to
 * a constant object (disallowing write access when dereferenced), while
 * #SmartPointer<ABC># is a mutable pointer.
 */
template<typename T>
class SmartPointer
{
  public:
				     /**
				      * Standard constructor for null pointer.
				      */
    SmartPointer();

				     /*
				      * Copy constructor for
				      * #SmartPointer#. We copy the object
				      * subscribed to from #tt#, but subscribe
				      * ourselves to it again.
				      */
    SmartPointer (const SmartPointer<T> &tt);

				     /**
				      * Constructor taking a normal pointer.
				      * If possible, i.e. if the pointer
				      * is not a null pointer, the constructor
				      * subscribes to the given object to
				      * lock it, i.e. to prevent its
				      * destruction before the end of its use.
				      */
    SmartPointer (T *t);


/**
				      * Destructor, removing the subscription.
				      */
    ~SmartPointer();
    
				     /**
				      * Assignment operator for normal
				      * pointers. The pointer
				      * subscribes to the new object
				      * automatically and unsubscribes
				      * to an old one if it exists.
				      */
    SmartPointer<T> & operator= (T *tt);

				     /**
				      *Assignment operator for
				      * #SmartPointer#.  The pointer
				      * subscribes to the new object
				      * automatically and unsubscribes
				      * to an old one if it exists.
				      */
    SmartPointer<T> & operator= (const SmartPointer<T>& tt);

				     /**
				      * Conversion to normal pointer.
				      */
    operator T* () const;
    
				     /**
				      * Dereferencing operator.
				      */
    T& operator * () const;
    
				     /**
				      * Dereferencing operator.
				      */
    T * operator -> () const;
    
  private:
				     /**
				      * Pointer to the object we want
				      * to subscribt to.
				      */
  T * t;
};


/* --------------------------- inline Template functions ------------------------------*/


template <typename T>
SmartPointer<T>::SmartPointer () :
		t (0)
{};


template <typename T>
SmartPointer<T>::SmartPointer (T *t) :
		t (t)
{
  if (t)
    t->subscribe();
};


template <typename T>
SmartPointer<T>::SmartPointer (const SmartPointer<T> &tt) :
		t (tt.t)
{
  if (t)
    t->subscribe();
};


template <typename T>
SmartPointer<T>::~SmartPointer ()
{
  if (t)
    t->unsubscribe();
};


template <typename T>
SmartPointer<T> & SmartPointer<T>::operator = (T *tt)
{
				   // optimize if no real action is
				   // requested
  if (t == tt)
    return *this;
  
  if (t)
    t->unsubscribe();
  t = tt;
  if (tt)
    tt->subscribe();
  return *this;
};


template <typename T>
SmartPointer<T> & SmartPointer<T>::operator = (const SmartPointer<T>& tt)
{
				   // if objects on the left and right
				   // hand side of the operator= are
				   // the same, then this is a no-op
  if (&tt == this)
    return *this;
  
  if (t)
    t->unsubscribe();
  t = static_cast<T*>(tt);
  if (tt)
    tt->subscribe();
  return *this;
};


template <typename T>
inline
SmartPointer<T>::operator T* () const
{
  return t;
};


template <typename T>
inline
T & SmartPointer<T>::operator * () const
{
  return *t;
};


template <typename T>
inline
T * SmartPointer<T>::operator -> () const
{
  return t;
};


#endif
