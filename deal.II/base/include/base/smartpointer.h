//----------------------------  smartpointer.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  smartpointer.h  ---------------------------
#ifndef __deal2__smartpointer_h
#define __deal2__smartpointer_h


#include <base/config.h>
#include <base/subscriptor.h>


/**
 * Smart pointers avoid destruction of an object in use. They can be used just
 * like a pointer (i.e. using the @p{*} and @p{->} operators and through casting)
 * but make sure that the object pointed to is not deleted in the course of
 * use of the pointer by signalling the pointee its use.
 *
 * Objects pointed to should inherit @p{Subscriptor} or must implement
 * the same functionality. Null pointers are an exception from this
 * rule and are allowed, too.
 *
 * @p{SmartPointer} does NOT implement any memory handling! Especially,
 * deleting a @p{SmartPointer} does not delete the object. Writing
 * @begin{verbatim}
 * SmartPointer<T> t = new T;
 * @end{verbatim}
 * is a sure way to program a memory leak! The secure version is
 * @begin{verbatim}
 * T* p = new T;
 * {
 *   SmartPointer<T> t = p;
 *   ...
 * }
 * delete p;
 * @end{verbatim}
 *
 * Note that a smart pointer can handle @p{const}ness of an object, i.e.
 * a @p{SmartPointer<const ABC>} really behaves as if it were a pointer to
 * a constant object (disallowing write access when dereferenced), while
 * @p{SmartPointer<ABC>} is a mutable pointer.
 *
 * @author Guido Kanschat, Wolfgang Bangerth, 1998, 1999, 2000
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
				      * @p{SmartPointer}. We copy the object
				      * subscribed to from @p{tt}, but subscribe
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
				      * to an old one if it exists. It
				      * will not try to subscribe to a
				      * null-pointer, but stilll
				      * delete the old subscription.
				      */
    SmartPointer<T> & operator= (T *tt);

				     /**
				      *Assignment operator for
				      * @p{SmartPointer}.  The pointer
				      * subscribes to the new object
				      * automatically and unsubscribes
				      * to an old one if it exists.
				      */
    SmartPointer<T> & operator= (const SmartPointer<T> &tt);

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
    void swap (SmartPointer<T> &tt);

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
				      * used by @em{this} object, not
				      * by the object pointed to.
				      */
    unsigned int memory_consumption () const;
    
  private:
				     /**
				      * Pointer to the object we want
				      * to subscribt to. Since it is
				      * often necessary to follow this
				      * pointer when debugging, we
				      * have deliberately chosen a
				      * short name.
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
  if (t != 0)
    t->subscribe();
};



template <typename T>
SmartPointer<T>::SmartPointer (const SmartPointer<T> &tt) :
		t (tt.t)
{
  if (t != 0)
    t->subscribe();
};



template <typename T>
SmartPointer<T>::~SmartPointer ()
{
  if (t != 0)
    t->unsubscribe();
};



template <typename T>
SmartPointer<T> & SmartPointer<T>::operator = (T *tt)
{
				   // optimize if no real action is
				   // requested
  if (t == tt)
    return *this;
  
  if (t != 0)
    t->unsubscribe();
  t = tt;
  if (tt != 0)
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
  
  if (t != 0)
    t->unsubscribe();
  t = static_cast<T*>(tt);
  if (tt != 0)
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



template <typename T>
inline
void SmartPointer<T>::swap (SmartPointer<T> &tt)
{
  swap (t, tt.t);
};



template <typename T>
inline
void SmartPointer<T>::swap (T *&tt)
{
  if (t != 0)
    t->unsubscribe ();
  
  std::swap (t, tt);

  if (t != 0)
    t->subscribe ();
};



template <typename T>
inline
unsigned int SmartPointer<T>::memory_consumption () const
{
  return sizeof(SmartPointer<T>);
};




/**
 * Global function to swap the contents of two smart pointers. As both
 * objects to which the pointers point retain to be subscribed to, we
 * do not have to change their subscription count.
 */
template <typename T>
inline
void swap (SmartPointer<T> &t1, SmartPointer<T> &t2)
{
  t1.swap (t2);
};



/**
 * Global function to swap the contents of a smart pointer and a
 * C-style pointer.
 *
 * Note that we indeed need a reference of a pointer, as we want to
 * change the pointer variable which we are given.
 */
template <typename T>
inline
void swap (SmartPointer<T> &t1, T *&t2)
{
  t1.swap (t2);
};



/**
 * Global function to swap the contents of a C-style pointer and a
 * smart pointer.
 *
 * Note that we indeed need a reference of a pointer, as we want to
 * change the pointer variable which we are given.
 */
template <typename T>
inline
void swap (T *&t1, SmartPointer<T> &t2)
{
  t2.swap (t1);
};


#endif
