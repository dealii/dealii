/*----------------------------   smartpointer.h     ---------------------------*/
/*      $Id$                 */
#ifndef __smartpointer_H
#define __smartpointer_H
/*----------------------------   smartpointer.h     ---------------------------*/

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
				      * Assignment operator for
				      * normal pointers. Change of
				      * subscription is necessary.
				      */
    SmartPointer<T> & operator= (T *tt);

				     /**
				      *Assignment operator for
				      * #SmartPointer#. Change of
				      * subscription is necessary.
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
SmartPointer<T>::~SmartPointer () {
  if (t)
    t->unsubscribe();
};



template <typename T>
SmartPointer<T> & SmartPointer<T>::operator = (T *tt) {
  if (t)
    t->unsubscribe();
  t = tt;
  if (tt)
    tt->subscribe();
  return *this;
};


template <typename T>
SmartPointer<T> & SmartPointer<T>::operator = (const SmartPointer<T>& tt) {
  if (t)
    t->unsubscribe();
  t = (T*) tt;
  if (tt)
    tt->subscribe();
  return *this;
};


template <typename T>
inline
SmartPointer<T>::operator T* () const {
  return t;
};



template <typename T>
inline
T & SmartPointer<T>::operator * () const {
  return *t;
};



template <typename T>
inline
T * SmartPointer<T>::operator -> () const {
  return t;
};





/*----------------------------   smartpointer.h     ---------------------------*/
/* end of #ifndef __smartpointer_H */
#endif
/*----------------------------   smartpointer.h     ---------------------------*/
