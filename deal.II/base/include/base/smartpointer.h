/*----------------------------   smartpointer.h     ---------------------------*/
/*      $Id$                 */
#ifndef __smartpointer_H
#define __smartpointer_H
/*----------------------------   smartpointer.h     ---------------------------*/

#ifndef #ifndef __subscriptor_H
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
 * <PRE>
 * SmartPointer<T> t = new T;
 * </PRE>
 * is a sure way to program a memory leak! The secure version is
 * <PRE>
 * T* p = new T;
 * {
 *   SmartPointer<T> t = p;
 *   ...
 * }
 * delete p;
 * </PRE> */
template<class T>
class SmartPointer
{
    T* t;

  public:
				     /**
				      * Constructor taking a normal pointer.  */
    SmartPointer(T* tt) :
		    t(tt) 
      {
	t->subscribe();
      }

				   /**
				    * Standard constructor for null pointer.
				    */
  SmartPointer() :
		  t(0)
      {}

				     /**
				      * Destructor, removing the subscription.
				      */
    ~SmartPointer()
      {
	if (t)
	  t->unsubscribe();
      }
				   /**
				    * Assignment operator. Change of
				    * subscription is necessary.
				    */
  SmartPointer<T>& operator=(T* tt)
      {
	if (t)
	  t->unsubscribe();
	t = tt;
	if (tt)
	  tt->subscribe();
	return *this;
      }
  

				     /**
				      * Conversion to normal pointer.
				      */
    operator T* () const
      {
	return t;
      }

				     /**
				      * Dereferencing operator.
				      */
    T& operator* () const
      {
	return *t;
      }

				     /**
				      * Dereferencing operator.
				      */
    T* operator -> () const
      {
	return t;
      }
};



/*----------------------------   smartpointer.h     ---------------------------*/
/* end of #ifndef __smartpointer_H */
#endif
/*----------------------------   smartpointer.h     ---------------------------*/
