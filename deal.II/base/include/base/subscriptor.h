/*----------------------------   subscriptor.h     ---------------------------*/
/*      $Id$                 */
#ifndef __subscriptor_H
#define __subscriptor_H
/*----------------------------   subscriptor.h     ---------------------------*/

/**
 * Handling of subscriptions.
 *
 * This class, as a base class, allows to keep track of other objects
 * using a specific object. It is used, when an object, given to a
 * constructor by reference, is stored. Then, the original object may
 * not be deleted before the dependent object is deleted.
 */
class Subscriptor
{
  unsigned counter;
public:
				   /// Object may not be deleted, since it is used.
  DeclException0(InUse);//, unsigned, "Object is used " << counter << " times");
  DeclException0(NotUsed);
  
				   /// Constructor setting the counter to zero.
  Subscriptor()
		  : counter(0)
      {}
				   /// Destructor, asserting that the counter is zero.
  ~Subscriptor()
      {
	Assert(counter==0, InUse());
      }
				   /**
				    * Copy-constructor.
				    *
				    * The counter of the copy is zero,
				    * since references point to the
				    * original object. */

  Subscriptor(const Subscriptor&)
		  : counter(0)
      {}

				   /**
				    * Assignment operator.
				    *
				    * This has to be handled with
				    * care, too, because the counter
				    * has to remain the same. */
  Subscriptor& operator=(const Subscriptor&)
      {
	return *this;
      }

				   /// Subscribes a user of the object.
  void subscribe() //mutable
      {
	++counter;
      }
  
				   /// Unsubscribes a user from the object.
  void unsubscribe() //mutable
      {
	Assert(counter>0, NotUsed());
	--counter;
      }
};

/**
 * Smart references avoid destruction of a referenced object.
 */
template<class T>
class SmartReference
{
  T& t;
public:
				   /// Constructor taking a normal reference.
  SmartReference(const T& tt)
		  : t(/*const_cast<T&>*/ tt) 
      {
	t.subscribe();
      }
				   /// Destructor, removing the subscription.
      ~SmartReference()
      {
	t.unsubscribe();
      }
  
				   /// Conversion to normal reference
      operator T& ()
      {
	return t;
      }
				   /// Conversion to normal const reference.
      operator const T& () const
      {
	return t;
      }
};

/**
 * Smart pointers avoid destruction of an object in use.
 */
template<class T>
class SmartPointer
{
  T* t;
public:
				   /// Constructor taking a normal pointer.
  SmartPointer(const T* tt)
		  : t(/*const_cast<T&>*/ tt) 
      {
	t->subscribe();
      }
				   /// Destructor, removing the subscription.
      ~SmartPointer()
      {
	t->unsubscribe();
      }
  
				   /// Conversion to normal pointer.
      operator T* ()
      {
	return t;
      }
				   /// Conversion to normal const pointer.
      operator const T* () const
      {
	return t;
      }
      T& operator* ()
      {
	return *t;
      }
      const T& operator* () const
      {
	return *t;
      }
      T* operator -> ()
      {
	return t;
      }
      const T* operator -> () const
      {
	return t;
      }

};




/*----------------------------   subscriptor.h     ---------------------------*/
/* end of #ifndef __subscriptor_H */
#endif
/*----------------------------   subscriptor.h     ---------------------------*/
