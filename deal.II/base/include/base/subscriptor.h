/*----------------------------   subscriptor.h     ---------------------------*/
/*      $Id$                 */
#ifndef __subscriptor_H
#define __subscriptor_H
/*----------------------------   subscriptor.h     ---------------------------*/


#ifndef __exceptions_H
#include <base/exceptions.h>
#endif


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
				     /**
				      * Store the number of objects which
				      * subscribed to this object. Initialally,
				      * this number is zero, and upon
				      * destruction it shall be zero again
				      * (i.e. all objects which subscribed
				      * should have unsubscribed again).
				      *
				      * The creator (and owner) of an object
				      * is not counted.
				      */
    mutable unsigned int counter;
  public:
				     /**
				      * Object may not be deleted, since
				      * it is used.
				      */
    DeclException0(InUse);//, unsigned, "Object is used " << counter << " times");
    DeclException0(NotUsed);
  
				     /**
				      * Constructor setting the counter to zero.
				      */
    Subscriptor()
		    : counter(0)
      {}
				     /**
				      * Destructor, asserting that the counter
				      * is zero.
				      */
    
    ~Subscriptor()
      {
	Assert(counter==0, InUse());
      }

				     /**
				      * Copy-constructor.
				      *
				      * The counter of the copy is zero,
				      * since references point to the
				      * original object.
				      */
    Subscriptor(const Subscriptor&)
		    : counter(0)
      {}

				     /**
				      * Assignment operator.
				      *
				      * This has to be handled with
				      * care, too, because the counter
				      * has to remain the same.
				      */
    Subscriptor& operator = (const Subscriptor&) 
      {
	return *this;
      }
    
				     /**
				      * Subscribes a user of the object.
				      */
    void subscribe () const
      {
	++counter;
      }
  
				     /**
				      * Unsubscribes a user from the object.
				      */
    void unsubscribe () const
      {
	Assert(counter>0, NotUsed());
	--counter;
      }
};





/**
 * Smart references avoid destruction of a referenced object.  This
 * class has not been fully developed, since the compiler could not
 * resolve the dot operator in a convenient manner. The use of
 * #SmartPointer# is recommended, instead.
 */
template<class T>
class SmartReference
{
    T& t;

  public:
				     /**
				      * Constructor taking a normal reference.
				      */
    SmartReference(const T& tt)
		    : t(tt) 
      {
	t.subscribe();
      }
    
				     /**
				      * Destructor, removing the subscription.
				      */
    ~SmartReference()
      {
	t.unsubscribe();
      }
  
				     /**
				      * Conversion to normal reference
				      */
    operator T& () const
      {
	return t;
      }
};







/*----------------------------   subscriptor.h     ---------------------------*/
/* end of #ifndef __subscriptor_H */
#endif
/*----------------------------   subscriptor.h     ---------------------------*/
