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
 * not be deleted before the dependent object is deleted. You can assert
 * this constraint by letting the object passed be derived from this class
 * and let the user #subscribe# to this object. The destructor the used
 * object inherits from the #Subscriptor# class then will lead to an error
 * when destruction is attempted while there are still subscriptions.
 */
class Subscriptor
{
  public:
				     /**
				      * Constructor setting the counter to zero.
				      */
    Subscriptor();

				     /**
				      * Destructor, asserting that the counter
				      * is zero.
				      */
    
    ~Subscriptor();
    
				     /**
				      * Copy-constructor.
				      *
				      * The counter of the copy is zero,
				      * since references point to the
				      * original object.
				      */
    Subscriptor(const Subscriptor&);
    
				     /**
				      * Assignment operator.
				      *
				      * This has to be handled with
				      * care, too, because the counter
				      * has to remain the same. It therefore
				      * does nothing more than returning
				      * #*this#.
				      */
    Subscriptor& operator = (const Subscriptor&);
        
				     /**
				      * Subscribes a user of the object.
				      */
    void subscribe () const;
      
				     /**
				      * Unsubscribes a user from the object.
				      */
    void unsubscribe () const;
    
				     /**
				      * Exception:
				      * Object may not be deleted, since
				      * it is used.
				      */
    DeclException0(InUse);
				     /**
				      * Exception: object should be used
				      * when #unsubscribe# is called.
				      */
    DeclException0(NotUsed);
  

  private:
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
				      *
				      * We use the #mutable# keyword in order
				      * to allow subscription to constant
				      * objects also.
				      */
    mutable unsigned int counter;
};





/*----------------------------   subscriptor.h     ---------------------------*/
/* end of #ifndef __subscriptor_H */
#endif
/*----------------------------   subscriptor.h     ---------------------------*/
