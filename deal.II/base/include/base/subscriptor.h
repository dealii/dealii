//----------------------------  subscriptor.h  ---------------------------
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
//----------------------------  subscriptor.h  ---------------------------
#ifndef __deal2__subscriptor_h
#define __deal2__subscriptor_h


#include <base/exceptions.h>
#include <typeinfo>



/**
 * Handling of subscriptions.
 *
 * This class, as a base class, allows to keep track of other objects
 * using a specific object. It is used, when an object, given to a
 * constructor by reference, is stored. Then, the original object may
 * not be deleted before the dependent object is deleted. You can assert
 * this constraint by letting the object passed be derived from this class
 * and let the user @p{subscribe} to this object. The destructor the used
 * object inherits from the @p{Subscriptor} class then will lead to an error
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
				      * Copy-constructor.
				      *
				      * The counter of the copy is zero,
				      * since references point to the
				      * original object.
				      */
    Subscriptor(const Subscriptor&);

				     /**
				      * Destructor, asserting that the counter
				      * is zero.
				      */
    virtual ~Subscriptor();
    
				     /**
				      * Assignment operator.
				      *
				      * This has to be handled with
				      * care, too, because the counter
				      * has to remain the same. It therefore
				      * does nothing more than returning
				      * @p{*this}.
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
				      * Return the present number of
				      * subscriptions to this object.
				      * This allows to use this class
				      * for reference counted lifetime
				      * determination where the last one
				      * to unsubscribe also deletes the
				      * object.
				      */
    unsigned int n_subscriptions () const;

				     /**
				      * Exception:
				      * Object may not be deleted, since
				      * it is used.
				      */
    DeclException2(ExcInUse,
		   int, char *,
		   << "Object of class " << arg2
		   << " is still used by " << arg1 << " other objects.");

				     /**
				      * Exception: object should be used
				      * when @p{unsubscribe} is called.
				      */
    DeclException0(ExcNotUsed);


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
				      * We use the @p{mutable} keyword in order
				      * to allow subscription to constant
				      * objects also.
				      */
    mutable unsigned int counter;

				     /**
				      * Pointer to the typeinfo object
				      * of this object, from which we
				      * can later deduce the class
				      * name. Since this information
				      * on the derived class is
				      * neither available in the
				      * destructor, nor in the
				      * constructor, we obtain it in
				      * between and store it here.
				      */
    mutable const type_info * object_info;
};


#endif
