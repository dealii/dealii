//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__subscriptor_h
#define __deal2__subscriptor_h


#include <base/config.h>
#include <base/exceptions.h>

#include <typeinfo>
#include <map>
#include <string>

/**
 * Handling of subscriptions.
 *
 * This class, as a base class, allows to keep track of other objects
 * using a specific object. It is used, when an object, given to a
 * constructor by reference, is stored. Then, the original object may
 * not be deleted before the dependent object is deleted. You can assert
 * this constraint by letting the object passed be derived from this class
 * and let the user subscribe() to this object. The destructor the used
 * object inherits from the Subscriptor class then will lead to an error
 * when destruction is attempted while there are still subscriptions.
 *
 * The utility of this class is even enhanced by providing identifying
 * strings to the functions subscribe() and unsubscribe(). In case of
 * a hanging subscription during destruction, this string will be
 * listed in the exception's message. For reasons of efficiency, these
 * strings are handled as <tt>const char*</tt>. Therefore, the
 * pointers provided to subscribe() and to unsubscribe() must be the
 * same. Strings with equal contents will not be recognized to be the
 * same. The handling in SmartPointer will take care of this.
 *
 * @note Due to a problem with <tt>volatile</tt> declarations, this
 * additional feature is switched of if multithreading is used.
 *
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
				      * <tt>*this</tt>.
				      */
    Subscriptor& operator = (const Subscriptor&);
        
				     /**
				      * Subscribes a user of the
				      * object. The subscriber may be
				      * identified by text supplied as
				      * <tt>identifier</tt>.
				      */
    void subscribe (const char* identifier = 0) const;
      
				     /**
				      * Unsubscribes a user from the
				      * object.
				      *
				      * @note The <tt>identifier</tt>
				      * must be the <b>same
				      * pointer</b> as the one
				      * supplied to subscribe(), not
				      * just the same text.
				      */
    void unsubscribe (const char* identifier = 0) const;

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

    				     /** @addtogroup Exceptions
				      * @{ */
    
				     /**
				      * Exception:
				      * Object may not be deleted, since
				      * it is used.
				      */
    DeclException3(ExcInUse,
		   int, char*, std::string&,
		   << "Object of class " << arg2
		   << " is still used by " << arg1 << " other objects."
		   << arg3);

				     /**
				      * A subscriber with the
				      * identification string given to
				      * Subscriptor::unsubscribe() did
				      * not subscribe to the object.
				      */
    DeclException2(ExcNoSubscriber, char*, char*,
		   << "No subscriber with identifier \"" << arg2
		   << "\" did subscribe to this object of class " << arg1);
    
				     /**
				      * Exception: object should be
				      * used when
				      * Subscriptor::unsubscribe() is
				      * called.
				      */
    DeclException0(ExcNotUsed);
				     //@}

  private:
				     /**
				      * The data type used in
				      * #counter_map.
				      */
    typedef std::map<const char*, unsigned int>::value_type
    map_value_type;
    
				     /**
				      * The iterator type used in
				      * #counter_map.
				      */
    typedef std::map<const char*, unsigned int>::iterator
    map_iterator;
    
    				     /**
				      * Store the number of objects
				      * which subscribed to this
				      * object. Initialally, this
				      * number is zero, and upon
				      * destruction it shall be zero
				      * again (i.e. all objects which
				      * subscribed should have
				      * unsubscribed again).
				      *
				      * The creator (and owner) of an
				      * object is counted in the map
				      * below if HE manages to supply
				      * identification.
				      *
				      * We use the <tt>mutable</tt> keyword
				      * in order to allow subscription
				      * to constant objects also.
				      *
				      * In multithreaded mode, this
				      * counter may be modified by
				      * different threads. We thus
				      * have to mark it
				      * <tt>volatile</tt>. However, this is
				      * counter-productive in non-MT
				      * mode since it may pessimize
				      * code. So use the macro
				      * defined above to selectively
				      * add volatility.
				      */
    mutable DEAL_VOLATILE unsigned int counter;
    
				     /**
				      * In this map, we count
				      * subscriptions for each
				      * different identification strig
				      * supplied to subscribe().
				      */
    mutable std::map<const char*, unsigned int> counter_map;
    
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
    mutable const std::type_info * object_info;
};

//---------------------------------------------------------------------------

// If we are in optimized mode, the subscription checking is turned
// off. Therefore, we provide inline definitions of subscribe and
// unsubscribe here. The definitions for debug mode are in
// subscriptor.cc.

#ifndef DEBUG

inline void
Subscriptor::subscribe(const char*) const
{}


inline void
Subscriptor::unsubscribe(const char*) const
{}

#endif

#endif
