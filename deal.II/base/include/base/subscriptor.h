//----------------------------  subscriptor.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  subscriptor.h  ---------------------------
#ifndef __deal2__subscriptor_h
#define __deal2__subscriptor_h


#include <base/config.h>
#include <base/exceptions.h>
#include <typeinfo>


namespace internal
{
				   /**
				    * @internal
				    * A namespace in which we
				    * implement helper classes for the
				    * subscriptor class.
				    */
  namespace Subscriptor
  {
                                     /**
				      * @internal
                                      * Template class that declares
                                      * in inner typedef with the
                                      * following semantics: if the
                                      * first template parameter is
                                      * <tt>true</tt>, then the inner
                                      * typedef is <tt>volatile T</tt>,
                                      * otherwise <tt>T</tt>. We achieve
                                      * this behavior by partial
                                      * specialization of the general
                                      * template for both values of
                                      * the boolean argument.
                                      *
                                      * This trick is used to declare
                                      * the type of the counter
                                      * variable to be eiter volatile
                                      * or not, depending on whether
                                      * we are in multithreaded mode
                                      * or not. (If we are in MT mode,
                                      * then we need the <tt>volatile</tt>
                                      * specifier since more than one
                                      * thread my modify the counter
                                      * at the same time.
                                      *
                                      * Since we only partially
                                      * specialize the case that the
                                      * boolean template argument is
                                      * <tt>false</tt>, this general
                                      * template catches the case that
                                      * it is <tt>true</tt>, i.e. in a
                                      * sense it is also a
                                      * specialization since a
                                      * <tt>bool</tt> can only have two
                                      * states.
                                      *
                                      * @author Wolfgang Bangerth, 2003
                                      */
    template <bool, typename T> struct PossiblyVolatile
    {
        typedef volatile T type;
    };

                                     /**
				      * @internal
                                      * Specialization of the template
                                      * for the case that the first
                                      * template argument is
                                      * <tt>false</tt>, i.e. the
                                      * non-volatile case.
                                      */
    template <typename T> struct PossiblyVolatile<false,T>
    {
        typedef T type;
    };
  }
}



/**
 * Handling of subscriptions.
 *
 * This class, as a base class, allows to keep track of other objects
 * using a specific object. It is used, when an object, given to a
 * constructor by reference, is stored. Then, the original object may
 * not be deleted before the dependent object is deleted. You can assert
 * this constraint by letting the object passed be derived from this class
 * and let the user <tt>subscribe</tt> to this object. The destructor the used
 * object inherits from the <tt>Subscriptor</tt> class then will lead to an error
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
				      * <tt>*this</tt>.
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
				      * when <tt>unsubscribe</tt> is called.
				      */
    DeclException0(ExcNotUsed);


  private:
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
				      * object is not counted.
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
				      * code. So use the above
				      * template class to selectively
				      * add volatility.
				      */
    mutable internal::Subscriptor::PossiblyVolatile<DEAL_II_USE_MT,unsigned int>::type counter;

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


#endif
