//----------------------------  template_constraints.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2003 by the deal authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  template_constraints.h  ---------------------------
#ifndef __deal2__template_constraints_h
#define __deal2__template_constraints_h


#include <base/config.h>


template <bool, typename> struct constraint_and_return_value;


/**
 * This specialization of the general template for the case of a
 * @p{true} first template argument declares a local typedef @p{type}
 * to the second template argument. It is used in order to construct
 * constraints on template arguments in template (and member template)
 * functions. The negative specialization is missing.
 * 
 * Here's how the trick works, called SFINAE (substitution failure is
 * not an error): The C++ standard prescribes that a template function
 * is only considered in a call, if all parts of its signature can be
 * instantiated with the template parameter replaced by the respective
 * types/values in this particular call. Example:
 * @begin{verbatim}
 *   template <typename T>
 *   typename T::type  foo(T) {...};
 *   ...
 *   foo(1);
 * @end{verbatim}
 * The compiler should detect that in this call, the template
 * parameter T must be identified with the type "int". However,
 * the return type T::type does not exist. The trick now is
 * that this is not considered an error: this template is simply
 * not considered, the compiler keeps on looking for another 
 * possible function foo.
 * 
 * The idea is then to make the return type un-instantiatable if
 * certain constraints on the template types are not satisfied:
 * @begin{verbatim}
 *   template <bool, typename> struct constraint_and_return_value;
 *   template <typename T> struct constraint_and_return_value<true,T> {
 *     typedef T type;
 *   };
 * @end{verbatim}
 * constraint_and_return_value<false,T> is not defined. Given something like
 * @begin{verbatim}
 *   template <typename>
 *   struct int_or_double         { static const bool value = false;};
 *   template <>
 *   struct int_or_double<int>    { static const bool value = true; };
 *   template <>
 *   struct int_or_double<double> { static const bool value = true; };
 * @end{verbatim}
 * we can write a template
 * @begin{verbatim}
 *   template <typename T>
 *   typename constraint_and_return_value<int_or_double<T>::value,void>::type
 *   f (T);
 * @end{verbatim}
 * which can only be instantiated if T=int or T=double. A call to
 * f('c') will just fail with a compiler error: "no instance of
 * f(char) found". On the other hand, if the predicate in the first
 * argument to the constraint_and_return_value template is true, then
 * the return type is just the second type in the template.
 *
 * @author Wolfgang Bangerth, 2003
 */ 
template <typename T> struct constraint_and_return_value<true,T>
{
    typedef T type;
};


#ifdef DEAL_II_SFINAE_BUG

/**
 * Closure class in case the compiler lacks support for the SFINAE
 * concept. If the compiler supports it, only the specialization for
 * the positive case is available.
 *
 * @author Wolfgang Bangerth, 2003
 */
template <typename T> struct constraint_and_return_value<false,T>
{
    typedef T type;
};

#endif


#endif
