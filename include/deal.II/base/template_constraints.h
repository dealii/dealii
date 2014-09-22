// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef __deal2__template_constraints_h
#define __deal2__template_constraints_h


#include <deal.II/base/config.h>

DEAL_II_NAMESPACE_OPEN

template <bool, typename> struct constraint_and_return_value;


/**
 * This specialization of the general template for the case of a
 * <tt>true</tt> first template argument declares a local typedef <tt>type</tt>
 * to the second template argument. It is used in order to construct
 * constraints on template arguments in template (and member template)
 * functions. The negative specialization is missing.
 *
 * Here's how the trick works, called SFINAE (substitution failure is
 * not an error): The C++ standard prescribes that a template function
 * is only considered in a call, if all parts of its signature can be
 * instantiated with the template parameter replaced by the respective
 * types/values in this particular call. Example:
 * @code
 *   template <typename T>
 *   typename T::type  foo(T) {...};
 *   ...
 *   foo(1);
 * @endcode
 * The compiler should detect that in this call, the template
 * parameter T must be identified with the type "int". However,
 * the return type T::type does not exist. The trick now is
 * that this is not considered an error: this template is simply
 * not considered, the compiler keeps on looking for another
 * possible function foo.
 *
 * The idea is then to make the return type un-instantiatable if
 * certain constraints on the template types are not satisfied:
 * @code
 *   template <bool, typename> struct constraint_and_return_value;
 *   template <typename T> struct constraint_and_return_value<true,T> {
 *     typedef T type;
 *   };
 * @endcode
 * constraint_and_return_value<false,T> is not defined. Given something like
 * @code
 *   template <typename>
 *   struct int_or_double         { static const bool value = false;};
 *   template <>
 *   struct int_or_double<int>    { static const bool value = true; };
 *   template <>
 *   struct int_or_double<double> { static const bool value = true; };
 * @endcode
 * we can write a template
 * @code
 *   template <typename T>
 *   typename constraint_and_return_value<int_or_double<T>::value,void>::type
 *   f (T);
 * @endcode
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



/**
 * A template class that simply exports its template argument as a local
 * typedef. This class, while at first appearing useless, makes sense in the
 * following context: if you have a function template as follows:
 * @code
 *   template <typename T> void f(T, T);
 * @endcode
 * then it can't be called in an expression like <code>f(1, 3.141)</code>
 * because the type <code>T</code> of the template can not be deduced
 * in a unique way from the types of the arguments. However, if the
 * template is written as
 * @code
 *   template <typename T> void f(T, typename identity<T>::type);
 * @endcode
 * then the call becomes valid: the type <code>T</code> is not deducible
 * from the second argument to the function, so only the first argument
 * participates in template type resolution.
 *
 * The context for this feature is as follows: consider
 * @code
 * template <typename RT, typename A>
 * void forward_call(RT (*p) (A), A a)  { p(a); }
 *
 * void h (double);
 *
 * void g()
 * {
 *   forward_call(&h, 1);
 * }
 * @endcode
 * This code fails to compile because the compiler can't decide whether the
 * template type <code>A</code> should be <code>double</code> (from the
 * signature of the function given as first argument to
 * <code>forward_call</code>, or <code>int</code> because the expression
 * <code>1</code> has that type. Of course, what we would like the compiler
 * to do is simply cast the <code>1</code> to <code>double</code>. We can
 * achieve this by writing the code as follows:
 * @code
 * template <typename RT, typename A>
 * void forward_call(RT (*p) (A), typename identity<A>::type a)  { p(a); }
 *
 * void h (double);
 *
 * void g()
 * {
 *   forward_call(&h, 1);
 * }
 * @endcode
 *
 * @author Wolfgang Bangerth, 2008
 */
template <typename T>
struct identity
{
  typedef T type;
};



/**
 * A class to perform comparisons of arbitrary pointers for equality. In some
 * circumstances, one would like to make sure that two arguments to a function
 * are not the same object. One would, in this case, make sure that their
 * addresses are not the same. However, sometimes the types of these two
 * arguments may be template types, and they may be the same type or not. In
 * this case, a simple comparison as in <tt>&object1 != &object2</tt> does
 * only work if the types of the two objects are equal, but the compiler will
 * barf if they are not. However, in the latter case, since the types of the
 * two objects are different, we can be sure that the two objects cannot be
 * the same.
 *
 * This class implements a comparison function that always returns @p false if
 * the types of its two arguments are different, and returns <tt>p1 == p2</tt>
 * otherwise.
 *
 * @author Wolfgang Bangerth, 2004
 */
struct PointerComparison
{
  /**
   * Comparison function for pointers of
   * the same type. Returns @p true if the
   * two pointers are equal.
   */
  template <typename T>
  static bool equal (const T *p1, const T *p2);

  /**
   * Comparison function for pointers of
   * different types. The C++ language does
   * not allow comparing these pointers
   * using <tt>operator==</tt>. However,
   * since the two pointers have different
   * types, we know that they can't be the
   * same, so we always return @p false.
   */
  template <typename T, typename U>
  static bool equal (const T *, const U *);
};



namespace internal
{
  /**
   * A type that is sometimes used for template tricks. For example, in
   * some situations one would like to do this:
   *
   * @code
   *   template <int dim>
   *   class X {
   *     // do something on subdim-dimensional sub-objects of the big
   *     // dim-dimensional thing (for example on vertices/lines/quads of
   *     // cells):
   *     template <int subdim> void f();
   *   };
   *
   *   template <int dim>
   *   template <>
   *   void X<dim>::f<0> () { ...operate on the vertices of a cell... }
   *
   *   template <int dim, int subdim> void g(X<dim> &x) {
   *     x.f<subdim> ();
   *   }
   * @endcode
   *
   * The problem is: the language doesn't allow us to specialize
   * <code>X::f()</code> without specializing the outer class first. One
   * of the common tricks is therefore to use something like this:
   *
   * @code
   *   template <int N> struct int2type {};
   *
   *   template <int dim>
   *   class X {
   *     // do something on subdim-dimensional sub-objects of the big
   *     // dim-dimensional thing (for example on vertices/lines/quads of
   *     // cells):
   *     void f(int2type<0>);
   *     void f(int2type<1>);
   *     void f(int2type<2>);
   *     void f(int2type<3>);
   *   };
   *
   *   template <int dim>
   *   void X<dim>::f (int2type<0>) { ...operate on the vertices of a cell... }
   *
   *   template <int dim>
   *   void X<dim>::f (int2type<1>) { ...operate on the lines of a cell... }
   *
   *   template <int dim, int subdim> void g(X<dim> &x) {
   *     x.f (int2type<subdim>());
   *   }
   * @endcode
   *
   * Note that we have replaced specialization of <code>X::f()</code> by
   * overloading, but that from inside the function <code>g()</code>, we
   * can still select which of the different <code>X::f()</code> we want
   * based on the <code>subdim</code> template argument.
   *
   * @author Wolfgang Bangerth, 2006
   */
  template <int N>
  struct int2type
  {};


  /**
   * The equivalent of the int2type class for boolean arguments.
   *
   * @author Wolfgang Bangerth, 2009
   */
  template <bool B>
  struct bool2type
  {};
}



/**
 * A type that can be used to determine whether two types are equal.
 * It allows to write code like
 * @code
 *   template <typename T>
 *   void Vector<T>::some_operation () {
 *     if (types_are_equal<T,double>::value == true)
 *       call_some_blas_function_for_doubles;
 *     else
 *       do_it_by_hand;
 *   }
 * @endcode
 *
 * This construct is made possible through the existence of a partial
 * specialization of the class for template arguments that are equal.
 */
template <typename T, typename U>
struct types_are_equal
{
  static const bool value = false;
};


/**
 * Partial specialization of the general template for the case that
 * both template arguments are equal. See the documentation of the
 * general template for more information.
 */
template <typename T>
struct types_are_equal<T,T>
{
  static const bool value = true;
};



// --------------- inline functions -----------------


template <typename T, typename U>
inline
bool
PointerComparison::equal (const T *, const U *)
{
  return false;
}



template <typename T>
inline
bool
PointerComparison::equal (const T *p1, const T *p2)
{
  return (p1==p2);
}



DEAL_II_NAMESPACE_CLOSE

#endif
