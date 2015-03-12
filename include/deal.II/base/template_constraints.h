// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2015 by the deal.II authors
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

#include <complex>


DEAL_II_NAMESPACE_OPEN

template <bool, typename> struct constraint_and_return_value;


/**
 * This specialization of the general template for the case of a <tt>true</tt>
 * first template argument declares a local typedef <tt>type</tt> to the
 * second template argument. It is used in order to construct constraints on
 * template arguments in template (and member template) functions. The
 * negative specialization is missing.
 *
 * Here's how the trick works, called SFINAE (substitution failure is not an
 * error): The C++ standard prescribes that a template function is only
 * considered in a call, if all parts of its signature can be instantiated
 * with the template parameter replaced by the respective types/values in this
 * particular call. Example:
 * @code
 *   template <typename T>
 *   typename T::type  foo(T) {...};
 *   ...
 *   foo(1);
 * @endcode
 * The compiler should detect that in this call, the template parameter T must
 * be identified with the type "int". However, the return type T::type does
 * not exist. The trick now is that this is not considered an error: this
 * template is simply not considered, the compiler keeps on looking for
 * another possible function foo.
 *
 * The idea is then to make the return type un-instantiatable if certain
 * constraints on the template types are not satisfied:
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
 * which can only be instantiated if T=int or T=double. A call to f('c') will
 * just fail with a compiler error: "no instance of f(char) found". On the
 * other hand, if the predicate in the first argument to the
 * constraint_and_return_value template is true, then the return type is just
 * the second type in the template.
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
 * because the type <code>T</code> of the template can not be deduced in a
 * unique way from the types of the arguments. However, if the template is
 * written as
 * @code
 *   template <typename T> void f(T, typename identity<T>::type);
 * @endcode
 * then the call becomes valid: the type <code>T</code> is not deducible from
 * the second argument to the function, so only the first argument
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
 * <code>1</code> has that type. Of course, what we would like the compiler to
 * do is simply cast the <code>1</code> to <code>double</code>. We can achieve
 * this by writing the code as follows:
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
   * Comparison function for pointers of the same type. Returns @p true if the
   * two pointers are equal.
   */
  template <typename T>
  static bool equal (const T *p1, const T *p2);

  /**
   * Comparison function for pointers of different types. The C++ language
   * does not allow comparing these pointers using <tt>operator==</tt>.
   * However, since the two pointers have different types, we know that they
   * can't be the same, so we always return @p false.
   */
  template <typename T, typename U>
  static bool equal (const T *, const U *);
};



namespace internal
{
  /**
   * A type that is sometimes used for template tricks. For example, in some
   * situations one would like to do this:
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
   * <code>X::f()</code> without specializing the outer class first. One of
   * the common tricks is therefore to use something like this:
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
   * overloading, but that from inside the function <code>g()</code>, we can
   * still select which of the different <code>X::f()</code> we want based on
   * the <code>subdim</code> template argument.
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
 * A type that can be used to determine whether two types are equal. It allows
 * to write code like
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
 * Partial specialization of the general template for the case that both
 * template arguments are equal. See the documentation of the general template
 * for more information.
 */
template <typename T>
struct types_are_equal<T,T>
{
  static const bool value = true;
};



/**
 * A class with a local typedef that represents the type that results from the
 * product of two variables of type @p T and @p U. In other words, we would
 * like to infer the type of the <code>product</code> variable in code like
 * this:
 * @code
 *   T t;
 *   U u;
 *   auto product = t*u;
 * @endcode
 * The local typedef of this structure represents the type the variable
 * <code>product</code> would have.
 *
 *
 * <h3>Where is this useful</h3>
 *
 * The purpose of this class is principally to represent the type one needs to
 * use to represent the values or gradients of finite element fields at
 * quadrature points. For example, assume you are storing the values $U_j$ of
 * unknowns in a Vector<float>, then evaluating $u_h(x_q) = \sum_j U_j
 * \varphi_j(x_q)$ at quadrature points results in values $u_h(x_q)$ that need
 * to be stored as @p double variables because the $U_j$ are @p float values
 * and the $\varphi_j(x_q)$ are computed as @p double values, and the product
 * are then @p double values. On the other hand, if you store your unknowns
 * $U_j$ as <code>std::complex@<double@></code> values and you try to evaluate
 * $\nabla u_h(x_q) = \sum_j U_j \nabla\varphi_j(x_q)$ at quadrature points,
 * then the gradients $\nabla u_h(x_q)$ need to be stored as objects of type
 * <code>Tensor@<1,dim,std::complex@<double@>@></code> because that's what you
 * get when you multiply a complex number by a <code>Tensor@<1,dim@></code>
 * (the type used to represent the gradient of shape functions of scalar
 * finite elements).
 *
 * Likewise, if you are using a vector valued element (with dim components)
 * and the $U_j$ are stored as @p double variables, then $u_h(x_q) = \sum_j
 * U_j \varphi_j(x_q)$ needs to have type <code>Tensor@<1,dim@></code>
 * (because the shape functions have type <code>Tensor@<1,dim@></code>).
 * Finally, if you store the $U_j$ as objects of type
 * <code>std::complex@<double@></code> and you have a vector valued element,
 * then the gradients $\nabla u_h(x_q) = \sum_j U_j \nabla\varphi_j(x_q)$ will
 * result in objects of type <code>Tensor@<2,dim,std::complex@<double@>
 * @></code>.
 *
 * In all of these cases, this type is used to identify which type needs to be
 * used for the result of computing the product of unknowns and the values,
 * gradients, or other properties of shape functions.
 *
 * @author Wolfgang Bangerth, 2015
 */
template <typename T, typename U>
struct ProductType
{
#ifdef DEAL_II_WITH_CXX11
  typedef decltype(T() * U()) type;
#endif
};

#ifndef DEAL_II_WITH_CXX11

template <typename T>
struct ProductType<T,T>
{
  typedef T type;
};

template <typename T>
struct ProductType<T,bool>
{
  typedef T type;
};

template <typename T>
struct ProductType<bool, T>
{
  typedef T type;
};

template <>
struct ProductType<bool,double>
{
  typedef double type;
};

template <>
struct ProductType<double,bool>
{
  typedef double type;
};

template <>
struct ProductType<double,float>
{
  typedef double type;
};

template <>
struct ProductType<float,double>
{
  typedef double type;
};

template <>
struct ProductType<double,long double>
{
  typedef long double type;
};

template <>
struct ProductType<long double,double>
{
  typedef long double type;
};

template <>
struct ProductType<double,int>
{
  typedef double type;
};

template <>
struct ProductType<int,double>
{
  typedef double type;
};

template <>
struct ProductType<float,int>
{
  typedef float type;
};

template <>
struct ProductType<int,float>
{
  typedef float type;
};


#endif


// Annoyingly, there is no std::complex<T>::operator*(U) for scalars U
// other than T. Consequently, even with C++11, we need the following
// specializations:
template <typename T>
struct ProductType<std::complex<T>,std::complex<T> >
{
  typedef std::complex<T> type;
};

template <typename T, typename U>
struct ProductType<std::complex<T>,std::complex<U> >
{
  typedef std::complex<typename ProductType<T,U>::type> type;
};

template <typename U>
struct ProductType<double,std::complex<U> >
{
  typedef std::complex<typename ProductType<double,U>::type> type;
};

template <typename T>
struct ProductType<std::complex<T>,double>
{
  typedef std::complex<typename ProductType<T,double>::type> type;
};


template <typename U>
struct ProductType<float,std::complex<U> >
{
  typedef std::complex<typename ProductType<float,U>::type> type;
};

template <typename T>
struct ProductType<std::complex<T>,float>
{
  typedef std::complex<typename ProductType<T,float>::type> type;
};



/**
 * This class provides a local typedef @p type that is equal to the template
 * argument but only if the template argument corresponds to a scalar type
 * (i.e., one of the floating point types, signed or unsigned integer, or a
 * complex number). If the template type @p T is not a scalar, then no class
 * <code>EnableIfScalar@<T@></code> is declared and, consequently, no local
 * typedef is available.
 *
 * The purpose of the class is to disable certain template functions if one of
 * the arguments is not a scalar number. By way of (nonsensical) example,
 * consider the following function:
 * @code
 *   template <typename T>
 *   T multiply (const T t1, const T t2) { return t1*t2; }
 * @endcode
 * This function can be called with any two arguments of the same type @p T.
 * This includes arguments for which this clearly makes no sense.
 * Consequently, one may want to restrict the function to only scalars, and
 * this can be written as
 * @code
 *   template <typename T>
 *   typename EnableIfScalar<T>::type
 *   multiply (const T t1, const T t2) { return t1*t2; }
 * @endcode
 * At a place where you call the function, the compiler will deduce the type
 * @p T from the arguments. For example, in
 * @code
 *   multiply(1.234, 2.345);
 * @endcode
 * it will deduce @p T to be @p double, and because
 * <code>EnableIfScalar@<double@>::type</code> equals @p double, the compiler
 * will instantiate a function <code>double multiply(const double, const
 * double)</code> from the template above. On the other hand, in a context
 * like
 * @code
 *   std::vector<char> v1, v2;
 *   multiply(v1, v2);
 * @endcode
 * the compiler will deduce @p T to be <code>std::vector@<char@></code> but
 * because <code>EnableIfScalar@<std::vector@<char@>@>::type</code> does not
 * exist the compiler does not consider the template for instantiation. This
 * technique is called "Substitution Failure is not an Error (SFINAE)". It
 * makes sure that the template function can not even be called, rather than
 * leading to a later error about the fact that the operation
 * <code>t1*t2</code> is not defined (or may lead to some nonsensical result).
 * It also allows the declaration of overloads of a function such as @p
 * multiply for different types of arguments, without resulting in ambiguous
 * call errors by the compiler.
 *
 * @author Wolfgang Bangerth, 2015
 */
template <typename T>
struct EnableIfScalar;


template <> struct EnableIfScalar<double>
{
  typedef double type;
};


template <> struct EnableIfScalar<float>
{
  typedef float type;
};


template <> struct EnableIfScalar<long double>
{
  typedef long double type;
};


template <> struct EnableIfScalar<int>
{
  typedef int type;
};


template <> struct EnableIfScalar<unsigned int>
{
  typedef unsigned int type;
};



template <typename T> struct EnableIfScalar<std::complex<T> >
{
  typedef std::complex<T> type;
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
