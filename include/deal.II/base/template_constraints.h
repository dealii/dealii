// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#ifndef dealii_template_constraints_h
#define dealii_template_constraints_h


#include <deal.II/base/config.h>

#include <deal.II/base/complex_overloads.h>

#include <complex>
#include <iterator>
#include <utility>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace TemplateConstraints
  {
    // helper struct for is_base_of_all and all_same_as
    template <bool... Values>
    struct BoolStorage;


    /**
     * A helper class whose `value` member is true or false depending on
     * whether all of the given boolean template arguments are true.
     */
    template <bool... Values>
    struct all_true
    {
      static constexpr bool value =
        std::is_same<BoolStorage<Values..., true>,
                     BoolStorage<true, Values...>>::value;
    };
  } // namespace TemplateConstraints
} // namespace internal

/**
 * This struct is a generalization of std::is_base_of<Base, Derived>
 * to template parameter packs and tests if all of the Derived...
 * classes have Base as base class or are Base itself. The result
 * is stored in the member variable value.
 */
template <class Base, class... Derived>
struct is_base_of_all
{
  static constexpr bool value = internal::TemplateConstraints::all_true<
    std::is_base_of<Base, Derived>::value...>::value;
};



/**
 * This struct is a generalization of std::is_same to template
 * parameter packs and tests if all of the types in the `Types...`
 * parameter pack are equal to the `Type` given as first template
 * argument. The result is stored in the member variable value.
 */
template <class Type, class... Types>
struct all_same_as
{
  static constexpr bool value = internal::TemplateConstraints::all_true<
    std::is_same<Type, Types>::value...>::value;
};



/*
 * A generalization of `std::enable_if` that only works if
 * <i>all</i> of the given boolean template parameters are
 * true.
 */
template <bool... Values>
struct enable_if_all
  : std::enable_if<internal::TemplateConstraints::all_true<Values...>::value>
{};



/**
 * A type trait that checks to see if a class behaves as an iterable container
 * that has a beginning and an end. This implies that the class either defines
 * the `begin()` and `end()` functions, or is a C-style array.
 */
template <typename T>
class has_begin_and_end
{
  template <typename C>
  static std::false_type
  test(...);

  template <typename C>
  static auto
  test(int) -> decltype(std::begin(std::declval<C>()),
                        std::end(std::declval<C>()),
                        std::true_type());

public:
  using type = decltype(test<T>(0));

  static const bool value = type::value;
};



template <bool, typename>
struct constraint_and_return_value;


/**
 * This specialization of the general template for the case of a <tt>true</tt>
 * first template argument declares a local alias <tt>type</tt> to the
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
 *   typename T::type  foo(T)
 *   {
 *     ...
 *   };
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
 *   template <bool, typename>
 *   struct constraint_and_return_value;
 *
 *   template <typename T>
 *   struct constraint_and_return_value<true,T>
 *   {
 *     using type = T;
 *   };
 * @endcode
 * constraint_and_return_value<false,T> is not defined. Given something like
 * @code
 *   template <typename>
 *   struct int_or_double
 *   {
 *     static const bool value = false;
 *   };
 *
 *   template <>
 *   struct int_or_double<int>
 *   {
 *     static const bool value = true;
 *   };
 *
 *   template <>
 *   struct int_or_double<double>
 *   {
 *     static const bool value = true;
 *   };
 * @endcode
 * we can write a template
 * @code
 *   template <typename T>
 *   typename constraint_and_return_value<int_or_double<T>::value,void>::type
 *     f (T);
 * @endcode
 * which can only be instantiated if T=int or T=double. A call to f('c') will
 * just fail with a compiler error: "no instance of f(char) found". On the
 * other hand, if the predicate in the first argument to the
 * constraint_and_return_value template is true, then the return type is just
 * the second type in the template.
 *
 * @deprecated Use std::enable_if instead.
 *
 * @author Wolfgang Bangerth, 2003
 */
template <typename T>
struct DEAL_II_DEPRECATED constraint_and_return_value<true, T>
{
  using type = T;
};



/**
 * A template class that simply exports its template argument as a local
 * alias. This class, while at first appearing useless, makes sense in the
 * following context: if you have a function template as follows:
 * @code
 * template <typename T>
 * void f(T, T);
 * @endcode
 * then it can't be called in an expression like <code>f(1, 3.141)</code>
 * because the type <code>T</code> of the template can not be deduced in a
 * unique way from the types of the arguments. However, if the template is
 * written as
 * @code
 * template <typename T>
 * void f(T, typename identity<T>::type);
 * @endcode
 * then the call becomes valid: the type <code>T</code> is not deducible from
 * the second argument to the function, so only the first argument
 * participates in template type resolution.
 *
 * The context for this feature is as follows: consider
 * @code
 * template <typename RT, typename A>
 * void forward_call(RT (*p) (A), A a)
 * {
 *   p(a);
 * }
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
 * void forward_call(RT (*p) (A), typename identity<A>::type a)
 * {
 *   p(a);
 * }
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
  using type = T;
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
  static bool
  equal(const T *p1, const T *p2)
  {
    return (p1 == p2);
  }


  /**
   * Comparison function for pointers of different types. The C++ language
   * does not allow comparing these pointers using <tt>operator==</tt>.
   * However, since the two pointers have different types, we know that they
   * can't be the same, so we always return @p false.
   */
  template <typename T, typename U>
  static bool
  equal(const T *, const U *)
  {
    return false;
  }
};



/**
 * A type that can be used to determine whether two types are equal. It allows
 * to write code like
 * @code
 *   template <typename T>
 *   void Vector<T>::some_operation ()
 *   {
 *     if (std::is_same<T,double>::value == true)
 *       call_some_blas_function_for_doubles;
 *     else
 *       do_it_by_hand;
 *   }
 * @endcode
 *
 * This construct is made possible through the existence of a partial
 * specialization of the class for template arguments that are equal.
 *
 * @deprecated Use the standard library type trait <code>std::is_same</code>
 * instead of this class.
 */
template <typename T, typename U>
struct DEAL_II_DEPRECATED types_are_equal : std::is_same<T, U>
{};



namespace internal
{
  /**
   * A struct that implements the default product type resulting from the
   * multiplication of two types.
   *
   * @note Care should be taken when @p T or @p U have qualifiers (@p const or
   * @p volatile) or are @p lvalue or @p rvalue references! It is recommended
   * that specialization of this class is only made for unqualified (fully
   * stripped) types and that the ProductType class be used to determine the
   * result of operating with (potentially) qualified types.
   *
   * @author Wolfgang Bangerth, Jean-Paul Pelteret, 2017
   */
  template <typename T, typename U>
  struct ProductTypeImpl
  {
    using type = decltype(std::declval<T>() * std::declval<U>());
  };

} // namespace internal



/**
 * A class with a local alias that represents the type that results from the
 * product of two variables of type @p T and @p U. In other words, we would
 * like to infer the type of the <code>product</code> variable in code like
 * this:
 * @code
 *   T t;
 *   U u;
 *   auto product = t*u;
 * @endcode
 * The local alias of this structure represents the type the variable
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
 * @author Wolfgang Bangerth, 2015, 2017
 */
template <typename T, typename U>
struct ProductType
{
  using type =
    typename internal::ProductTypeImpl<typename std::decay<T>::type,
                                       typename std::decay<U>::type>::type;
};

namespace internal
{
  // Annoyingly, there is no std::complex<T>::operator*(U) for scalars U
  // other than T (not even in C++11, or C++14). We provide our own overloads
  // in base/complex_overloads.h, but in order for them to work, we have to
  // manually specify all products we want to allow:

  template <typename T>
  struct ProductTypeImpl<std::complex<T>, std::complex<T>>
  {
    using type = std::complex<T>;
  };

  template <typename T, typename U>
  struct ProductTypeImpl<std::complex<T>, std::complex<U>>
  {
    using type = std::complex<typename ProductType<T, U>::type>;
  };

  template <typename U>
  struct ProductTypeImpl<double, std::complex<U>>
  {
    using type = std::complex<typename ProductType<double, U>::type>;
  };

  template <typename T>
  struct ProductTypeImpl<std::complex<T>, double>
  {
    using type = std::complex<typename ProductType<T, double>::type>;
  };

  template <typename U>
  struct ProductTypeImpl<float, std::complex<U>>
  {
    using type = std::complex<typename ProductType<float, U>::type>;
  };

  template <typename T>
  struct ProductTypeImpl<std::complex<T>, float>
  {
    using type = std::complex<typename ProductType<T, float>::type>;
  };

} // namespace internal



/**
 * This class provides a local alias @p type that is equal to the template
 * argument but only if the template argument corresponds to a scalar type
 * (i.e., one of the floating point types, signed or unsigned integer, or a
 * complex number). If the template type @p T is not a scalar, then no class
 * <code>EnableIfScalar@<T@></code> is declared and, consequently, no local
 * alias is available.
 *
 * The purpose of the class is to disable certain template functions if one of
 * the arguments is not a scalar number. By way of (nonsensical) example,
 * consider the following function:
 * @code
 *   template <typename T>
 *   T multiply (const T t1, const T t2)
 *   {
 *     return t1*t2;
 *   }
 * @endcode
 * This function can be called with any two arguments of the same type @p T.
 * This includes arguments for which this clearly makes no sense.
 * Consequently, one may want to restrict the function to only scalars, and
 * this can be written as
 * @code
 *   template <typename T>
 *   typename EnableIfScalar<T>::type
 *   multiply (const T t1, const T t2)
 *   {
 *     return t1*t2;
 *   }
 * @endcode
 * At a place where you call the function, the compiler will deduce the type
 * @p T from the arguments. For example, in
 * @code
 *   multiply(1.234, 2.345);
 * @endcode
 * it will deduce @p T to be @p double, and because
 * <code>EnableIfScalar@<double@>::%type</code> equals @p double, the compiler
 * will instantiate a function <code>double multiply(const double, const
 * double)</code> from the template above. On the other hand, in a context
 * like
 * @code
 *   std::vector<char> v1, v2;
 *   multiply(v1, v2);
 * @endcode
 * the compiler will deduce @p T to be <code>std::vector@<char@></code> but
 * because <code>EnableIfScalar@<std::vector@<char@>@>::%type</code> does not
 * exist the compiler does not consider the template for instantiation. This
 * technique is called "Substitution Failure is not an Error (SFINAE)". It
 * makes sure that the template function can not even be called, rather than
 * leading to a later error about the fact that the operation
 * <code>t1*t2</code> is not defined (or may lead to some nonsensical result).
 * It also allows the declaration of overloads of a function such as @p
 * multiply for different types of arguments, without resulting in ambiguous
 * call errors by the compiler.
 *
 * @author Wolfgang Bangerth, Matthias Maier, 2015 - 2017
 */
template <typename T>
struct EnableIfScalar;


template <>
struct EnableIfScalar<double>
{
  using type = double;
};

template <>
struct EnableIfScalar<float>
{
  using type = float;
};

template <>
struct EnableIfScalar<long double>
{
  using type = long double;
};

template <>
struct EnableIfScalar<int>
{
  using type = int;
};

template <>
struct EnableIfScalar<unsigned int>
{
  using type = unsigned int;
};

template <typename T>
struct EnableIfScalar<std::complex<T>>
{
  using type = std::complex<T>;
};


DEAL_II_NAMESPACE_CLOSE

#endif
