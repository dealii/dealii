// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2022 by the deal.II authors
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
#include <utility>

DEAL_II_NAMESPACE_OPEN

// Detection idiom adapted from Version 2 of the C++ Extensions for Library
// Fundamentals, ISO/IEC TS 19568:2017
namespace internal
{
  /**
   * A namespace used to declare the machinery for detecting whether a specific
   * class supports an operation. This approach simulates C++20-style
   * concepts with language standards before C++20.
   */
  namespace SupportsOperation
  {
    template <class...>
    using void_t = void;

    /**
     * The primary template class used in detecting operations. If the
     * compiler does not choose the specialization, then the fall-back
     * case is this general template, which then declares member variables
     * and types according to the failed detection.
     */
    template <class Default,
              class AlwaysVoid,
              template <class...>
              class Op,
              class... Args>
    struct detector
    {
      using value_t = std::false_type;
      using type    = Default;
    };

    /**
     * A specialization of the general template.
     *
     * The trick this class uses is that, just like the general template,
     * its second argument is always `void`, but here it is written as
     * `void_t<Op<Args...>>` and consequently the compiler will only select this
     * specialization if `Op<Args...>` is in fact a valid type. This means that
     * the operation we seek to understand is indeed supported.
     *
     * This specialization then declares member variables and types according
     * to the successful detection.
     */
    template <class Default, template <class...> class Op, class... Args>
    struct detector<Default, void_t<Op<Args...>>, Op, Args...>
    {
      using value_t = std::true_type;
      using type    = Op<Args...>;
    };


    /**
     * A base class for the `nonesuch` type to inherit from so it is not an
     * aggregate.
     */
    struct nonesuch_base
    {};

    /**
     * A type that can not be used in any reasonable way and consequently
     * can be used to indicate a failed detection in template metaprogramming.
     */
    struct nonesuch : private nonesuch_base
    {
      ~nonesuch()                = delete;
      nonesuch(nonesuch const &) = delete;
      void
      operator=(nonesuch const &) = delete;
    };

    template <class Default, template <class...> class Op, class... Args>
    using detected_or = detector<Default, void, Op, Args...>;

    template <template <class...> class Op, class... Args>
    using is_detected = typename detected_or<nonesuch, Op, Args...>::value_t;

    template <template <class...> class Op, class... Args>
    using detected_t = typename detected_or<nonesuch, Op, Args...>::type;

    template <class Default, template <class...> class Op, class... Args>
    using detected_or_t = typename detected_or<Default, Op, Args...>::type;

    template <class Expected, template <class...> class Op, class... Args>
    using is_detected_exact = std::is_same<Expected, detected_t<Op, Args...>>;

    template <class To, template <class...> class Op, class... Args>
    using is_detected_convertible =
      std::is_convertible<detected_t<Op, Args...>, To>;
  } // namespace SupportsOperation


  /**
   * A `constexpr` variable that describes whether or not `Op<Args...>` is a
   * valid expression.
   *
   * The way this is used is to define an `Op` operation template that
   * describes the operation we want to perform, and `Args` is a template
   * pack that describes the arguments to the operation. This variable
   * then states whether the operation, with these arguments, leads to
   * a valid C++ expression.
   *
   * An example is if one wanted to find out whether a type `T` has
   * a `get_mpi_communicator()` member function. In that case, one would write
   * the operation as
   * @code
   * template <typename T>
   * using get_mpi_communicator_op
   *   = decltype(std::declval<T>().get_mpi_communicator());
   * @endcode
   * and could define a variable like
   * @code
   * template <typename T>
   * constexpr bool has_get_mpi_communicator =
   * is_supported_operation<get_mpi_communicator_op, T>;
   * @endcode
   *
   * The trick used here is that `get_mpi_communicator_op` is a general
   * template, but when used with a type that does *not* have a
   * `get_mpi_communicator()` member variable, the `decltype(...)` operation
   * will fail because its argument does not represent a valid expression for
   * such a type. In other words, for such types `T` that do not have such a
   * member function, the general template `get_mpi_communicator_op` represents
   * a valid declaration, but the instantiation `get_mpi_communicator_op<T>`
   * is not, and the variable declared here detects and reports this.
   */
  template <template <class...> class Op, class... Args>
  constexpr bool is_supported_operation =
    SupportsOperation::is_detected<Op, Args...>::value;
} // namespace internal



namespace internal
{
  namespace TemplateConstraints
  {
    // TODO: Once we are able to use DEAL_II_HAVE_CXX17, the following classes
    // can be made much simpler with the help of fold expressions, see
    // https://en.cppreference.com/w/cpp/language/fold

    // helper struct for is_base_of_all and all_same_as
    template <bool... Values>
    struct BoolStorage;


    /**
     * A helper class whose `value` member is true or false depending on
     * whether all of the given boolean template arguments are `true`.
     * The class works by comparing the list of boolean values
     * `true, Values...` with the list `Values..., true` (i.e., with
     * its rotated self). The two are only the same if `Values...` is
     * a list of only `true` values.
     */
    template <bool... Values>
    struct all_true
    {
      static constexpr bool value =
        std::is_same<BoolStorage<Values..., true>,
                     BoolStorage<true, Values...>>::value;
    };


    /**
     * A class whose `value` member is set to `true` if any of the
     * boolean template arguments are true.
     */
    template <bool... Values>
    struct any_true;


    template <bool V1, bool... Values>
    struct any_true<V1, Values...>
    {
      static constexpr bool value = V1 || any_true<Values...>::value;
    };


    template <>
    struct any_true<>
    {
      static constexpr bool value = false;
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
 * argument. The result is stored in the member variable `value`.
 */
template <class Type, class... Types>
struct all_same_as
{
  static constexpr bool value = internal::TemplateConstraints::all_true<
    std::is_same<Type, Types>::value...>::value;
};



/**
 * This struct is a generalization of std::is_same to template
 * parameter packs and tests if any of the types in the `Types...`
 * parameter pack are equal to the `Type` given as first template
 * argument. The result is stored in the member variable `value`.
 */
template <class Type, class... Types>
struct is_same_as_any_of
{
  static constexpr bool value = internal::TemplateConstraints::any_true<
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
using begin_and_end_t =
  decltype(std::begin(std::declval<T>()), std::end(std::declval<T>()));

template <typename T>
constexpr bool has_begin_and_end =
  internal::is_supported_operation<begin_and_end_t, T>;



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
 */
template <typename T>
struct identity
{
  using type = T;
};



/**
 * A class that always returns a given value.
 * This is needed as a workaround for lambdas used as default parameters
 * some compilers struggle to deal with.
 */
template <typename ArgType, typename ValueType>
struct always_return
{
  ValueType value;
  ValueType
  operator()(const ArgType &)
  {
    return value;
  }
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
