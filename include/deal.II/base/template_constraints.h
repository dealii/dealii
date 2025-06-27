// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2003 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_template_constraints_h
#define dealii_template_constraints_h


#include <deal.II/base/config.h>

#include <deal.II/base/complex_overloads.h>
#include <deal.II/base/mpi_stub.h>
#include <deal.II/base/std_cxx20/type_traits.h>

#include <complex>
#include <type_traits>
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
      nonesuch(const nonesuch &) = delete;
      void
      operator=(const nonesuch &) = delete;
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
      static constexpr bool value = (Values && ...);
    };


    /**
     * A class whose `value` member is set to `true` if any of the
     * boolean template arguments are true.
     */
    template <bool... Values>
    struct any_true
    {
      static constexpr bool value = (Values || ...);
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
    std::is_base_of_v<Base, Derived>...>::value;
};



/**
 * This struct is a generalization of std::is_same to template
 * parameter packs and tests if all of the types in the `Types...`
 * parameter pack are equal to the `Type` given as first template
 * argument. The result is stored in the member variable `value`.
 */
template <typename Type, class... Types>
struct all_same_as
{
  static constexpr bool value = internal::TemplateConstraints::all_true<
    std::is_same_v<Type, Types>...>::value;
};



/**
 * This struct is a generalization of std::is_same to template
 * parameter packs and tests if any of the types in the `Types...`
 * parameter pack are equal to the `Type` given as first template
 * argument. The result is stored in the member variable `value`.
 */
template <typename Type, class... Types>
struct is_same_as_any_of
{
  static constexpr bool value = internal::TemplateConstraints::any_true<
    std::is_same_v<Type, Types>...>::value;
};



/*
 * A generalization of `std::enable_if` that only works if
 * <i>all</i> of the given boolean template parameters are
 * true. See [here](https://en.cppreference.com/w/cpp/types/enable_if)
 * for what `std::enable_if` does.
 *
 * @note
 * In contrast to `std::enable_if`, this template has no additional
 * template type (which for `std::enable_if` is defaulted to `void`).
 * As a consequence, this structure cannot be used for anything other
 * than enabling or disabling a template declaration; in particular
 * it cannot be used to set the return type of a function as one might
 * do with something like
 * @code
 *   template <typename T>
 *   typename std::enable_if<std::is_floating_point_v<T>, T>::type
 *   abs (const T t);
 * @endcode
 * which declares a function template `abs()` that can only be
 * instantiated if `T` is a floating point type; this function then
 * returns an object of type `T` as indicated by the last argument to
 * `std::enable_if`. The reason `enable_if_all` does not allow providing
 * this additional type is that variadic templates (here, the list of
 * `bool` arguments) must be the last template argument.
 */
template <bool... Values>
struct enable_if_all
  : std::enable_if<internal::TemplateConstraints::all_true<Values...>::value>
{};



/*
 * A generalization of `std::enable_if_t` that only works if
 * <i>all</i> of the given boolean template parameters are
 * true. See [here](https://en.cppreference.com/w/cpp/types/enable_if)
 * for what `std::enable_if_t` does.
 *
 * @note
 * In contrast to `std::enable_if_t`, this template has no additional
 * template type (which for `std::enable_if` is defaulted to `void`).
 * As a consequence, this structure cannot be used for anything other
 * than enabling or disabling a template declaration; in particular
 * it cannot be used to set the return type of a function as one might
 * do with something like
 * @code
 *   template <typename T>
 *   std::enable_if_t<std::is_floating_point_v<T>, T>
 *   abs (const T t);
 * @endcode
 * which declares a function template `abs()` that can only be
 * instantiated if `T` is a floating point type; this function then
 * returns an object of type `T` as indicated by the last argument to
 * `std::enable_if`. The reason `enable_if_all` does not allow providing
 * this additional type is that variadic templates (here, the list of
 * `bool` arguments) must be the last template argument.
 */
template <bool... Values>
using enable_if_all_t = typename enable_if_all<Values...>::type;


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
 * A `using` declaration to make the
 * [std::identity_type](https://en.cppreference.com/w/cpp/types/type_identity)
 * class available under the name that deal.II has used for a long time.
 *
 * @deprecated Use `std_cxx20::type_identity` instead.
 */
template <typename T>
using identity DEAL_II_DEPRECATED = std_cxx20::type_identity<T>;


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
    typename internal::ProductTypeImpl<std::decay_t<T>, std::decay_t<U>>::type;
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


// Forward declarations of vector types
template <typename Number>
class Vector;

template <typename Number>
class BlockVector;

namespace LinearAlgebra
{
  template <typename Number>
  class Vector;

  template <typename Number>
  class BlockVector;

  namespace distributed
  {
    template <typename Number, typename MemorySpace>
    class Vector;

    template <typename Number, typename MemorySpace>
    class BlockVector;
  } // namespace distributed
} // namespace LinearAlgebra

#ifdef DEAL_II_WITH_PETSC
namespace PETScWrappers
{
  class VectorBase;
  class Vector;
  class BlockVector;

  namespace MPI
  {
    class Vector;
    class BlockVector;

    class SparseMatrix;
    class BlockSparseMatrix;
  } // namespace MPI
} // namespace PETScWrappers
#endif

#ifdef DEAL_II_WITH_TRILINOS
namespace TrilinosWrappers
{
  namespace MPI
  {
    class Vector;
    class BlockVector;
  } // namespace MPI
} // namespace TrilinosWrappers

namespace LinearAlgebra
{
  namespace EpetraWrappers
  {
    class Vector;
  }

#  ifdef DEAL_II_TRILINOS_WITH_TPETRA
  namespace TpetraWrappers
  {
    template <typename Number, typename MemorySpace>
    class Vector;

    template <typename Number, typename MemorySpace>
    class BlockVector;

    template <typename Number, typename MemorySpace>
    class SparseMatrix;

    template <typename Number, typename MemorySpace>
    class BlockSparseMatrix;
  } // namespace TpetraWrappers
#  endif
} // namespace LinearAlgebra
#endif



/**
 * A namespace that is used to declare concepts used in C++20-style
 * `requires` clauses.
 */
namespace concepts
{
#if defined(DEAL_II_HAVE_CXX20) || defined(DOXYGEN)
  /**
   * A concept that identifies whether a template argument `C`
   * represents a [contiguous
   * container](https://en.cppreference.com/w/cpp/named_req/ContiguousContainer).
   * A contiguous container is a container object (such as `std::vector`,
   * `std::array`, or `boost::container::small_vector` that stores its elements
   * in one contiguous array in which we access all elements via a pointer to
   * the first element plus an offset. In contrast, linked lists, maps, and
   * similar objects are typically not stored as contiguous containers.
   */
  template <typename C>
  concept is_contiguous_container = requires(C &c) {
    {
      std::data(c)
    };
    {
      std::size(c)
    };
  };


  /**
   * A concept that tests that a combination of `dim` and `spacedim`
   * template arguments is valid. Specifically, we must have that
   * - `dim>=1`
   * - `spacedim<=3`
   * - `dim<=spacedim`.
   * These are the kinds of requirements that are imposed, for
   * example, on class Triangulation.
   */
  template <int dim, int spacedim>
  concept is_valid_dim_spacedim =
    (dim >= 1 && spacedim <= 3 && dim <= spacedim);

  namespace internal
  {
    /**
     * A template variable that returns whether the template argument is
     * a valid deal.II vector type. Its general definition is `false`, with
     * specializations dealing with actual vector types for which the
     * predicate is `true`.
     */
    template <typename T>
    inline constexpr bool is_dealii_vector_type = false;

    template <typename Number>
    inline constexpr bool is_dealii_vector_type<dealii::Vector<Number>> = true;

    template <typename Number>
    inline constexpr bool is_dealii_vector_type<dealii::BlockVector<Number>> =
      true;

    template <typename Number>
    inline constexpr bool
      is_dealii_vector_type<dealii::LinearAlgebra::BlockVector<Number>> = true;

    template <typename Number, typename MemorySpace>
    inline constexpr bool is_dealii_vector_type<
      dealii::LinearAlgebra::distributed::Vector<Number, MemorySpace>> = true;

    template <typename Number, typename MemorySpace>
    inline constexpr bool is_dealii_vector_type<
      dealii::LinearAlgebra::distributed::BlockVector<Number, MemorySpace>> =
      true;

#  ifdef DEAL_II_WITH_PETSC
    template <>
    inline constexpr bool is_dealii_vector_type<dealii::PETScWrappers::Vector> =
      true;

    template <>
    inline constexpr bool
      is_dealii_vector_type<dealii::PETScWrappers::BlockVector> = true;

    template <>
    inline constexpr bool
      is_dealii_vector_type<dealii::PETScWrappers::MPI::Vector> = true;

    template <>
    inline constexpr bool
      is_dealii_vector_type<dealii::PETScWrappers::MPI::BlockVector> = true;
#  endif

#  ifdef DEAL_II_WITH_TRILINOS
    template <>
    inline constexpr bool
      is_dealii_vector_type<dealii::TrilinosWrappers::MPI::Vector> = true;

    template <>
    inline constexpr bool
      is_dealii_vector_type<dealii::TrilinosWrappers::MPI::BlockVector> = true;

    template <>
    inline constexpr bool
      is_dealii_vector_type<dealii::LinearAlgebra::EpetraWrappers::Vector> =
        true;

#    ifdef DEAL_II_TRILINOS_WITH_TPETRA
    template <typename Number, typename MemorySpace>
    inline constexpr bool is_dealii_vector_type<
      dealii::LinearAlgebra::TpetraWrappers::Vector<Number, MemorySpace>> =
      true;

    template <typename Number, typename MemorySpace>
    inline constexpr bool is_dealii_vector_type<
      dealii::LinearAlgebra::TpetraWrappers::BlockVector<Number, MemorySpace>> =
      true;
#    endif
#  endif


    /**
     * A template variable that returns whether the template argument is
     * a valid deal.II vector type that is internally built on PETSc
     * functionality. Its general definition is `false`, with
     * specializations dealing with actual vector types for which the
     * predicate is `true`.
     */
    template <typename T>
    inline constexpr bool is_dealii_petsc_vector_type = false;

#  ifdef DEAL_II_WITH_PETSC
    template <>
    inline constexpr bool
      is_dealii_petsc_vector_type<dealii::PETScWrappers::VectorBase> = true;

    template <>
    inline constexpr bool
      is_dealii_petsc_vector_type<dealii::PETScWrappers::Vector> = true;

    template <>
    inline constexpr bool
      is_dealii_petsc_vector_type<dealii::PETScWrappers::BlockVector> = true;

    template <>
    inline constexpr bool
      is_dealii_petsc_vector_type<dealii::PETScWrappers::MPI::Vector> = true;

    template <>
    inline constexpr bool
      is_dealii_petsc_vector_type<dealii::PETScWrappers::MPI::BlockVector> =
        true;
#  endif


    /**
     * A template variable that returns whether the template argument is
     * a valid deal.II matrix type that is internally built on PETSc
     * functionality. Its general definition is `false`, with
     * specializations dealing with actual matrix types for which the
     * predicate is `true`.
     */
    template <typename T>
    inline constexpr bool is_dealii_petsc_matrix_type = false;

#  ifdef DEAL_II_WITH_PETSC
    template <>
    inline constexpr bool
      is_dealii_petsc_matrix_type<dealii::PETScWrappers::MPI::SparseMatrix> =
        true;

    template <>
    inline constexpr bool is_dealii_petsc_matrix_type<
      dealii::PETScWrappers::MPI::BlockSparseMatrix> = true;
#  endif
  } // namespace internal


  /**
   * A concept that tests whether a given template argument is a deal.II
   * vector type. This concept is used in many places, such as for the
   * functions in namespace VectorTools, where functions take a vector
   * as argument, but the type of the vector is a template argument. The
   * concept ensures that these functions cannot be called with an `int`
   * argument, for example, for which the declaration without the concept
   * would be perfectly valid but for which one would later get a linker
   * error because the function is only instantiated for deal.II vector
   * types.
   */
  template <typename VectorType>
  concept is_dealii_vector_type =
    internal::is_dealii_vector_type<std::remove_cv_t<VectorType>>;

  /**
   * A concept that tests whether a given template argument is a deal.II
   * vector type into which one can write. It is defined by asking
   * whether the is_dealii_vector_type concept is satisfied, and
   * that the template argument is not a `const`-qualified type. For
   * example, `is_writable_dealii_vector_type<dealii::Vector>` is true,
   * whereas `is_writable_dealii_vector_type<const dealii::Vector>` is
   * not.
   */
  template <typename VectorType>
  concept is_writable_dealii_vector_type =
    is_dealii_vector_type<VectorType> && (std::is_const_v<VectorType> == false);

  /**
   * A concept that tests whether a given template argument is a deal.II
   * vector type that internally builds on PETSc functionality. This
   * concept is used to constrain some classes that implement advanced
   * functionality based on PETSc and that requires that the vector
   * it works on are PETSc vectors. This includes, for example, the
   * time stepping and nonlinear solver classes in namespace PETScWrappers.
   */
  template <typename VectorType>
  concept is_dealii_petsc_vector_type =
    internal::is_dealii_petsc_vector_type<VectorType>;

  /**
   * A concept that tests whether a given template argument is a deal.II
   * matrix type that internally builds on PETSc functionality. This
   * concept is used to constrain some classes that implement advanced
   * functionality based on PETSc and that requires that the matrix
   * it works on are PETSc matrices. This includes, for example, the
   * time stepping and nonlinear solver classes in namespace PETScWrappers.
   */
  template <typename VectorType>
  concept is_dealii_petsc_matrix_type =
    internal::is_dealii_petsc_matrix_type<VectorType>;
#endif
} // namespace concepts


template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
class Triangulation;

template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
class DoFHandler;

namespace parallel
{
  namespace distributed
  {
    template <int dim, int spacedim>
    DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
    class Triangulation;
  }
  namespace shared
  {
    template <int dim, int spacedim>
    DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
    class Triangulation;
  }
  namespace fullydistributed
  {
    template <int dim, int spacedim>
    DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
    class Triangulation;
  }
} // namespace parallel

namespace concepts
{
#if defined(DEAL_II_HAVE_CXX20) || defined(DOXYGEN)
  namespace internal
  {
    template <typename T>
    inline constexpr bool is_triangulation_or_dof_handler = false;

    template <int dim, int spacedim>
    inline constexpr bool
      is_triangulation_or_dof_handler<Triangulation<dim, spacedim>> = true;

    template <int dim, int spacedim>
    inline constexpr bool is_triangulation_or_dof_handler<
      parallel::distributed::Triangulation<dim, spacedim>> = true;

    template <int dim, int spacedim>
    inline constexpr bool is_triangulation_or_dof_handler<
      parallel::shared::Triangulation<dim, spacedim>> = true;

    template <int dim, int spacedim>
    inline constexpr bool is_triangulation_or_dof_handler<
      parallel::fullydistributed::Triangulation<dim, spacedim>> = true;

    template <int dim, int spacedim>
    inline constexpr bool
      is_triangulation_or_dof_handler<DoFHandler<dim, spacedim>> = true;
  } // namespace internal


  /**
   * A concept that is used to check whether the `MeshType` template
   * type used in many functions in namespace GridTools and
   * VectorTools is in fact a "mesh" in the sense expected by these
   * functions. Specifically, this means that the type is either a
   * Triangulation or a DoFHandler type.
   */
  template <typename MeshType>
  concept is_triangulation_or_dof_handler =
    internal::is_triangulation_or_dof_handler<MeshType>;

  /**
   * A concept that tests whether a class `VectorType` has the required
   * interface to serve as a vector in vector-space operations -- principally
   * what is required to run iterative solvers: things such as norms, dot
   * products, etc.
   */
  template <typename VectorType>
  concept is_vector_space_vector = requires(VectorType                      U,
                                            VectorType                      V,
                                            VectorType                      W,
                                            typename VectorType::value_type a,
                                            typename VectorType::value_type b,
                                            typename VectorType::value_type s) {
    // Check local type requirements:
    typename VectorType::value_type;
    typename VectorType::size_type;
    typename VectorType::real_type;

    // Check some assignment and reinit operations;
    U.reinit(V);
    U.reinit(V, /* omit_zeroing_entries= */ true);

    U = V;
    U = a; // assignment of scalar

    U.equ(a, V);

    // Check scaling operations
    U *= a;
    U /= a;

    U.scale(V);

    // Vector additions:
    U += V;
    U -= V;

    U.add(a);
    U.add(a, V);
    U.add(a, V, b, W);

    U.sadd(s, a, V);

    // Norms and similar stuff:
    {
      U.mean_value()
    } -> std::convertible_to<typename VectorType::value_type>;

    {
      U.l1_norm()
    } -> std::convertible_to<typename VectorType::real_type>;

    {
      U.l2_norm()
    } -> std::convertible_to<typename VectorType::real_type>;

    {
      U.linfty_norm()
    } -> std::convertible_to<typename VectorType::real_type>;

    // Dot products:
    {
      U *V
    } -> std::convertible_to<typename VectorType::value_type>;

    {
      U.add_and_dot(a, V, W)
    } -> std::convertible_to<typename VectorType::value_type>;

    // Some queries:
    {
      U.size()
    } -> std::convertible_to<typename VectorType::size_type>;

    {
      U.all_zero()
    } -> std::same_as<bool>;

    {
      U.get_mpi_communicator()
    } -> std::same_as<MPI_Comm>;
  };


  /**
   * A concept that tests whether objects of type `MatrixType` can act
   * as linear operators on `VectorType`. In practice, that means that
   * `MatrixType` must have a `vmult()` member function that can take
   * a `VectorType` object as input and produce another `VectorType`
   * as output (both objects being taken as arguments to the `vmult()`
   * function).
   */
  template <typename MatrixType, typename VectorType>
  concept is_linear_operator_on =
    requires(const MatrixType &A, VectorType &dst, const VectorType &src) {
      A.vmult(dst, src);
    };


  /**
   * A concept that tests whether objects of type `MatrixType` can act
   * as the transposes of linear operators on `VectorType`. In practice, that
   * means that `MatrixType` must have a `Tvmult()` member function that can
   * take a `VectorType` object as input and produce another `VectorType`
   * as output (both objects being taken as arguments to the `vmult()`
   * function).
   */
  template <typename MatrixType, typename VectorType>
  concept is_transpose_linear_operator_on =
    requires(const MatrixType &A, VectorType &dst, const VectorType &src) {
      A.Tvmult(dst, src);
    };

#endif

} // namespace concepts



DEAL_II_NAMESPACE_CLOSE

#endif
