// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_base_mutable_bind_h
#define dealii_base_mutable_bind_h

#include <deal.II/base/config.h>

#include <deal.II/base/patterns.h>

#include <tuple>
#include <utility>

DEAL_II_NAMESPACE_OPEN

namespace Utilities
{
  /**
   * A mutable version of std::bind, that binds all arguments of a function
   * pointer to a stored tuple, and allows you to update the tuple between
   * calls.
   *
   * An example usage of this class is through the helper function
   * mutable_bind() that creates a MutableBind object on the fly, based on its
   * arguments:
   *
   * @code
   * void my_function(const int &a, const double &b);
   *
   * auto bound = mutable_bind(my_function, 1, 2.0);
   *
   * bound(); // will execute my_function(1, 2.0);
   *
   * bound.set_arguments(2, 3.0);
   * bound(); // will execute my_function(2, 3.0);
   *
   * bound.parse_arguments("3: 4.0");
   * bound(); // will execute my_function(3, 4.0);
   * @endcode
   *
   * The arguments are copied to the tuple, with their reference and const
   * attributes removed. Only copy constructible objects are allowed as
   * function arguments. If you need to keep some references around, you may
   * wrap your function into a lambda function:
   *
   * @code
   *  void
   *  example_function(const Point<2> &p,
   *                   const double &d,
   *                   const unsigned int i = 3) {
   *  ...
   *  };
   *
   *  const Point<2> p(1, 2);
   *
   *  Utilities::MutableBind<void, double, unsigned int> exp = {
   *    [&p](const double &d,
   *         const unsigned int i)
   *    {
   *      example_function(p, d, i);
   *    },
   *    {}};
   *
   *  exp.parse_arguments("3.0 : 4");
   *  exp(); // calls example_function(p, 3.0, 4);
   * @endcode
   */
  template <typename ReturnType, class... FunctionArgs>
  class MutableBind
  {
  public:
    /**
     * An alias to the stored std::tuple type. Only copy constructible
     * objects are allowed as tuple members.
     */
    using TupleType =
      std::tuple<std::remove_cv_t<std::remove_reference_t<FunctionArgs>>...>;

    /**
     * Construct a MutableBind object specifying the function, and
     * each arguments separately.
     */
    template <typename FunctionType>
    MutableBind(FunctionType function, FunctionArgs &&...arguments);

    /**
     * Construct a MutableBind object specifying the function, and
     * the arguments as a tuple.
     */
    template <typename FunctionType>
    MutableBind(FunctionType function, TupleType &&arguments);

    /**
     * Construct a MutableBind object specifying only the function. By default,
     * the arguments are left to their default constructor values.
     */
    template <typename FunctionType>
    MutableBind(FunctionType function);

    /**
     * Call the original function, passing as arguments the elements of the
     * tuple of bound arguments.
     */
    ReturnType
    operator()() const;

    /**
     * Set the arguments to use in @p function, for next time
     * operator()() is called, using move semantic.
     */
    void
    set_arguments(TupleType &&arguments);

    /**
     * Set the arguments to use in @p function, for next time
     * operator()() is called, using move semantic.
     */
    void
    set_arguments(FunctionArgs &&...arguments);

    /**
     * Parse the arguments to use in @p function from a string, for next time
     * operator()() is called.
     *
     * The conversion is performed using a user supplied Patterns::PatternBase
     * object. By default, Patterns::Tools::Convert<TupleType>::to_pattern() is
     * used to determine how to convert from @p value_string to a TupleType
     * object.
     *
     * @param value_string The string to convert from
     * @param pattern A unique pointer to the pattern to use when performing
     * the conversion
     */
    void
    parse_arguments(const std::string           &value_string,
                    const Patterns::PatternBase &pattern =
                      *Patterns::Tools::Convert<TupleType>::to_pattern());

  private:
    /**
     * An std::function that stores the original function.
     */
    const std::function<ReturnType(FunctionArgs...)> function;

    /**
     * Currently stored arguments. These are forwarded to the function object
     * above, when calling operator()().
     */
    TupleType arguments;
  };



  /**
   * Create a MutableBind object from a function pointer and a list of
   * arguments.
   *
   * An example usage is given by:
   * @code
   * void my_function(const int &a, const double &b);
   *
   * auto bound = mutable_bind(my_function, 1, 2.0);
   *
   * bound(); // will execute my_function(1, 2.0);
   *
   * bound.set_arguments(2, 3.0);
   * bound(); // will execute my_function(2, 3.0);
   *
   * bound.parse_arguments("3: 4.0");
   * bound(); // will execute my_function(3, 4.0);
   * @endcode
   */
  template <typename ReturnType, class... FunctionArgs>
  MutableBind<ReturnType, FunctionArgs...>
  mutable_bind(ReturnType (*function)(FunctionArgs...),
               std_cxx20::type_identity_t<FunctionArgs> &&...arguments);

  /**
   * Same as above, using a std::function object.
   */
  template <typename ReturnType, class... FunctionArgs>
  MutableBind<ReturnType, FunctionArgs...>
  mutable_bind(std::function<ReturnType(FunctionArgs...)>,
               std_cxx20::type_identity_t<FunctionArgs> &&...arguments);

  /**
   * Create a MutableBind object from a function pointer, with uninitialized
   * arguments.
   *
   * Notice that if you do not call one of the MutableBind::set_arguments()
   * methods, or the MutableBind::parse_arguments() function on the returned
   * object, then the arguments passed to the function object will be
   * initialized with the values coming from each of the arguments' default
   * constructors.
   */
  template <typename ReturnType, class... FunctionArgs>
  MutableBind<ReturnType, FunctionArgs...>
    mutable_bind(ReturnType (*function)(FunctionArgs...));

  /**
   * Same as above, using a std::function object.
   */
  template <typename ReturnType, class... FunctionArgs>
  MutableBind<ReturnType, FunctionArgs...>
    mutable_bind(std::function<ReturnType(FunctionArgs...)>);



#ifndef DOXYGEN
  template <typename ReturnType, class... FunctionArgs>
  template <typename FunctionType>
  MutableBind<ReturnType, FunctionArgs...>::MutableBind(
    FunctionType function,
    FunctionArgs &&...arguments)
    : function(function)
    , arguments(std::make_tuple(std::move(arguments)...))
  {}



  template <typename ReturnType, class... FunctionArgs>
  template <typename FunctionType>
  MutableBind<ReturnType, FunctionArgs...>::MutableBind(FunctionType function,
                                                        TupleType  &&arguments)
    : function(function)
    , arguments(std::move(arguments))
  {}



  template <typename ReturnType, class... FunctionArgs>
  template <typename FunctionType>
  MutableBind<ReturnType, FunctionArgs...>::MutableBind(FunctionType function)
    : function(function)
  {}



  template <typename ReturnType, class... FunctionArgs>
  ReturnType
  MutableBind<ReturnType, FunctionArgs...>::operator()() const
  {
    return std::apply(function, arguments);
  }



  template <typename ReturnType, class... FunctionArgs>
  void
  MutableBind<ReturnType, FunctionArgs...>::set_arguments(
    FunctionArgs &&...args)
  {
    arguments = std::make_tuple(std::move(args)...);
  }



  template <typename ReturnType, class... FunctionArgs>
  void
  MutableBind<ReturnType, FunctionArgs...>::set_arguments(TupleType &&args)
  {
    arguments = std::move(args);
  }



  template <typename ReturnType, class... FunctionArgs>
  void
  MutableBind<ReturnType, FunctionArgs...>::parse_arguments(
    const std::string           &value_string,
    const Patterns::PatternBase &pattern)
  {
    arguments =
      Patterns::Tools::Convert<TupleType>::to_value(value_string, pattern);
  }



  template <typename ReturnType, class... FunctionArgs>
  MutableBind<ReturnType, FunctionArgs...>
  mutable_bind(ReturnType (*function)(FunctionArgs...),
               std_cxx20::type_identity_t<FunctionArgs> &&...arguments)
  {
    return MutableBind<ReturnType, FunctionArgs...>(function,
                                                    std::move(arguments)...);
  }



  template <typename ReturnType, class... FunctionArgs>
  MutableBind<ReturnType, FunctionArgs...>
  mutable_bind(ReturnType (*function)(FunctionArgs...))
  {
    return MutableBind<ReturnType, FunctionArgs...>(function);
  }



  template <typename ReturnType, class... FunctionArgs>
  MutableBind<ReturnType, FunctionArgs...>
  mutable_bind(std::function<ReturnType(FunctionArgs...)> function,
               std_cxx20::type_identity_t<FunctionArgs> &&...arguments)
  {
    return MutableBind<ReturnType, FunctionArgs...>(function,
                                                    std::move(arguments)...);
  }



  template <typename ReturnType, class... FunctionArgs>
  MutableBind<ReturnType, FunctionArgs...>
  mutable_bind(std::function<ReturnType(FunctionArgs...)> function)
  {
    return MutableBind<ReturnType, FunctionArgs...>(function);
  }
#endif
} // namespace Utilities

DEAL_II_NAMESPACE_CLOSE

#endif
