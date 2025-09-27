// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------
#ifndef dealii_cxx20_functional_h
#define dealii_cxx20_functional_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_HAVE_CXX20
#  include <functional>
#else
#  include <functional>
#  include <tuple>
#  include <utility>
#endif

DEAL_II_NAMESPACE_OPEN

namespace std_cxx20
{
#ifndef DEAL_II_HAVE_CXX20
#  ifndef DOXYGEN

  namespace internal
  {
    /**
     * This call is used to invoke the lambda with stored bound arguments.
     *
     * @note The extra call is required to pattern-match the index_sequence.
     */
    template <typename F,
              size_t... Ind,
              typename BoundArgsTuple,
              typename... CallArgs>
    constexpr decltype(auto)
    call_bind(F &&function,
              std::index_sequence<Ind...>,
              BoundArgsTuple &&bound_args,
              CallArgs &&...call_args)
    {
      return std::invoke(std::forward<F>(function),
                         std::get<Ind>(
                           std::forward<BoundArgsTuple>(bound_args))...,
                         std::forward<CallArgs>(call_args)...);
    }

    /**
     * Return a callable which closely approximates what std::bind_front() is
     * doing.
     */
    template <typename F, typename... BoundArgs>
    decltype(auto)
    make_bind_front(F &&f, BoundArgs &&...bound_args)
    {
      return [f          = std::forward<F>(f),
              bound_args = std::make_tuple(
                std::forward<BoundArgs>(bound_args)...)](auto &&...call_args) {
        // Perform actual call inside a helper functions which allows to use
        // pattern-matching to the index sequence.
        return call_bind(f,
                         std::index_sequence_for<BoundArgs...>{},
                         bound_args,
                         std::forward<decltype(call_args)>(call_args)...);
      };
    }

  } // namespace internal

#  endif

  /**
   * This function generates a forwarding call wrapper to the function @p f
   * which has its first `n` arguments bound to the `n` arguments passed in
   * @p bound_args. The implementation is  an approximation to the functionality
   * provided by `std::bind_front`, see
   * https://en.cppreference.com/w/cpp/utility/functional/bind_front.
   *
   * This function allows to remove boilerplate in user code. Often functions in
   * the deal.II library take a `std::function` object as an argument. While you
   * can pass many different types of callable objects to a `std::function`
   * directly, you cannot pass a pointer to a member function. Instead you can
   * use the present function to create a compatible callable on-the-fly which
   * can be assigned to `std::function`. This is demonstrated in the following
   * example:
   *
   * @code
   *
   * // An exemplary function that takes an std::function as argument.
   * void example_function(const std::function<int(double)>& fn);
   *
   * class MyClass
   * {
   * public:
   *  void do_something()
   *  {
   *    // Pass the member function
   *    example_function(
   *      std_cxx20::bind_front(&MyClass::my_function, this));
   *  }
   *
   * private:
   *  // This function has a signature which is compatible with the function
   *  // that is expected by example_function().
   *  int my_function(double);
   * };
   *
   * @endcode
   */
  template <typename F, typename... BoundArgs>
  decltype(auto)
  bind_front(F &&f, BoundArgs &&...bound_args)
  {
    return internal::make_bind_front(std::forward<F>(f),
                                     std::forward<BoundArgs>(bound_args)...);
  }
#else
  using std::bind_front;
#endif

} // namespace std_cxx20

DEAL_II_NAMESPACE_CLOSE

#endif
