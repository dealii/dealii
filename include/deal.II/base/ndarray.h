// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_ndarray_h
#define dealii_ndarray_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

#include <array>

DEAL_II_NAMESPACE_OPEN

#ifndef DOXYGEN
namespace internal
{
  namespace ndarray
  {
    // clang-format off
    /**
     * A variadic template helper class to recursively "unroll" the size
     * information of the ndarray. This is best explained on an example:
     * @code
     *    HelperArray<double, 1, 2, 3, 4>::type
     * == std::array<HelperArray<double, 2, 3, 4>::type, 1>
     * == std::array<std::array<HelperArray<double, 3, 4>::type, 2>, 1>
     * == std::array<std::array<std::array<HelperArray<double, 4>::type, 3>, 2>, 1>
     * == std::array<std::array<std::array<std::array<HelperArray<double>::type, 4>, 3>, 2>, 1>
     * == std::array<std::array<std::array<std::array<double, 4>, 3>, 2>, 1>
     * @endcode
     */
    template <typename T, std::size_t... Ns>
    struct HelperArray;
    // clang-format on

    /**
     * Recursively define the type alias "type" of HelperArray<T, N, ...Ns>
     * by wrapping a std::array around HelperArray<T, Ns...>::type
     */
    template <typename T, std::size_t N, std::size_t... Ns>
    struct HelperArray<T, N, Ns...>
    {
      using type = std::array<typename HelperArray<T, Ns...>::type, N>;
    };

    /**
     * End recursion once no std::size_t template parameters are left and
     * simply set the type alias to type T
     */
    template <typename T>
    struct HelperArray<T>
    {
      using type = T;
    };
  } // namespace ndarray
} // namespace internal
#endif // DOXYGEN

/**
 * A (variadic template) type alias for conveniently defining multidimensional
 * <a href="https://en.cppreference.com/w/cpp/container/array">std::array</a>s.
 *
 * The problem we try to address with the type alias is the following.
 * Suppose you want to create a multdimensional array of doubles of, for
 * example, rank 3, with sizes 2, 3, 4 for the first, middle, and last
 * index. Then using C-style arrays you could simply write
 * @code
 *   double my_array[2][3][4] = { ... };
 * @endcode
 * There are a number of good reasons why using C-style arrays is
 * generally discouraged (ranging from incompatibilities with STL functions
 * requiring awkward wrappers, surprises when comparing for equality, etc.).
 * If you want to do the same, however, using the more modern (and
 * encouraged) `std::array` class, then you would have to declare
 * @code
 *   std::array<std::array<std::array<double, 4>, 3>, 2> = { ... };
 * @endcode
 * The repetitions of `std::array` look awkward and, worse, the index
 * ranges have reversed: the leftmost index has range [0,2), the middle
 * index has range [0,3) and the rightmost index has range [0,4).
 * We address this issue by providing a class ndarray that allows to you
 * declare the above stacked `std::array` type by simply writing:
 * @code
 *   dealii::ndarray<double, 2, 3, 4> my_array = { ... };
 * @endcode
 *
 * @note dealii::ndarray is merely syntactic sugar in the form of a
 * <a href="https://en.cppreference.com/w/cpp/language/type_alias">type
 * alias</a>
 * (`using` declaration). It is not a deal.II specific class, but merely a
 * helper to cleanly define multidimensional arrays realized by "stacked"
 * `std::array` classes.
 */
template <typename T, std::size_t... Ns>
using ndarray = typename internal::ndarray::HelperArray<T, Ns...>::type;

DEAL_II_NAMESPACE_CLOSE

#endif
