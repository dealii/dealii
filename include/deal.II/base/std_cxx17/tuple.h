// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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
#ifndef dealii_cxx17_tuple_h
#define dealii_cxx17_tuple_h

#include <deal.II/base/config.h>

#include <tuple>

DEAL_II_NAMESPACE_OPEN
namespace std_cxx17
{
#ifndef DEAL_II_HAVE_CXX17
  template <typename F, typename Tuple, size_t... S>
  auto
  apply_impl(F &&fn, Tuple &&t, std::index_sequence<S...>)
    -> decltype(std::forward<F>(fn)(std::get<S>(std::forward<Tuple>(t))...))
  {
    return std::forward<F>(fn)(std::get<S>(std::forward<Tuple>(t))...);
  }

  template <typename F, typename Tuple>
  auto
  apply(F &&fn, Tuple &&t) -> decltype(apply_impl(
    std::forward<F>(fn),
    std::forward<Tuple>(t),
    std::make_index_sequence<
      std::tuple_size<typename std::remove_reference<Tuple>::type>::value>()))
  {
    std::size_t constexpr tSize =
      std::tuple_size<typename std::remove_reference<Tuple>::type>::value;
    return apply_impl(std::forward<F>(fn),
                      std::forward<Tuple>(t),
                      std::make_index_sequence<tSize>());
  }
#else
  using std::apply;
#endif
} // namespace std_cxx17
DEAL_II_NAMESPACE_CLOSE

#endif
