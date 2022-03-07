// ---------------------------------------------------------------------
//
// Copyright (C) 2022 by the deal.II authors
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
#ifndef dealii_cxx17_variant_h
#define dealii_cxx17_variant_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_HAVE_CXX17
#  include <variant>
#else
DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#  include <boost/variant.hpp>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS
#endif

DEAL_II_NAMESPACE_OPEN
namespace std_cxx17
{
#ifndef DEAL_II_HAVE_CXX17
  using boost::get;
  using boost::variant;
#else

  template <std::size_t I, class... Types>
  constexpr std::add_pointer_t<
    std::variant_alternative_t<I, std::variant<Types...>>>
  get(std::variant<Types...> *pv) noexcept
  {
    return std::get_if<I>(pv);
  }

  template <std::size_t I, class... Types>
  constexpr std::add_pointer_t<
    const std::variant_alternative_t<I, std::variant<Types...>>>
  get(const std::variant<Types...> *pv) noexcept
  {
    return std::get_if<I>(pv);
  }

  template <class T, class... Types>
  constexpr std::add_pointer_t<T>
  get(std::variant<Types...> *pv) noexcept
  {
    return std::get_if<T>(pv);
  }

  template <class T, class... Types>
  constexpr std::add_pointer_t<const T>
  get(const std::variant<Types...> *pv) noexcept
  {
    return std::get_if<T>(pv);
  }


  using std::variant;
#endif
} // namespace std_cxx17
DEAL_II_NAMESPACE_CLOSE

#endif // dealii_cxx17_variant_h
