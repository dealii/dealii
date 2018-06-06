// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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
#ifndef dealii_cxx14_utility_h
#define dealii_cxx14_utility_h

#include <deal.II/base/config.h>

#include <utility>

#ifndef DEAL_II_WITH_CXX14
// needed for array type check
#  include <type_traits>
#endif

DEAL_II_NAMESPACE_OPEN
namespace std_cxx14
{
#ifdef DEAL_II_WITH_CXX14
  using std::index_sequence;
  using std::index_sequence_for;
  using std::make_index_sequence;
#else
  template <size_t... Ints>
  struct index_sequence
  {
    using type       = index_sequence;
    using value_type = size_t;
    static constexpr std::size_t
    size() noexcept
    {
      return sizeof...(Ints);
    }
  };

  // --------------------------------------------------------------
  namespace internal
  {
    template <class Sequence1, class Sequence2>
    struct merge_and_renumber;

    template <size_t... I1, size_t... I2>
    struct merge_and_renumber<index_sequence<I1...>, index_sequence<I2...>>
      : index_sequence<I1..., (sizeof...(I1) + I2)...>
    {};
  } // namespace internal
  // --------------------------------------------------------------

  template <size_t N>
  struct make_index_sequence : internal::merge_and_renumber<
                                 typename make_index_sequence<N / 2>::type,
                                 typename make_index_sequence<N - N / 2>::type>
  {};

  template <>
  struct make_index_sequence<0> : index_sequence<>
  {};
  template <>
  struct make_index_sequence<1> : index_sequence<0>
  {};

  template <class... T>
  using index_sequence_for = make_index_sequence<sizeof...(T)>;
#endif
} // namespace std_cxx14
DEAL_II_NAMESPACE_CLOSE

#endif // dealii_cxx14_memory_h
