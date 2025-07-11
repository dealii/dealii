// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_ADOLC

#  include <deal.II/differentiation/ad/ad_number_traits.h>
#  include <deal.II/differentiation/ad/adolc_number_types.h>

#  include <functional>
#  include <utility>

#  ifdef DEAL_II_WITH_ADOLC
#    include <adolc/adouble.h> // Taped double
#    include <adolc/adtl.h>    // Tapeless double
#  endif

DEAL_II_NAMESPACE_OPEN


#  ifdef DEAL_II_ADOLC_WITH_ADVANCED_BRANCHING

namespace numbers
{
  namespace internal
  {
    namespace
    {
      // Apply some comparator and extract the boolean result of the operation,
      // instead of the "adub" return type typically returned by ADOL-C for
      // such a comparison. This is implemented as a general function so that
      // the list of implemented comparative operations can be easily extended.
      bool
      adouble_boolean_comparator(
        const adouble &value_1,
        const adouble &value_2,
        const std::function<adouble(const adouble &, const adouble &)>
          &comparator)
      {
        using ad_type = typename Differentiation::AD::NumberTraits<
          double,
          Differentiation::AD::NumberTypes::adolc_taped>::ad_type;
        static_assert(
          std::is_same_v<adouble, ad_type>,
          "The type of the AD number is not that which was expected.");
        const ad_type result = comparator(value_1, value_2);
        return !(Differentiation::AD::ADNumberTraits<ad_type>::get_scalar_value(
                   result) == 0.0);
      }
    } // namespace
  }   // namespace internal

  bool
  values_are_equal(const adouble &value_1, const adouble &value_2)
  {
    return internal::adouble_boolean_comparator(
      value_1, value_2, [](const adouble &a, const adouble &b) -> adouble {
        return dealii::internal::NumberType<adouble>::value(a == b);
      });
  }

  bool
  value_is_less_than(const adouble &value_1, const adouble &value_2)
  {
    return internal::adouble_boolean_comparator(
      value_1, value_2, [](const adouble &a, const adouble &b) -> adouble {
        return dealii::internal::NumberType<adouble>::value(a < b);
      });
  }
} // namespace numbers

#  endif


/*---------------------- Explicit Instantiations ----------------------*/

#  include "differentiation/ad/adolc_number_types.inst"

DEAL_II_NAMESPACE_CLOSE

#endif
