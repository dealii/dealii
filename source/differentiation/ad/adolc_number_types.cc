// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2019 by the deal.II authors
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


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_ADOLC

#  include <deal.II/differentiation/ad/ad_number_traits.h>
#  include <deal.II/differentiation/ad/adolc_number_types.h>

#  include <functional>
#  include <utility>

DEAL_II_NAMESPACE_OPEN


#  ifdef DEAL_II_ADOLC_WITH_ADVANCED_BRANCHING

namespace numbers
{
  namespace internal
  {
    namespace
    {
      // Apply some comparator and extract the boolean result of the operation,
      // instead of the "adub" return type tpyically returned by ADOL-C for
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
          std::is_same<adouble, ad_type>::value,
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

#  include "adolc_number_types.inst"

DEAL_II_NAMESPACE_CLOSE

#endif
