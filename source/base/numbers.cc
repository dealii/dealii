// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/numbers.h>

#include <cmath>
#include <complex>
#include <limits>

DEAL_II_NAMESPACE_OPEN

namespace numbers
{
  template <typename number>
  constexpr bool NumberTraits<number>::is_complex;

#ifndef DOXYGEN
  template <typename number>
  constexpr bool NumberTraits<std::complex<number>>::is_complex;
#endif
  // explicit instantiations
  template struct NumberTraits<double>;
  template struct NumberTraits<float>;
  template struct NumberTraits<long double>;

  template struct NumberTraits<std::complex<double>>;
  template struct NumberTraits<std::complex<float>>;
  template struct NumberTraits<std::complex<long double>>;
} // namespace numbers

DEAL_II_NAMESPACE_CLOSE
