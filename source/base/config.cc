// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2015 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#include <deal.II/base/config.h>

#include <cmath>
#include <complex>
#include <limits>

DEAL_II_NAMESPACE_OPEN

namespace numbers
{
  template <typename number>
  const bool NumberTraits<number>::is_complex;

  template <typename number>
  const bool NumberTraits<std::complex<number>>::is_complex;

  // explicit instantiations
  template struct NumberTraits<double>;
  template struct NumberTraits<float>;
  template struct NumberTraits<long double>;

  template struct NumberTraits<std::complex<double>>;
  template struct NumberTraits<std::complex<float>>;
  template struct NumberTraits<std::complex<long double>>;
} // namespace numbers

DEAL_II_NAMESPACE_CLOSE
