// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2020 by the deal.II authors
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


// Test that all of the ADOL-C number traits combinations are constructible

#include <deal.II/differentiation/ad/ad_number_types.h>
#include <deal.II/differentiation/ad/adolc_number_types.h>

#include <complex>
#include <fstream>
#include <iomanip>

#include "../tests.h"

namespace AD = Differentiation::AD;

int
main()
{
  initlog();

  // --- Taped ---
  AD::NumberTraits<float, AD::NumberTypes::adolc_taped>();
  AD::NumberTraits<std::complex<float>, AD::NumberTypes::adolc_taped>();
  AD::NumberTraits<double, AD::NumberTypes::adolc_taped>();
  AD::NumberTraits<std::complex<double>, AD::NumberTypes::adolc_taped>();

  AD::ADNumberTraits<adouble>();
  AD::ADNumberTraits<std::complex<adouble>>();
  AD::NumberTraits<adouble, AD::NumberTypes::adolc_taped>();
  AD::NumberTraits<std::complex<adouble>, AD::NumberTypes::adolc_taped>();
  deallog << "Taped OK" << std::endl;

  // --- Tapeless ---
  AD::NumberTraits<float, AD::NumberTypes::adolc_tapeless>();
  AD::NumberTraits<std::complex<float>, AD::NumberTypes::adolc_tapeless>();
  AD::NumberTraits<double, AD::NumberTypes::adolc_tapeless>();
  AD::NumberTraits<std::complex<double>, AD::NumberTypes::adolc_tapeless>();

  AD::ADNumberTraits<adtl::adouble>();
  AD::ADNumberTraits<std::complex<adtl::adouble>>();
  AD::NumberTraits<adtl::adouble, AD::NumberTypes::adolc_tapeless>();
  AD::NumberTraits<std::complex<adtl::adouble>,
                   AD::NumberTypes::adolc_tapeless>();
  deallog << "Tapeless OK" << std::endl;

  deallog << "OK" << std::endl;
}
