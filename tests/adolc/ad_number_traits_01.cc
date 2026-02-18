// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2017 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


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
