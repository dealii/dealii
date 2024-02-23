// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test that all of the Sacado number traits combinations are constructible

#include <deal.II/differentiation/ad/ad_number_types.h>
#include <deal.II/differentiation/ad/sacado_number_types.h>

#include "../tests.h"

// #  include <Sacado.hpp>

#include <complex>
#include <fstream>
#include <iomanip>

namespace AD = Differentiation::AD;

int
main()
{
  initlog();

  // --- Sacado::Fad::DFad ---
  AD::NumberTraits<float, AD::NumberTypes::sacado_dfad>();
  AD::NumberTraits<std::complex<float>, AD::NumberTypes::sacado_dfad>();
  AD::NumberTraits<double, AD::NumberTypes::sacado_dfad>();
  AD::NumberTraits<std::complex<double>, AD::NumberTypes::sacado_dfad>();

  AD::ADNumberTraits<Sacado::Fad::DFad<float>>();
  AD::ADNumberTraits<std::complex<Sacado::Fad::DFad<float>>>();
  AD::ADNumberTraits<Sacado::Fad::DFad<double>>();
  AD::ADNumberTraits<std::complex<Sacado::Fad::DFad<double>>>();
  AD::NumberTraits<Sacado::Fad::DFad<float>, AD::NumberTypes::sacado_dfad>();
  AD::NumberTraits<std::complex<Sacado::Fad::DFad<float>>,
                   AD::NumberTypes::sacado_dfad>();
  AD::NumberTraits<Sacado::Fad::DFad<double>, AD::NumberTypes::sacado_dfad>();
  AD::NumberTraits<std::complex<Sacado::Fad::DFad<double>>,
                   AD::NumberTypes::sacado_dfad>();
  deallog << "Sacado::Fad::DFad OK" << std::endl;

  // --- Nested Sacado::Fad::DFad ---
  AD::NumberTraits<float, AD::NumberTypes::sacado_dfad_dfad>();
  AD::NumberTraits<std::complex<float>, AD::NumberTypes::sacado_dfad_dfad>();
  AD::NumberTraits<double, AD::NumberTypes::sacado_dfad_dfad>();
  AD::NumberTraits<std::complex<double>, AD::NumberTypes::sacado_dfad_dfad>();

  AD::ADNumberTraits<Sacado::Fad::DFad<Sacado::Fad::DFad<float>>>();
  AD::ADNumberTraits<
    std::complex<Sacado::Fad::DFad<Sacado::Fad::DFad<float>>>>();
  AD::ADNumberTraits<Sacado::Fad::DFad<Sacado::Fad::DFad<double>>>();
  AD::ADNumberTraits<
    std::complex<Sacado::Fad::DFad<Sacado::Fad::DFad<double>>>>();
  AD::NumberTraits<Sacado::Fad::DFad<Sacado::Fad::DFad<float>>,
                   AD::NumberTypes::sacado_dfad_dfad>();
  AD::NumberTraits<std::complex<Sacado::Fad::DFad<Sacado::Fad::DFad<float>>>,
                   AD::NumberTypes::sacado_dfad_dfad>();
  AD::NumberTraits<Sacado::Fad::DFad<Sacado::Fad::DFad<double>>,
                   AD::NumberTypes::sacado_dfad_dfad>();
  AD::NumberTraits<std::complex<Sacado::Fad::DFad<Sacado::Fad::DFad<double>>>,
                   AD::NumberTypes::sacado_dfad_dfad>();
  deallog << "Nested Sacado::Fad::DFad OK" << std::endl;

  // --- Sacado::Rad::ADvar ---
  AD::NumberTraits<float, AD::NumberTypes::sacado_rad>();
  AD::NumberTraits<double, AD::NumberTypes::sacado_rad>();

  AD::ADNumberTraits<Sacado::Rad::ADvar<float>>();
  AD::ADNumberTraits<Sacado::Rad::ADvar<double>>();
  AD::NumberTraits<Sacado::Rad::ADvar<float>, AD::NumberTypes::sacado_rad>();
  AD::NumberTraits<Sacado::Rad::ADvar<double>, AD::NumberTypes::sacado_rad>();

#ifdef DEAL_II_WITH_TRILINOS_CXX_SUPPORTS_SACADO_COMPLEX_RAD
  AD::NumberTraits<std::complex<float>, AD::NumberTypes::sacado_rad>();
  AD::NumberTraits<std::complex<double>, AD::NumberTypes::sacado_rad>();

  AD::ADNumberTraits<std::complex<Sacado::Rad::ADvar<float>>>();
  AD::ADNumberTraits<std::complex<Sacado::Rad::ADvar<double>>>();
  AD::NumberTraits<std::complex<Sacado::Rad::ADvar<float>>,
                   AD::NumberTypes::sacado_rad>();
  AD::NumberTraits<std::complex<Sacado::Rad::ADvar<double>>,
                   AD::NumberTypes::sacado_rad>();
#endif
  deallog << "Sacado::Rad::ADvar OK" << std::endl;

  // --- Nested Sacado::Rad::ADvar<Sacado::Fad::DFad> ---
  AD::NumberTraits<float, AD::NumberTypes::sacado_rad_dfad>();
  AD::NumberTraits<double, AD::NumberTypes::sacado_rad_dfad>();

  AD::ADNumberTraits<Sacado::Rad::ADvar<Sacado::Fad::DFad<float>>>();
  AD::ADNumberTraits<Sacado::Rad::ADvar<Sacado::Fad::DFad<double>>>();
  AD::NumberTraits<Sacado::Rad::ADvar<Sacado::Fad::DFad<float>>,
                   AD::NumberTypes::sacado_rad_dfad>();
  AD::NumberTraits<Sacado::Rad::ADvar<Sacado::Fad::DFad<double>>,
                   AD::NumberTypes::sacado_rad_dfad>();

#ifdef DEAL_II_WITH_TRILINOS_CXX_SUPPORTS_SACADO_COMPLEX_RAD
  AD::NumberTraits<std::complex<Sacado::Fad::DFad<float>>,
                   AD::NumberTypes::sacado_rad_dfad>();
  AD::NumberTraits<std::complex<Sacado::Fad::DFad<double>>,
                   AD::NumberTypes::sacado_rad_dfad>();

  AD::ADNumberTraits<
    std::complex<Sacado::Rad::ADvar<Sacado::Fad::DFad<float>>>>();
  AD::ADNumberTraits<
    std::complex<Sacado::Rad::ADvar<Sacado::Fad::DFad<double>>>>();
  AD::NumberTraits<std::complex<Sacado::Rad::ADvar<Sacado::Fad::DFad<float>>>,
                   AD::NumberTypes::sacado_rad_dfad>();
  AD::NumberTraits<std::complex<Sacado::Rad::ADvar<Sacado::Fad::DFad<double>>>,
                   AD::NumberTypes::sacado_rad_dfad>();
#endif
  deallog << "Sacado::Rad::ADvar<Sacado::Fad::DFad> OK" << std::endl;
  deallog << "OK" << std::endl;
}
