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


// Test that all of the ADOL-C number traits return the correct information
// ::adouble (Taped)

#include <deal.II/differentiation/ad/ad_number_types.h>
#include <deal.II/differentiation/ad/adolc_number_types.h>

#include <complex>
#include <fstream>
#include <iomanip>
#include <type_traits>

#include "../tests.h"

namespace AD = Differentiation::AD;

template <typename NumberTraitsType>
void
print_info()
{
  deallog << "type_code: "
          << static_cast<std::underlying_type<AD::NumberTypes>::type>(
               NumberTraitsType::type_code)
          << std::endl;
  deallog << "is_taped: " << NumberTraitsType::is_taped << std::endl;
  deallog << "is_tapeless: " << NumberTraitsType::is_tapeless << std::endl;
  deallog << "is_real_valued: " << NumberTraitsType::is_real_valued
          << std::endl;
  deallog << "is_complex_valued: " << NumberTraitsType::is_complex_valued
          << std::endl;
  deallog << "n_supported_derivative_levels: "
          << NumberTraitsType::n_supported_derivative_levels << std::endl;

  deallog << "is_ad_number: "
          << AD::is_ad_number<typename NumberTraitsType::ad_type>::value
          << std::endl;
  deallog << "is_adolc_number: "
          << AD::is_adolc_number<typename NumberTraitsType::ad_type>::value
          << std::endl;
  deallog << "is_taped_ad_number: "
          << AD::is_taped_ad_number<typename NumberTraitsType::ad_type>::value
          << std::endl;
  deallog
    << "is_tapeless_ad_number: "
    << AD::is_tapeless_ad_number<typename NumberTraitsType::ad_type>::value
    << std::endl;
}

int
main()
{
  initlog();

  deallog.push("NumberTraits");

  deallog.push("float");
  print_info<AD::NumberTraits<float, AD::NumberTypes::adolc_taped>>();
  deallog.pop();

  deallog.push("std::complex<float>");
  print_info<
    AD::NumberTraits<std::complex<float>, AD::NumberTypes::adolc_taped>>();
  deallog.pop();

  deallog.push("double");
  print_info<AD::NumberTraits<double, AD::NumberTypes::adolc_taped>>();
  deallog.pop();

  deallog.push("std::complex<double>");
  print_info<
    AD::NumberTraits<std::complex<double>, AD::NumberTypes::adolc_taped>>();
  deallog.pop();

  deallog.push("adouble");
  print_info<AD::NumberTraits<adouble, AD::NumberTypes::adolc_taped>>();
  deallog.pop();

  deallog.push("std::complex<adouble>");
  print_info<
    AD::NumberTraits<std::complex<adouble>, AD::NumberTypes::adolc_taped>>();
  deallog.pop();

  deallog.pop(); // NumberTraits

  deallog << std::endl;

  deallog.push("ADNumberTraits");

  deallog.push("adouble");
  print_info<AD::ADNumberTraits<adouble>>();
  deallog.pop();

  deallog.push("std::complex<adouble>");
  print_info<AD::ADNumberTraits<std::complex<adouble>>>();
  deallog.pop();

  deallog.pop(); // ADNumberTraits

  deallog << "OK" << std::endl;
}
