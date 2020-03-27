// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2018 by the deal.II authors
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


// Test that all of the ADOL-C number traits return the correct information
// adtl::adouble (Tapeless)

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
  print_info<AD::NumberTraits<float, AD::NumberTypes::adolc_tapeless>>();
  deallog.pop();

  deallog.push("std::complex<float>");
  print_info<
    AD::NumberTraits<std::complex<float>, AD::NumberTypes::adolc_tapeless>>();
  deallog.pop();

  deallog.push("double");
  print_info<AD::NumberTraits<double, AD::NumberTypes::adolc_tapeless>>();
  deallog.pop();

  deallog.push("std::complex<double>");
  print_info<
    AD::NumberTraits<std::complex<double>, AD::NumberTypes::adolc_tapeless>>();
  deallog.pop();

  deallog.push("adouble");
  print_info<
    AD::NumberTraits<adtl::adouble, AD::NumberTypes::adolc_tapeless>>();
  deallog.pop();

  deallog.push("std::complex<adouble>");
  print_info<AD::NumberTraits<std::complex<adtl::adouble>,
                              AD::NumberTypes::adolc_tapeless>>();
  deallog.pop();

  deallog.pop(); // NumberTraits

  deallog << std::endl;

  deallog.push("ADNumberTraits");

  deallog.push("adtl::adouble");
  print_info<AD::ADNumberTraits<adtl::adouble>>();
  deallog.pop();

  deallog.push("std::complex<adtl::adouble>");
  print_info<AD::ADNumberTraits<std::complex<adtl::adouble>>>();
  deallog.pop();

  deallog.pop(); // ADNumberTraits

  deallog << "OK" << std::endl;
}
