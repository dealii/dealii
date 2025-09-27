// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test that all of the ADOL-C number traits return the correct information

#include <deal.II/differentiation/ad.h>

#include <complex>
#include <fstream>
#include <iomanip>
#include <type_traits>

#include "../tests.h"

namespace AD = Differentiation::AD;

template <typename ADNumber>
void
print()
{
  deallog << "is_ad_number: " << AD::is_ad_number<ADNumber>::value << std::endl;
  deallog << "is_adolc_number: " << AD::is_adolc_number<ADNumber>::value
          << std::endl;
  deallog << "is_taped_ad_number: " << AD::is_taped_ad_number<ADNumber>::value
          << std::endl;
  deallog << "is_tapeless_ad_number: "
          << AD::is_tapeless_ad_number<ADNumber>::value << std::endl;

  Assert(AD::is_ad_number<const ADNumber>::value ==
           AD::is_ad_number<ADNumber>::value,
         ExcMessage("Error const"));
  Assert(AD::is_adolc_number<const ADNumber>::value ==
           AD::is_adolc_number<ADNumber>::value,
         ExcMessage("Error const"));
  Assert(AD::is_taped_ad_number<const ADNumber>::value ==
           AD::is_taped_ad_number<ADNumber>::value,
         ExcMessage("Error const"));
  Assert(AD::is_tapeless_ad_number<const ADNumber>::value ==
           AD::is_tapeless_ad_number<ADNumber>::value,
         ExcMessage("Error const"));

  Assert(AD::is_ad_number<ADNumber &>::value ==
           AD::is_ad_number<ADNumber>::value,
         ExcMessage("Error reference"));
  Assert(AD::is_adolc_number<ADNumber &>::value ==
           AD::is_adolc_number<ADNumber>::value,
         ExcMessage("Error reference"));
  Assert(AD::is_taped_ad_number<ADNumber &>::value ==
           AD::is_taped_ad_number<ADNumber>::value,
         ExcMessage("Error reference"));
  Assert(AD::is_tapeless_ad_number<ADNumber &>::value ==
           AD::is_tapeless_ad_number<ADNumber>::value,
         ExcMessage("Error reference"));

  Assert(AD::is_ad_number<const ADNumber &>::value ==
           AD::is_ad_number<ADNumber>::value,
         ExcMessage("Error const reference"));
  Assert(AD::is_adolc_number<const ADNumber &>::value ==
           AD::is_adolc_number<ADNumber>::value,
         ExcMessage("Error const reference"));
  Assert(AD::is_taped_ad_number<const ADNumber &>::value ==
           AD::is_taped_ad_number<ADNumber>::value,
         ExcMessage("Error const reference"));
  Assert(AD::is_tapeless_ad_number<const ADNumber &>::value ==
           AD::is_tapeless_ad_number<ADNumber>::value,
         ExcMessage("Error const reference"));

  Assert(AD::is_ad_number<ADNumber &&>::value ==
           AD::is_ad_number<ADNumber>::value,
         ExcMessage("Error rvalue"));
  Assert(AD::is_adolc_number<ADNumber &&>::value ==
           AD::is_adolc_number<ADNumber>::value,
         ExcMessage("Error rvalue"));
  Assert(AD::is_taped_ad_number<ADNumber &&>::value ==
           AD::is_taped_ad_number<ADNumber>::value,
         ExcMessage("Error rvalue"));
  Assert(AD::is_tapeless_ad_number<ADNumber &&>::value ==
           AD::is_tapeless_ad_number<ADNumber>::value,
         ExcMessage("Error rvalue"));
}

int
main()
{
  initlog();

  deallog.push("adouble");
  print<::adouble>();
  deallog.pop();

  deallog.push("std::complex<::adouble>");
  print<std::complex<::adouble>>();
  deallog.pop();

  deallog.push("adtl::adouble");
  print<adtl::adouble>();
  deallog.pop();

  deallog.push("std::complex<adtl::adouble>");
  print<std::complex<adtl::adouble>>();
  deallog.pop();

  deallog << "OK" << std::endl;
}
