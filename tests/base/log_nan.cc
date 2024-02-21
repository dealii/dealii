// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check for a bug reported by Luca Heltai 2006-03-07 on the mailing
// list. the test should actually output "nan", but prints "0"

#include <cfenv>
#include <limits>

#include "../tests.h"

int
main()
{
  initlog();
  deallog << std::setprecision(3);

  // the isnan() function (which we call in is_finite()) helpfully
  // produces a floating point exception when called with a signalling
  // NaN if FP exceptions are on. this of course makes it completely
  // unusable since we can no longer detect whether something is a
  // NaN. that said, to make the test work in these cases, simply
  // switch off floating point exceptions for invalid arguments
#ifdef DEAL_II_HAVE_FP_EXCEPTIONS
  fedisableexcept(FE_INVALID);
#endif

  deallog << std::numeric_limits<double>::quiet_NaN() << std::endl;
  deallog << std::numeric_limits<double>::signaling_NaN() << std::endl;

  return 0;
}
