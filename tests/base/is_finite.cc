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


// check numbers::is_finite

#include <cfenv>
#include <limits>

#include "../tests.h"


template <typename T>
void
check()
{
  using namespace numbers;

  deallog << std::numeric_limits<T>::quiet_NaN() << "   -->   ";
  deallog << is_finite(std::numeric_limits<T>::quiet_NaN()) << std::endl;

  deallog << std::numeric_limits<T>::signaling_NaN() << "   -->   ";
  deallog << is_finite(std::numeric_limits<T>::signaling_NaN()) << std::endl;

  deallog << std::numeric_limits<T>::infinity() << "   -->   ";
  deallog << is_finite(std::numeric_limits<T>::infinity()) << std::endl;

  deallog << -std::numeric_limits<T>::infinity() << "   -->   ";
  deallog << is_finite(-std::numeric_limits<T>::infinity()) << std::endl;

  deallog << std::numeric_limits<T>::min() << "   -->   ";
  deallog << is_finite(std::numeric_limits<T>::min()) << std::endl;

  deallog << std::numeric_limits<T>::max() << "   -->   ";
  deallog << is_finite(std::numeric_limits<T>::max()) << std::endl;

  deallog << static_cast<T>(1) << "   -->   ";
  deallog << is_finite(static_cast<T>(1)) << std::endl;

  deallog << static_cast<T>(-1) << "   -->   ";
  deallog << is_finite(static_cast<T>(-1)) << std::endl;
}


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


  check<double>();
  check<long double>();

  check<int>();
  check<unsigned int>();

  return 0;
}
