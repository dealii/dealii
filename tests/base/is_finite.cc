// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2021 by the deal.II authors
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
