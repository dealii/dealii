// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2017 by the deal.II authors
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


// check numbers::is_finite for complex arguments

#include <cfenv>
#include <limits>

#include "../tests.h"


template <typename T>
void
check()
{
  using namespace numbers;

  deallog << std::numeric_limits<
               typename numbers::NumberTraits<T>::real_type>::quiet_NaN()
          << "   -->   ";
  deallog << is_finite(
               T(std::numeric_limits<
                   typename numbers::NumberTraits<T>::real_type>::quiet_NaN(),
                 0))
          << ' ';
  deallog << is_finite(
               T(0,
                 std::numeric_limits<
                   typename numbers::NumberTraits<T>::real_type>::quiet_NaN()))
          << std::endl;

  deallog << std::numeric_limits<
               typename numbers::NumberTraits<T>::real_type>::signaling_NaN()
          << "   -->   ";
  deallog << is_finite(T(
               std::numeric_limits<
                 typename numbers::NumberTraits<T>::real_type>::signaling_NaN(),
               0))
          << ' ';
  deallog
    << is_finite(
         T(0,
           std::numeric_limits<
             typename numbers::NumberTraits<T>::real_type>::signaling_NaN()))
    << std::endl;

  deallog << std::numeric_limits<
               typename numbers::NumberTraits<T>::real_type>::infinity()
          << "   -->   ";
  deallog << is_finite(
               T(std::numeric_limits<
                   typename numbers::NumberTraits<T>::real_type>::infinity(),
                 0))
          << ' ';
  deallog << is_finite(
               T(0,
                 std::numeric_limits<
                   typename numbers::NumberTraits<T>::real_type>::infinity()))
          << std::endl;

  deallog << -std::numeric_limits<
               typename numbers::NumberTraits<T>::real_type>::infinity()
          << "   -->   ";
  deallog << is_finite(
               T(-std::numeric_limits<
                   typename numbers::NumberTraits<T>::real_type>::infinity(),
                 0))
          << ' ';
  deallog << is_finite(
               T(0,
                 -std::numeric_limits<
                   typename numbers::NumberTraits<T>::real_type>::infinity()))
          << std::endl;

  deallog
    << std::numeric_limits<typename numbers::NumberTraits<T>::real_type>::min()
    << "   -->   ";
  deallog << is_finite(T(std::numeric_limits<
                           typename numbers::NumberTraits<T>::real_type>::min(),
                         0))
          << ' ';
  deallog << is_finite(
               T(0,
                 std::numeric_limits<
                   typename numbers::NumberTraits<T>::real_type>::min()))
          << std::endl;

  deallog
    << std::numeric_limits<typename numbers::NumberTraits<T>::real_type>::max()
    << "   -->   ";
  deallog << is_finite(T(std::numeric_limits<
                           typename numbers::NumberTraits<T>::real_type>::max(),
                         0))
          << ' ';
  deallog << is_finite(
               T(0,
                 std::numeric_limits<
                   typename numbers::NumberTraits<T>::real_type>::max()))
          << std::endl;

  deallog << static_cast<typename numbers::NumberTraits<T>::real_type>(1)
          << "   -->   ";
  deallog << is_finite(T(
               static_cast<typename numbers::NumberTraits<T>::real_type>(1), 0))
          << ' ';
  deallog << is_finite(T(
               0, static_cast<typename numbers::NumberTraits<T>::real_type>(1)))
          << std::endl;

  deallog << static_cast<typename numbers::NumberTraits<T>::real_type>(-1)
          << "   -->   ";
  deallog << is_finite(
               T(static_cast<typename numbers::NumberTraits<T>::real_type>(-1),
                 0))
          << ' ';
  deallog << is_finite(
               T(0,
                 static_cast<typename numbers::NumberTraits<T>::real_type>(-1)))
          << std::endl;
}


int
main()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);

  // the isnan() function (which we call in is_finite()) helpfully
  // produces a floating point exception when called with a signalling
  // NaN if FP exceptions are on. this of course makes it completely
  // unusable since we can no longer detect whether something is a
  // NaN. that said, to make the test work in these cases, simply
  // switch off floating point exceptions for invalid arguments
#if defined(DEAL_II_HAVE_FP_EXCEPTIONS)
  fedisableexcept(FE_INVALID);
#endif

  check<std::complex<double>>();
  check<std::complex<long double>>();

  return 0;
}
