// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2013 by the deal.II authors
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

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <fstream>
#include <iomanip>
#include <limits>


template <typename T>
void check ()
{
  using namespace numbers;

  deallog << std::numeric_limits<typename numbers::NumberTraits<T>::real_type>::quiet_NaN() << "   -->   ";
  deallog << is_finite(T(std::numeric_limits<typename numbers::NumberTraits<T>::real_type>::quiet_NaN(), 0)) << ' ';
  deallog << is_finite(T(0,std::numeric_limits<typename numbers::NumberTraits<T>::real_type>::quiet_NaN())) << std::endl;

  deallog << std::numeric_limits<typename numbers::NumberTraits<T>::real_type>::signaling_NaN() << "   -->   ";
  deallog << is_finite(T(std::numeric_limits<typename numbers::NumberTraits<T>::real_type>::signaling_NaN(), 0)) << ' ';
  deallog << is_finite(T(0,std::numeric_limits<typename numbers::NumberTraits<T>::real_type>::signaling_NaN())) << std::endl;

  deallog << std::numeric_limits<typename numbers::NumberTraits<T>::real_type>::infinity() << "   -->   ";
  deallog << is_finite(T(std::numeric_limits<typename numbers::NumberTraits<T>::real_type>::infinity(), 0)) << ' ';
  deallog << is_finite(T(0, std::numeric_limits<typename numbers::NumberTraits<T>::real_type>::infinity())) << std::endl;

  deallog << -std::numeric_limits<typename numbers::NumberTraits<T>::real_type>::infinity() << "   -->   ";
  deallog << is_finite(T(-std::numeric_limits<typename numbers::NumberTraits<T>::real_type>::infinity(), 0)) << ' ';
  deallog << is_finite(T(0, -std::numeric_limits<typename numbers::NumberTraits<T>::real_type>::infinity())) << std::endl;

  deallog << std::numeric_limits<typename numbers::NumberTraits<T>::real_type>::min() << "   -->   ";
  deallog << is_finite(T(std::numeric_limits<typename numbers::NumberTraits<T>::real_type>::min(), 0)) << ' ';
  deallog << is_finite(T(0,std::numeric_limits<typename numbers::NumberTraits<T>::real_type>::min())) << std::endl;

  deallog << std::numeric_limits<typename numbers::NumberTraits<T>::real_type>::max() << "   -->   ";
  deallog << is_finite(T(std::numeric_limits<typename numbers::NumberTraits<T>::real_type>::max(), 0)) << ' ';
  deallog << is_finite(T(0, std::numeric_limits<typename numbers::NumberTraits<T>::real_type>::max())) << std::endl;

  deallog << static_cast<typename numbers::NumberTraits<T>::real_type> (1) << "   -->   ";
  deallog << is_finite(T(static_cast<typename numbers::NumberTraits<T>::real_type> (1), 0)) << ' ';
  deallog << is_finite(T(0, static_cast<typename numbers::NumberTraits<T>::real_type> (1))) << std::endl;

  deallog << static_cast<typename numbers::NumberTraits<T>::real_type> (-1) << "   -->   ";
  deallog << is_finite(T(static_cast<typename numbers::NumberTraits<T>::real_type> (-1), 0)) << ' ';
  deallog << is_finite(T(0, static_cast<typename numbers::NumberTraits<T>::real_type> (-1))) << std::endl;
}


int main ()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check<std::complex<double> > ();
  check<std::complex<long double> > ();

  return 0;
}

