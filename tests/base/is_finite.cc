//----------------------------  is_finite.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005, 2006, 2007 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  is_finite.cc  ---------------------------

// check numbers::is_finite

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <fstream>
#include <iomanip>
#include <limits>


template <typename T>
void check ()
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

  deallog << static_cast<T> (1) << "   -->   ";
  deallog << is_finite(static_cast<T> (1)) << std::endl;

  deallog << static_cast<T> (-1) << "   -->   ";
  deallog << is_finite(static_cast<T> (-1)) << std::endl;
}


int main ()
{
  std::ofstream logfile("is_finite/output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check<double> ();
  check<long double> ();

  check<int> ();
  check<unsigned int> ();
  
  return 0;
}

