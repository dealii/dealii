//----------------------------  number_traits_complex.cc  ---------------------------
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
//----------------------------  number_traits_complex.cc  ---------------------------

// check numbers::NumberTraits for real data types

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <fstream>
#include <iomanip>
#include <limits>
#include <typeinfo>


template <typename number>
void check (const number &x)
{
  deallog << "typeid(x).name() = " << typeid(x).name()
	  << std::endl;

  deallog << "typeid(NumberTraits<number>::real_type).name() = "
	  << typeid(typename numbers::NumberTraits<number>::real_type).name()
	  << std::endl;

  deallog << numbers::NumberTraits<number>::conjugate (x)
	  << std::endl;

  deallog << numbers::NumberTraits<number>::abs_square (x)
	  << std::endl;

  deallog << numbers::NumberTraits<number>::abs (x)
	  << std::endl;
} 



int main ()
{
  std::ofstream logfile("number_traits_complex/output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check (std::complex<float>(1.5, 2.5));
  check (std::complex<float>(-1.5, -2.5));

  check (std::complex<double>(1.5, 2.5));
  check (std::complex<double>(-1.5, -2.5));

  check (std::complex<long double>(1.5, 2.5));
  check (std::complex<long double>(-1.5, -2.5));

  return 0;
}

