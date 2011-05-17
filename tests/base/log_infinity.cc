//----------------------------  log_infinity.cc  ---------------------------
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
//----------------------------  log_infinity.cc  ---------------------------

// similar to log_nan, but test for infinities

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <fstream>
#include <iomanip>
#include <limits>

int main ()
{
  std::ofstream logfile("log_infinity/output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  deallog << std::numeric_limits<double>::infinity() << std::endl;
  deallog << -std::numeric_limits<double>::infinity() << std::endl;

  return 0;
}

