//----------------------------  log_nan.cc  ---------------------------
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
//----------------------------  log_nan.cc  ---------------------------

// check for a bug reported by Luca Heltai 2006-03-07 on the mailing
// list. the test should actually output "nan", but prints "0"

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <fstream>
#include <iomanip>
#include <limits>

int main ()
{
  std::ofstream logfile("log_nan/output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  deallog << std::numeric_limits<double>::quiet_NaN() << std::endl;
  deallog << std::numeric_limits<double>::signaling_NaN() << std::endl;

  return 0;
}

