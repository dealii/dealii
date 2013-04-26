//----------------------------  full_matrix_print.cc,v  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2007, 2008, 2013 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  full_matrix_print.cc,v  ---------------------------

// verify that FullMatrix::print saves the stream state around the output

#include "../tests.h"
#include <cmath>
#include <fstream>
#include <iomanip>

#include <deal.II/base/logstream.h>
#include <deal.II/lac/full_matrix.h>

const double entries[9] = { 11.1,12.2,13.3,21.456,22.12345678901,23,31,32,33 };


int
main ()
{
  std::ofstream logfile("full_matrix_print/output");
  deallog << std::fixed;
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  logfile << numbers::PI << std::endl;
  
  FullMatrix<double> T(3,3,entries);

  // try writing to a regular stream
  T.print(logfile, 15, 8);

  // print something normal -- should use the old settings, 6 digits of
  // precision
  logfile << numbers::PI << std::endl;

  // now try the same with a LogStream
  T.print(deallog, 15, 8);

  logfile << numbers::PI << std::endl;
}

      
