//----------------------------  print_formatted_ez_01.cc,v  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005, 2007, 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  print_formatted_ez_01.cc,v  ---------------------------

// Check SparseMatrixEZ::print_formatted

#include "../tests.h"

#include <deal.II/base/logstream.h>

#include <deal.II/lac/sparse_matrix_ez.h>

#include <iomanip>
#include <fstream>

int main()
{
  std::ofstream logfile("print_formatted_ez_01/output");
  logfile.setf(std::ios::fixed);
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  
  SparseMatrixEZ<double> ez(5,4);
  ez.set(0,0,2.);
  ez.set(0,2,3.);
  ez.set(0,3,4.);
  ez.set(1,0,5.);
  ez.set(1,1,6.);
  ez.set(1,3,7.);
  ez.set(2,0,8.);
  ez.set(2,1,9.);
  ez.set(2,2,10.);
  ez.set(2,3,11.);
  ez.set(4,0,12.);
  ez.set(4,2,13.);
  ez.set(4,3,14.);

  ez.print_formatted(logfile, 0, false, 5, "~");
}
