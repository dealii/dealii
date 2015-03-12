// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2014 by the deal.II authors
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


// Test SparseMatrixEZ::copy_from

#include "../tests.h"

#include <deal.II/base/logstream.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_matrix_ez.h>

#include <fstream>

int main()
{
  std::ofstream logfile("output");
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

  deallog << "SparseMatrixEZ<float>::copy_from  SparseMatrixEZ<double>"
          << std::endl;
  SparseMatrixEZ<float> ff;
  ff.copy_from(ez);
  ff.print_formatted(logfile, 0, false, 5, "~");
}
