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


// Tests Householder class for QR-decomposition

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/householder.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <fstream>
#include <iomanip>

const double rect[] =
{
  4., 3., 2., 1.,
  5., 8., 1., -2.,
  11., 13., -4., -5
};


int main()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);

  FullMatrix<double> A(4,3,rect);
  Householder<double> H(A);

  Vector<double> u(4);
  Vector<double> v1(3);
  Vector<double> v2(3);

  for (unsigned int i=0; i<u.size(); ++i)
    u(i) = i*i;
  deallog << "Distance " << H.least_squares(v1,u) << std::endl;
}
