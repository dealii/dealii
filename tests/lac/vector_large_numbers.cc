// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2013 by the deal.II authors
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



#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <cmath>
#include <fstream>
#include <iomanip>



void check_large_numbers()
{
  Vector<float> v(10);
  v(0) = 1e13;
  v(1) = 1e19;
  v(2) = 1e21;
  v(4) = 1.;
  v(7) = 1e20;

  const double correct = std::sqrt(1e26 + 1e38 + 1e42 + 1e40 + 1.);
  Assert (std::abs(v.l2_norm() - correct) < 1e-6*correct, ExcInternalError());

  v = 0.;
  v(5) = 1e-30;
  v(9) = 1e-32;
  const double correct2 = std::sqrt(1e-60 + 1e-64);
  Assert (std::abs(v.l2_norm() - correct2) < 1e-6*correct2, ExcInternalError());

  Vector<double> w(7);
  w(0) = 1e232;
  w(6) = 1e231;
  w(3) = 1e200;
  const double correct3 = std::sqrt(100. + 1.) * 1e231;
  Assert (std::abs(w.l2_norm() - correct3) < 1e-13*correct3, ExcInternalError());

  w = 0;
  w(1) = 1e-302;
  w(2) = 1e-303;
  w(3) = 2e-303;
  w(4) = 3e-303;
  w(5) = -1e-303;
  const double correct4 = std::sqrt(100. + 1. + 4. + 9 + 1.) * 1e-303;
  Assert (std::abs(w.l2_norm() - correct4) < 1e-13*correct4, ExcInternalError());

  const double correct5 = std::pow(1000. + 1. + 8. + 27 + 1., 1./3.) * 1e-303;
  Assert (std::abs(w.lp_norm(3.) - correct5) < 1e-13*correct5, ExcInternalError());
  deallog << "OK" << std::endl;
}



int main()
{
  std::ofstream logfile("output");
  deallog << std::fixed;
  deallog << std::setprecision(2);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check_large_numbers();
}


