// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------



#include <deal.II/lac/la_vector.h>

#include <limits>

#include "../tests.h"



void
check_large_numbers()
{
  LinearAlgebra::Vector<float> v(10);
  v(0) = 1e13;
  v(1) = 1e19;
  v(2) = 1e21;
  v(4) = 1.;
  v(7) = 1e20;

  const double correct = std::sqrt(1e26 + 1e38 + 1e42 + 1e40 + 1.);
  AssertThrow(std::abs(v.l2_norm() - correct) < 1e-6 * correct,
              ExcInternalError());

  for (unsigned int i = 0; i < v.size(); ++i)
    v[i] = (float)0.;
  v(5)                  = 1e-30;
  v(9)                  = 1e-32;
  const double correct2 = std::sqrt(1e-60 + 1e-64);
  AssertThrow(std::abs(v.l2_norm() - correct2) < 1e-6 * correct2,
              ExcInternalError());

  LinearAlgebra::Vector<double> w(7);
  w(0)                  = 1e232;
  w(6)                  = 1e231;
  w(3)                  = 1e200;
  const double correct3 = std::sqrt(100. + 1.) * 1e231;
  AssertThrow(std::abs(w.l2_norm() - correct3) < 1e-13 * correct3,
              ExcInternalError());

  for (unsigned int i = 0; i < w.size(); ++i)
    w[i] = (float)0.;
  w(1)                  = 1e-302;
  w(2)                  = 1e-303;
  w(3)                  = 2e-303;
  w(4)                  = 3e-303;
  w(5)                  = -1e-303;
  const double correct4 = std::sqrt(100. + 1. + 4. + 9 + 1.) * 1e-303;
  AssertThrow(std::abs(w.l2_norm() - correct4) <= 1e-13 * correct4,
              ExcInternalError());

  deallog << "OK" << std::endl;
}



int
main()
{
  std::ofstream logfile("output");
  deallog << std::fixed;
  deallog << std::setprecision(2);
  deallog.attach(logfile);

  check_large_numbers();
}
