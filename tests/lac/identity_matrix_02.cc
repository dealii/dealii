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


// check conversion constructor from IdentityMatrix to FullMatrix


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/identity_matrix.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <fstream>
#include <iomanip>
#include <cmath>

template<typename number>
void
check_vmult()
{
  FullMatrix<number> M(IdentityMatrix(4));
  Vector<number> u(4);
  Vector<number> v(4);

  for (unsigned int i=0; i<4; ++i)
    u(i) = i+1;

  M.vmult(v,u);
  Assert (v == u, ExcInternalError());
  for (unsigned int i=0; i<v.size(); ++i)
    deallog << ' ' << v(i);
  deallog << std::endl;

  M.vmult_add(v,u);
  v /= 2;
  Assert (v == u, ExcInternalError());
  for (unsigned int i=0; i<v.size(); ++i)
    deallog << ' ' << v(i);
  deallog << std::endl;

  M.Tvmult(v,u);
  Assert (v == u, ExcInternalError());
  for (unsigned int i=0; i<v.size(); ++i)
    deallog << ' ' << v(i);
  deallog << std::endl;

  M.Tvmult_add(v,u);
  v /= 2;
  Assert (v == u, ExcInternalError());
  for (unsigned int i=0; i<v.size(); ++i)
    deallog << ' ' << v(i);
  deallog << std::endl;
}


int main()
{
  std::ofstream logfile("output");
  logfile.setf(std::ios::fixed);
  deallog << std::setprecision(0);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check_vmult<double>();
  check_vmult<float>();
}
