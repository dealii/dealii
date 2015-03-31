// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2015 by the deal.II authors
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


// check TableIndices in various ways

#include "../tests.h"
#include <deal.II/base/tensor.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <fstream>
#include <iomanip>

int main ()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  const TableIndices<2> t1(84,42);
  TableIndices<2> t2;
  t2[0] = 84;
  t2[1] = 42;

  AssertThrow (t1 == t2, ExcInternalError());
  AssertThrow (t1[0] == t2[0], ExcInternalError());
  AssertThrow (t1[1] == t2[1], ExcInternalError());

  AssertThrow (! (t1 != t2), ExcInternalError());

  t2.sort();
  AssertThrow (t1 != t2, ExcInternalError());
  AssertThrow (t1[0] == t2[1], ExcInternalError());
  AssertThrow (t1[1] == t2[0], ExcInternalError());

  deallog << "OK" << std::endl;
}
