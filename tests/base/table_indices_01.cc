//----------------------------  table_indices_01.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2006, 2013 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  table_indices_01.cc  ---------------------------

// check TableIndices in various ways

#include "../tests.h"
#include <deal.II/base/tensor.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <fstream>
#include <iomanip>

int main ()
{
  std::ofstream logfile("table_indices_01/output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  const TableIndices<2> t1(84,42);
  TableIndices<2> t2;
  t2[0] = 84;
  t2[1] = 42;

  Assert (t1 == t2, ExcInternalError());
  Assert (t1[0] == t2[0], ExcInternalError());
  Assert (t1[1] == t2[1], ExcInternalError());

  Assert (! (t1 != t2), ExcInternalError());

  t2.sort();
  Assert (t1 != t2, ExcInternalError());
  Assert (t1[0] == t2[1], ExcInternalError());
  Assert (t1[1] == t2[0], ExcInternalError());

  deallog << "OK" << std::endl;
}
