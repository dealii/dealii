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


// test for AlignedVector<bool>

#include "../tests.h"
#include <iomanip>
#include <fstream>
#include <cmath>

#include <deal.II/base/aligned_vector.h>


void test ()
{
  typedef AlignedVector<bool> VEC;
  VEC a(4);
  deallog << "Constructor: ";
  for (unsigned int i=0; i<a.size(); ++i)
    deallog << a[i] << " ";
  deallog << std::endl;

  a[2] = true;
  a.push_back (true);
  a.push_back (false);

  VEC b (a);
  b.push_back (true);
  a.insert_back (b.begin(), b.end());

  deallog << "Insertion: ";
  for (unsigned int i=0; i<a.size(); ++i)
    deallog << a[i] << " ";
  deallog << std::endl;

  a.resize(4);
  deallog << "Shrinking: ";
  for (unsigned int i=0; i<a.size(); ++i)
    deallog << a[i] << " ";
  deallog << std::endl;

  a.reserve(100);
  deallog << "Reserve: ";
  for (unsigned int i=0; i<a.size(); ++i)
    deallog << a[i] << " ";
  deallog << std::endl;

  a = b;
  deallog << "Assignment: ";
  for (unsigned int i=0; i<a.size(); ++i)
    deallog << a[i] << " ";
  deallog << std::endl;

  // check setting elements for large vectors
  a.resize (0);
  a.resize (100000, true);
  deallog << "Check large initialization: ";
  for (unsigned int i=0; i<100000; ++i)
    AssertDimension (a[i], true);
  deallog << "OK" << std::endl;

  // check resize for large vectors
  deallog << "Check large resize: ";
  a.resize (200000, false);
  a.resize (400000, true);
  for (unsigned int i=0; i<100000; ++i)
    AssertDimension (a[i], true);
  for (unsigned int i=100000; i<200000; ++i)
    AssertDimension (a[i], false);
  for (unsigned int i=200000; i<400000; ++i)
    AssertDimension (a[i], true);
  deallog << "OK" << std::endl;

  deallog << "Memory consumption: " << a.memory_consumption()
          << " for a.size()=" << a.size() << std::endl;
}




int main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();
}
