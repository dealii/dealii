//-----------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2012 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------------------------

// test for AlignedVector<unsigned int> which tests the basic stuff in the
// aligned vector

#include "../tests.h"
#include <iomanip>
#include <fstream>
#include <cmath>

#include <deal.II/base/vectorization.h>


void test ()
{
  typedef AlignedVector<unsigned int> VEC;
  VEC a(4);
  deallog << "Constructor: ";
  for (unsigned int i=0; i<a.size(); ++i)
    deallog << a[i] << " ";
  deallog << std::endl;

  a[2] = 1;
  a.push_back (5);
  a.push_back (42);

  VEC b (a);
  b.push_back (27);
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
  a.resize (100000, 1);
  deallog << "Check large initialization: ";
  for (unsigned int i=0; i<100000; ++i)
    AssertDimension (a[i], 1);
  deallog << "OK" << std::endl;

				// check resize for large vectors
  deallog << "Check large resize: ";
  a.resize (200000, 2);
  a.resize (400000);
  for (unsigned int i=0; i<100000; ++i)
    AssertDimension (a[i], 1);
  for (unsigned int i=100000; i<200000; ++i)
    AssertDimension (a[i], 2);
  for (unsigned int i=200000; i<400000; ++i)
    AssertDimension (a[i], 0);
  deallog << "OK" << std::endl;
}




int main()
{
  std::ofstream logfile("vectorization_01/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();
}
