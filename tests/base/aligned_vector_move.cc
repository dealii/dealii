// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2016 by the deal.II authors
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


// test for C++11 move operations on AlignedVector<unsigned int>

#include "../tests.h"
#include <iomanip>
#include <fstream>

#include <deal.II/base/aligned_vector.h>


void test()
{
  typedef AlignedVector<unsigned int> VEC;

  VEC a(4, 2);

  deallog << "Size: " << a.size() << std::endl;

  VEC b = std::move(a);

  deallog << "Size of new VEC: " << b.size() << std::endl;
  deallog << "Size of moved VEC: " << a.size() << std::endl;

  deallog << "Contents of new VEC: ";
  for (unsigned int i = 0; i < b.size(); ++i)
    deallog << b[i] << " ";
  deallog << std::endl;

  a.resize(6, 42);

  deallog << "Size of moved VEC after resize: " << a.size() << std::endl;

  deallog << "Contents of new VEC: ";
  for (unsigned int i = 0; i < a.size(); ++i)
    deallog << a[i] << " ";
  deallog << std::endl;

  a = std::move(b);

  deallog << "Size of move-assigned VEC: " << a.size() << std::endl;
  deallog << "Size of moved VEC: " << b.size() << std::endl;
}


int main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.threshold_double(1.e-10);

  test ();

  return 0;
}
