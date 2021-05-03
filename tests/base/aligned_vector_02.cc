// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2018 by the deal.II authors
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


// test for AlignedVector<AlignedVector<unsigned int> >

#include <deal.II/base/aligned_vector.h>

#include "../tests.h"


using VEC    = AlignedVector<unsigned int>;
using VECVEC = AlignedVector<VEC>;
void
print_vec(VECVEC &v)
{
  for (unsigned int i = 0; i < v.size(); ++i)
    {
      deallog << "[";
      for (unsigned int j = 0; j < v[i].size(); ++j)
        deallog << v[i][j] << " ";
      deallog << "]";
    }
  deallog << std::endl;
}

void
test()
{
  using VEC = AlignedVector<unsigned int>;
  VEC a(4);
  a[0] = 2;
  a[1] = 1;
  a[2] = 42;
  VECVEC v(2);
  deallog << "Constructor: ";
  print_vec(v);

  v[0] = a;
  v[1] = a;

  deallog << "Assignment: ";
  print_vec(v);

  VECVEC w(v);
  deallog << "Assignment vector: ";
  print_vec(w);
  deallog << "Data consistency after assignment: ";
  print_vec(v);

  a[1] = 41;
  a.push_back(100);
  v.push_back(a);
  deallog << "Insertion: ";
  print_vec(v);

  v.resize(1);
  deallog << "Shrinking: ";
  print_vec(v);

  v.reserve(100);
  deallog << "Reserve: ";
  print_vec(v);

  v.resize(10);
  deallog << "Resize: ";
  print_vec(v);

  v.resize(0);
  deallog << "Clear: ";
  print_vec(v);
}



int
main()
{
  initlog();

  test();
}
