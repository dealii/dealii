// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test for AlignedVector<VectorizedArray<double> >

#include <deal.II/base/aligned_vector.h>
#include <deal.II/base/vectorization.h>

#include "../tests.h"

void
test()
{
  using vector_t               = VectorizedArray<double>;
  const unsigned int n_vectors = VectorizedArray<double>::size();
  using VEC                    = AlignedVector<vector_t>;
  std::vector<double> a_ref(4), b_ref;
  VEC                 a(4);
  deallog << "Constructor: ";
  for (unsigned int i = 0; i < a.size(); ++i)
    for (unsigned int d = 0; d < n_vectors; ++d)
      AssertThrow(a[i][d] == 0, ExcInternalError());
  deallog << "OK" << std::endl;

  a[2] = 1;
  a.push_back(make_vectorized_array(5.));
  a.push_back(make_vectorized_array(42.));
  a_ref[2] = 1;
  a_ref.push_back(5);
  a_ref.push_back(42);

  VEC b(a);
  b.push_back(make_vectorized_array(27.));
  a.insert_back(b.begin(), b.end());

  // keep reference in sync
  b_ref = a_ref;
  b_ref.push_back(27);
  a_ref.insert(a_ref.end(), b_ref.begin(), b_ref.end());

  deallog << "Insertion: ";
  for (unsigned int i = 0; i < a.size(); ++i)
    for (unsigned int d = 0; d < n_vectors; ++d)
      AssertThrow(a[i][d] == a_ref[i], ExcInternalError());
  deallog << "OK" << std::endl;

  a.resize(4);
  deallog << "Shrinking: ";
  for (unsigned int i = 0; i < a.size(); ++i)
    for (unsigned int d = 0; d < n_vectors; ++d)
      AssertThrow(a[i][d] == a_ref[i], ExcInternalError());
  deallog << "OK" << std::endl;

  a.reserve(100);
  deallog << "Reserve: ";
  for (unsigned int i = 0; i < a.size(); ++i)
    for (unsigned int d = 0; d < n_vectors; ++d)
      AssertThrow(a[i][d] == a_ref[i], ExcInternalError());
  deallog << "OK" << std::endl;

  a     = b;
  a_ref = b_ref;
  deallog << "Assignment: ";
  for (unsigned int i = 0; i < a.size(); ++i)
    for (unsigned int d = 0; d < n_vectors; ++d)
      AssertThrow(a[i][d] == a_ref[i], ExcInternalError());
  deallog << "OK" << std::endl;

  // check setting elements for large vectors
  a.resize(0);
  a.resize(100000, make_vectorized_array(1.));
  deallog << "Check large initialization: ";
  for (unsigned int i = 0; i < 100000; ++i)
    for (unsigned int d = 0; d < n_vectors; ++d)
      AssertThrow(a[i][d] == 1., ExcInternalError());
  deallog << "OK" << std::endl;

  // check resize for large vectors
  deallog << "Check large resize: ";
  a.resize(200000, make_vectorized_array(2.));
  a.resize(400000);
  for (unsigned int i = 0; i < 100000; ++i)
    for (unsigned int d = 0; d < n_vectors; ++d)
      AssertThrow(a[i][d] == 1., ExcInternalError());
  for (unsigned int i = 100000; i < 200000; ++i)
    for (unsigned int d = 0; d < n_vectors; ++d)
      AssertThrow(a[i][d] == 2., ExcInternalError());
  for (unsigned int i = 200000; i < 400000; ++i)
    for (unsigned int d = 0; d < n_vectors; ++d)
      AssertThrow(a[i][d] == 0., ExcInternalError());
  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  test();
}
