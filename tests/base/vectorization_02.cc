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


// test for AlignedVector<VectorizedArray<double> >

#include "../tests.h"
#include <iomanip>
#include <fstream>
#include <cmath>

#include <deal.II/base/vectorization.h>
#include <deal.II/base/aligned_vector.h>

void test ()
{
  typedef VectorizedArray<double> vector_t;
  const unsigned int n_vectors = VectorizedArray<double>::n_array_elements;
  typedef AlignedVector<vector_t> VEC;
  std::vector<double> a_ref(4), b_ref;
  VEC a(4);
  deallog << "Constructor: ";
  for (unsigned int i=0; i<a.size(); ++i)
    for (unsigned int d=0; d<n_vectors; ++d)
      AssertThrow (a[i][d] == 0, ExcInternalError());
  deallog << "OK" << std::endl;

  a[2] = 1;
  a.push_back (make_vectorized_array(5.));
  a.push_back (make_vectorized_array(42.));
  a_ref[2] = 1;
  a_ref.push_back(5);
  a_ref.push_back(42);

  VEC b (a);
  b.push_back (make_vectorized_array(27.));
  a.insert_back (b.begin(), b.end());

  // keep reference in sync
  b_ref = a_ref;
  b_ref.push_back (27);
  a_ref.insert (a_ref.end(), b_ref.begin(), b_ref.end());

  deallog << "Insertion: ";
  for (unsigned int i=0; i<a.size(); ++i)
    for (unsigned int d=0; d<n_vectors; ++d)
      AssertThrow (a[i][d] == a_ref[i], ExcInternalError());
  deallog << "OK" << std::endl;

  a.resize(4);
  deallog << "Shrinking: ";
  for (unsigned int i=0; i<a.size(); ++i)
    for (unsigned int d=0; d<n_vectors; ++d)
      AssertThrow (a[i][d] == a_ref[i], ExcInternalError());
  deallog << "OK" << std::endl;

  a.reserve(100);
  deallog << "Reserve: ";
  for (unsigned int i=0; i<a.size(); ++i)
    for (unsigned int d=0; d<n_vectors; ++d)
      AssertThrow (a[i][d] == a_ref[i], ExcInternalError());
  deallog << "OK" << std::endl;

  a = b;
  a_ref = b_ref;
  deallog << "Assignment: ";
  for (unsigned int i=0; i<a.size(); ++i)
    for (unsigned int d=0; d<n_vectors; ++d)
      AssertThrow (a[i][d] == a_ref[i], ExcInternalError());
  deallog << "OK" << std::endl;

  // check setting elements for large vectors
  a.resize (0);
  a.resize (100000, make_vectorized_array(1.));
  deallog << "Check large initialization: ";
  for (unsigned int i=0; i<100000; ++i)
    for (unsigned int d=0; d<n_vectors; ++d)
      AssertThrow (a[i][d] == 1., ExcInternalError());
  deallog << "OK" << std::endl;

  // check resize for large vectors
  deallog << "Check large resize: ";
  a.resize (200000, make_vectorized_array(2.));
  a.resize (400000);
  for (unsigned int i=0; i<100000; ++i)
    for (unsigned int d=0; d<n_vectors; ++d)
      AssertThrow (a[i][d] == 1., ExcInternalError());
  for (unsigned int i=100000; i<200000; ++i)
    for (unsigned int d=0; d<n_vectors; ++d)
      AssertThrow (a[i][d] == 2., ExcInternalError());
  for (unsigned int i=200000; i<400000; ++i)
    for (unsigned int d=0; d<n_vectors; ++d)
      AssertThrow (a[i][d] == 0., ExcInternalError());
  deallog << "OK" << std::endl;
}




int main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();
}
