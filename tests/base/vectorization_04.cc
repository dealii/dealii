// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2014 by the deal.II authors
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


// test for VectorizedArray::load and VectorizedArray::store

#include "../tests.h"
#include <iomanip>
#include <fstream>
#include <cmath>

#include <deal.II/base/vectorization.h>
#include <deal.II/base/aligned_vector.h>

template <typename Number>
void test ()
{
  const unsigned int n_vectors = VectorizedArray<Number>::n_array_elements;
  std::vector<Number> values(n_vectors*5);
  for (unsigned int i=0; i<values.size(); ++i)
    values[i] = i;
  AlignedVector<VectorizedArray<Number> > copied(4);

  // test load operation for all possible values of alignment
  for (unsigned int shift=0; shift < n_vectors; ++shift)
    {
      for (unsigned int i=0; i<4; ++i)
        copied[i].load(&values[i*n_vectors+shift]);
      for (unsigned int i=0; i<4; ++i)
        for (unsigned int v=0; v<n_vectors; ++v)
          AssertThrow(copied[i][v] == values[i*n_vectors+v+shift],
                      ExcInternalError());
    }
  deallog << "load OK" << std::endl;

  // test store operation
  std::vector<Number> stored(n_vectors*5);
  for (unsigned int shift=0; shift < n_vectors; ++shift)
    {
      for (unsigned int i=0; i<4; ++i)
        {
          VectorizedArray<Number> tmp;
          tmp.load(&values[i*n_vectors]);
          tmp.store(&stored[i*n_vectors+shift]);
        }
      for (unsigned int i=0; i<4*n_vectors; ++i)
        AssertThrow(stored[i+shift] == i, ExcInternalError());
    }
  deallog << "store OK" << std::endl;
}




int main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  deallog.push("double");
  test<double> ();
  deallog.pop();
  deallog.push("float");
  test<float> ();
  deallog.pop();

  // test long double and unsigned int: in these cases, the default path of
  // VectorizedArray is taken no matter what was done for double or float
  deallog.push("long double");
  test<long double> ();
  deallog.pop();

  deallog.push("unsigned int");
  test<unsigned int> ();
  deallog.pop();
}
