// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2020 by the deal.II authors
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


// test VectorizedArray::streaming_store

#include <deal.II/base/aligned_vector.h>
#include <deal.II/base/vectorization.h>

#include "../tests.h"

template <typename Number>
void
test()
{
  const unsigned int  n_chunks  = 50000;
  const unsigned int  n_vectors = VectorizedArray<Number>::size();
  std::vector<Number> values(n_vectors * n_chunks);
  for (unsigned int i = 0; i < values.size(); ++i)
    values[i] = i;

  // test store operation
  AlignedVector<Number> stored(n_vectors * n_chunks);
  for (unsigned int i = 0; i < n_chunks; ++i)
    {
      VectorizedArray<Number> tmp;
      tmp.load(&values[i * n_vectors]);
      tmp.streaming_store(&stored[i * n_vectors]);
    }
  for (unsigned int i = 0; i < n_chunks * n_vectors; ++i)
    AssertThrow(stored[i] == i, ExcInternalError());
  deallog << "streaming store OK" << std::endl;
}



int
main()
{
  initlog();

  deallog.push("double");
  test<double>();
  deallog.pop();
  deallog.push("float");
  test<float>();
  deallog.pop();

  // test long double and unsigned int: in these cases, the default path of
  // VectorizedArray is taken no matter what was done for double or float
  deallog.push("long double");
  test<long double>();
  deallog.pop();

  deallog.push("unsigned int");
  test<unsigned int>();
  deallog.pop();
}
