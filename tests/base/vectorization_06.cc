// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test gather and scatter operations of vectorized array

#include <deal.II/base/vectorization.h>

#include <deal.II/lac/vector.h>

#include <limits>

#include "../tests.h"


template <typename Number>
void
test()
{
  // since the number of array elements is system dependent, it is not a good
  // idea to print them to an output file. Instead, check the values manually
  std::vector<Number> vec(200);
  for (unsigned int i = 0; i < vec.size(); ++i)
    vec[i] = i + 1;

  const unsigned int n_vectors = VectorizedArray<Number>::size();
  unsigned int       indices[n_vectors];
  for (unsigned int i = 0; i < n_vectors; ++i)
    indices[i] = i;
  VectorizedArray<Number> arr;
  arr = Number(-1000);

  arr.gather(&vec[0], indices);
  unsigned int n_errors = 0;
  for (unsigned int i = 0; i < n_vectors; ++i)
    if (arr[i] != i + 1)
      ++n_errors;
  deallog << "gather contiguous #errors:      " << n_errors << std::endl;
  if (n_errors > 0)
    {
      for (unsigned int j = 0; j < n_vectors; ++j)
        deallog << arr[j] << " vs " << j + 1 << "   ";
      deallog << std::endl;
    }

  arr.gather(&vec[2], indices);
  n_errors = 0;
  for (unsigned int i = 0; i < n_vectors; ++i)
    if (arr[i] != i + 3)
      ++n_errors;
  deallog << "gather contiguous #errors:      " << n_errors << std::endl;
  if (n_errors > 0)
    {
      for (unsigned int j = 0; j < n_vectors; ++j)
        deallog << arr[j] << " vs " << j + 3 << "   ";
      deallog << std::endl;
    }

  for (unsigned int i = 0; i < n_vectors; ++i)
    indices[i] = 3 * i + 1;
  arr.gather(&vec[0], indices);
  n_errors = 0;
  for (unsigned int i = 0; i < n_vectors; ++i)
    if (arr[i] != 3 * i + 2)
      ++n_errors;
  deallog << "gather non-contiguous #errors:  " << n_errors << std::endl;
  if (n_errors > 0)
    {
      for (unsigned int j = 0; j < n_vectors; ++j)
        deallog << arr[j] << " vs " << 3 * j + 2 << "   ";
      deallog << std::endl;
    }

  arr.gather(&vec[3], indices);
  n_errors = 0;
  for (unsigned int i = 0; i < n_vectors; ++i)
    if (arr[i] != 3 * i + 5)
      ++n_errors;
  deallog << "gather non-contiguous #errors:  " << n_errors << std::endl;
  if (n_errors > 0)
    {
      for (unsigned int j = 0; j < n_vectors; ++j)
        deallog << arr[j] << " vs " << 3 * j + 5 << "   ";
      deallog << std::endl;
    }

  arr = Number(-1);
  arr.scatter(indices, &vec[0]);
  n_errors = 0;
  for (unsigned int i = 0; i < n_vectors; ++i)
    if (vec[indices[i]] != -1)
      ++n_errors;
  unsigned int *start = &indices[0], *end = start + n_vectors;
  for (unsigned int i = 0; i < vec.size(); ++i)
    if (std::find(start, end, i) == end && vec[i] < 0)
      ++n_errors;
  deallog << "scatter non-contiguous #errors: " << n_errors << std::endl;
  if (n_errors > 0)
    {
      for (unsigned int i = 0; i < vec.size(); ++i)
        deallog << vec[i] << ' ';
      deallog << std::endl;
    }
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

  // test long double: in that case, the default
  // path of VectorizedArray is taken no matter
  // what was done for double or float
  deallog.push("long double");
  test<long double>();
  deallog.pop();
}
