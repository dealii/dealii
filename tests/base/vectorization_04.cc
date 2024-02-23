// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test for VectorizedArray::load and VectorizedArray::store

#include <deal.II/base/aligned_vector.h>
#include <deal.II/base/vectorization.h>

#include "../tests.h"

template <typename Number, typename VectorizedArrayType>
void
test()
{
  const unsigned int  n_vectors = VectorizedArrayType::size();
  std::vector<Number> values(n_vectors * 5);
  for (unsigned int i = 0; i < values.size(); ++i)
    values[i] = i;
  AlignedVector<VectorizedArrayType> copied(4);

  // test load operation for all possible values of alignment
  for (unsigned int shift = 0; shift < n_vectors; ++shift)
    {
      for (unsigned int i = 0; i < 4; ++i)
        copied[i].load(&values[i * n_vectors + shift]);
      for (unsigned int i = 0; i < 4; ++i)
        for (unsigned int v = 0; v < n_vectors; ++v)
          AssertThrow(copied[i][v] == values[i * n_vectors + v + shift],
                      ExcInternalError());
    }
  deallog << "load OK" << std::endl;

  // test store operation
  std::vector<Number> stored(n_vectors * 5);
  for (unsigned int shift = 0; shift < n_vectors; ++shift)
    {
      for (unsigned int i = 0; i < 4; ++i)
        {
          VectorizedArrayType tmp;
          tmp.load(&values[i * n_vectors]);
          tmp.store(&stored[i * n_vectors + shift]);
        }
      for (unsigned int i = 0; i < 4 * n_vectors; ++i)
        AssertThrow(stored[i + shift] == i, ExcInternalError());
    }
  deallog << "store OK" << std::endl;
}



int
main()
{
  initlog();

  deallog.push("double <-> VectorizedArray<double, 1>");
  test<double, VectorizedArray<double, 1>>();
  deallog.pop();
  deallog.push("float <-> VectorizedArray<float, 1>");
  test<float, VectorizedArray<float, 1>>();
  deallog.pop();
  deallog.push("float <-> VectorizedArray<double, 1>");
  test<float, VectorizedArray<double, 1>>();
  deallog.pop();

#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 128
  deallog.push("double <-> VectorizedArray<double, 2>");
  test<double, VectorizedArray<double, 2>>();
  deallog.pop();
  deallog.push("float <-> VectorizedArray<float, 4>");
  test<float, VectorizedArray<float, 4>>();
  deallog.pop();
  deallog.push("float <-> VectorizedArray<double, 2>");
  test<float, VectorizedArray<double, 2>>();
  deallog.pop();
#endif

#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 256
  deallog.push("double <-> VectorizedArray<double, 4>");
  test<double, VectorizedArray<double, 4>>();
  deallog.pop();
  deallog.push("float <-> VectorizedArray<float, 8>");
  test<float, VectorizedArray<float, 8>>();
  deallog.pop();
  deallog.push("float <-> VectorizedArray<double, 4>");
  test<float, VectorizedArray<double, 4>>();
  deallog.pop();
#endif

#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 512
  deallog.push("double <-> VectorizedArray<double, 8>");
  test<double, VectorizedArray<double, 8>>();
  deallog.pop();
  deallog.push("float <-> VectorizedArray<float, 16>");
  test<float, VectorizedArray<float, 16>>();
  deallog.pop();
  deallog.push("float <-> VectorizedArray<double, 8>");
  test<float, VectorizedArray<double, 8>>();
  deallog.pop();
#endif

  // test long double and unsigned int: in these cases, the default path of
  // VectorizedArray is taken no matter what was done for double or float
  deallog.push("long double");
  test<long double, VectorizedArray<long double>>();
  deallog.pop();

  deallog.push("unsigned int");
  test<unsigned int, VectorizedArray<unsigned int>>();
  deallog.pop();
}
