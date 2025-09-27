// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test that single contraction of SymmetricTensors and Tensors work
// as expected.

#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>

#include "../tests.h"

template <int rank, int dim, typename NumberType>
void
fill_tensor(SymmetricTensor<rank, dim, NumberType> &t)
{
  for (unsigned int i = 0; i != t.n_independent_components; ++i)
    {
      t.access_raw_entry(i) = i + 1;
    }
}


template <int rank, int dim, typename NumberType>
void
fill_tensor(Tensor<rank, dim, NumberType> &t)
{
  for (unsigned int i = 0; i != t.n_independent_components; ++i)
    {
      t[t.unrolled_to_component_indices(i)] = i + 1;
    }
}


template <typename NumberType1, int rank, int dim, typename NumberType2>
void
test_ST(const SymmetricTensor<rank, dim, NumberType2> &S)
{
  Tensor<1, dim, NumberType1> T1;
  Tensor<2, dim, NumberType1> T2;
  Tensor<3, dim, NumberType1> T3;
  Tensor<4, dim, NumberType1> T4;

  fill_tensor(T1);
  fill_tensor(T2);
  fill_tensor(T3);
  fill_tensor(T4);

  deallog << "T1*S: " << T1 * S << std::endl;
  deallog << "T2*S: " << T2 * S << std::endl;
  deallog << "T3*S: " << T3 * S << std::endl;
  deallog << "T4*S: " << T4 * S << std::endl;

  deallog << "S*T1: " << S * T1 << std::endl;
  deallog << "S*T2: " << S * T2 << std::endl;
  deallog << "S*T3: " << S * T3 << std::endl;
  deallog << "S*T4: " << S * T4 << std::endl;
}


template <typename NumberType1, typename NumberType2>
void
test_ST()
{
  const unsigned int dim = 2;

  deallog.push("rank 2");
  {
    SymmetricTensor<2, dim, NumberType2> S;
    fill_tensor(S);
    test_ST<NumberType1>(S);
  }
  deallog.pop();

  deallog.push("rank 4");
  {
    SymmetricTensor<4, dim, NumberType2> S;
    fill_tensor(S);
    test_ST<NumberType1>(S);
  }
  deallog.pop();
}

int
main()
{
  initlog();

  // Benchmark results for same types
  deallog.push("d,d");
  test_ST<double, double>();
  deallog.pop();

  // Differing number types:
  deallog.push("d,f");
  test_ST<double, float>();
  deallog.pop();
  deallog.push("f,d");
  test_ST<float, double>();
  deallog.pop();

  deallog << "OK" << std::endl;
}
