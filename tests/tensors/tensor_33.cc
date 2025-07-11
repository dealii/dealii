// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test that we can convert a Tensor of a given rank and dimension to a Tensor
// of Tensors of rank-1 and dimension dim.

#include <deal.II/base/table_indices.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>

#include "../tests.h"


template <int rank, int dim>
void
test()
{
  std::vector<double> values(Utilities::fixed_power<rank>(dim), 0.0);
  std::iota(values.begin(), values.end(), 0.0);

  Tensor<rank, dim>  t(make_array_view(values));
  TableIndices<rank> indices;

  Tensor<1, dim, Tensor<rank - 1, dim>> t2(t);

  // Check for Tensor<rank,dim>
  deallog << "Rank: " << std::to_string(rank)
          << ". Dim: " << std::to_string(dim) << ". Tensor1: " << t
          << std::endl;

  // Check for Tensor<1,dim,Tensor<rank-1,dim>>
  deallog << "Rank: " << std::to_string(rank)
          << ". Dim: " << std::to_string(dim) << ". Tensor2: " << t2
          << std::endl;
}

int
main(int argc, char *argv[])
{
  initlog();

  test<1, 1>();
  test<2, 1>();
  test<3, 1>();
  test<1, 2>();
  test<2, 2>();
  test<3, 2>();
  test<1, 3>();
  test<2, 3>();
  test<3, 3>();

  deallog << "OK" << std::endl;

  return 0;
}
