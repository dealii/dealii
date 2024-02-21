// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test operator[], norm and norm_square of cuda_tensor.

#include <deal.II/base/tensor.h>

#include "../tests.h"

void
test_cpu()
{
  double             a[3][3] = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  const unsigned int dim     = 3;
  Tensor<2, dim>     t;
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      t[i][j] = a[i][j];

  deallog.push("values");
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      deallog << t[i][j] << std::endl;
  deallog.pop();

  deallog << "norm: " << t.norm() << std::endl;
  deallog << "norm_square: " << t.norm_square() << std::endl;
}

void
test_gpu()
{
  const unsigned int                                       dim = 3;
  Kokkos::View<double, MemorySpace::Default::kokkos_space> norm_dev("norm_dev");
  double                                                   norm_host;
  Kokkos::View<double, MemorySpace::Default::kokkos_space> norm_square_dev(
    "norm_square_dev");
  double norm_square_host;
  Kokkos::View<Tensor<2, dim>, MemorySpace::Default::kokkos_space> t_dev(
    "t_dev");

  using ExecutionSpace = MemorySpace::Default::kokkos_space::execution_space;
  ExecutionSpace exec;

  // Launch the kernels.
  Kokkos::parallel_for(
    Kokkos::MDRangePolicy<ExecutionSpace, Kokkos::Rank<2>>({{0, 0}},
                                                           {{dim, dim}}),
    KOKKOS_LAMBDA(int i, int j) { t_dev()[i][j] = j + i * dim + 1.; });
  Kokkos::parallel_for(
    Kokkos::RangePolicy<ExecutionSpace>(exec, 0, 1), KOKKOS_LAMBDA(int) {
      norm_dev()        = t_dev().norm();
      norm_square_dev() = t_dev().norm_square();
    });

  // Copy the result to the host
  Kokkos::deep_copy(norm_host, norm_dev);
  Kokkos::deep_copy(norm_square_host, norm_square_dev);

  // Output result
  deallog << "norm GPU: " << norm_host << std::endl;
  deallog << "norm_square GPU: " << norm_square_host << std::endl;
}

int
main()
{
  initlog();

  Kokkos::initialize();

  test_cpu();

  test_gpu();

  Kokkos::finalize();
}
