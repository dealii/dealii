/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2021 - 2022 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------
 */


// Compute the inner product of two vectors on the GPU using Kokkos.

#include <Kokkos_Core.hpp>

#include "../tests.h"

int
main(int argc, char *argv[])
{
  initlog();
  const int N = 100;

  Kokkos::initialize(argc, argv);
  {
    using MemSpace = Kokkos::CudaSpace;

    using ExecSpace = MemSpace::execution_space;

    using range_policy = Kokkos::RangePolicy<ExecSpace>;

    using ViewVectorType = Kokkos::View<double *, Kokkos::LayoutLeft, MemSpace>;

    ViewVectorType y("y", N);
    ViewVectorType x("x", N);

    // Create host mirrors of device views.
    ViewVectorType::HostMirror h_y = Kokkos::create_mirror_view(y);
    ViewVectorType::HostMirror h_x = Kokkos::create_mirror_view(x);

    // Initialize y vector on host.
    for (int i = 0; i < N; ++i)
      {
        h_y(i) = 1;
      }

    // Initialize x vector on host.
    for (int i = 0; i < N; ++i)
      {
        h_x(i) = 1;
      }

    // Deep copy host views to device views.
    Kokkos::deep_copy(y, h_y);
    Kokkos::deep_copy(x, h_x);

    double result = 0;

    Kokkos::parallel_reduce(
      range_policy(0, N),
      KOKKOS_LAMBDA(int j, double &update) { update += x(j) * y(j); },
      result);

    double h_result = 0.;
    for (int i = 0; i < N; ++i)
      h_result += h_x(i) * h_y(i);

    deallog << " Error: " << result - h_result << std::endl;
  }
  Kokkos::finalize();

  return 0;
}
