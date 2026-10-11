// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2022 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#include <Kokkos_Core.hpp>

struct FillFunctor
{
  KOKKOS_FUNCTION void
  operator()(const int i) const
  {
    m_values(i) = i;
  }
  Kokkos::View<int *> m_values;
};

struct SumFunctor
{
  KOKKOS_FUNCTION void
  operator()(const int i, int &sum) const
  {
    sum += m_values(i);
  }
  Kokkos::View<int *> m_values;
};

int
main(int argc, char *argv[])
{
  const long n = 100;
  int        sum;
  Kokkos::initialize(argc, argv);
  {
    Kokkos::View<int *> values("values", n);
    Kokkos::parallel_for(n, FillFunctor{values});
    Kokkos::parallel_reduce(n, SumFunctor{values}, sum);
  }
  Kokkos::finalize();
  return (sum != (n * (n - 1)) / 2);
}
