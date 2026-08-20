// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2009 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


#include <deal.II/base/config.h>

#include <deal.II/base/exception_macros.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/memory_consumption.h>

#include <deal.II/lac/block_indices.h>

#include <deal.II/meshworker/local_results.h>

#include <Kokkos_Macros.hpp>

#include <cstddef>
#include <vector>

DEAL_II_NAMESPACE_OPEN

namespace MeshWorker
{
  template <typename number>
  void
  LocalResults<number>::reinit(const BlockIndices &bi)
  {
    for (unsigned int i = 0; i < J.size(); ++i)
      J[i] = 0.;
    for (unsigned int i = 0; i < R.size(); ++i)
      R[i].reinit(bi);
    for (unsigned int i = 0; i < M1.size(); ++i)
      M1[i].matrix.reinit(bi.block_size(M1[i].row),
                          bi.block_size(M1[i].column));
    for (unsigned int i = 0; i < M2.size(); ++i)
      M2[i].matrix.reinit(bi.block_size(M2[i].row),
                          bi.block_size(M2[i].column));
    quadrature_data.reset_values();
  }


  template <typename number>
  std::size_t
  LocalResults<number>::memory_consumption() const
  {
    std::size_t mem = sizeof(*this) + MemoryConsumption::memory_consumption(J) +
                      MemoryConsumption::memory_consumption(R) +
                      MemoryConsumption::memory_consumption(M1) +
                      MemoryConsumption::memory_consumption(M2) +
                      MemoryConsumption::memory_consumption(quadrature_data);
    return mem;
  }


  template class LocalResults<float>;
  template class LocalResults<double>;
} // namespace MeshWorker


DEAL_II_NAMESPACE_CLOSE
