// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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

#include <deal.II/lac/affine_constraints.templates.h>

DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace AffineConstraintsImplementation
  {
    template void
    set_zero_all(
      const std::vector<types::global_dof_index> &                  cm,
      LinearAlgebra::distributed::Vector<float, MemorySpace::CUDA> &vec);

    template void
    set_zero_all(
      const std::vector<types::global_dof_index> &                   cm,
      LinearAlgebra::distributed::Vector<double, MemorySpace::CUDA> &vec);
  } // namespace AffineConstraintsImplementation
} // namespace internal

template void
AffineConstraints<float>::set_zero<
  LinearAlgebra::distributed::Vector<float, MemorySpace::CUDA>>(
  LinearAlgebra::distributed::Vector<float, MemorySpace::CUDA> &) const;

template void
AffineConstraints<double>::set_zero<
  LinearAlgebra::distributed::Vector<double, MemorySpace::CUDA>>(
  LinearAlgebra::distributed::Vector<double, MemorySpace::CUDA> &) const;

DEAL_II_NAMESPACE_CLOSE
