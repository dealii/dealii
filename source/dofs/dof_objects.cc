// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2006 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#include <deal.II/base/exceptions.h>
#include <deal.II/base/memory_consumption.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_objects.h>

#include <deal.II/fe/fe.h>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace DoFHandlerImplementation
  {
    template <int dim>
    std::size_t
    DoFObjects<dim>::memory_consumption() const
    {
      return (MemoryConsumption::memory_consumption(dofs));
    }
  } // namespace DoFHandlerImplementation
} // namespace internal


// explicit instantiations
namespace internal
{
  namespace DoFHandlerImplementation
  {
#include "dofs/dof_objects.inst"
  }
} // namespace internal

DEAL_II_NAMESPACE_CLOSE
