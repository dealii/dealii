// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2020 by the deal.II authors
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
#include "dof_objects.inst"
  }
} // namespace internal

DEAL_II_NAMESPACE_CLOSE
