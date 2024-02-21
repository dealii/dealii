// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/exceptions.h>
#include <deal.II/base/memory_consumption.h>

#include <deal.II/dofs/dof_faces.h>

DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace DoFHandlerImplementation
  {
    std::size_t
    DoFFaces<1>::memory_consumption() const
    {
      return 0;
    }


    std::size_t
    DoFFaces<2>::memory_consumption() const
    {
      return MemoryConsumption::memory_consumption(lines);
    }


    std::size_t
    DoFFaces<3>::memory_consumption() const
    {
      return (MemoryConsumption::memory_consumption(quads) +
              MemoryConsumption::memory_consumption(lines));
    }
  } // namespace DoFHandlerImplementation
} // namespace internal

DEAL_II_NAMESPACE_CLOSE
