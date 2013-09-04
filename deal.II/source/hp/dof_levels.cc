// ---------------------------------------------------------------------
// $Id$
//
// Copyright (C) 2003 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#include <deal.II/base/memory_consumption.h>
#include <deal.II/hp/dof_levels.h>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace hp
  {
    template <int dim>
    std::size_t
    DoFLevel<dim>::memory_consumption () const
    {
      return (MemoryConsumption::memory_consumption (active_fe_indices) +
	      MemoryConsumption::memory_consumption (dof_indices) +
              MemoryConsumption::memory_consumption (dof_offsets));
    }


    // explicit instantiations
    template std::size_t DoFLevel<1>::memory_consumption () const;
    template std::size_t DoFLevel<2>::memory_consumption () const;
    template std::size_t DoFLevel<3>::memory_consumption () const;
  }
}

DEAL_II_NAMESPACE_CLOSE
