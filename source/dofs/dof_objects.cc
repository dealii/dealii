// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2013 by the deal.II authors
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

#include <deal.II/base/exceptions.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/dofs/dof_objects.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe.h>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace DoFHandler
  {
    template <int dim>
    std::size_t
    DoFObjects<dim>::memory_consumption () const
    {
      return (MemoryConsumption::memory_consumption (dofs));
    }



    template <int dim>
    template <int dh_dim, int spacedim>
    void
    DoFObjects<dim>::
    set_dof_index (const dealii::DoFHandler<dh_dim, spacedim> &dof_handler,
                   const unsigned int       obj_index,
                   const unsigned int       fe_index,
                   const unsigned int       local_index,
                   const types::global_dof_index       global_index)
    {
      Assert ((fe_index == dealii::DoFHandler<dh_dim, spacedim>::default_fe_index),
              ExcMessage ("Only the default FE index is allowed for non-hp DoFHandler objects"));
      Assert (local_index<dof_handler.get_fe().template n_dofs_per_object<dim>(),
              ExcIndexRange (local_index, 0, dof_handler.get_fe().template n_dofs_per_object<dim>()));
      Assert (obj_index * dof_handler.get_fe().template n_dofs_per_object<dim>()+local_index
              <
              dofs.size(),
              ExcInternalError());

      dofs[obj_index * dof_handler.get_fe()
           .template n_dofs_per_object<dim>() + local_index] = global_index;
    }
  }
}


// explicit instantiations
namespace internal
{
  namespace DoFHandler
  {
#include "dof_objects.inst"
  }
}

DEAL_II_NAMESPACE_CLOSE
