// ---------------------------------------------------------------------
//
// Copyright (C) 2016 by the deal.II authors
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


#ifndef dealii__mg_transfer_internal_h
#define dealii__mg_transfer_internal_h

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/multigrid/mg_constrained_dofs.h>


DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace MGTransfer
  {

    /**
     * Internal function for filling the copy indices from global to level
     * indices
     */
    template <int dim, int spacedim>
    void fill_copy_indices(const dealii::DoFHandler<dim,spacedim>                                                  &mg_dof,
                           const MGConstrainedDoFs                                                                 *mg_constrained_dofs,
                           std::vector<std::vector<std::pair<types::global_dof_index, types::global_dof_index> > > &copy_indices,
                           std::vector<std::vector<std::pair<types::global_dof_index, types::global_dof_index> > > &copy_indices_global_mine,
                           std::vector<std::vector<std::pair<types::global_dof_index, types::global_dof_index> > > &copy_indices_level_mine);
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif
