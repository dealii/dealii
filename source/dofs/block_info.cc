// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2013 by the deal.II authors
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

#include <deal.II/dofs/block_info.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/multigrid/mg_dof_handler.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_tools.h>

DEAL_II_NAMESPACE_OPEN


template <int dim, int spacedim>
void
BlockInfo::initialize(const DoFHandler<dim, spacedim> &dof, bool levels_only, bool active_only)
{
  if (!levels_only && dof.has_active_dofs())
    {
      const FiniteElement<dim, spacedim> &fe = dof.get_fe();
      std::vector<types::global_dof_index> sizes(fe.n_blocks());
      DoFTools::count_dofs_per_block(dof, sizes);
      bi_global.reinit(sizes);
    }

  if (!active_only && dof.has_level_dofs())
    {
      std::vector<std::vector<types::global_dof_index> > sizes (dof.get_tria ().n_levels ());

      for (unsigned int i = 0; i < sizes.size (); ++i)
        sizes[i].resize (dof.get_fe ().n_blocks ());

      MGTools::count_dofs_per_block (dof, sizes);
      levels.resize (sizes.size ());

      for (unsigned int i = 0; i < sizes.size (); ++i)
        levels[i].reinit (sizes[i]);
    }
}


template <int dim, int spacedim>
void
BlockInfo::initialize_local(const DoFHandler<dim, spacedim> &dof)
{
  const FiniteElement<dim, spacedim> &fe = dof.get_fe();
  std::vector<types::global_dof_index> sizes(fe.n_blocks());

  base_elements.resize(fe.n_blocks());

  for (unsigned int i=0; i<base_elements.size(); ++i)
    base_elements[i] = fe.block_to_base_index(i).first;

  local_renumbering.resize(fe.n_dofs_per_cell());
  FETools::compute_block_renumbering(fe,
                                     local_renumbering,
                                     sizes, false);
  bi_local.reinit(sizes);
}


// explicit instantiations
#include "block_info.inst"

DEAL_II_NAMESPACE_CLOSE
