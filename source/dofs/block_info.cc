// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2009 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/dofs/block_info.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_tools.h>

#include <deal.II/multigrid/mg_tools.h>

DEAL_II_NAMESPACE_OPEN


template <int dim, int spacedim>
void
BlockInfo::initialize(const DoFHandler<dim, spacedim> &dof,
                      bool                             levels_only,
                      bool                             active_only)
{
  Assert(dof.has_hp_capabilities() == false,
         (typename DoFHandler<dim, spacedim>::ExcNotImplementedWithHP()));

  if (!levels_only && dof.has_active_dofs())
    {
      const std::vector<types::global_dof_index> sizes =
        DoFTools::count_dofs_per_fe_block(dof);
      bi_global.reinit(sizes);
    }

  if (!active_only && dof.has_level_dofs())
    {
      std::vector<std::vector<types::global_dof_index>> sizes(
        dof.get_triangulation().n_levels(),
        std::vector<types::global_dof_index>(dof.get_fe().n_blocks()));

      MGTools::count_dofs_per_block(dof, sizes);
      levels.resize(sizes.size());

      for (unsigned int i = 0; i < sizes.size(); ++i)
        levels[i].reinit(sizes[i]);
    }
}



template <int dim, int spacedim>
void
BlockInfo::initialize_local(const DoFHandler<dim, spacedim> &dof)
{
  Assert(dof.has_hp_capabilities() == false,
         (typename DoFHandler<dim, spacedim>::ExcNotImplementedWithHP()));

  const FiniteElement<dim, spacedim>  &fe = dof.get_fe();
  std::vector<types::global_dof_index> sizes(fe.n_blocks());

  base_elements.resize(fe.n_blocks());

  for (unsigned int i = 0; i < base_elements.size(); ++i)
    base_elements[i] = fe.block_to_base_index(i).first;

  local_renumbering.resize(fe.n_dofs_per_cell());
  FETools::compute_block_renumbering(fe, local_renumbering, sizes, false);
  bi_local.reinit(sizes);
}

// explicit instantiations
#include "dofs/block_info.inst"

DEAL_II_NAMESPACE_CLOSE
