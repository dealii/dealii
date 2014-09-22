// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2014 by the deal.II authors
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

#include <deal.II/meshworker/dof_info.h>
#include <deal.II/base/quadrature_lib.h>

DEAL_II_NAMESPACE_OPEN


namespace MeshWorker
{
  template <int dim, int spacedim, typename number>
  DoFInfo<dim,spacedim,number>::DoFInfo(const BlockInfo &info)
    :
    block_info(&info, typeid(*this).name()),
    level_cell (false)
  {
    indices_by_block.resize(info.local().size());
    for (unsigned int i=0; i<indices_by_block.size(); ++i)
      indices_by_block[i].resize(info.local().block_size(i));
  }


  template <int dim, int spacedim, typename number>
  DoFInfo<dim,spacedim,number>::DoFInfo()
  {}


  template <int dim, int spacedim, typename number>
  void
  DoFInfo<dim,spacedim,number>::set_block_indices()
  {
    for (unsigned int i=0; i<indices.size(); ++i)
      {
        const std::pair<unsigned int, unsigned int>
        bi = block_info->local().global_to_local(this->block_info->renumber(i));
        indices_by_block[bi.first][bi.second] = indices_org[i];
      }
    // Remove this after
    // changing block codes
    for (unsigned int i=0; i<indices.size(); ++i)
      indices[this->block_info->renumber(i)] = indices_org[i];
  }
}


DEAL_II_NAMESPACE_CLOSE

