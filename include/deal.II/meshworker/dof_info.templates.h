// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2009 - 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_dof_info_templates_h
#define dealii_dof_info_templates_h


#include <deal.II/base/config.h>

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/meshworker/dof_info.h>

DEAL_II_NAMESPACE_OPEN


namespace MeshWorker
{
  template <int dim, int spacedim, typename number>
  DoFInfo<dim, spacedim, number>::DoFInfo(const BlockInfo &info)
    : face_number(numbers::invalid_unsigned_int)
    , sub_number(numbers::invalid_unsigned_int)
    , block_info(&info, typeid(*this).name())
    , level_cell(false)
  {
    indices_by_block.resize(info.local().size());
    for (unsigned int i = 0; i < indices_by_block.size(); ++i)
      indices_by_block[i].resize(info.local().block_size(i));
  }



  template <int dim, int spacedim, typename number>
  void
  DoFInfo<dim, spacedim, number>::set_block_indices()
  {
    for (unsigned int i = 0; i < indices.size(); ++i)
      {
        const std::pair<unsigned int, unsigned int> bi =
          block_info->local().global_to_local(this->block_info->renumber(i));
        indices_by_block[bi.first][bi.second] = indices_org[i];
      }
    // Remove this after
    // changing block codes
    for (unsigned int i = 0; i < indices.size(); ++i)
      indices[this->block_info->renumber(i)] = indices_org[i];
  }
} // namespace MeshWorker


DEAL_II_NAMESPACE_CLOSE

#endif
