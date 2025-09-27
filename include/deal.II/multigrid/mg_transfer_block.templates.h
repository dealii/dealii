// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2007 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#ifndef dealii_mg_transfer_block_templates_h
#define dealii_mg_transfer_block_templates_h

#include <deal.II/base/config.h>

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe.h>

#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/sparse_matrix.h>

#include <deal.II/multigrid/mg_base.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_transfer_block.h>

#include <algorithm>

DEAL_II_NAMESPACE_OPEN

/* --------------------- MGTransferBlockSelect -------------- */

// Simplify some things below
using IT = std::vector<std::pair<unsigned int, unsigned int>>::const_iterator;


template <typename number>
template <int dim, typename number2, int spacedim>
void
MGTransferBlockSelect<number>::copy_from_mg(
  const DoFHandler<dim, spacedim>     &dof_handler,
  BlockVector<number2>                &dst,
  const MGLevelObject<Vector<number>> &src) const
{
  for (unsigned int level = 0;
       level < dof_handler.get_triangulation().n_levels();
       ++level)
    for (IT i = copy_indices[selected_block][level].begin();
         i != copy_indices[selected_block][level].end();
         ++i)
      dst.block(selected_block)(i->first) = src[level](i->second);
}



template <typename number>
template <int dim, typename number2, int spacedim>
void
MGTransferBlockSelect<number>::copy_from_mg(
  const DoFHandler<dim, spacedim>     &dof_handler,
  Vector<number2>                     &dst,
  const MGLevelObject<Vector<number>> &src) const
{
  for (unsigned int level = 0;
       level < dof_handler.get_triangulation().n_levels();
       ++level)
    for (IT i = copy_indices[selected_block][level].begin();
         i != copy_indices[selected_block][level].end();
         ++i)
      dst(i->first) = src[level](i->second);
}



template <typename number>
template <int dim, typename number2, int spacedim>
void
MGTransferBlockSelect<number>::copy_from_mg_add(
  const DoFHandler<dim, spacedim>     &dof_handler,
  BlockVector<number2>                &dst,
  const MGLevelObject<Vector<number>> &src) const
{
  for (unsigned int level = 0;
       level < dof_handler.get_triangulation().n_levels();
       ++level)
    for (IT i = copy_indices[selected_block][level].begin();
         i != copy_indices[selected_block][level].end();
         ++i)
      dst.block(selected_block)(i->first) += src[level](i->second);
}



template <typename number>
template <int dim, typename number2, int spacedim>
void
MGTransferBlockSelect<number>::copy_from_mg_add(
  const DoFHandler<dim, spacedim>     &dof_handler,
  Vector<number2>                     &dst,
  const MGLevelObject<Vector<number>> &src) const
{
  for (unsigned int level = 0;
       level < dof_handler.get_triangulation().n_levels();
       ++level)
    for (IT i = copy_indices[selected_block][level].begin();
         i != copy_indices[selected_block][level].end();
         ++i)
      dst(i->first) += src[level](i->second);
}



template <typename number>
std::size_t
MGTransferBlockSelect<number>::memory_consumption() const
{
  return sizeof(int) + MGTransferBlockBase::memory_consumption();
}


/* --------------------- MGTransferBlock -------------- */



template <typename number>
template <int dim, typename number2, int spacedim>
void
MGTransferBlock<number>::copy_from_mg(
  const DoFHandler<dim, spacedim>          &dof_handler,
  BlockVector<number2>                     &dst,
  const MGLevelObject<BlockVector<number>> &src) const
{
  for (unsigned int block = 0; block < selected.size(); ++block)
    if (selected[block])
      for (unsigned int level = 0;
           level < dof_handler.get_triangulation().n_levels();
           ++level)
        for (IT i = copy_indices[block][level].begin();
             i != copy_indices[block][level].end();
             ++i)
          dst.block(block)(i->first) =
            src[level].block(mg_block[block])(i->second);
}



template <typename number>
template <int dim, typename number2, int spacedim>
void
MGTransferBlock<number>::copy_from_mg_add(
  const DoFHandler<dim, spacedim>          &dof_handler,
  BlockVector<number2>                     &dst,
  const MGLevelObject<BlockVector<number>> &src) const
{
  for (unsigned int block = 0; block < selected.size(); ++block)
    if (selected[block])
      for (unsigned int level = 0;
           level < dof_handler.get_triangulation().n_levels();
           ++level)
        for (IT i = copy_indices[block][level].begin();
             i != copy_indices[block][level].end();
             ++i)
          dst.block(block)(i->first) +=
            src[level].block(mg_block[block])(i->second);
}

DEAL_II_NAMESPACE_CLOSE

#endif
