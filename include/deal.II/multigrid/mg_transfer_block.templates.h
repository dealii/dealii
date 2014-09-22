// ---------------------------------------------------------------------
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


#ifndef __deal2__mg_transfer_block_templates_h
#define __deal2__mg_transfer_block_templates_h

#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/fe/fe.h>
#include <deal.II/multigrid/mg_base.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_transfer_block.h>

#include <algorithm>

DEAL_II_NAMESPACE_OPEN

/* --------------------- MGTransferBlockSelect -------------- */

// Simplify some things below
typedef std::vector<std::pair<unsigned int, unsigned int> >::const_iterator IT;


template <typename number>
template <int dim, typename number2, int spacedim>
void
MGTransferBlockSelect<number>::copy_from_mg (
  const DoFHandler<dim,spacedim>              &mg_dof_handler,
  BlockVector<number2>                 &dst,
  const MGLevelObject<Vector<number> > &src) const
{
  for (unsigned int level=0; level<mg_dof_handler.get_tria().n_levels(); ++level)
    for (IT i= copy_indices[selected_block][level].begin();
         i != copy_indices[selected_block][level].end(); ++i)
      dst.block(selected_block)(i->first) = src[level](i->second);
}



template <typename number>
template <int dim, typename number2, int spacedim>
void
MGTransferBlockSelect<number>::copy_from_mg (
  const DoFHandler<dim,spacedim>              &mg_dof_handler,
  Vector<number2>                      &dst,
  const MGLevelObject<Vector<number> > &src) const
{
  for (unsigned int level=0; level<mg_dof_handler.get_tria().n_levels(); ++level)
    for (IT i= copy_indices[selected_block][level].begin();
         i != copy_indices[selected_block][level].end(); ++i)
      dst(i->first) = src[level](i->second);
}



template <typename number>
template <int dim, typename number2, int spacedim>
void
MGTransferBlockSelect<number>::copy_from_mg_add (
  const DoFHandler<dim,spacedim>              &mg_dof_handler,
  BlockVector<number2>                 &dst,
  const MGLevelObject<Vector<number> > &src) const
{
  for (unsigned int level=0; level<mg_dof_handler.get_tria().n_levels(); ++level)
    for (IT i= copy_indices[selected_block][level].begin();
         i != copy_indices[selected_block][level].end(); ++i)
      dst.block(selected_block)(i->first) += src[level](i->second);
}



template <typename number>
template <int dim, typename number2, int spacedim>
void
MGTransferBlockSelect<number>::copy_from_mg_add (
  const DoFHandler<dim,spacedim>              &mg_dof_handler,
  Vector<number2>                      &dst,
  const MGLevelObject<Vector<number> > &src) const
{
  for (unsigned int level=0; level<mg_dof_handler.get_tria().n_levels(); ++level)
    for (IT i= copy_indices[selected_block][level].begin();
         i != copy_indices[selected_block][level].end(); ++i)
      dst(i->first) += src[level](i->second);
}



template <typename number>
std::size_t
MGTransferBlockSelect<number>::memory_consumption () const
{
  return sizeof(int) + MGTransferBlockBase::memory_consumption();
}


/* --------------------- MGTransferBlock -------------- */



template <typename number>
template <int dim, typename number2, int spacedim>
void
MGTransferBlock<number>::copy_from_mg (
  const DoFHandler<dim,spacedim> &mg_dof_handler,
  BlockVector<number2> &dst,
  const MGLevelObject<BlockVector<number> > &src) const
{
  for (unsigned int block=0; block<selected.size(); ++block)
    if (selected[block])
      for (unsigned int level=0; level<mg_dof_handler.get_tria().n_levels(); ++level)
        for (IT i= copy_indices[block][level].begin();
             i != copy_indices[block][level].end(); ++i)
          dst.block(block)(i->first) = src[level].block(mg_block[block])(i->second);
}



template <typename number>
template <int dim, typename number2, int spacedim>
void
MGTransferBlock<number>::copy_from_mg_add (
  const DoFHandler<dim,spacedim> &mg_dof_handler,
  BlockVector<number2> &dst,
  const MGLevelObject<BlockVector<number> > &src) const
{
  for (unsigned int block=0; block<selected.size(); ++block)
    if (selected[block])
      for (unsigned int level=0; level<mg_dof_handler.get_tria().n_levels(); ++level)
        for (IT i= copy_indices[block][level].begin();
             i != copy_indices[block][level].end(); ++i)
          dst.block(block)(i->first) += src[level].block(mg_block[block])(i->second);
}

DEAL_II_NAMESPACE_CLOSE

#endif
