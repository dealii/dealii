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


#ifndef __deal2__mg_transfer_component_templates_h
#define __deal2__mg_transfer_component_templates_h

#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/fe/fe.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/multigrid/mg_base.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_transfer_component.h>
#include <deal.II/numerics/data_out.h>

#include <algorithm>
#include <sstream>
#include <fstream>

DEAL_II_NAMESPACE_OPEN

/* --------------------- MGTransferSelect -------------- */



template <typename number>
template <int dim, typename number2, int spacedim>
void
MGTransferSelect<number>::copy_to_mg (
  const DoFHandler<dim,spacedim>        &mg_dof_handler,
  MGLevelObject<Vector<number> > &dst,
  const BlockVector<number2>     &src) const
{
  do_copy_to_mg (mg_dof_handler, dst, src.block(target_component[selected_component]));
}



template <typename number>
template <int dim, typename number2, int spacedim>
void
MGTransferSelect<number>::copy_to_mg (
  const DoFHandler<dim,spacedim>        &mg_dof_handler,
  MGLevelObject<Vector<number> > &dst,
  const Vector<number2>          &src) const
{
  do_copy_to_mg (mg_dof_handler, dst, src);
}



template <typename number>
template <int dim, typename number2, int spacedim>
void
MGTransferSelect<number>::copy_from_mg (
  const DoFHandler<dim,spacedim>              &mg_dof_handler,
  BlockVector<number2>                 &dst,
  const MGLevelObject<Vector<number> > &src) const
{
  dst = 0;
  do_copy_from_mg (mg_dof_handler,
                   dst.block(target_component[selected_component]), src);
  if (constraints != 0)
    constraints->condense(dst);
}



template <typename number>
template <int dim, typename number2, int spacedim>
void
MGTransferSelect<number>::copy_from_mg (
  const DoFHandler<dim,spacedim>              &mg_dof_handler,
  Vector<number2>                      &dst,
  const MGLevelObject<Vector<number> > &src) const
{
  dst = 0;
  do_copy_from_mg (mg_dof_handler, dst, src);
  if (constraints != 0)
    {
      //If we were given constraints
      //apply them to the dst that goes
      //back now to the linear solver.
      //Since constraints are globally
      //defined create a global vector here
      //and copy dst to the right component,
      //apply the constraints then and copy
      //the block back to dst.
      const unsigned int n_blocks =
        *std::max_element(target_component.begin(), target_component.end()) + 1;
      std::vector<types::global_dof_index> dofs_per_block (n_blocks);
      DoFTools::count_dofs_per_block (mg_dof_handler, dofs_per_block, target_component);
      BlockVector<number> tmp;
      tmp.reinit(n_blocks);
      for (unsigned int b=0; b<n_blocks; ++b)
        tmp.block(b).reinit(dofs_per_block[b]);
      tmp.collect_sizes ();
      tmp.block(target_component[selected_component]) = dst;
      constraints->condense(tmp);
      dst = tmp.block(target_component[selected_component]);
    }
}



template <typename number>
template <int dim, typename number2, int spacedim>
void
MGTransferSelect<number>::copy_from_mg_add (
  const DoFHandler<dim,spacedim>              &mg_dof_handler,
  BlockVector<number2>                 &dst,
  const MGLevelObject<Vector<number> > &src) const
{
  do_copy_from_mg_add (mg_dof_handler, dst, src);
}



template <typename number>
template <int dim, typename number2, int spacedim>
void
MGTransferSelect<number>::copy_from_mg_add (
  const DoFHandler<dim,spacedim>              &mg_dof_handler,
  Vector<number2>                      &dst,
  const MGLevelObject<Vector<number> > &src) const
{
  do_copy_from_mg_add (mg_dof_handler, dst, src);
}



template <typename number>
template <int dim, class OutVector, int spacedim>
void
MGTransferSelect<number>::do_copy_from_mg (
  const DoFHandler<dim,spacedim>              &mg_dof_handler,
  OutVector                            &dst,
  const MGLevelObject<Vector<number> > &src) const
{
//  const FiniteElement<dim>& fe = mg_dof_handler.get_fe();

//  const unsigned int dofs_per_cell = fe.dofs_per_cell;
//  std::vector<unsigned int> global_dof_indices (dofs_per_cell);
//  std::vector<unsigned int> level_dof_indices (dofs_per_cell);

  typename DoFHandler<dim,spacedim>::active_cell_iterator
  level_cell = mg_dof_handler.begin_active();
  const typename DoFHandler<dim,spacedim>::active_cell_iterator
  endc = mg_dof_handler.end();

  // traverse all cells and copy the
  // data appropriately to the output
  // vector

  // Note that the level is
  // monotonuosly increasing
  dst = 0;
  for (; level_cell != endc; ++level_cell)
    {
      const unsigned int level = level_cell->level();
      typedef std::vector<std::pair<types::global_dof_index, unsigned int> >::const_iterator IT;
      for (IT i=copy_to_and_from_indices[level].begin();
           i != copy_to_and_from_indices[level].end(); ++i)
        dst(i->first) = src[level](i->second);
    }
}


template <typename number>
template <int dim, class OutVector, int spacedim>
void
MGTransferSelect<number>::do_copy_from_mg_add (
  const DoFHandler<dim,spacedim>              &mg_dof_handler,
  OutVector                            &dst,
  const MGLevelObject<Vector<number> > &src) const
{
  const FiniteElement<dim> &fe = mg_dof_handler.get_fe();
  const unsigned int dofs_per_cell = fe.dofs_per_cell;

  std::vector<types::global_dof_index> global_dof_indices (dofs_per_cell);
  std::vector<types::global_dof_index> level_dof_indices (dofs_per_cell);

  typename DoFHandler<dim,spacedim>::active_cell_iterator
  level_cell = mg_dof_handler.begin_active();
  const typename DoFHandler<dim,spacedim>::active_cell_iterator
  endc = mg_dof_handler.end();

  // traverse all cells and copy the
  // data appropriately to the output
  // vector

  // Note that the level is
  // monotonuosly increasing
  dst = 0;
  for (; level_cell != endc; ++level_cell)
    {
      const unsigned int level = level_cell->level();
      typedef std::vector<std::pair<types::global_dof_index, unsigned int> >::const_iterator IT;
      for (IT i=copy_to_and_from_indices[level].begin();
           i != copy_to_and_from_indices[level].end(); ++i)
        dst(i->first) += src[level](i->second);
    }
}


template <typename number>
std::size_t
MGTransferSelect<number>::memory_consumption () const
{
  return sizeof(int) + MGTransferComponentBase::memory_consumption();
}



DEAL_II_NAMESPACE_CLOSE

#endif
