//---------------------------------------------------------------------------
//    mg_transfer.templates.h,v 1.22 2006/01/29 15:03:55 guido Exp
//    Version:
//
//    Copyright (C) 2003, 2004, 2005, 2006, 2007, 2009, 2010, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#ifndef __deal2__mg_transfer_component_templates_h
#define __deal2__mg_transfer_component_templates_h

#include <lac/sparse_matrix.h>
#include <grid/tria_iterator.h>
#include <fe/fe.h>
#include <lac/constraint_matrix.h>
#include <multigrid/mg_base.h>
#include <multigrid/mg_dof_accessor.h>
#include <multigrid/mg_tools.h>
#include <multigrid/mg_transfer_component.h>
#include <numerics/data_out.h>

#include <algorithm>
#include <sstream>
#include <fstream>

DEAL_II_NAMESPACE_OPEN

/* --------------------- MGTransferSelect -------------- */



template <typename number>
template <int dim, typename number2, int spacedim>
void
MGTransferSelect<number>::copy_to_mg (
  const MGDoFHandler<dim,spacedim>        &mg_dof_handler,
  MGLevelObject<Vector<number> > &dst,
  const BlockVector<number2>     &src) const
{
  do_copy_to_mg (mg_dof_handler, dst, src.block(target_component[selected_component]));
}



template <typename number>
template <int dim, typename number2, int spacedim>
void
MGTransferSelect<number>::copy_to_mg (
  const MGDoFHandler<dim,spacedim>        &mg_dof_handler,
  MGLevelObject<Vector<number> > &dst,
  const Vector<number2>          &src) const
{
  do_copy_to_mg (mg_dof_handler, dst, src);
}



template <typename number>
template <int dim, typename number2, int spacedim>
void
MGTransferSelect<number>::copy_from_mg (
  const MGDoFHandler<dim,spacedim>&              mg_dof_handler,
  BlockVector<number2>&                 dst,
  const MGLevelObject<Vector<number> >& src) const
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
  const MGDoFHandler<dim,spacedim>&              mg_dof_handler,
  Vector<number2>&                      dst,
  const MGLevelObject<Vector<number> >& src) const
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
    std::vector<unsigned int> dofs_per_block (n_blocks);
    DoFTools::count_dofs_per_block (mg_dof_handler, dofs_per_block, target_component);
    BlockVector<number> tmp;
    tmp.reinit(n_blocks);
    for(unsigned int b=0; b<n_blocks; ++b)
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
  const MGDoFHandler<dim,spacedim>&              mg_dof_handler,
  BlockVector<number2>&                 dst,
  const MGLevelObject<Vector<number> >& src) const
{
  do_copy_from_mg_add (mg_dof_handler, dst, src);
}



template <typename number>
template <int dim, typename number2, int spacedim>
void
MGTransferSelect<number>::copy_from_mg_add (
  const MGDoFHandler<dim,spacedim>&              mg_dof_handler,
  Vector<number2>&                      dst,
  const MGLevelObject<Vector<number> >& src) const
{
  do_copy_from_mg_add (mg_dof_handler, dst, src);
}



template <typename number>
template <int dim, class OutVector, int spacedim>
void
MGTransferSelect<number>::do_copy_from_mg (
  const MGDoFHandler<dim,spacedim>              &mg_dof_handler,
  OutVector                            &dst,
  const MGLevelObject<Vector<number> > &src) const
{
//  const FiniteElement<dim>& fe = mg_dof_handler.get_fe();

//  const unsigned int dofs_per_cell = fe.dofs_per_cell;
//  std::vector<unsigned int> global_dof_indices (dofs_per_cell);
//  std::vector<unsigned int> level_dof_indices (dofs_per_cell);

  typename MGDoFHandler<dim,spacedim>::active_cell_iterator
    level_cell = mg_dof_handler.begin_active();
  const typename MGDoFHandler<dim,spacedim>::active_cell_iterator
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
    typedef std::vector<std::pair<unsigned int, unsigned int> >::const_iterator IT;
    for (IT i=copy_to_and_from_indices[level].begin();
        i != copy_to_and_from_indices[level].end(); ++i)
      dst(i->first) = src[level](i->second);
  }
}


template <typename number>
template <int dim, class OutVector, int spacedim>
void
MGTransferSelect<number>::do_copy_from_mg_add (
  const MGDoFHandler<dim,spacedim>              &mg_dof_handler,
  OutVector                            &dst,
  const MGLevelObject<Vector<number> > &src) const
{
  const FiniteElement<dim>& fe = mg_dof_handler.get_fe();
  const unsigned int dofs_per_cell = fe.dofs_per_cell;

  std::vector<unsigned int> global_dof_indices (dofs_per_cell);
  std::vector<unsigned int> level_dof_indices (dofs_per_cell);

  typename MGDoFHandler<dim,spacedim>::active_cell_iterator
    level_cell = mg_dof_handler.begin_active();
  const typename MGDoFHandler<dim,spacedim>::active_cell_iterator
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
      typedef std::vector<std::pair<unsigned int, unsigned int> >::const_iterator IT;
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
