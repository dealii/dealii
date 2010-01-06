//---------------------------------------------------------------------------
//    mg_transfer.templates.h,v 1.22 2006/01/29 15:03:55 guido Exp
//    Version:
//
//    Copyright (C) 2003, 2004, 2005, 2006, 2007, 2009, 2010 by the deal.II authors
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

#include <algorithm>

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
  do_copy_to_mg (mg_dof_handler, dst, src, 0);
}



template <typename number>
template <int dim, typename number2, int spacedim>
void
MGTransferSelect<number>::copy_to_mg (
  const MGDoFHandler<dim,spacedim>        &mg_dof_handler,
  MGLevelObject<Vector<number> > &dst,
  const Vector<number2>          &src) const
{
  do_copy_to_mg (mg_dof_handler, dst, src, component_start[selected_component]);
}



template <typename number>
template <int dim, typename number2, int spacedim>
void
MGTransferSelect<number>::copy_from_mg (
  const MGDoFHandler<dim,spacedim>&              mg_dof_handler,
  BlockVector<number2>&                 dst,
  const MGLevelObject<Vector<number> >& src) const
{
  do_copy_from_mg (mg_dof_handler, dst, src, 0);
}



template <typename number>
template <int dim, typename number2, int spacedim>
void
MGTransferSelect<number>::copy_from_mg (
  const MGDoFHandler<dim,spacedim>&              mg_dof_handler,
  Vector<number2>&                      dst,
  const MGLevelObject<Vector<number> >& src) const
{
  do_copy_from_mg (mg_dof_handler, dst, src,
		   component_start[selected_component]);
}



template <typename number>
template <int dim, typename number2, int spacedim>
void
MGTransferSelect<number>::copy_from_mg_add (
  const MGDoFHandler<dim,spacedim>&              mg_dof_handler,
  BlockVector<number2>&                 dst,
  const MGLevelObject<Vector<number> >& src) const
{
  do_copy_from_mg_add (mg_dof_handler, dst, src, 0);
}



template <typename number>
template <int dim, typename number2, int spacedim>
void
MGTransferSelect<number>::copy_from_mg_add (
  const MGDoFHandler<dim,spacedim>&              mg_dof_handler,
  Vector<number2>&                      dst,
  const MGLevelObject<Vector<number> >& src) const
{
  do_copy_from_mg_add (mg_dof_handler, dst, src,
		       component_start[selected_component]);
}



template <typename number>
template <int dim, class OutVector, int spacedim>
void
MGTransferSelect<number>::do_copy_from_mg (
  const MGDoFHandler<dim,spacedim>              &mg_dof_handler,
  OutVector                            &dst,
  const MGLevelObject<Vector<number> > &src,
  const unsigned int offset) const
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
  for (; level_cell != endc; ++level_cell)
    {
      const unsigned int level = level_cell->level();

				       // get the dof numbers of
				       // this cell for the global
				       // and the level-wise
				       // numbering
      level_cell->get_dof_indices (global_dof_indices);
      level_cell->get_mg_dof_indices(level_dof_indices);

				       // copy level-wise data to
				       // global vector
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	{
	  const unsigned int component
	    = mg_target_component[fe.system_to_component_index(i).first];
	if (mg_selected[component])
	  {
	    const unsigned int level_start
	      = mg_component_start[level][component];
	    dst(global_dof_indices[i] - offset)
	      = src[level](level_dof_indices[i]-level_start);
	  }
	}
    }
}



template <typename number>
template <int dim, class OutVector, int spacedim>
void
MGTransferSelect<number>::do_copy_from_mg_add (
  const MGDoFHandler<dim,spacedim>              &mg_dof_handler,
  OutVector                            &dst,
  const MGLevelObject<Vector<number> > &src,
  const unsigned int offset) const
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
  for (; level_cell != endc; ++level_cell)
    {
      const unsigned int level = level_cell->level();

				       // get the dof numbers of
				       // this cell for the global
				       // and the level-wise
				       // numbering
      level_cell->get_dof_indices (global_dof_indices);
      level_cell->get_mg_dof_indices(level_dof_indices);
				       // copy level-wise data to
				       // global vector
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	{
	  const unsigned int component
	    = mg_target_component[fe.system_to_component_index(i).first];
	  if (mg_selected[component])
	    {
	      const unsigned int level_start
		= mg_component_start[level][component];
	      dst(global_dof_indices[i] - offset)
		+= src[level](level_dof_indices[i] - level_start);
	    }
	}
    }
}


template <typename number>
unsigned int
MGTransferSelect<number>::memory_consumption () const
{
  return sizeof(int) + MGTransferComponentBase::memory_consumption();
}



DEAL_II_NAMESPACE_CLOSE

#endif
