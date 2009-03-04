//---------------------------------------------------------------------------
//    mg_transfer.templates.h,v 1.22 2006/01/29 15:03:55 guido Exp
//    Version: 
//
//    Copyright (C) 2003, 2004, 2005, 2006, 2007, 2009 by the deal.II authors
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
template <int dim, class InVector, int spacedim>
void
MGTransferSelect<number>::do_copy_to_mg (
  const MGDoFHandler<dim,spacedim>&        mg_dof_handler,
  MGLevelObject<Vector<number> >& dst,
  const InVector&                 src,
  const unsigned int              offset) const
{
  const FiniteElement<dim>& fe = mg_dof_handler.get_fe();
  
  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  
				   // set the elements of the vectors
				   // on all levels to zero
  unsigned int minlevel = dst.get_minlevel();
  unsigned int maxlevel = dst.get_maxlevel();
  
  dst=0;
  
  Assert(sizes.size()==mg_dof_handler.get_tria().n_levels(),
	 ExcMatricesNotBuilt());

  MGTools::reinit_vector_by_components(
    mg_dof_handler, dst, mg_selected, mg_target_component, sizes);
  
  std::vector<unsigned int> global_dof_indices (dofs_per_cell);
  std::vector<unsigned int> level_dof_indices  (dofs_per_cell);
  
				   // Build a vector of the selected
				   // indices, since traversing all
				   // indices on each cell is too
				   // slow.
  std::vector<unsigned int> selected_indices;
  selected_indices.reserve(dofs_per_cell);
  for (unsigned int i=0; i<dofs_per_cell; ++i)
    if (mg_selected[mg_target_component[fe.system_to_component_index(i).first]])
      selected_indices.push_back(i);
  unsigned int selected_component
    = mg_target_component[fe.system_to_component_index(selected_indices[0]).first];

				   // traverse the grid top-down
				   // (i.e. starting with the most
				   // refined grid). this way, we can
				   // always get that part of one
				   // level of the output vector which
				   // corresponds to a region which is
				   // more refined, by restriction of
				   // the respective vector on the
				   // next finer level, which we then
				   // already have built.
  for (int level=maxlevel; level>=static_cast<signed int>(minlevel); --level)
    {

  typename MGDoFHandler<dim,spacedim>::active_cell_iterator
	level_cell = mg_dof_handler.begin_active(level);
      const typename MGDoFHandler<dim,spacedim>::active_cell_iterator
	level_end  = mg_dof_handler.end_active(level);
				       // Compute coarse level right hand side
				       // by restricting from fine level.
      for (; level_cell!=level_end; ++level_cell)
	{
					   // get the dof numbers of
					   // this cell for the global
					   // and the level-wise
					   // numbering
	  level_cell->get_dof_indices(global_dof_indices);
	  level_cell->get_mg_dof_indices (level_dof_indices);

					   // transfer the global
					   // defect in the vector
					   // into the level-wise one
	  const unsigned int level_start
	    = mg_component_start[level][selected_component];
	  const typename std::vector<unsigned int>::const_iterator
	    end = selected_indices.end();
	  
	  for (typename std::vector<unsigned int>::const_iterator
		 i=selected_indices.begin();
	       i != end; ++i)
	    {
	      dst[level](level_dof_indices[*i] - level_start)
		= src(global_dof_indices[*i] - offset);
	    }
	}

				       // for that part of the level
				       // which is further refined:
				       // get the defect by
				       // restriction of the defect on
				       // one level higher
      if (static_cast<unsigned int>(level) < maxlevel)
	{
	  restrict_and_add (level+1, dst[level], dst[level+1]);
	}
    };
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
