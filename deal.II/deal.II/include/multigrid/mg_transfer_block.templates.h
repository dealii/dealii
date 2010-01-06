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

#ifndef __deal2__mg_transfer_block_templates_h
#define __deal2__mg_transfer_block_templates_h

#include <lac/sparse_matrix.h>
#include <lac/constraint_matrix.h>
#include <grid/tria_iterator.h>
#include <fe/fe.h>
#include <multigrid/mg_base.h>
#include <multigrid/mg_dof_accessor.h>
#include <multigrid/mg_tools.h>
#include <multigrid/mg_transfer_block.h>

#include <algorithm>

DEAL_II_NAMESPACE_OPEN

/* --------------------- MGTransferBlockSelect -------------- */

// Simplify some things below
typedef std::map<unsigned int, unsigned int>::const_iterator IT;



template <typename number>
template <int dim, typename number2, int spacedim>
void
MGTransferBlockSelect<number>::copy_from_mg (
  const MGDoFHandler<dim,spacedim>&              mg_dof_handler,
  BlockVector<number2>&                 dst,
  const MGLevelObject<Vector<number> >& src) const
{
  for (unsigned int level=0;level<mg_dof_handler.get_tria().n_levels();++level)
    for (IT i= copy_indices[selected_block][level].begin();
	 i != copy_indices[selected_block][level].end();++i)
      dst.block(selected_block)(i->first) = src[level](i->second);
}



template <typename number>
template <int dim, typename number2, int spacedim>
void
MGTransferBlockSelect<number>::copy_from_mg (
  const MGDoFHandler<dim,spacedim>&              mg_dof_handler,
  Vector<number2>&                      dst,
  const MGLevelObject<Vector<number> >& src) const
{
  for (unsigned int level=0;level<mg_dof_handler.get_tria().n_levels();++level)
    for (IT i= copy_indices[selected_block][level].begin();
	 i != copy_indices[selected_block][level].end();++i)
      dst(i->first) = src[level](i->second);
}



template <typename number>
template <int dim, typename number2, int spacedim>
void
MGTransferBlockSelect<number>::copy_from_mg_add (
  const MGDoFHandler<dim,spacedim>&              mg_dof_handler,
  BlockVector<number2>&                 dst,
  const MGLevelObject<Vector<number> >& src) const
{
  for (unsigned int level=0;level<mg_dof_handler.get_tria().n_levels();++level)
    for (IT i= copy_indices[selected_block][level].begin();
	 i != copy_indices[selected_block][level].end();++i)
      dst.block(selected_block)(i->first) += src[level](i->second);
}



template <typename number>
template <int dim, typename number2, int spacedim>
void
MGTransferBlockSelect<number>::copy_from_mg_add (
  const MGDoFHandler<dim,spacedim>&              mg_dof_handler,
  Vector<number2>&                      dst,
  const MGLevelObject<Vector<number> >& src) const
{
  for (unsigned int level=0;level<mg_dof_handler.get_tria().n_levels();++level)
    for (IT i= copy_indices[selected_block][level].begin();
	 i != copy_indices[selected_block][level].end();++i)
      dst(i->first) += src[level](i->second);
}



template <typename number>
unsigned int
MGTransferBlockSelect<number>::memory_consumption () const
{
  return sizeof(int) + MGTransferBlockBase::memory_consumption();
}


/* --------------------- MGTransferBlock -------------- */



template <typename number>
template <int dim, typename number2, int spacedim>
void
MGTransferBlock<number>::copy_from_mg (
  const MGDoFHandler<dim,spacedim>& mg_dof_handler,
  BlockVector<number2>& dst,
  const MGLevelObject<BlockVector<number> >& src) const
{
  for (unsigned int block=0;block<selected.size();++block)
    if (selected[block])
      for (unsigned int level=0;level<mg_dof_handler.get_tria().n_levels();++level)
	for (IT i= copy_indices[block][level].begin();
	     i != copy_indices[block][level].end();++i)
	  dst.block(block)(i->first) = src[level].block(mg_block[block])(i->second);
}



template <typename number>
template <int dim, typename number2, int spacedim>
void
MGTransferBlock<number>::copy_from_mg_add (
  const MGDoFHandler<dim,spacedim>& mg_dof_handler,
  BlockVector<number2>& dst,
  const MGLevelObject<BlockVector<number> >& src) const
{
  for (unsigned int block=0;block<selected.size();++block)
    if (selected[block])
      for (unsigned int level=0;level<mg_dof_handler.get_tria().n_levels();++level)
	for (IT i= copy_indices[block][level].begin();
	     i != copy_indices[block][level].end();++i)
	  dst.block(block)(i->first) += src[level].block(mg_block[block])(i->second);
}

DEAL_II_NAMESPACE_CLOSE

#endif
