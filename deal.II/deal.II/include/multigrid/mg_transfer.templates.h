//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2003, 2004, 2005, 2006, 2007 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#ifndef __deal2__mg_transfer_templates_h
#define __deal2__mg_transfer_templates_h

#include <lac/sparse_matrix.h>
#include <grid/tria_iterator.h>
#include <fe/fe.h>
#include <dofs/dof_constraints.h>
#include <multigrid/mg_base.h>
#include <multigrid/mg_dof_accessor.h>
#include <multigrid/mg_tools.h>
#include <multigrid/mg_transfer.h>

#include <algorithm>

DEAL_II_NAMESPACE_OPEN

/* --------------------- MGTransferPrebuilt -------------- */

// Simplify some things below
typedef std::map<unsigned int, unsigned int>::const_iterator IT;


template <class VECTOR>
template <int dim, class InVector, int spacedim>
void
MGTransferPrebuilt<VECTOR>::copy_to_mg (
  const MGDoFHandler<dim,spacedim>& mg_dof_handler,
  MGLevelObject<VECTOR>& dst,
  const InVector& src) const
{
  MGTools::reinit_vector(mg_dof_handler, dst);
  bool first = true;
  for (unsigned int level=mg_dof_handler.get_tria().n_levels();level != 0;)
    {
      --level;
      for (IT i= copy_indices[level].begin();
	   i != copy_indices[level].end();++i)
	dst[level](i->second) = src(i->first);
      if (!first)
	restrict_and_add (level+1, dst[level], dst[level+1]);
      first = false;
    }
}



template <class VECTOR>
template <int dim, class OutVector, int spacedim>
void
MGTransferPrebuilt<VECTOR>::copy_from_mg(
  const MGDoFHandler<dim,spacedim>&       mg_dof_handler,
  OutVector&                     dst,
  const MGLevelObject<VECTOR>& src) const
{
  for (unsigned int level=0;level<mg_dof_handler.get_tria().n_levels();++level)
    for (IT i= copy_indices[level].begin();
	 i != copy_indices[level].end();++i)
      dst(i->first) = src[level](i->second);
}



template <class VECTOR>
template <int dim, class OutVector, int spacedim>
void
MGTransferPrebuilt<VECTOR>::copy_from_mg_add (
  const MGDoFHandler<dim,spacedim>              &mg_dof_handler,
  OutVector                            &dst,
  const MGLevelObject<VECTOR> &src) const
{
  for (unsigned int level=0;level<mg_dof_handler.get_tria().n_levels();++level)
    for (IT i= copy_indices[level].begin();
	 i != copy_indices[level].end();++i)
      dst(i->first) += src[level](i->second);
}


template <class VECTOR>
unsigned int
MGTransferPrebuilt<VECTOR>::memory_consumption () const
{
  unsigned int result = sizeof(*this);
  result += sizeof(unsigned int) * sizes.size();
#ifdef DEAL_PREFER_MATRIX_EZ
  std::vector<std_cxx0x::shared_ptr<SparseMatrixEZ<double> > >::const_iterator m;
  const std::vector<std_cxx0x::shared_ptr<SparseMatrixEZ<double> > >::const_iterator end = prolongation_matrices.end();
  for (m = prolongation_matrices.begin(); m != end ; ++m)
    result += *m->memory_consumption();
#else
  for (unsigned int i=0;i<prolongation_matrices.size();++i)
    result += prolongation_matrices[i]->memory_consumption()
	      + prolongation_sparsities[i]->memory_consumption();
#endif
  return result;
}


DEAL_II_NAMESPACE_CLOSE

#endif
