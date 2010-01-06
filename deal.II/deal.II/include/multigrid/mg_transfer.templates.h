//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2003, 2004, 2005, 2006, 2007, 2009, 2010 by the deal.II authors
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
#include <multigrid/mg_base.h>
#include <multigrid/mg_dof_accessor.h>
#include <multigrid/mg_tools.h>
#include <multigrid/mg_transfer.h>

#include <algorithm>

DEAL_II_NAMESPACE_OPEN

/* --------------------- MGTransferPrebuilt -------------- */




template <class VECTOR>
template <int dim, class OutVector, int spacedim>
void
MGTransferPrebuilt<VECTOR>::copy_from_mg(
  const MGDoFHandler<dim,spacedim>&       mg_dof_handler,
  OutVector&                     dst,
  const MGLevelObject<VECTOR>& src) const
{
				       // For non-DG: degrees of
				       // freedom in the refinement
				       // face may need special
				       // attention, since they belong
				       // to the coarse level, but
				       // have fine level basis
				       // functions
  dst = 0;
  for (unsigned int level=0;level<mg_dof_handler.get_tria().n_levels();++level)
  {
    typedef std::map<unsigned int, unsigned int>::const_iterator IT;
    for (IT i= copy_indices[level].begin();
	 i != copy_indices[level].end();++i)
      dst(i->first) = src[level](i->second);
  }
  if (constraints != 0)
    constraints->condense(dst);
}



template <class VECTOR>
template <int dim, class OutVector, int spacedim>
void
MGTransferPrebuilt<VECTOR>::copy_from_mg_add (
  const MGDoFHandler<dim,spacedim>& mg_dof_handler,
  OutVector                            &dst,
  const MGLevelObject<VECTOR> &src) const
{
				       // For non-DG: degrees of
				       // freedom in the refinement
				       // face may need special
				       // attention, since they belong
				       // to the coarse level, but
				       // have fine level basis
				       // functions
  typedef std::map<unsigned int, unsigned int>::const_iterator IT;
  for (unsigned int level=0;level<mg_dof_handler.get_tria().n_levels();++level)
    for (IT i= copy_indices[level].begin();
	 i != copy_indices[level].end();++i)
      dst(i->first) += src[level](i->second);
}



template <class VECTOR>
void
MGTransferPrebuilt<VECTOR>::
set_component_to_block_map (const std::vector<unsigned int> &map)
{
  component_to_block_map = map;
}



template <class VECTOR>
template <int dim, int spacedim>
void MGTransferPrebuilt<VECTOR>::find_dofs_on_refinement_edges (
    const MGDoFHandler<dim,spacedim>& mg_dof)
{
  for (unsigned int level = 0; level<mg_dof.get_tria().n_levels(); ++level)
  {
    std::vector<bool> tmp (mg_dof.n_dofs(level));
    dofs_on_refinement_edge.push_back (tmp);
    dofs_on_refinement_boundary.push_back (tmp);
  }

  const unsigned int   dofs_per_cell   = mg_dof.get_fe().dofs_per_cell;
  const unsigned int   dofs_per_face   = mg_dof.get_fe().dofs_per_face;

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);
  std::vector<unsigned int> face_dof_indices (dofs_per_face);

  typename MGDoFHandler<dim>::cell_iterator cell = mg_dof.begin(),
           endc = mg_dof.end();

  for (; cell!=endc; ++cell)
  {
    std::vector<bool> cell_dofs(dofs_per_cell);
    std::vector<bool> boundary_cell_dofs(dofs_per_cell);
    const unsigned int level = cell->level();
    cell->get_mg_dof_indices (local_dof_indices);
    for (unsigned int face_nr=0;
        face_nr<GeometryInfo<dim>::faces_per_cell; ++face_nr)
    {
      typename DoFHandler<dim>::face_iterator face = cell->face(face_nr);
      if(!cell->at_boundary(face_nr))
      {
        //interior face
        typename MGDoFHandler<dim>::cell_iterator neighbor
          = cell->neighbor(face_nr);
        // Do refinement face
        // from the coarse side
        if (neighbor->level() < cell->level())
        {
          for(unsigned int j=0; j<dofs_per_face; ++j)
            cell_dofs[mg_dof.get_fe().face_to_cell_index(j,face_nr)] = true;
        }
      }
      //boundary face
      else
        for(unsigned int j=0; j<dofs_per_face; ++j)
          boundary_cell_dofs[mg_dof.get_fe().face_to_cell_index(j,face_nr)] = true;
    }//faces
    for(unsigned int i=0; i<dofs_per_cell; ++i)
    {
      if(cell_dofs[i])
      {
        dofs_on_refinement_edge[level][local_dof_indices[i]] = true;
        dofs_on_refinement_boundary[level][local_dof_indices[i]] = true;
      }
      if(boundary_cell_dofs[i])
        dofs_on_refinement_boundary[level][local_dof_indices[i]] = true;
    }
  }
}

template <class VECTOR>
unsigned int
MGTransferPrebuilt<VECTOR>::memory_consumption () const
{
  unsigned int result = sizeof(*this);
  result += sizeof(unsigned int) * sizes.size();
#ifdef DEAL_PREFER_MATRIX_EZ
  std::vector<std_cxx1x::shared_ptr<SparseMatrixEZ<double> > >::const_iterator m;
  const std::vector<std_cxx1x::shared_ptr<SparseMatrixEZ<double> > >::const_iterator end = prolongation_matrices.end();
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
