//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2003, 2004, 2005, 2006, 2007, 2009, 2010, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#ifndef __deal2__mg_transfer_templates_h
#define __deal2__mg_transfer_templates_h

#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/fe/fe.h>
#include <deal.II/multigrid/mg_base.h>
#include <deal.II/multigrid/mg_dof_accessor.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_transfer.h>

#include <algorithm>

DEAL_II_NAMESPACE_OPEN


namespace
{
				     /**
				      * Adjust vectors on all levels to
				      * correct size.  Here, we just
				      * count the numbers of degrees
				      * of freedom on each level and
				      * @p reinit each level vector
				      * to this length.
                                      * For compatibility reasons with
                                      * the next function
                                      * the target_component is added
                                      * here but is not used.
				      */
  template <int dim, typename number, int spacedim>
  void
  reinit_vector (const dealii::MGDoFHandler<dim,spacedim> &mg_dof,
		 std::vector<unsigned int> ,
		 MGLevelObject<dealii::Vector<number> > &v)
  {
    for (unsigned int level=v.get_minlevel();
	 level<=v.get_maxlevel();++level)
      {
	unsigned int n = mg_dof.n_dofs (level);
	v[level].reinit(n);
      }

  }


				     /**
				      * Adjust vectors on all levels to
				      * correct size.  Here, we just
				      * count the numbers of degrees
				      * of freedom on each level and
				      * @p reinit each level vector
				      * to this length. The target_component
                                      * is handed to MGTools::count_dofs_per_block.
                                      * See for documentation there.
				      */
  template <int dim, typename number, int spacedim>
  void
  reinit_vector (const dealii::MGDoFHandler<dim,spacedim> &mg_dof,
		 std::vector<unsigned int> target_component,
		 MGLevelObject<BlockVector<number> > &v)
  {
    const unsigned int n_blocks = mg_dof.get_fe().n_blocks();
    if (target_component.size()==0)
      {
        target_component.resize(n_blocks);
        for (unsigned int i=0;i<n_blocks;++i)
          target_component[i] = i;
      }
    Assert(target_component.size()==n_blocks,
	   ExcDimensionMismatch(target_component.size(),n_blocks));
    const unsigned int max_block
      = *std::max_element (target_component.begin(),
			   target_component.end());
    const unsigned int n_target_blocks = max_block + 1;

    std::vector<std::vector<unsigned int> >
      ndofs(mg_dof.get_tria().n_levels(),
	    std::vector<unsigned int>(n_target_blocks));
    MGTools::count_dofs_per_block (mg_dof, ndofs, target_component);

    for (unsigned int level=v.get_minlevel();
	 level<=v.get_maxlevel();++level)
      {
	v[level].reinit(n_target_blocks);
	for (unsigned int b=0; b<n_target_blocks; ++b)
	  v[level].block(b).reinit(ndofs[level][b]);
	v[level].collect_sizes();
      }
  }
}



/* --------------------- MGTransferPrebuilt -------------- */




template <class VECTOR>
template <int dim, class InVector, int spacedim>
void
MGTransferPrebuilt<VECTOR>::copy_to_mg (
  const MGDoFHandler<dim,spacedim>& mg_dof_handler,
  MGLevelObject<VECTOR>& dst,
  const InVector& src) const
{
  reinit_vector(mg_dof_handler, component_to_block_map, dst);
  bool first = true;
  for (unsigned int level=mg_dof_handler.get_tria().n_levels();level != 0;)
    {
      --level;
      VECTOR& dst_level = dst[level];

      typedef std::vector<std::pair<unsigned int, unsigned int> >::const_iterator IT;
      for (IT i= copy_indices[level].begin();
	   i != copy_indices[level].end();++i)
	dst_level(i->second) = src(i->first);

				       // For non-DG: degrees of
				       // freedom in the refinement
				       // face may need special
				       // attention, since they belong
				       // to the coarse level, but
				       // have fine level basis
				       // functions
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
    typedef std::vector<std::pair<unsigned int, unsigned int> >::const_iterator IT;

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
  for (unsigned int level=0;level<mg_dof_handler.get_tria().n_levels();++level)
  {
    typedef std::vector<std::pair<unsigned int, unsigned int> >::const_iterator IT;
    for (IT i= copy_indices[level].begin();
	 i != copy_indices[level].end();++i)
      dst(i->first) += src[level](i->second);
  }
}



template <class VECTOR>
void
MGTransferPrebuilt<VECTOR>::
set_component_to_block_map (const std::vector<unsigned int> &map)
{
  component_to_block_map = map;
}

template <class VECTOR>
std::size_t
MGTransferPrebuilt<VECTOR>::memory_consumption () const
{
  std::size_t result = sizeof(*this);
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
