//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2005, 2006, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <lac/vector.h>
#include <lac/block_vector.h>
#include <fe/fe.h>
#include <multigrid/mg_dof_handler.h>
#include <multigrid/mg_base.h>
#include <multigrid/mg_level_object.h>
#include <multigrid/mg_tools.h>

#include <numeric>

DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace mg
  {
    template <int dim, typename number, int spacedim>
    void
    reinit_vector (const dealii::MGDoFHandler<dim,spacedim> &mg_dof,
		   MGLevelObject<dealii::Vector<number> > &v)
    {
      for (unsigned int level=v.get_minlevel();
	   level<=v.get_maxlevel();++level)
	{
	  unsigned int n = mg_dof.n_dofs (level);
	  v[level].reinit(n);
	}

    }


    template <int dim, typename number, int spacedim>
    void
    reinit_vector (const dealii::MGDoFHandler<dim,spacedim> &mg_dof,
		   MGLevelObject<BlockVector<number> > &v)
    {
      const unsigned int n_blocks = mg_dof.get_fe().n_blocks();
      std::vector<std::vector<unsigned int> >
	ndofs(mg_dof.get_tria().n_levels(),
	      std::vector<unsigned int>(n_blocks));
      MGTools::count_dofs_per_block (mg_dof, ndofs);

      for (unsigned int level=v.get_minlevel();
	   level<=v.get_maxlevel();++level)
	{
	  v[level].reinit(ndofs[level]);
	}
    }


    template <int dim, typename number, int spacedim>
    void
    reinit_vector_by_components (
      const dealii::MGDoFHandler<dim,spacedim> &mg_dof,
      MGLevelObject<BlockVector<number> > &v,
      const std::vector<bool> &sel,
      const std::vector<unsigned int> &target_comp,
      std::vector<std::vector<unsigned int> >& ndofs)
    {
      std::vector<bool> selected=sel;
      std::vector<unsigned int> target_component=target_comp;
      const unsigned int ncomp = mg_dof.get_fe().n_components();

				       // If the selected and
				       // target_component have size 0,
				       // they must be replaced by default
				       // values.
				       //
				       // Since we already made copies
				       // directly after this function was
				       // called, we use the arguments
				       // directly.
      if (target_component.size() == 0)
	{
	  target_component.resize(ncomp);
	  for (unsigned int i=0;i<ncomp;++i)
	    target_component[i] = i;
	}

				       // If selected is an empty vector,
				       // all components are selected.
      if (selected.size() == 0)
	{
	  selected.resize(target_component.size());
	  std::fill_n (selected.begin(), ncomp, false);
	  for (unsigned int i=0;i<target_component.size();++i)
	    selected[target_component[i]] = true;
	}

      Assert (selected.size() == target_component.size(),
	      ExcDimensionMismatch(selected.size(), target_component.size()));

				       // Compute the number of blocks needed
      const unsigned int n_selected
	= std::accumulate(selected.begin(),
			  selected.end(),
			  0U);

      if (ndofs.size() == 0)
	{
	  std::vector<std::vector<unsigned int> >
	    new_dofs(mg_dof.get_tria().n_levels(),
		     std::vector<unsigned int>(target_component.size()));
	  std::swap(ndofs, new_dofs);
	  MGTools::count_dofs_per_component (mg_dof, ndofs,
					     true, target_component);
	}

      for (unsigned int level=v.get_minlevel();
	   level<=v.get_maxlevel();++level)
	{
	  v[level].reinit(n_selected, 0);
	  unsigned int k=0;
	  for (unsigned int i=0;i<selected.size() && (k<v[level].n_blocks());++i)
	    {
	      if (selected[i])
		{
		  v[level].block(k++).reinit(ndofs[level][i]);
		}
	      v[level].collect_sizes();
	    }
	}
    }


    template <int dim, typename number, int spacedim>
    void
    reinit_vector_by_components (
      const dealii::MGDoFHandler<dim,spacedim> &mg_dof,
      MGLevelObject<dealii::Vector<number> > &v,
      const std::vector<bool> &selected,
      const std::vector<unsigned int> &target_component,
      std::vector<std::vector<unsigned int> >& ndofs)
    {
      Assert (selected.size() == target_component.size(),
	      ExcDimensionMismatch(selected.size(), target_component.size()));

				       // Compute the number of blocks needed
#ifdef DEBUG
      const unsigned int n_selected
	= std::accumulate(selected.begin(),
			  selected.end(),
			  0U);
      Assert(n_selected == 1, ExcDimensionMismatch(n_selected, 1));
#endif

      unsigned int selected_block = 0;
      while (!selected[selected_block])
	++selected_block;

      if (ndofs.size() == 0)
	{
	  std::vector<std::vector<unsigned int> >
	    new_dofs(mg_dof.get_tria().n_levels(),
		     std::vector<unsigned int>(target_component.size()));
	  std::swap(ndofs, new_dofs);
	  MGTools::count_dofs_per_component (mg_dof, ndofs,
					     true, target_component);
	}

      for (unsigned int level=v.get_minlevel();
	   level<=v.get_maxlevel();++level)
	{
	  v[level].reinit(ndofs[level][selected_block]);
	}
    }


    template <int dim, typename number, int spacedim>
    void
    reinit_vector_by_blocks (
      const dealii::MGDoFHandler<dim,spacedim> &mg_dof,
      MGLevelObject<BlockVector<number> > &v,
      const std::vector<bool> &sel,
      std::vector<std::vector<unsigned int> >& ndofs)
    {
      std::vector<bool> selected=sel;
				       // Compute the number of blocks needed
      const unsigned int n_selected
	= std::accumulate(selected.begin(),
			  selected.end(),
			  0U);

      if (ndofs.size() == 0)
	{
	  std::vector<std::vector<unsigned int> >
	    new_dofs(mg_dof.get_tria().n_levels(),
		     std::vector<unsigned int>(selected.size()));
	  std::swap(ndofs, new_dofs);
	  MGTools::count_dofs_per_block (mg_dof, ndofs);
	}

      for (unsigned int level=v.get_minlevel();
	   level<=v.get_maxlevel();++level)
	{
	  v[level].reinit(n_selected, 0);
	  unsigned int k=0;
	  for (unsigned int i=0;i<selected.size() && (k<v[level].n_blocks());++i)
	    {
	      if (selected[i])
		{
		  v[level].block(k++).reinit(ndofs[level][i]);
		}
	      v[level].collect_sizes();
	    }
	}
    }


    template <int dim, typename number, int spacedim>
    void
    reinit_vector_by_blocks (
      const dealii::MGDoFHandler<dim,spacedim> &mg_dof,
      MGLevelObject<dealii::Vector<number> > &v,
      const unsigned int selected_block,
      std::vector<std::vector<unsigned int> >& ndofs)
    {
      const unsigned int n_blocks = mg_dof.get_fe().n_blocks();
      Assert(selected_block < n_blocks, ExcIndexRange(selected_block, 0, n_blocks));

      std::vector<bool> selected(n_blocks, false);
      selected[selected_block] = true;

      if (ndofs.size() == 0)
	{
	  std::vector<std::vector<unsigned int> >
	    new_dofs(mg_dof.get_tria().n_levels(),
		     std::vector<unsigned int>(selected.size()));
	  std::swap(ndofs, new_dofs);
	  MGTools::count_dofs_per_block (mg_dof, ndofs);
	}

      for (unsigned int level=v.get_minlevel();
	   level<=v.get_maxlevel();++level)
	{
	  v[level].reinit(ndofs[level][selected_block]);
	}
    }
  }
}



template <class VECTOR>
MGTransferBase<VECTOR>::~MGTransferBase()
{}


template <class VECTOR>
MGMatrixBase<VECTOR>::~MGMatrixBase()
{}


template <class VECTOR>
MGSmootherBase<VECTOR>::~MGSmootherBase()
{}


template <class VECTOR>
MGCoarseGridBase<VECTOR>::~MGCoarseGridBase()
{}


// Explicit instantiations

//TODO: Use the template expander script for this
namespace internal
{
  namespace mg
  {
    template void reinit_vector<deal_II_dimension> (
      const dealii::MGDoFHandler<deal_II_dimension>&,
      MGLevelObject<dealii::Vector<double> >&);
    template void reinit_vector<deal_II_dimension> (
      const dealii::MGDoFHandler<deal_II_dimension>&,
      MGLevelObject<dealii::Vector<float> >&);

    template void reinit_vector<deal_II_dimension> (
      const dealii::MGDoFHandler<deal_II_dimension>&,
      MGLevelObject<BlockVector<double> >&);
    template void reinit_vector<deal_II_dimension> (
      const dealii::MGDoFHandler<deal_II_dimension>&,
      MGLevelObject<BlockVector<float> >&);

    template void reinit_vector_by_components<deal_II_dimension> (
      const dealii::MGDoFHandler<deal_II_dimension>&,
      MGLevelObject<BlockVector<double> >&,
      const std::vector<bool> &,
      const std::vector<unsigned int> &,
      std::vector<std::vector<unsigned int> >&);
    template void reinit_vector_by_components<deal_II_dimension> (
      const dealii::MGDoFHandler<deal_II_dimension>&,
      MGLevelObject<BlockVector<float> >&,
      const std::vector<bool> &,
      const std::vector<unsigned int> &,
      std::vector<std::vector<unsigned int> >&);

    template void reinit_vector_by_components<deal_II_dimension> (
      const dealii::MGDoFHandler<deal_II_dimension>&, MGLevelObject<dealii::Vector<double> >&,
      const std::vector<bool>&, const std::vector<unsigned int>&,
      std::vector<std::vector<unsigned int> >&);
    template void reinit_vector_by_components<deal_II_dimension> (
      const dealii::MGDoFHandler<deal_II_dimension>&, MGLevelObject<dealii::Vector<float> >&,
      const std::vector<bool>&, const std::vector<unsigned int>&,
      std::vector<std::vector<unsigned int> >&);

    template void reinit_vector_by_blocks<deal_II_dimension> (
      const dealii::MGDoFHandler<deal_II_dimension>&, MGLevelObject<BlockVector<double> >&,
      const std::vector<bool> &, std::vector<std::vector<unsigned int> >&);
    template void reinit_vector_by_blocks<deal_II_dimension> (
      const dealii::MGDoFHandler<deal_II_dimension>&, MGLevelObject<BlockVector<float> >&,
      const std::vector<bool> &, std::vector<std::vector<unsigned int> >&);

    template void reinit_vector_by_blocks<deal_II_dimension> (
      const dealii::MGDoFHandler<deal_II_dimension>&, MGLevelObject<dealii::Vector<double> >&,
      unsigned int, std::vector<std::vector<unsigned int> >&);
    template void reinit_vector_by_blocks<deal_II_dimension> (
      const dealii::MGDoFHandler<deal_II_dimension>&, MGLevelObject<dealii::Vector<float> >&,
      unsigned int, std::vector<std::vector<unsigned int> >&);
  }
}


template class MGTransferBase<dealii::Vector<double> >;
template class MGTransferBase<dealii::Vector<float> >;
template class MGTransferBase<BlockVector<double> >;
template class MGTransferBase<BlockVector<float> >;

template class MGMatrixBase<dealii::Vector<double> >;
template class MGMatrixBase<dealii::Vector<float> >;
template class MGMatrixBase<BlockVector<double> >;
template class MGMatrixBase<BlockVector<float> >;

template class MGSmootherBase<dealii::Vector<float> >;
template class MGSmootherBase<dealii::Vector<double> >;
template class MGSmootherBase<BlockVector<float> >;
template class MGSmootherBase<BlockVector<double> >;

template class MGCoarseGridBase<dealii::Vector<double> >;
template class MGCoarseGridBase<dealii::Vector<float> >;
template class MGCoarseGridBase<BlockVector<double> >;
template class MGCoarseGridBase<BlockVector<float> >;

DEAL_II_NAMESPACE_CLOSE
