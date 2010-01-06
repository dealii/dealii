//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <base/logstream.h>

#include <lac/vector.h>
#include <lac/block_vector.h>
#include <lac/sparse_matrix.h>
#include <lac/block_sparse_matrix.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_tools.h>
#include <fe/fe.h>
#include <multigrid/mg_dof_handler.h>
#include <multigrid/mg_dof_accessor.h>
#include <multigrid/mg_transfer_component.h>
#include <multigrid/mg_transfer_component.templates.h>
#include <multigrid/mg_tools.h>

#include <algorithm>
#include <numeric>

DEAL_II_NAMESPACE_OPEN


namespace
{
				     /**
				      * Adjust block-vectors on all
				      * levels to correct size.  Count
				      * the numbers of degrees of
				      * freedom on each level
				      * component-wise. Then, assign
				      * each block of @p vector the
				      * corresponding size.
				      *
				      * The boolean field @p selected
				      * allows restricting this
				      * operation to certain
				      * components. In this case, @p
				      * vector will only have as many
				      * blocks as there are true
				      * values in @p selected (no
				      * blocks of length zero are
				      * padded in). If this argument
				      * is omitted, all blocks will be
				      * considered.
				      *
				      * Degrees of freedom must be
				      * sorted by component in order
				      * to obtain reasonable results
				      * from this function.
				      *
				      * The argument
				      * @p target_component allows to
				      * re-sort and group components
				      * as in
				      * DoFRenumbering::component_wise.
				      *
				      *
				      */
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


				     /**
				      * Adjust vectors on all levels
				      * to correct size.  Count the
				      * numbers of degrees of freedom
				      * on each level component-wise
				      * in a single component. Then,
				      * assign @p vector the
				      * corresponding size.
				      *
				      * The boolean field @p selected
				      * may be nonzero in a single
				      * component, indicating the
				      * block of a block vector the
				      * argument @p v corresponds to.
				      *
				      * Degrees of freedom must be
				      * sorted by component in order
				      * to obtain reasonable results
				      * from this function.
				      *
				      * The argument
				      * @p target_component allows to
				      * re-sort and group components
				      * as in
				      * DoFRenumbering::component_wise.
				      */
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

  reinit_vector_by_components(mg_dof_handler, dst, mg_selected,
			      mg_target_component, sizes);

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
    }
}


template <int dim, int spacedim>
void MGTransferComponentBase::build_matrices (
  const DoFHandler<dim,spacedim>&,
  const MGDoFHandler<dim,spacedim>& mg_dof)
{
				   // Fill target component with
				   // standard values (identity) if it
				   // is empty
  if (target_component.size() == 0)
    {
      target_component.resize(mg_dof.get_fe().n_components());
      for (unsigned int i=0;i<target_component.size();++i)
	target_component[i] = i;
    } else {
				       // otherwise, check it for consistency
      Assert (target_component.size() == mg_dof.get_fe().n_components(),
	      ExcDimensionMismatch(target_component.size(),
				   mg_dof.get_fe().n_components()));

      for (unsigned int i=0;i<target_component.size();++i)
	{
	  Assert(i<target_component.size(),
		 ExcIndexRange(i,0,target_component.size()));
	}
    }
				   // Do the same for the multilevel
				   // components. These may be
				   // different.
  if (mg_target_component.size() == 0)
    {
      mg_target_component.resize(mg_dof.get_fe().n_components());
      for (unsigned int i=0;i<mg_target_component.size();++i)
	mg_target_component[i] = target_component[i];
    } else {
      Assert (mg_target_component.size() == mg_dof.get_fe().n_components(),
	      ExcDimensionMismatch(mg_target_component.size(),
				   mg_dof.get_fe().n_components()));

      for (unsigned int i=0;i<mg_target_component.size();++i)
	{
	  Assert(i<mg_target_component.size(),
		 ExcIndexRange(i,0,mg_target_component.size()));
	}
    }

  const FiniteElement<dim>& fe = mg_dof.get_fe();

				   // Effective number of components
				   // is the maximum entry in
				   // mg_target_component. This
				   // assumes that the values in that
				   // vector don't have holes.
  const unsigned int n_components  =
    *std::max_element(mg_target_component.begin(), mg_target_component.end()) + 1;
  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const unsigned int n_levels      = mg_dof.get_tria().n_levels();

  Assert (mg_selected.size() == fe.n_components(),
	  ExcDimensionMismatch(mg_selected.size(), fe.n_components()));

				   // Compute the lengths of all blocks
  sizes.resize(n_levels);
  MGTools::count_dofs_per_component(mg_dof, sizes, true, mg_target_component);

				   // Fill some index vectors
				   // for later use.
  mg_component_start = sizes;
				   // Compute start indices from sizes
  for (unsigned int l=0;l<mg_component_start.size();++l)
    {
      unsigned int k=0;
      for (unsigned int i=0;i<mg_component_start[l].size();++i)
	{
	  const unsigned int t=mg_component_start[l][i];
	  mg_component_start[l][i] = k;
	  k += t;
	}
    }

  component_start.resize(*std::max_element (target_component.begin(),
					    target_component.end()) + 1);
  DoFTools::
    count_dofs_per_component (static_cast<const DoFHandler<dim,spacedim>&>(mg_dof),
                              component_start, true, target_component);

  unsigned int k=0;
  for (unsigned int i=0;i<component_start.size();++i)
    {
      const unsigned int t=component_start[i];
      component_start[i] = k;
      k += t;
    }

				   // Build index vectors for
				   // copy_to_mg and
				   // copy_from_mg. These vectors must
				   // be prebuilt, since the
				   // get_dof_indices functions are
				   // too slow

  copy_to_and_from_indices.resize(n_levels);

// Building the prolongation matrices starts here!

				   // reset the size of the array of
				   // matrices. call resize(0) first,
				   // in order to delete all elements
				   // and clear their memory. then
				   // repopulate these arrays
				   //
				   // note that on resize(0), the
				   // shared_ptr class takes care of
				   // deleting the object it points to
				   // by itself
  prolongation_matrices.resize (0);
  prolongation_sparsities.resize (0);

  for (unsigned int i=0; i<n_levels-1; ++i)
    {
      prolongation_sparsities
	.push_back (std_cxx1x::shared_ptr<BlockSparsityPattern> (new BlockSparsityPattern));
      prolongation_matrices
	.push_back (std_cxx1x::shared_ptr<BlockSparseMatrix<double> > (new BlockSparseMatrix<double>));
    }

				   // two fields which will store the
				   // indices of the multigrid dofs
				   // for a cell and one of its children
  std::vector<unsigned int> dof_indices_parent (dofs_per_cell);
  std::vector<unsigned int> dof_indices_child (dofs_per_cell);

				   // for each level: first build the
				   // sparsity pattern of the matrices
				   // and then build the matrices
				   // themselves. note that we only
				   // need to take care of cells on
				   // the coarser level which have
				   // children
  for (unsigned int level=0; level<n_levels-1; ++level)
    {
				       // reset the dimension of the
				       // structure.  note that for
				       // the number of entries per
				       // row, the number of parent
				       // dofs coupling to a child dof
				       // is necessary. this, is the
				       // number of degrees of freedom
				       // per cell
      prolongation_sparsities[level]->reinit (n_components, n_components);
      for (unsigned int i=0; i<n_components; ++i)
	for (unsigned int j=0; j<n_components; ++j)
	  if (i==j)
	    prolongation_sparsities[level]->block(i,j)
	      .reinit(sizes[level+1][i],
		      sizes[level][j],
		      dofs_per_cell+1, false);
	  else
	    prolongation_sparsities[level]->block(i,j)
	      .reinit(sizes[level+1][i],
		      sizes[level][j],
		      0, false);

      prolongation_sparsities[level]->collect_sizes();

      for (typename MGDoFHandler<dim,spacedim>::cell_iterator cell=mg_dof.begin(level);
	   cell != mg_dof.end(level); ++cell)
	if (cell->has_children())
	  {
	    cell->get_mg_dof_indices (dof_indices_parent);

	    for (unsigned int child=0; child<cell->n_children(); ++child)
	      {
						 // set an alias to the
						 // prolongation matrix for
						 // this child
		const FullMatrix<double> &prolongation
		  = mg_dof.get_fe().get_prolongation_matrix (child, cell->refinement_case());

		cell->child(child)->get_mg_dof_indices (dof_indices_child);

						 // now tag the entries in the
						 // matrix which will be used
						 // for this pair of parent/child
		for (unsigned int i=0; i<dofs_per_cell; ++i)
		  for (unsigned int j=0; j<dofs_per_cell; ++j)
		    if (prolongation(i,j) != 0)
		      {
			const unsigned int icomp
			  = fe.system_to_component_index(i).first;
			const unsigned int jcomp
			  = fe.system_to_component_index(j).first;
			if ((icomp==jcomp) && mg_selected[mg_target_component[icomp]])
			  prolongation_sparsities[level]->add(dof_indices_child[i],
							      dof_indices_parent[j]);
		      };
	      };
	  };
      prolongation_sparsities[level]->compress ();

      prolongation_matrices[level]->reinit (*prolongation_sparsities[level]);
				       // now actually build the matrices
      for (typename MGDoFHandler<dim,spacedim>::cell_iterator cell=mg_dof.begin(level);
	   cell != mg_dof.end(level); ++cell)
	if (cell->has_children())
	  {
	    cell->get_mg_dof_indices (dof_indices_parent);

	    for (unsigned int child=0; child<cell->n_children(); ++child)
	      {
						 // set an alias to the
						 // prolongation matrix for
						 // this child
		const FullMatrix<double> &prolongation
		  = mg_dof.get_fe().get_prolongation_matrix (child, cell->refinement_case());

		cell->child(child)->get_mg_dof_indices (dof_indices_child);

						 // now set the entries in the
						 // matrix
		for (unsigned int i=0; i<dofs_per_cell; ++i)
		  for (unsigned int j=0; j<dofs_per_cell; ++j)
		    if (prolongation(i,j) != 0)
		      {
			const unsigned int icomp = fe.system_to_component_index(i).first;
			const unsigned int jcomp = fe.system_to_component_index(j).first;
			if ((icomp==jcomp) && mg_selected[mg_target_component[icomp]])
			  prolongation_matrices[level]->set(dof_indices_child[i],
							    dof_indices_parent[j],
							    prolongation(i,j));
		      }
	      }
	  }
    }
}


template <typename number>
template <int dim, int spacedim>
void MGTransferSelect<number>::build_matrices (
  const DoFHandler<dim,spacedim> &dof,
  const MGDoFHandler<dim,spacedim> &mg_dof,
  unsigned int select,
  unsigned int mg_select,
  const std::vector<unsigned int>& t_component,
  const std::vector<unsigned int>& mg_t_component)
{
  const FiniteElement<dim>& fe = mg_dof.get_fe();
  unsigned int ncomp = mg_dof.get_fe().n_components();

  target_component = t_component;
  mg_target_component = mg_t_component;

  selected_component = select;
  mg_selected_component = mg_select;
  selected.resize(ncomp, false);
  selected[select] = true;
  mg_selected.resize(ncomp, false);
  mg_selected[mg_select] = true;
				   // If components are renumbered,
				   // find the first original
				   // component corresponding to the
				   // target component.
  for (unsigned int i=0;i<target_component.size();++i)
    {
      if (target_component[i] == select)
	{
	  selected_component = i;
	  break;
	}
    }

  for (unsigned int i=0;i<mg_target_component.size();++i)
    {
      if (mg_target_component[i] == mg_select)
	{
	  mg_selected_component = i;
	  break;
	}
    }
  MGTransferComponentBase::build_matrices (dof, mg_dof);

  std::vector<unsigned int> global_dof_indices (fe.dofs_per_cell);
  std::vector<unsigned int> level_dof_indices  (fe.dofs_per_cell);
  for (int level=dof.get_tria().n_levels()-1; level>=0; --level)
    {
      typename MGDoFHandler<dim,spacedim>::active_cell_iterator
	level_cell = mg_dof.begin_active(level);
      const typename MGDoFHandler<dim,spacedim>::active_cell_iterator
	level_end  = mg_dof.end_active(level);

				       // Compute coarse level right hand side
				       // by restricting from fine level.
      for (; level_cell!=level_end; ++level_cell)
	{
	  DoFAccessor<dim, DoFHandler<dim,spacedim> >& global_cell = *level_cell;
					   // get the dof numbers of
					   // this cell for the global
					   // and the level-wise
					   // numbering
	  global_cell.get_dof_indices(global_dof_indices);
	  level_cell->get_mg_dof_indices (level_dof_indices);

	  for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
	    {
	      const unsigned int component
		= mg_target_component[fe.system_to_component_index(i).first];
	      if (mg_selected[component])
		{
		  const unsigned int level_start
		    = mg_component_start[level][component];
		  copy_to_and_from_indices[level].insert(
		    std::make_pair(global_dof_indices[i]
				   - component_start[selected_component],
				   level_dof_indices[i]-level_start));
		}
	    }
	}
    }
}



// explicit instantiations

template
void MGTransferSelect<float>::build_matrices<deal_II_dimension>
(const DoFHandler<deal_II_dimension> &d,
 const MGDoFHandler<deal_II_dimension> &,
 unsigned int, unsigned int,
 const std::vector<unsigned int>&,
 const std::vector<unsigned int>&);

template
void MGTransferSelect<double>::build_matrices<deal_II_dimension>
(const DoFHandler<deal_II_dimension> &d,
 const MGDoFHandler<deal_II_dimension> &,
 unsigned int, unsigned int,
 const std::vector<unsigned int>&,
 const std::vector<unsigned int>&);

template void
MGTransferSelect<float>::copy_to_mg (
  const MGDoFHandler<deal_II_dimension>&,
  MGLevelObject<Vector<float> >&,
  const Vector<double>&) const;
template void
MGTransferSelect<float>::copy_to_mg (
  const MGDoFHandler<deal_II_dimension>&,
  MGLevelObject<Vector<float> >&,
  const BlockVector<double>&) const;
template void
MGTransferSelect<float>::copy_from_mg (
  const MGDoFHandler<deal_II_dimension>&,
  Vector<double>&,
  const MGLevelObject<Vector<float> >&) const;
template void
MGTransferSelect<float>::copy_from_mg (
  const MGDoFHandler<deal_II_dimension>&,
  BlockVector<double>&,
  const MGLevelObject<Vector<float> >&) const;
template void
MGTransferSelect<float>::copy_from_mg_add (
  const MGDoFHandler<deal_II_dimension>&,
  Vector<double>&,
  const MGLevelObject<Vector<float> >&) const;
template void
MGTransferSelect<float>::copy_from_mg_add (
  const MGDoFHandler<deal_II_dimension>&,
  BlockVector<double>&,
  const MGLevelObject<Vector<float> >&) const;

template void
MGTransferSelect<double>::copy_to_mg (
  const MGDoFHandler<deal_II_dimension>&,
  MGLevelObject<Vector<double> >&,
  const Vector<double>&) const;
template void
MGTransferSelect<double>::copy_to_mg (
  const MGDoFHandler<deal_II_dimension>&,
  MGLevelObject<Vector<double> >&,
  const BlockVector<double>&) const;
template void
MGTransferSelect<double>::copy_from_mg (
  const MGDoFHandler<deal_II_dimension>&,
  Vector<double>&,
  const MGLevelObject<Vector<double> >&) const;
template void
MGTransferSelect<double>::copy_from_mg (
  const MGDoFHandler<deal_II_dimension>&,
  BlockVector<double>&,
  const MGLevelObject<Vector<double> >&) const;
template void
MGTransferSelect<double>::copy_from_mg_add (
  const MGDoFHandler<deal_II_dimension>&,
  Vector<double>&,
  const MGLevelObject<Vector<double> >&) const;
template void
MGTransferSelect<double>::copy_from_mg_add (
  const MGDoFHandler<deal_II_dimension>&,
  BlockVector<double>&,
  const MGLevelObject<Vector<double> >&) const;

DEAL_II_NAMESPACE_CLOSE
