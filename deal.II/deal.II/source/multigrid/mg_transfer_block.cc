//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2003, 2004, 2005, 2006, 2007, 2008, 2010 by the deal.II authors
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
#include <multigrid/mg_transfer_block.h>
#include <multigrid/mg_transfer_block.templates.h>
#include <multigrid/mg_tools.h>

#include <algorithm>
#include <numeric>
#include <iostream>

DEAL_II_NAMESPACE_OPEN

namespace
{
				     /**
				      * Adjust vectors on all levels
				      * to correct size. The degrees
				      * of freedom on each level are
				      * counted by block and only the
				      * block selected is used.
				      */
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


				     /**
				      * Adjust block vectors on all
				      * levels to correct size. The
				      * degrees of freedom on each
				      * level are counted by block.
				      */
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


template <typename number>
template <int dim, typename number2, int spacedim>
void
MGTransferBlockSelect<number>::copy_to_mg (
  const MGDoFHandler<dim,spacedim>        &mg_dof_handler,
  MGLevelObject<Vector<number> > &dst,
  const BlockVector<number2>     &src) const
{
  reinit_vector_by_blocks(mg_dof_handler, dst, selected_block, sizes);
				   // For MGTransferBlockSelect, the
				   // multilevel block is always the
				   // first, since only one block is
				   // selected.
  bool first = true;
  for (unsigned int level=mg_dof_handler.get_tria().n_levels();level != 0;)
    {
      --level;
      for (IT i= copy_indices[selected_block][level].begin();
	   i != copy_indices[selected_block][level].end();++i)
	dst[level](i->second) = src.block(selected_block)(i->first);
      if (!first)
	restrict_and_add (level+1, dst[level], dst[level+1]);
      first = false;
    }
}



template <typename number>
template <int dim, typename number2, int spacedim>
void
MGTransferBlockSelect<number>::copy_to_mg (
  const MGDoFHandler<dim,spacedim>        &mg_dof_handler,
  MGLevelObject<Vector<number> > &dst,
  const Vector<number2>          &src) const
{
  reinit_vector_by_blocks(mg_dof_handler, dst, selected_block, sizes);
				   // For MGTransferBlockSelect, the
				   // multilevel block is always the
				   // first, since only one block is selected.
  bool first = true;
  for (unsigned int level=mg_dof_handler.get_tria().n_levels();level != 0;)
    {
      --level;
      for (IT i= copy_indices[selected_block][level].begin();
	   i != copy_indices[selected_block][level].end();++i)
	dst[level](i->second) = src(i->first);
      if (!first)
	restrict_and_add (level+1, dst[level], dst[level+1]);
      first = false;
    }
}



template <typename number>
template <int dim, typename number2, int spacedim>
void
MGTransferBlock<number>::copy_to_mg (
  const MGDoFHandler<dim,spacedim>& mg_dof_handler,
  MGLevelObject<BlockVector<number> >& dst,
  const BlockVector<number2>& src) const
{
  reinit_vector_by_blocks(mg_dof_handler, dst, selected, sizes);
  bool first = true;
  for (unsigned int level=mg_dof_handler.get_tria().n_levels();level != 0;)
    {
      --level;
      for (unsigned int block=0;block<selected.size();++block)
	if (selected[block])
	  for (IT i= copy_indices[block][level].begin();
	       i != copy_indices[block][level].end();++i)
	    dst[level].block(mg_block[block])(i->second) = src.block(block)(i->first);
      if (!first)
	restrict_and_add (level+1, dst[level], dst[level+1]);
      first = false;
    }
}



template <int dim, int spacedim>
void MGTransferBlockBase::build_matrices (
  const DoFHandler<dim,spacedim>&,
  const MGDoFHandler<dim,spacedim>& mg_dof)
{
  const FiniteElement<dim>& fe = mg_dof.get_fe();
  const unsigned int n_blocks  = fe.n_blocks();
  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const unsigned int n_levels      = mg_dof.get_tria().n_levels();

  Assert (selected.size() == n_blocks,
	  ExcDimensionMismatch(selected.size(), n_blocks));

				   // Compute the mapping between real
				   // blocks and blocks used for
				   // multigrid computations.
  mg_block.resize(n_blocks);
  n_mg_blocks = 0;
  for (unsigned int i=0;i<n_blocks;++i)
    if (selected[i])
      mg_block[i] = n_mg_blocks++;
    else
      mg_block[i] = numbers::invalid_unsigned_int;

				   // Compute the lengths of all blocks
  sizes.clear ();
  sizes.resize(n_levels, std::vector<unsigned int>(fe.n_blocks()));
  MGTools::count_dofs_per_block(mg_dof, sizes);

				   // Fill some index vectors
				   // for later use.
  mg_block_start = sizes;
				   // Compute start indices from sizes
  for (unsigned int l=0;l<mg_block_start.size();++l)
    {
      unsigned int k=0;
      for (unsigned int i=0;i<mg_block_start[l].size();++i)
	{
	  const unsigned int t=mg_block_start[l][i];
	  mg_block_start[l][i] = k;
	  k += t;
	}
    }

  block_start.resize(n_blocks);
  DoFTools::count_dofs_per_block (static_cast<const DoFHandler<dim,spacedim>&>(mg_dof),
				  block_start);

  unsigned int k=0;
  for (unsigned int i=0;i<block_start.size();++i)
    {
      const unsigned int t=block_start[i];
      block_start[i] = k;
      k += t;
    }
				   // Build index vectors for
				   // copy_to_mg and
				   // copy_from_mg. These vectors must
				   // be prebuilt, since the
				   // get_dof_indices functions are
				   // too slow
  copy_indices.resize(n_blocks);
  for (unsigned int block=0;block<n_blocks;++block)
    if (selected[block])
      copy_indices[block].resize(n_levels);

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
      prolongation_sparsities[level]->reinit (n_blocks, n_blocks);
      for (unsigned int i=0; i<n_blocks; ++i)
	for (unsigned int j=0; j<n_blocks; ++j)
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

	    Assert(cell->n_children()==GeometryInfo<dim>::max_children_per_cell,
		   ExcNotImplemented());
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
			  = fe.system_to_block_index(i).first;
			const unsigned int jcomp
			  = fe.system_to_block_index(j).first;
			if ((icomp==jcomp) && selected[icomp])
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

	    Assert(cell->n_children()==GeometryInfo<dim>::max_children_per_cell,
		   ExcNotImplemented());
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
			const unsigned int icomp = fe.system_to_block_index(i).first;
			const unsigned int jcomp = fe.system_to_block_index(j).first;
			if ((icomp==jcomp) && selected[icomp])
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
void MGTransferBlockSelect<number>::build_matrices (
  const DoFHandler<dim,spacedim> &dof,
  const MGDoFHandler<dim,spacedim> &mg_dof,
  unsigned int select)
{
  const FiniteElement<dim>& fe = mg_dof.get_fe();
  unsigned int n_blocks = mg_dof.get_fe().n_blocks();

  selected_block = select;
  selected.resize(n_blocks, false);
  selected[select] = true;

  MGTransferBlockBase::build_matrices (dof, mg_dof);

  std::vector<unsigned int> temp_copy_indices;
  std::vector<unsigned int> global_dof_indices (fe.dofs_per_cell);
  std::vector<unsigned int> level_dof_indices  (fe.dofs_per_cell);
  for (int level=dof.get_tria().n_levels()-1; level>=0; --level)
    {
      typename MGDoFHandler<dim,spacedim>::active_cell_iterator
	level_cell = mg_dof.begin_active(level);
      const typename MGDoFHandler<dim,spacedim>::active_cell_iterator
	level_end  = mg_dof.end_active(level);

      temp_copy_indices.resize (0);
      temp_copy_indices.resize (sizes[level][selected_block],
				numbers::invalid_unsigned_int);

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
	      const unsigned int block = fe.system_to_block_index(i).first;
	      if (selected[block])
		temp_copy_indices[level_dof_indices[i] - mg_block_start[level][block]]
		  = global_dof_indices[i] - block_start[block];
	    }
	}

				// now all the active dofs got a valid entry,
				// the other ones have an invalid entry. Count
				// the invalid entries and then resize the
				// copy_indices object. Then, insert the pairs
				// of global index and level index into
				// copy_indices.
      const unsigned int n_active_dofs =
	std::count_if (temp_copy_indices.begin(), temp_copy_indices.end(),
		       std::bind2nd(std::not_equal_to<unsigned int>(),
				    numbers::invalid_unsigned_int));
      copy_indices[selected_block][level].resize (n_active_dofs);
      unsigned int counter = 0;
      for (unsigned int i=0; i<temp_copy_indices.size(); ++i)
	if (temp_copy_indices[i] != numbers::invalid_unsigned_int)
	  copy_indices[selected_block][level][counter++] =
	    std::make_pair<unsigned int,unsigned int> (temp_copy_indices[i], i);
      Assert (counter == n_active_dofs, ExcInternalError());
    }
}




template <typename number>
template <int dim, int spacedim>
void MGTransferBlock<number>::build_matrices (
  const DoFHandler<dim,spacedim> &dof,
  const MGDoFHandler<dim,spacedim> &mg_dof,
  const std::vector<bool>& sel)
{
  const FiniteElement<dim>& fe = mg_dof.get_fe();
  unsigned int n_blocks = mg_dof.get_fe().n_blocks();

  if (sel.size() != 0)
    {
      Assert(sel.size() == n_blocks,
	     ExcDimensionMismatch(sel.size(), n_blocks));
      selected = sel;
    }
  if (selected.size() == 0)
    selected = std::vector<bool> (n_blocks, true);

  MGTransferBlockBase::build_matrices (dof, mg_dof);

  std::vector<std::vector<unsigned int> > temp_copy_indices (n_blocks);
  std::vector<unsigned int> global_dof_indices (fe.dofs_per_cell);
  std::vector<unsigned int> level_dof_indices  (fe.dofs_per_cell);
  for (int level=dof.get_tria().n_levels()-1; level>=0; --level)
    {
      typename MGDoFHandler<dim,spacedim>::active_cell_iterator
	level_cell = mg_dof.begin_active(level);
      const typename MGDoFHandler<dim,spacedim>::active_cell_iterator
	level_end  = mg_dof.end_active(level);

      for (unsigned int block=0; block<n_blocks; ++block)
	if (selected[block])
	  {
	    temp_copy_indices[block].resize (0);
	    temp_copy_indices[block].resize (sizes[level][block],
					     numbers::invalid_unsigned_int);
	  }

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
	      const unsigned int block = fe.system_to_block_index(i).first;
	      if (selected[block])
		temp_copy_indices[block][level_dof_indices[i] - mg_block_start[level][block]]
		  = global_dof_indices[i] - block_start[block];
	    }
	}

      for (unsigned int block=0; block<n_blocks; ++block)
	if (selected[block])
	  {
	    const unsigned int n_active_dofs =
	      std::count_if (temp_copy_indices[block].begin(),
			     temp_copy_indices[block].end(),
			     std::bind2nd(std::not_equal_to<unsigned int>(),
					  numbers::invalid_unsigned_int));
	    copy_indices[block][level].resize (n_active_dofs);
	    unsigned int counter = 0;
	    for (unsigned int i=0; i<temp_copy_indices[block].size(); ++i)
	      if (temp_copy_indices[block][i] != numbers::invalid_unsigned_int)
		copy_indices[block][level][counter++] =
		  std::make_pair<unsigned int,unsigned int>
		  (temp_copy_indices[block][i], i);
	    Assert (counter == n_active_dofs, ExcInternalError());
	  }
    }
}



// explicit instantiations

template
void MGTransferBlock<float>::build_matrices<deal_II_dimension>
(const DoFHandler<deal_II_dimension>&, const MGDoFHandler<deal_II_dimension>&,
 const std::vector<bool>&);

template
void MGTransferBlock<double>::build_matrices<deal_II_dimension>
(const DoFHandler<deal_II_dimension>&, const MGDoFHandler<deal_II_dimension>&,
 const std::vector<bool>&);

template
void MGTransferBlockSelect<float>::build_matrices<deal_II_dimension>
(const DoFHandler<deal_II_dimension>&, const MGDoFHandler<deal_II_dimension>&,
 const unsigned int);

template
void MGTransferBlockSelect<double>::build_matrices<deal_II_dimension>
(const DoFHandler<deal_II_dimension>&, const MGDoFHandler<deal_II_dimension>&,
 const unsigned int);

template void
MGTransferBlock<float>::copy_to_mg (
  const MGDoFHandler<deal_II_dimension>&,
  MGLevelObject<BlockVector<float> >&,
  const BlockVector<double>&) const;
template void
MGTransferBlock<float>::copy_from_mg (
  const MGDoFHandler<deal_II_dimension>&,
  BlockVector<double>&,
  const MGLevelObject<BlockVector<float> >&) const;
template void
MGTransferBlock<float>::copy_from_mg_add (
  const MGDoFHandler<deal_II_dimension>&,
  BlockVector<double>&,
  const MGLevelObject<BlockVector<float> >&) const;

template void
MGTransferBlock<double>::copy_to_mg (
  const MGDoFHandler<deal_II_dimension>&,
  MGLevelObject<BlockVector<double> >&,
  const BlockVector<double>&) const;
template void
MGTransferBlock<double>::copy_from_mg (
  const MGDoFHandler<deal_II_dimension>&,
  BlockVector<double>&,
  const MGLevelObject<BlockVector<double> >&) const;
template void
MGTransferBlock<double>::copy_from_mg_add (
  const MGDoFHandler<deal_II_dimension>&,
  BlockVector<double>&,
  const MGLevelObject<BlockVector<double> >&) const;

template void
MGTransferBlockSelect<float>::copy_to_mg (
  const MGDoFHandler<deal_II_dimension>&,
  MGLevelObject<Vector<float> >&,
  const Vector<double>&) const;
template void
MGTransferBlockSelect<float>::copy_to_mg (
  const MGDoFHandler<deal_II_dimension>&,
  MGLevelObject<Vector<float> >&,
  const BlockVector<double>&) const;
template void
MGTransferBlockSelect<float>::copy_from_mg (
  const MGDoFHandler<deal_II_dimension>&,
  Vector<double>&,
  const MGLevelObject<Vector<float> >&) const;
template void
MGTransferBlockSelect<float>::copy_from_mg (
  const MGDoFHandler<deal_II_dimension>&,
  BlockVector<double>&,
  const MGLevelObject<Vector<float> >&) const;
template void
MGTransferBlockSelect<float>::copy_from_mg_add (
  const MGDoFHandler<deal_II_dimension>&,
  Vector<double>&,
  const MGLevelObject<Vector<float> >&) const;
template void
MGTransferBlockSelect<float>::copy_from_mg_add (
  const MGDoFHandler<deal_II_dimension>&,
  BlockVector<double>&,
  const MGLevelObject<Vector<float> >&) const;

template void
MGTransferBlockSelect<double>::copy_to_mg (
  const MGDoFHandler<deal_II_dimension>&,
  MGLevelObject<Vector<double> >&,
  const Vector<double>&) const;
template void
MGTransferBlockSelect<double>::copy_to_mg (
  const MGDoFHandler<deal_II_dimension>&,
  MGLevelObject<Vector<double> >&,
  const BlockVector<double>&) const;
template void
MGTransferBlockSelect<double>::copy_from_mg (
  const MGDoFHandler<deal_II_dimension>&,
  Vector<double>&,
  const MGLevelObject<Vector<double> >&) const;
template void
MGTransferBlockSelect<double>::copy_from_mg (
  const MGDoFHandler<deal_II_dimension>&,
  BlockVector<double>&,
  const MGLevelObject<Vector<double> >&) const;
template void
MGTransferBlockSelect<double>::copy_from_mg_add (
  const MGDoFHandler<deal_II_dimension>&,
  Vector<double>&,
  const MGLevelObject<Vector<double> >&) const;
template void
MGTransferBlockSelect<double>::copy_from_mg_add (
  const MGDoFHandler<deal_II_dimension>&,
  BlockVector<double>&,
  const MGLevelObject<Vector<double> >&) const;

DEAL_II_NAMESPACE_CLOSE
