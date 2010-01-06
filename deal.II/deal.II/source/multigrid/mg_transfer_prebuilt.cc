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
#include <base/function.h>

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
#include <multigrid/mg_tools.h>
#include <multigrid/mg_transfer.h>
#include <multigrid/mg_transfer.templates.h>

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

      typedef std::map<unsigned int, unsigned int>::const_iterator IT;
      for (IT i= copy_indices[level].begin();
	   i != copy_indices[level].end();++i)
	dst[level](i->second) = src(i->first);
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



template <typename VECTOR>
template <int dim, int spacedim>
void MGTransferPrebuilt<VECTOR>::build_matrices (
  const MGDoFHandler<dim,spacedim> &mg_dof)
{
  //start building the matrices here
  const unsigned int n_levels      = mg_dof.get_tria().n_levels();
  const unsigned int dofs_per_cell = mg_dof.get_fe().dofs_per_cell;

  sizes.resize(n_levels);
  for (unsigned int l=0;l<n_levels;++l)
    sizes[l] = mg_dof.n_dofs(l);

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
	.push_back (std_cxx1x::shared_ptr<SparsityPattern> (new SparsityPattern));
      prolongation_matrices
	.push_back (std_cxx1x::shared_ptr<SparseMatrix<double> > (new SparseMatrix<double>));
    }

				   // two fields which will store the
				   // indices of the multigrid dofs
				   // for a cell and one of its children
  std::vector<unsigned int> dof_indices_parent (dofs_per_cell);
  std::vector<unsigned int> dof_indices_child (dofs_per_cell);

				   // for each level: first build the sparsity
				   // pattern of the matrices and then build the
				   // matrices themselves. note that we only
				   // need to take care of cells on the coarser
				   // level which have children
  for (unsigned int level=0; level<n_levels-1; ++level)
    {
				       // reset the dimension of the structure.
				       // note that for the number of entries
				       // per row, the number of parent dofs
				       // coupling to a child dof is
				       // necessary. this, of course, is the
				       // number of degrees of freedom per
				       // cell
				       // increment dofs_per_cell
				       // since a useless diagonal
				       // element will be stored
      prolongation_sparsities[level]->reinit (sizes[level+1],
					      sizes[level],
					      dofs_per_cell+1);

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

		Assert (prolongation.n() != 0, ExcNoProlongation());

		cell->child(child)->get_mg_dof_indices (dof_indices_child);

						 // now tag the entries in the
						 // matrix which will be used
						 // for this pair of parent/child
		for (unsigned int i=0; i<dofs_per_cell; ++i)
		  for (unsigned int j=0; j<dofs_per_cell; ++j)
		    if (prolongation(i,j) != 0)
		      {
			prolongation_sparsities[level]->add (dof_indices_child[i],
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
		      prolongation_matrices[level]->set (dof_indices_child[i],
							 dof_indices_parent[j],
							 prolongation(i,j));
	      }
	  }
    }
  // impose boundary conditions
  // but only in the column of
  // the prolongation matrix
  //TODO: this way is not very efficient
  std::vector<std::set<unsigned int> > boundary_indices (mg_dof.get_tria().n_levels());
  typename FunctionMap<dim>::type boundary;
  ZeroFunction<dim>                    homogeneous_dirichlet_bc (1);
  boundary[0] = &homogeneous_dirichlet_bc;
  MGTools::make_boundary_list(mg_dof, boundary, boundary_indices);

  for (unsigned int level=0; level<n_levels-1; ++level)
  {
    if (boundary_indices[level].size() == 0)
      continue;

    const unsigned int n_dofs = prolongation_matrices[level]->m();
    for (unsigned int i=0; i<n_dofs; ++i)
    {
      SparseMatrix<double>::iterator anfang = prolongation_matrices[level]->begin(i),
        ende = prolongation_matrices[level]->end(i);
      for(; anfang != ende; ++anfang)
      {
        std::set<unsigned int>::const_iterator dof = boundary_indices[level].begin(),
          endd = boundary_indices[level].end();
        for (; dof != endd; ++dof)
        {
          const unsigned int column_number = anfang->column();
          const unsigned int dof_number = *dof;
          if(dof_number == column_number)
            prolongation_matrices[level]->set(i,dof_number,0);
        }
      }
    }
  }

  find_dofs_on_refinement_edges (mg_dof);
				   // Prepare indices to copy from a
				   // global vector to the parts of an
				   // MGVector corresponding to active
				   // cells.
  copy_indices.resize(mg_dof.get_tria().n_levels());
  std::vector<unsigned int> global_dof_indices (dofs_per_cell);
  std::vector<unsigned int> level_dof_indices  (dofs_per_cell);
  for (int level=mg_dof.get_tria().n_levels()-1; level>=0; --level)
    {
      copy_indices[level].clear();
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

	  for (unsigned int i=0; i<dofs_per_cell; ++i)
	    {
              if(!dofs_on_refinement_edge[level][level_dof_indices[i]])
	      copy_indices[level].insert(
		std::make_pair(global_dof_indices[i], level_dof_indices[i]));
	    }
	}
    }
}


//TODO: Use template expander script

template
void MGTransferPrebuilt<Vector<float> >::build_matrices<deal_II_dimension>
(const MGDoFHandler<deal_II_dimension> &mg_dof);

template
void MGTransferPrebuilt<Vector<double> >::build_matrices<deal_II_dimension>
(const MGDoFHandler<deal_II_dimension> &mg_dof);

template
void MGTransferPrebuilt<BlockVector<float> >::build_matrices<deal_II_dimension>
(const MGDoFHandler<deal_II_dimension> &mg_dof);

template
void MGTransferPrebuilt<BlockVector<double> >::build_matrices<deal_II_dimension>
(const MGDoFHandler<deal_II_dimension> &mg_dof);

template void
MGTransferPrebuilt<Vector<float> >::copy_to_mg (
  const MGDoFHandler<deal_II_dimension>&,
  MGLevelObject<Vector<float> >&,
  const Vector<double>&) const;
template void
MGTransferPrebuilt<Vector<float> >::copy_to_mg (
  const MGDoFHandler<deal_II_dimension>&,
  MGLevelObject<Vector<float> >&,
  const BlockVector<double>&) const;
template void
MGTransferPrebuilt<Vector<float> >::copy_from_mg (
  const MGDoFHandler<deal_II_dimension>&,
  Vector<double>&,
  const MGLevelObject<Vector<float> >&) const;
template void
MGTransferPrebuilt<Vector<float> >::copy_from_mg (
  const MGDoFHandler<deal_II_dimension>&,
  BlockVector<double>&,
  const MGLevelObject<Vector<float> >&) const;
template void
MGTransferPrebuilt<Vector<float> >::copy_from_mg_add (
  const MGDoFHandler<deal_II_dimension>&,
  Vector<double>&,
  const MGLevelObject<Vector<float> >&) const;
template void
MGTransferPrebuilt<Vector<float> >::copy_from_mg_add (
  const MGDoFHandler<deal_II_dimension>&,
  BlockVector<double>&,
  const MGLevelObject<Vector<float> >&) const;

template void
MGTransferPrebuilt<BlockVector<float> >::copy_to_mg (
  const MGDoFHandler<deal_II_dimension>&,
  MGLevelObject<BlockVector<float> >&,
  const Vector<double>&) const;
template void
MGTransferPrebuilt<BlockVector<float> >::copy_to_mg (
  const MGDoFHandler<deal_II_dimension>&,
  MGLevelObject<BlockVector<float> >&,
  const BlockVector<double>&) const;
template void
MGTransferPrebuilt<BlockVector<float> >::copy_from_mg (
  const MGDoFHandler<deal_II_dimension>&,
  Vector<double>&,
  const MGLevelObject<BlockVector<float> >&) const;
template void
MGTransferPrebuilt<BlockVector<float> >::copy_from_mg (
  const MGDoFHandler<deal_II_dimension>&,
  BlockVector<double>&,
  const MGLevelObject<BlockVector<float> >&) const;
template void
MGTransferPrebuilt<BlockVector<float> >::copy_from_mg_add (
  const MGDoFHandler<deal_II_dimension>&,
  Vector<double>&,
  const MGLevelObject<BlockVector<float> >&) const;
template void
MGTransferPrebuilt<BlockVector<float> >::copy_from_mg_add (
  const MGDoFHandler<deal_II_dimension>&,
  BlockVector<double>&,
  const MGLevelObject<BlockVector<float> >&) const;

template void
MGTransferPrebuilt<Vector<double> >::copy_to_mg (
  const MGDoFHandler<deal_II_dimension>&,
  MGLevelObject<Vector<double> >&,
  const Vector<double>&) const;
template void
MGTransferPrebuilt<Vector<double> >::copy_to_mg (
  const MGDoFHandler<deal_II_dimension>&,
  MGLevelObject<Vector<double> >&,
  const BlockVector<double>&) const;
template void
MGTransferPrebuilt<Vector<double> >::copy_from_mg (
  const MGDoFHandler<deal_II_dimension>&,
  Vector<double>&,
  const MGLevelObject<Vector<double> >&) const;
template void
MGTransferPrebuilt<Vector<double> >::copy_from_mg (
  const MGDoFHandler<deal_II_dimension>&,
  BlockVector<double>&,
  const MGLevelObject<Vector<double> >&) const;
template void
MGTransferPrebuilt<Vector<double> >::copy_from_mg_add (
  const MGDoFHandler<deal_II_dimension>&,
  Vector<double>&,
  const MGLevelObject<Vector<double> >&) const;
template void
MGTransferPrebuilt<Vector<double> >::copy_from_mg_add (
  const MGDoFHandler<deal_II_dimension>&,
  BlockVector<double>&,
  const MGLevelObject<Vector<double> >&) const;

template void
MGTransferPrebuilt<BlockVector<double> >::copy_to_mg (
  const MGDoFHandler<deal_II_dimension>&,
  MGLevelObject<BlockVector<double> >&,
  const Vector<double>&) const;
template void
MGTransferPrebuilt<BlockVector<double> >::copy_to_mg (
  const MGDoFHandler<deal_II_dimension>&,
  MGLevelObject<BlockVector<double> >&,
  const BlockVector<double>&) const;
template void
MGTransferPrebuilt<BlockVector<double> >::copy_from_mg (
  const MGDoFHandler<deal_II_dimension>&,
  Vector<double>&,
  const MGLevelObject<BlockVector<double> >&) const;
template void
MGTransferPrebuilt<BlockVector<double> >::copy_from_mg (
  const MGDoFHandler<deal_II_dimension>&,
  BlockVector<double>&,
  const MGLevelObject<BlockVector<double> >&) const;
template void
MGTransferPrebuilt<BlockVector<double> >::copy_from_mg_add (
  const MGDoFHandler<deal_II_dimension>&,
  Vector<double>&,
  const MGLevelObject<BlockVector<double> >&) const;
template void
MGTransferPrebuilt<BlockVector<double> >::copy_from_mg_add (
  const MGDoFHandler<deal_II_dimension>&,
  BlockVector<double>&,
  const MGLevelObject<BlockVector<double> >&) const;

DEAL_II_NAMESPACE_CLOSE
