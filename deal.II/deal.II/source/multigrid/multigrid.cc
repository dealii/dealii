//-----------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------------------


#include <lac/vector.h>
#include <lac/sparse_matrix.h>
#include <lac/block_sparse_matrix.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_tools.h>
#include <fe/fe.h>
#include <multigrid/mg_dof_handler.h>
#include <multigrid/mg_dof_accessor.h>
#include <multigrid/mg_transfer.h>
#include <multigrid/mg_dof_tools.h>


/* ----------------------------- MGTransferPrebuilt ------------------------ */


template <int dim>
void MGTransferPrebuilt::build_matrices (const MGDoFHandler<dim> &mg_dof) 
{
  const unsigned int n_levels      = mg_dof.get_tria().n_levels();
  const unsigned int dofs_per_cell = mg_dof.get_fe().dofs_per_cell;
  
  sizes.resize(n_levels);
  for (unsigned int l=0;l<n_levels;++l)
    sizes[l] = mg_dof.n_dofs(l);
  
				   // reset the size of the array of
				   // matrices
  prolongation_sparsities.resize (n_levels-1);
  prolongation_matrices.resize (n_levels-1);


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
      prolongation_sparsities[level].reinit (sizes[level+1],
					      sizes[level],
//TODO:[WB,GK] evil hack, must be corrected!
					      dofs_per_cell+1);
      
      for (typename MGDoFHandler<dim>::cell_iterator cell=mg_dof.begin(level);
	   cell != mg_dof.end(level); ++cell)
	if (cell->has_children())
	  {
	    cell->get_mg_dof_indices (dof_indices_parent);

	    for (unsigned int child=0;
		 child<GeometryInfo<dim>::children_per_cell; ++child)
	      {
						 // set an alias to the
						 // prolongation matrix for
						 // this child
		const FullMatrix<double> &prolongation
		  = mg_dof.get_fe().prolongate(child);
	    
		cell->child(child)->get_mg_dof_indices (dof_indices_child);

						 // now tag the entries in the
						 // matrix which will be used
						 // for this pair of parent/child
		for (unsigned int i=0; i<dofs_per_cell; ++i)
		  for (unsigned int j=0; j<dofs_per_cell; ++j)
		    if (prolongation(i,j) != 0)
		      {
			prolongation_sparsities[level].add
			  (dof_indices_child[i],
			   dof_indices_parent[j]);
		      };
	      };
	  };
      prolongation_sparsities[level].compress ();

      prolongation_matrices[level].reinit (prolongation_sparsities[level]);

				       // now actually build the matrices
      for (typename MGDoFHandler<dim>::cell_iterator cell=mg_dof.begin(level);
	   cell != mg_dof.end(level); ++cell)
	if (cell->has_children())
	  {
	    cell->get_mg_dof_indices (dof_indices_parent);

	    for (unsigned int child=0;
		 child<GeometryInfo<dim>::children_per_cell; ++child)
	      {
						 // set an alias to the
						 // prolongation matrix for
						 // this child
		const FullMatrix<double> &prolongation
		  = mg_dof.get_fe().prolongate(child);
	    
		cell->child(child)->get_mg_dof_indices (dof_indices_child);

						 // now set the entries in the
						 // matrix
		for (unsigned int i=0; i<dofs_per_cell; ++i)
		  for (unsigned int j=0; j<dofs_per_cell; ++j)
		    if (prolongation(i,j) != 0)
		      {
			prolongation_matrices[level].set
			  (dof_indices_child[i],
			   dof_indices_parent[j],
			   prolongation(i,j));
		      };
	      };
	  };
    };
};


/* ----------------------------- MGTransferBlock ------------------------ */


template <int dim>
void MGTransferBlockBase::build_matrices (const MGDoFHandler<dim> &mg_dof,
				      std::vector<bool> select) 
{
  const FiniteElement<dim>& fe = mg_dof.get_fe();
  const unsigned int n_components  = fe.n_components();
  const unsigned int dofs_per_cell = fe.dofs_per_cell;  
  const unsigned int n_levels      = mg_dof.get_tria().n_levels();
  
  selected = select;
  
  Assert (selected.size() == n_components,
	  ExcDimensionMismatch(selected.size(), n_components))

  sizes.resize(n_levels);
  MGTools::count_dofs_per_component(mg_dof, sizes);
  
				   // reset the size of the array of
				   // matrices
  prolongation_sparsities.resize (n_levels-1);
  prolongation_matrices.resize (n_levels-1);

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
      prolongation_sparsities[level].reinit (n_components, n_components);
      for (unsigned int i=0;i<n_components;++i)
	for (unsigned int j=0;j<n_components;++j)
	  if (i==j)
	    prolongation_sparsities[level].block(i,j).reinit(
	      sizes[level+1][i],
	      sizes[level][j],
//TODO:[GK] Split this by component to save memory
	      dofs_per_cell);
	  else
	    prolongation_sparsities[level].block(i,j).reinit(
	      sizes[level+1][i],
	      sizes[level][j],
	      0);

      prolongation_sparsities[level].collect_sizes();
      
      for (typename MGDoFHandler<dim>::cell_iterator cell=mg_dof.begin(level);
	   cell != mg_dof.end(level); ++cell)
	if (cell->has_children())
	  {
	    cell->get_mg_dof_indices (dof_indices_parent);

	    for (unsigned int child=0;
		 child<GeometryInfo<dim>::children_per_cell; ++child)
	      {
						 // set an alias to the
						 // prolongation matrix for
						 // this child
		const FullMatrix<double> &prolongation
		  = mg_dof.get_fe().prolongate(child);
	    
		cell->child(child)->get_mg_dof_indices (dof_indices_child);

						 // now tag the entries in the
						 // matrix which will be used
						 // for this pair of parent/child
		for (unsigned int i=0; i<dofs_per_cell; ++i)
		  for (unsigned int j=0; j<dofs_per_cell; ++j)
		    if (prolongation(i,j) != 0)
		      {
			const unsigned int icomp = fe.system_to_component_index(i).first;
			const unsigned int jcomp = fe.system_to_component_index(j).first;
			if ((icomp==jcomp) && selected[icomp])
			  prolongation_sparsities[level].add
			    (dof_indices_child[i],
			     dof_indices_parent[j]);
		      };
	      };
	  };
      prolongation_sparsities[level].compress ();

      prolongation_matrices[level].reinit (prolongation_sparsities[level]);

				       // now actually build the matrices
      for (typename MGDoFHandler<dim>::cell_iterator cell=mg_dof.begin(level);
	   cell != mg_dof.end(level); ++cell)
	if (cell->has_children())
	  {
	    cell->get_mg_dof_indices (dof_indices_parent);

	    for (unsigned int child=0;
		 child<GeometryInfo<dim>::children_per_cell; ++child)
	      {
						 // set an alias to the
						 // prolongation matrix for
						 // this child
		const FullMatrix<double> &prolongation
		  = mg_dof.get_fe().prolongate(child);
	    
		cell->child(child)->get_mg_dof_indices (dof_indices_child);

						 // now set the entries in the
						 // matrix
		for (unsigned int i=0; i<dofs_per_cell; ++i)
		  for (unsigned int j=0; j<dofs_per_cell; ++j)
		    if (prolongation(i,j) != 0)
		      {
			const unsigned int icomp = fe.system_to_component_index(i).first;
			const unsigned int jcomp = fe.system_to_component_index(j).first;
			if ((icomp==jcomp) && selected[icomp])
			  prolongation_matrices[level].set
			    (dof_indices_child[i],
			     dof_indices_parent[j],
			     prolongation(i,j));
		      };
	      };
	  };
    };
				   // Finally, fill some index vectors
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

  component_start.resize(n_components);
  DoFTools::count_dofs_per_component(mg_dof, component_start);

  unsigned int k=0;
  for (unsigned int i=0;i<component_start.size();++i)
    {
      const unsigned int t=component_start[i];
      component_start[i] = k;
      k += t;
    }
  
};


template <int dim>
void MGTransferBlock::build_matrices (const MGDoFHandler<dim> &mg_dof,
					  std::vector<bool> select) 
{
  if (select.size() == 0)
    select = std::vector<bool> (mg_dof.get_fe().n_components(), true);

  MGTransferBlockBase::build_matrices (mg_dof, select);
}


template <int dim>
void MGTransferSelect::build_matrices (const MGDoFHandler<dim> &mg_dof,
					 unsigned int select)
{
  selected = select;
  std::vector<bool> s(mg_dof.get_fe().n_components(), false);
  s[select] = true;

  MGTransferBlockBase::build_matrices (mg_dof, s);
}



// explicit instatations

template
void MGTransferPrebuilt::build_matrices (const MGDoFHandler<deal_II_dimension> &mg_dof);

template
void MGTransferBlock::build_matrices (const MGDoFHandler<deal_II_dimension> &mg_dof,
				      std::vector<bool>);

template
void MGTransferSelect::build_matrices (const MGDoFHandler<deal_II_dimension> &mg_dof,
				       unsigned int);
