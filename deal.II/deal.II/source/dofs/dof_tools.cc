// $Id$
// Copyright Guido Kanschat, 1999

#include <grid/dof.h>
#include <grid/mg_dof.h>
#include <grid/dof_accessor.h>
#include <grid/mg_dof_accessor.h>
#include <grid/tria_iterator.h>
#include <fe/fe.h>
#include <fe/fe_system.h>
#include <basic/dof_tools.h>
#include <lac/sparsematrix.h>

#if deal_II_dimension == 1

//TODO: Implement make_flux_sparsity_pattern in 1d

template <>
void
DoFTools::make_flux_sparsity_pattern (const DoFHandler<1>&,
				      SparseMatrixStruct &)
{
  Assert(false, ExcNotImplemented());
}

#else

//TODO: Check this function for potential of optimization. (G)

template<int dim>
void
DoFTools::make_flux_sparsity_pattern (const DoFHandler<dim>& dof,
				      SparseMatrixStruct &sparsity)
{
  const unsigned int n_dofs = dof.n_dofs();
  
  Assert (sparsity.n_rows() == n_dofs,
	  ExcDimensionMismatch (sparsity.n_rows(), n_dofs));
  Assert (sparsity.n_cols() == n_dofs,
	  ExcDimensionMismatch (sparsity.n_cols(), n_dofs));

  const unsigned int total_dofs = dof.get_fe().dofs_per_cell;
  vector<int> dofs_on_this_cell(total_dofs);
  vector<int> dofs_on_other_cell(total_dofs);
  DoFHandler<dim>::active_cell_iterator cell = dof.begin_active(),
					endc = dof.end();
  for (; cell!=endc; ++cell) 
    {
      cell->get_dof_indices (dofs_on_this_cell);
				       // make sparsity pattern for this cell
      for (unsigned int i=0; i<total_dofs; ++i)
	for (unsigned int j=0; j<total_dofs; ++j)
	  sparsity.add (dofs_on_this_cell[i],
			dofs_on_this_cell[j]);

				       // Loop over all interior neighbors
      for (unsigned int face = 0;
	   face < GeometryInfo<dim>::faces_per_cell;
	   ++face)
	{
	  if (cell->face(face)->boundary_indicator() == 255)
	    {
	      DoFHandler<dim>::active_cell_iterator neighbor = cell->neighbor(face);
	      neighbor->get_dof_indices (dofs_on_other_cell);
					       // only add one direction
					       // The other is taken care of
					       // by neighbor.
	      for (unsigned int i=0; i<total_dofs; ++i)
		{
		  for (unsigned int j=0; j<total_dofs; ++j)
		    {
		      sparsity.add (dofs_on_this_cell[i],
				    dofs_on_other_cell[j]);
		    }
		}
	    }
	} 
    }
}

#endif

template<int dim>
void
DoFTools::extract_dofs(const DoFHandler<dim> &dof,
		       const vector<bool>    &local_select,
		       vector<bool>          &selected_dofs)
{
  const FiniteElement<dim> &fe = dof.get_fe();
  Assert(local_select.size() == fe.n_components(),
	 ExcDimensionMismatch(local_select.size(), fe.n_components()));
  Assert(selected_dofs.size() == dof.n_dofs(),
	 ExcDimensionMismatch(selected_dofs.size(), dof.n_dofs()));

				   // preset all values by false
  fill_n (selected_dofs.begin(), dof.n_dofs(), false);
  
  vector<int> indices(fe.dofs_per_cell);
  
  DoFHandler<dim>::active_cell_iterator c;
  for (c = dof.begin_active() ; c != dof.end() ; ++ c)
    {
      c->get_dof_indices(indices);
      for (unsigned int i=0;i<fe.dofs_per_cell;++i)
	{
	  const unsigned int component = fe.system_to_component_index(i).first;

	  if (local_select[component] == true)
	    selected_dofs[indices[i]] = true;
	}
    }
}



template<int dim>
void
DoFTools::extract_level_dofs(const unsigned int       level,
			     const MGDoFHandler<dim> &dof,
			     const vector<bool>      &local_select,
			     vector<bool>            &selected_dofs)
{
  const FiniteElement<dim>& fe = dof.get_fe();
  Assert(local_select.size() == fe.n_components(),
	 ExcDimensionMismatch(local_select.size(), fe.n_components()));
  Assert(selected_dofs.size() == dof.n_dofs(level),
	 ExcDimensionMismatch(selected_dofs.size(), dof.n_dofs(level)));

				   // preset all values by false
  fill_n (selected_dofs.begin(), dof.n_dofs(level), false);

  vector<int> indices(fe.dofs_per_cell);
  
  MGDoFHandler<dim>::cell_iterator c;
  for (c = dof.begin(level) ; c != dof.end(level) ; ++ c)
    {
      c->get_mg_dof_indices(indices);
      for (unsigned int i=0;i<fe.dofs_per_cell;++i)
	{
	  const unsigned int component = fe.system_to_component_index(i).first;
	  if (local_select[component]  == true)
	    selected_dofs[indices[i]] = true;
	}
    }
}



// explicit instantiations
#if deal_II_dimension > 1
template void DoFTools::make_flux_sparsity_pattern (const DoFHandler<deal_II_dimension>& dof,
						    SparseMatrixStruct &sparsity);
#endif
template void DoFTools::extract_dofs(const DoFHandler<deal_II_dimension>& dof,
				     const vector<bool>& local_select,
				     vector<bool>& selected_dofs);
template void DoFTools::extract_level_dofs(unsigned int level,
					   const MGDoFHandler<deal_II_dimension>& dof,
					   const vector<bool>& local_select,
					   vector<bool>& selected_dofs);

