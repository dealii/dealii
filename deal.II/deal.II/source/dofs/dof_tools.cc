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



template<int dim>
void
DoFTools::extract_dofs(const DoFHandler<dim> &dof,
		       const vector<bool>    &local_select,
		       vector<bool>          &selected_dofs)
{
  const FiniteElement<dim> &fe = dof.get_fe();
  Assert(local_select.size() == fe.n_components,
	 ExcDimensionMismatch(local_select.size(), fe.n_components));
  Assert(selected_dofs.size() == dof.n_dofs(),
	 ExcDimensionMismatch(selected_dofs.size(), dof.n_dofs()));

  vector<int> indices(fe.total_dofs);
  
  DoFHandler<dim>::active_cell_iterator c;
  for (c = dof.begin_active() ; c != dof.end() ; ++ c)
    {
      c->get_dof_indices(indices);
      for (unsigned int i=0;i<fe.total_dofs;++i)
	{
	  const pair<unsigned int, unsigned int> component
	    = fe.system_to_component_index(i).first;

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
  Assert(local_select.size() == fe.n_components,
	 ExcDimensionMismatch(local_select.size(), fe.n_components));
  Assert(selected_dofs.size() == dof.n_dofs(level),
	 ExcDimensionMismatch(selected_dofs.size(), dof.n_dofs(level)));

  vector<int> indices(fe.total_dofs);
  
  MGDoFHandler<dim>::cell_iterator c;
  for (c = dof.begin(level) ; c != dof.end(level) ; ++ c)
    {
      c->get_mg_dof_indices(indices);
      for (unsigned int i=0;i<fe.total_dofs;++i)
	{
	  const pair<unsigned int, unsigned int> component
	    = fe.system_to_component_index(i).first;
	  if (local_select[component]  == true)
	    selected_dofs[indices[i]] = true;
	}
    }
}



// explicit instantiations
template void DoFTools::extract_dofs(const DoFHandler<deal_II_dimension>& dof,
				     const vector<bool>& local_select,
				     vector<bool>& selected_dofs);
template void DoFTools::extract_level_dofs(unsigned int level,
					   const MGDoFHandler<deal_II_dimension>& dof,
					   const vector<bool>& local_select,
					   vector<bool>& selected_dofs);

