//----------------------------  mg_dof_tools.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  mg_dof_tools.cc  ---------------------------


#include <lac/sparsity_pattern.h>
#include <multigrid/mg_dof_handler.h>
#include <multigrid/mg_dof_accessor.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <multigrid/mg_dof_tools.h>
#include <fe/fe.h>


template <int dim>
void MGDoFTools::make_sparsity_pattern (const MGDoFHandler<dim> &dof,
					SparsityPattern         &sparsity,
					const unsigned int       level)
{
  const unsigned int n_dofs = dof.n_dofs(level);

  Assert (sparsity.n_rows() == n_dofs,
	  ExcDimensionMismatch (sparsity.n_rows(), n_dofs));
  Assert (sparsity.n_cols() == n_dofs,
	  ExcDimensionMismatch (sparsity.n_cols(), n_dofs));

  const unsigned int dofs_per_cell = dof.get_fe().dofs_per_cell;
  vector<unsigned int> dofs_on_this_cell(dofs_per_cell);
  MGDoFHandler<dim>::cell_iterator cell = dof.begin(level),
				   endc = dof.end(level);
  for (; cell!=endc; ++cell) 
    {
      cell->get_mg_dof_indices (dofs_on_this_cell);
				       // make sparsity pattern for this cell
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	for (unsigned int j=0; j<dofs_per_cell; ++j)
	  sparsity.add (dofs_on_this_cell[i],
			dofs_on_this_cell[j]);
    }
}

template<int dim>
void
MGDoFTools::make_flux_sparsity_pattern (const MGDoFHandler<dim> &dof,
					SparsityPattern       &sparsity,
					const unsigned int level)
{
  const unsigned int n_dofs = dof.n_dofs(level);
  
  Assert (sparsity.n_rows() == n_dofs,
	  ExcDimensionMismatch (sparsity.n_rows(), n_dofs));
  Assert (sparsity.n_cols() == n_dofs,
	  ExcDimensionMismatch (sparsity.n_cols(), n_dofs));

  const unsigned int dofs_per_cell = dof.get_fe().dofs_per_cell;
  vector<unsigned int> dofs_on_this_cell(dofs_per_cell);
  vector<unsigned int> dofs_on_other_cell(dofs_per_cell);
  MGDoFHandler<dim>::cell_iterator cell = dof.begin(level),
				   endc = dof.end(level);
  for (; cell!=endc; ++cell)
    {
      cell->get_mg_dof_indices (dofs_on_this_cell);
				       // make sparsity pattern for this cell
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	for (unsigned int j=0; j<dofs_per_cell; ++j)
	  sparsity.add (dofs_on_this_cell[i],
			dofs_on_this_cell[j]);

				       // Loop over all interior neighbors
      for (unsigned int face = 0;
	   face < GeometryInfo<dim>::faces_per_cell;
	   ++face)
	{
	  if ( (! cell->at_boundary(face)) &&
	       (static_cast<unsigned int>(cell->neighbor_level(face)) == level) )
	    {
	      MGDoFHandler<dim>::cell_iterator neighbor = cell->neighbor(face);
	      neighbor->get_mg_dof_indices (dofs_on_other_cell);
					       // only add one direction
					       // The other is taken care of
					       // by neighbor.
	      for (unsigned int i=0; i<dofs_per_cell; ++i)
		{
		  for (unsigned int j=0; j<dofs_per_cell; ++j)
		    {
		      sparsity.add (dofs_on_this_cell[i],
				    dofs_on_other_cell[j]);
		    }
		}
	    }
	} 
    }
}



template<int dim>
void
MGDoFTools::make_flux_sparsity_pattern (const MGDoFHandler<dim> &dof,
					SparsityPattern       &sparsity,
					const unsigned int level,
					const FullMatrix<double>& int_mask,
					const FullMatrix<double>& flux_mask)
{
  const unsigned int n_dofs = dof.n_dofs(level);
  const unsigned int n_comp = dof.get_fe().n_components();
  
  Assert (sparsity.n_rows() == n_dofs,
	  ExcDimensionMismatch (sparsity.n_rows(), n_dofs));
  Assert (sparsity.n_cols() == n_dofs,
	  ExcDimensionMismatch (sparsity.n_cols(), n_dofs));
  Assert (int_mask.m() == n_comp,
	  ExcDimensionMismatch (int_mask.m(), n_comp));
  Assert (int_mask.n() == n_comp,
	  ExcDimensionMismatch (int_mask.n(), n_comp));
  Assert (flux_mask.m() == n_comp,
	  ExcDimensionMismatch (flux_mask.m(), n_comp));
  Assert (flux_mask.n() == n_comp,
	  ExcDimensionMismatch (flux_mask.n(), n_comp));
  
  const unsigned int total_dofs = dof.get_fe().dofs_per_cell;
  vector<unsigned int> dofs_on_this_cell(total_dofs);
  vector<unsigned int> dofs_on_other_cell(total_dofs);
  MGDoFHandler<dim>::cell_iterator cell = dof.begin(level),
				   endc = dof.end(level);


  vector<vector<bool> > int_dof_mask(total_dofs,
				 vector<bool>(total_dofs, false));
  vector<vector<bool> > flux_dof_mask(total_dofs,
				 vector<bool>(total_dofs, false));
  for (unsigned int i=0; i<total_dofs; ++i)
    for (unsigned int j=0; j<total_dofs; ++j)
      {
	unsigned int ii = dof.get_fe().system_to_component_index(i).first;
	unsigned int jj = dof.get_fe().system_to_component_index(j).first;
	
	if (int_mask(ii,jj) != 0)
	  int_dof_mask[i][j] = true;
	if (flux_mask(ii,jj) != 0)
	  flux_dof_mask[i][j] = true;
      }
  
  (const_cast<Triangulation<dim>& > (dof.get_tria())).clear_user_flags();
  
  for (; cell!=endc; ++cell)
    {
      cell->get_mg_dof_indices (dofs_on_this_cell);
				       // make sparsity pattern for this cell
      for (unsigned int i=0; i<total_dofs; ++i)
	for (unsigned int j=0; j<total_dofs; ++j)
	  if (int_dof_mask[i][j])
	    sparsity.add (dofs_on_this_cell[i],
			  dofs_on_this_cell[j]);

				       // Loop over all interior neighbors
      for (unsigned int face = 0;
	   face < GeometryInfo<dim>::faces_per_cell;
	   ++face)
	{
	  MGDoFHandler<dim>::face_iterator cell_face = cell->face(face);
	  if (cell_face->user_flag_set ())
	    continue;

	  if (! cell->at_boundary (face) )
	    {
	      MGDoFHandler<dim>::cell_iterator neighbor = cell->neighbor(face);

	      if (neighbor->level() < cell->level())
		continue;

	      unsigned int neighbor_face = cell->neighbor_of_neighbor(face);

	      neighbor->get_mg_dof_indices (dofs_on_other_cell);
	      for (unsigned int i=0; i<total_dofs; ++i)
		{
		  for (unsigned int j=0; j<total_dofs; ++j)
		    if (flux_dof_mask[i][j])
		      {
			sparsity.add (dofs_on_this_cell[i],
				      dofs_on_other_cell[j]);
			sparsity.add (dofs_on_other_cell[i],
				      dofs_on_this_cell[j]);
		      }
		}
	      neighbor->face(neighbor_face)->set_user_flag (); 
	    }
	}
    }
}




// explicit instantiations
template void MGDoFTools::make_sparsity_pattern (const MGDoFHandler<deal_II_dimension> &,
						 SparsityPattern &,
						 const unsigned int);
template void MGDoFTools::make_flux_sparsity_pattern (const MGDoFHandler<deal_II_dimension> &,
						      SparsityPattern &,
						      const unsigned int);

#if deal_II_dimension > 1
template void MGDoFTools::make_flux_sparsity_pattern (const MGDoFHandler<deal_II_dimension> &,
						      SparsityPattern &,
						      const unsigned int,
						      const FullMatrix<double>&,
						      const FullMatrix<double>&);
#endif
