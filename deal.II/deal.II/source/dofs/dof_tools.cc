//----------------------------  dof_tools.cc  ---------------------------
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
//----------------------------  dof_tools.cc  ---------------------------


#include <base/multithread_info.h>
#include <base/thread_management.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <grid/intergrid_map.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_constraints.h>
#include <multigrid/mg_dof_handler.h>
#include <multigrid/mg_dof_accessor.h>
#include <fe/fe.h>
#include <fe/fe_system.h>
#include <dofs/dof_tools.h>
#include <lac/sparsity_pattern.h>
#include <lac/block_sparsity_pattern.h>
#include <lac/vector.h>

#include <algorithm>


template <int dim, class SparsityPattern>
void
DoFTools::make_sparsity_pattern (const DoFHandler<dim> &dof,
				 SparsityPattern       &sparsity)
{
  const unsigned int n_dofs = dof.n_dofs();

  Assert (sparsity.n_rows() == n_dofs,
	  ExcDimensionMismatch (sparsity.n_rows(), n_dofs));
  Assert (sparsity.n_cols() == n_dofs,
	  ExcDimensionMismatch (sparsity.n_cols(), n_dofs));

  const unsigned int dofs_per_cell = dof.get_fe().dofs_per_cell;
  vector<unsigned int> dofs_on_this_cell(dofs_per_cell);
  DoFHandler<dim>::active_cell_iterator cell = dof.begin_active(),
					endc = dof.end();
  for (; cell!=endc; ++cell) 
    {
      cell->get_dof_indices (dofs_on_this_cell);
				       // make sparsity pattern for this cell
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	for (unsigned int j=0; j<dofs_per_cell; ++j)
	  sparsity.add (dofs_on_this_cell[i],
			dofs_on_this_cell[j]);
    }
}


template <int dim, class SparsityPattern>
void
DoFTools::make_sparsity_pattern (const DoFHandler<dim>       &dof,
				 const vector<vector<bool> > &mask,
				 SparsityPattern             &sparsity)
{
  const unsigned int n_dofs = dof.n_dofs();
  const unsigned int dofs_per_cell = dof.get_fe().dofs_per_cell;

  Assert (sparsity.n_rows() == n_dofs,
	  ExcDimensionMismatch (sparsity.n_rows(), n_dofs));
  Assert (sparsity.n_cols() == n_dofs,
	  ExcDimensionMismatch (sparsity.n_cols(), n_dofs));
  Assert (mask.size() == dof.get_fe().n_components(),
	  ExcDimensionMismatch(mask.size(), dof.get_fe().n_components()));
  for (unsigned int i=0; i<mask.size(); ++i)
    Assert (mask[i].size() == dof.get_fe().n_components(),
	    ExcDimensionMismatch(mask[i].size(), dof.get_fe().n_components()));

				   // first build a mask for each dof,
				   // not like the one given which represents
				   // components
  vector<vector<bool> > dof_mask(dofs_per_cell,
				 vector<bool>(dofs_per_cell, false));
  for (unsigned int i=0; i<dofs_per_cell; ++i)
    for (unsigned int j=0; j<dofs_per_cell; ++j)
      dof_mask[i][j] = mask
		       [dof.get_fe().system_to_component_index(i).first]
		       [dof.get_fe().system_to_component_index(j).first];


  vector<unsigned int> dofs_on_this_cell(dofs_per_cell);
  DoFHandler<dim>::active_cell_iterator cell = dof.begin_active(),
		       endc = dof.end();
  for (; cell!=endc; ++cell) 
    {
      cell->get_dof_indices (dofs_on_this_cell);
				       // make sparsity pattern for this cell
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	for (unsigned int j=0; j<dofs_per_cell; ++j)
	  if (dof_mask[i][j] == true)
	    sparsity.add (dofs_on_this_cell[i],
			  dofs_on_this_cell[j]);
    }
}


#if deal_II_dimension == 1

template <>
void DoFTools::make_boundary_sparsity_pattern (const DoFHandler<1>&,
					       const vector<unsigned int>  &,
					       SparsityPattern    &)
{
    Assert (false, ExcInternalError());
};


template <>
void
DoFTools::make_boundary_sparsity_pattern (const DoFHandler<1>&,
					  const DoFHandler<1>::FunctionMap  &,
					  const vector<unsigned int>  &,
					  SparsityPattern    &)
{
  Assert (false, ExcInternalError());
}

#endif


template <int dim>
void
DoFTools::make_boundary_sparsity_pattern (const DoFHandler<dim>& dof,
					  const vector<unsigned int>  &dof_to_boundary_mapping,
					  SparsityPattern    &sparsity)
{
  const unsigned int n_dofs = dof.n_dofs();

  Assert (dof_to_boundary_mapping.size() == n_dofs, ExcInternalError());
  Assert (sparsity.n_rows() == dof.n_boundary_dofs(),
	  ExcDimensionMismatch (sparsity.n_rows(), dof.n_boundary_dofs()));
  Assert (sparsity.n_cols() == dof.n_boundary_dofs(),
	  ExcDimensionMismatch (sparsity.n_cols(), dof.n_boundary_dofs()));
  Assert (*max_element(dof_to_boundary_mapping.begin(),
		       dof_to_boundary_mapping.end()) == sparsity.n_rows()-1,
	  ExcInternalError());

  const unsigned int dofs_per_face = dof.get_fe().dofs_per_face;
  vector<unsigned int> dofs_on_this_face(dofs_per_face);
  DoFHandler<dim>::active_face_iterator face = dof.begin_active_face(),
		       endf = dof.end_face();
  for (; face!=endf; ++face)
    if (face->at_boundary())
      {
	face->get_dof_indices (dofs_on_this_face);

					 // make sure all dof indices have a
					 // boundary index
	Assert (*min_element(dofs_on_this_face.begin(),
			     dofs_on_this_face.end()) != DoFHandler<dim>::invalid_dof_index,
		ExcInternalError());
	
					 // make sparsity pattern for this cell
	for (unsigned int i=0; i<dofs_per_face; ++i)
	  for (unsigned int j=0; j<dofs_per_face; ++j) 
	    sparsity.add (dof_to_boundary_mapping[dofs_on_this_face[i]],
			  dof_to_boundary_mapping[dofs_on_this_face[j]]);
      };
};


template <int dim>
void DoFTools::make_boundary_sparsity_pattern (const DoFHandler<dim>& dof,
					       const DoFHandler<dim>::FunctionMap  &boundary_indicators,
					       const vector<unsigned int>  &dof_to_boundary_mapping,
					       SparsityPattern    &sparsity)
{
  const unsigned int n_dofs = dof.n_dofs();

  Assert (dof_to_boundary_mapping.size() == n_dofs, ExcInternalError());
  Assert (boundary_indicators.find(255) == boundary_indicators.end(),
	  DoFHandler<dim>::ExcInvalidBoundaryIndicator());
  Assert (sparsity.n_rows() == dof.n_boundary_dofs(boundary_indicators),
	  ExcDimensionMismatch (sparsity.n_rows(), dof.n_boundary_dofs(boundary_indicators)));
  Assert (sparsity.n_cols() == dof.n_boundary_dofs(boundary_indicators),
	  ExcDimensionMismatch (sparsity.n_cols(), dof.n_boundary_dofs(boundary_indicators)));
#ifdef DEBUG
  if (true)
    {
      unsigned int max_element = 0;
      for (vector<unsigned int>::const_iterator i=dof_to_boundary_mapping.begin();
	   i!=dof_to_boundary_mapping.end(); ++i)
	if ((*i != DoFHandler<dim>::invalid_dof_index) &&
	    (*i > max_element))
	  max_element = *i;
      Assert (max_element  == sparsity.n_rows()-1,
	      ExcInternalError());
    };
#endif

  const unsigned int dofs_per_face = dof.get_fe().dofs_per_face;
  vector<unsigned int> dofs_on_this_face(dofs_per_face);
  DoFHandler<dim>::active_face_iterator face = dof.begin_active_face(),
		       endf = dof.end_face();
  for (; face!=endf; ++face)
    if (boundary_indicators.find(face->boundary_indicator()) !=
	boundary_indicators.end())
      {
	face->get_dof_indices (dofs_on_this_face);

					 // make sure all dof indices have a
					 // boundary index
	Assert (*min_element(dofs_on_this_face.begin(),
			     dofs_on_this_face.end()) != DoFHandler<dim>::invalid_dof_index,
		ExcInternalError());
					 // make sparsity pattern for this cell
	for (unsigned int i=0; i<dofs_per_face; ++i)
	  for (unsigned int j=0; j<dofs_per_face; ++j)
	    sparsity.add (dof_to_boundary_mapping[dofs_on_this_face[i]],
			  dof_to_boundary_mapping[dofs_on_this_face[j]]);
      };
};


template<int dim>
void
DoFTools::make_flux_sparsity_pattern (const DoFHandler<dim> &dof,
				      SparsityPattern       &sparsity)
{
  const unsigned int n_dofs = dof.n_dofs();
  
  Assert (sparsity.n_rows() == n_dofs,
	  ExcDimensionMismatch (sparsity.n_rows(), n_dofs));
  Assert (sparsity.n_cols() == n_dofs,
	  ExcDimensionMismatch (sparsity.n_cols(), n_dofs));

  const unsigned int total_dofs = dof.get_fe().dofs_per_cell;
  vector<unsigned int> dofs_on_this_cell(total_dofs);
  vector<unsigned int> dofs_on_other_cell(total_dofs);
  DoFHandler<dim>::active_cell_iterator cell = dof.begin_active(),
					endc = dof.end();

  (const_cast<Triangulation<dim>& > (dof.get_tria())).clear_user_flags();
  
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
	  DoFHandler<dim>::face_iterator cell_face = cell->face(face);
	  if (cell_face->user_flag_set ())
	    continue;

	  if (! cell->at_boundary (face) )
	    {
	      DoFHandler<dim>::cell_iterator neighbor = cell->neighbor(face);
					       // Refinement edges are taken care of
					       // by coarser cells
	      if (neighbor->level() < cell->level())
		continue;

	      unsigned int neighbor_face = cell->neighbor_of_neighbor(face);

	      if (neighbor->has_children())
		{
		  for (unsigned int sub_nr = 0;
		       sub_nr != GeometryInfo<dim>::subfaces_per_face;
		       ++sub_nr)
		    {
		      DoFHandler<dim>::cell_iterator sub_neighbor
			= neighbor->
			child(GeometryInfo<dim>::child_cell_on_face(neighbor_face, sub_nr));

		      sub_neighbor->get_dof_indices (dofs_on_other_cell);
		      for (unsigned int i=0; i<total_dofs; ++i)
			{
			  for (unsigned int j=0; j<total_dofs; ++j)
			    {
			      sparsity.add (dofs_on_this_cell[i],
					    dofs_on_other_cell[j]);
			      sparsity.add (dofs_on_other_cell[i],
					    dofs_on_this_cell[j]);
			    }
			}
		      sub_neighbor->face(neighbor_face)->set_user_flag ();
		    }
		} else {
		  neighbor->get_dof_indices (dofs_on_other_cell);
		  for (unsigned int i=0; i<total_dofs; ++i)
		    {
		      for (unsigned int j=0; j<total_dofs; ++j)
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
}



template<int dim>
void
DoFTools::make_flux_sparsity_pattern (const DoFHandler<dim> &dof,
				      SparsityPattern       &sparsity,
				      const FullMatrix<double>& int_mask,
				      const FullMatrix<double>& flux_mask)
{
  const unsigned int n_dofs = dof.n_dofs();
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
  DoFHandler<dim>::active_cell_iterator cell = dof.begin_active(),
					endc = dof.end();


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
      cell->get_dof_indices (dofs_on_this_cell);
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
	  DoFHandler<dim>::face_iterator cell_face = cell->face(face);
	  if (cell_face->user_flag_set ())
	    continue;

	  if (! cell->at_boundary (face) )
	    {
	      DoFHandler<dim>::cell_iterator neighbor = cell->neighbor(face);
					       // Refinement edges are taken care of
					       // by coarser cells
	      if (neighbor->level() < cell->level())
		continue;

	      unsigned int neighbor_face = cell->neighbor_of_neighbor(face);

	      if (neighbor->has_children())
		{
		  for (unsigned int sub_nr = 0;
		       sub_nr != GeometryInfo<dim>::subfaces_per_face;
		       ++sub_nr)
		    {
		      DoFHandler<dim>::cell_iterator sub_neighbor
			= neighbor->
			child(GeometryInfo<dim>::child_cell_on_face(neighbor_face, sub_nr));

		      sub_neighbor->get_dof_indices (dofs_on_other_cell);
		      for (unsigned int i=0; i<total_dofs; ++i)
			{
			  for (unsigned int j=0; j<total_dofs; ++j)
			    {
			      if (flux_dof_mask[i][j])
				sparsity.add (dofs_on_this_cell[i],
					      dofs_on_other_cell[j]);
			      if (flux_dof_mask[j][i])
				sparsity.add (dofs_on_other_cell[i],
					      dofs_on_this_cell[j]);
			    }
			}
		      sub_neighbor->face(neighbor_face)->set_user_flag ();
		    }
		} else {
		  neighbor->get_dof_indices (dofs_on_other_cell);
		  for (unsigned int i=0; i<total_dofs; ++i)
		    {
		      for (unsigned int j=0; j<total_dofs; ++j)
			{
			  if (flux_dof_mask[i][j])
			    sparsity.add (dofs_on_this_cell[i],
					  dofs_on_other_cell[j]);
			  if (flux_dof_mask[j][i])
			    sparsity.add (dofs_on_other_cell[i],
					  dofs_on_this_cell[j]);
			}
		    }
		  neighbor->face(neighbor_face)->set_user_flag (); 
		}
	    } 
	}
    }
}



#if deal_II_dimension == 1

template <>
void DoFTools::make_hanging_node_constraints (const DoFHandler<1> &,
					      ConstraintMatrix &)
{
				   // nothing to be done here
};

#endif



#if deal_II_dimension == 2

template <>
void DoFTools::make_hanging_node_constraints (const DoFHandler<2> &dof_handler,
					      ConstraintMatrix    &constraints)
{
  const unsigned int dim = 2;
  
  const FiniteElement<dim> &fe   = dof_handler.get_fe();
  
				   // have space for the degrees of
				   // freedom on mother and child
				   // lines
  const unsigned int n_dofs_on_mother   = 2*fe.dofs_per_vertex + fe.dofs_per_line,
		     n_dofs_on_children = fe.dofs_per_vertex + 2*fe.dofs_per_line;

  vector<unsigned int> dofs_on_mother(n_dofs_on_mother);
  vector<unsigned int> dofs_on_children(n_dofs_on_children);

  Assert(n_dofs_on_mother == fe.constraints().n(),
	 ExcDimensionMismatch(n_dofs_on_mother,
			      fe.constraints().n()));
  Assert(n_dofs_on_children == fe.constraints().m(),
	 ExcDimensionMismatch(n_dofs_on_children,
			      fe.constraints().m()));

				   // loop over all lines; only on
				   // lines there can be constraints.
				   // We do so by looping over all
				   // active cells and checking
				   // whether any of the faces are
				   // refined which can only be from
				   // the neighboring cell because
				   // this one is active. In that
				   // case, the face is subject to
				   // constraints
				   //
				   // note that even though we may
				   // visit a face twice if the
				   // neighboring cells are equally
				   // refined, we can only visit each
				   // face with hanging nodes once
  DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
					endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
      if (cell->face(face)->has_children()) 
	{
	  const DoFHandler<dim>::line_iterator line = cell->face(face);
	  
					   // fill the dofs indices. Use same
					   // enumeration scheme as in
					   // @p{FiniteElement::constraints()}
	  unsigned int next_index = 0;
	  for (unsigned int vertex=0; vertex<2; ++vertex)
	    for (unsigned int dof=0; dof!=fe.dofs_per_vertex; ++dof)
	      dofs_on_mother[next_index++] = line->vertex_dof_index(vertex,dof);
	  for (unsigned int dof=0; dof!=fe.dofs_per_line; ++dof)
	    dofs_on_mother[next_index++] = line->dof_index(dof);
	  Assert (next_index == dofs_on_mother.size(),
		  ExcInternalError());
	  
	  next_index = 0;
	  for (unsigned int dof=0; dof!=fe.dofs_per_vertex; ++dof)
	    dofs_on_children[next_index++] = line->child(0)->vertex_dof_index(1,dof);
	  for (unsigned int child=0; child<2; ++child)
	    for (unsigned int dof=0; dof!=fe.dofs_per_line; ++dof)
	      dofs_on_children[next_index++] = line->child(child)->dof_index(dof);
	  Assert (next_index == dofs_on_children.size(),
		  ExcInternalError());
	  
					   // for each row in the constraint
					   // matrix for this line:
	  for (unsigned int row=0; row!=dofs_on_children.size(); ++row) 
	    {
	      constraints.add_line (dofs_on_children[row]);
	      for (unsigned int i=0; i!=dofs_on_mother.size(); ++i)
		constraints.add_entry (dofs_on_children[row],
				       dofs_on_mother[i],
				       fe.constraints()(row,i));
	    };
	};
};

#endif



#if deal_II_dimension == 3

template <>
void DoFTools::make_hanging_node_constraints (const DoFHandler<3> &dof_handler,
					      ConstraintMatrix    &constraints)
{
  const unsigned int dim = 3;
  
  const FiniteElement<dim> &fe   = dof_handler.get_fe();
  
				   // have space for the degrees of
				   // freedom on mother and child
				   // lines
  const unsigned int
    n_dofs_on_mother   = (4*fe.dofs_per_vertex+
			  4*fe.dofs_per_line+
			  fe.dofs_per_quad),
    n_dofs_on_children = (5*fe.dofs_per_vertex+
			  12*fe.dofs_per_line+
			  4*fe.dofs_per_quad);

  vector<unsigned int> dofs_on_mother(n_dofs_on_mother);
  vector<unsigned int> dofs_on_children(n_dofs_on_children);

  Assert(n_dofs_on_mother == fe.constraints().n(),
	 ExcDimensionMismatch(n_dofs_on_mother,
			      fe.constraints().n()));
  Assert(n_dofs_on_children == fe.constraints().m(),
	 ExcDimensionMismatch(n_dofs_on_children,
			      fe.constraints().m()));

				   // loop over all lines; only on
				   // lines there can be constraints.
				   // We do so by looping over all
				   // active cells and checking
				   // whether any of the faces are
				   // refined which can only be from
				   // the neighboring cell because
				   // this one is active. In that
				   // case, the face is subject to
				   // constraints
				   //
				   // note that even though we may
				   // visit a face twice if the
				   // neighboring cells are equally
				   // refined, we can only visit each
				   // face with hanging nodes once
  DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
					endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
      if (cell->face(f)->has_children()) 
	{
	  const DoFHandler<dim>::face_iterator face = cell->face(f);
	  
					   // fill the dofs indices. Use same
					   // enumeration scheme as in
					   // @p{FiniteElement::constraints()}
	  unsigned int next_index = 0;
	  for (unsigned int vertex=0; vertex<4; ++vertex)
	    for (unsigned int dof=0; dof!=fe.dofs_per_vertex; ++dof)
	      dofs_on_mother[next_index++] = face->vertex_dof_index(vertex,dof);
	  for (unsigned int line=0; line<4; ++line)
	    for (unsigned int dof=0; dof!=fe.dofs_per_line; ++dof)
	      dofs_on_mother[next_index++] = face->line(line)->dof_index(dof);
	  for (unsigned int dof=0; dof!=fe.dofs_per_quad; ++dof)
	    dofs_on_mother[next_index++] = face->dof_index(dof);
	  Assert (next_index == dofs_on_mother.size(),
		  ExcInternalError());
	  
	  next_index = 0;
	  Assert ((face->child(0)->vertex_dof_index(2,0) ==
		   face->child(1)->vertex_dof_index(3,0)) &&
		  (face->child(0)->vertex_dof_index(2,0) ==
		   face->child(2)->vertex_dof_index(0,0)) &&
		  (face->child(0)->vertex_dof_index(2,0) ==
		   face->child(3)->vertex_dof_index(1,0)),
		  ExcInternalError());
	  for (unsigned int dof=0; dof!=fe.dofs_per_vertex; ++dof)
	    dofs_on_children[next_index++]
	      = face->child(0)->vertex_dof_index(2,dof);
	  
					   // dof numbers on the centers of
					   // the lines bounding this face
	  for (unsigned int line=0; line<4; ++line)
	    for (unsigned int dof=0; dof!=fe.dofs_per_vertex; ++dof)
	      dofs_on_children[next_index++]
		= face->line(line)->child(0)->vertex_dof_index(1,dof);
	  
					   // next the dofs on the lines interior
					   // to the face; the order of these
					   // lines is laid down in the
					   // FiniteElement class documentation
	  for (unsigned int dof=0; dof<fe.dofs_per_line; ++dof)
	    dofs_on_children[next_index++]
	      = face->child(0)->line(1)->dof_index(dof);
	  for (unsigned int dof=0; dof<fe.dofs_per_line; ++dof)
	    dofs_on_children[next_index++]
	      = face->child(1)->line(2)->dof_index(dof);
	  for (unsigned int dof=0; dof<fe.dofs_per_line; ++dof)
	    dofs_on_children[next_index++]
	      = face->child(2)->line(3)->dof_index(dof);
	  for (unsigned int dof=0; dof<fe.dofs_per_line; ++dof)
	    dofs_on_children[next_index++]
	      = face->child(3)->line(0)->dof_index(dof);
	  
					   // dofs on the bordering lines
	  for (unsigned int line=0; line<4; ++line)
	    for (unsigned int child=0; child<2; ++child)
	      for (unsigned int dof=0; dof!=fe.dofs_per_line; ++dof)
		dofs_on_children[next_index++]
		  = face->line(line)->child(child)->dof_index(dof);
	  
					   // finally, for the dofs interior
					   // to the four child faces
	  for (unsigned int child=0; child<4; ++child)
	    for (unsigned int dof=0; dof!=fe.dofs_per_quad; ++dof)
	      dofs_on_children[next_index++]
		= face->child(child)->dof_index(dof);
	  Assert (next_index == dofs_on_children.size(),
		  ExcInternalError());
	  
					   // for each row in the constraint
					   // matrix for this line:
	  for (unsigned int row=0; row!=dofs_on_children.size(); ++row) 
	    {
	      constraints.add_line (dofs_on_children[row]);
	      for (unsigned int i=0; i!=dofs_on_mother.size(); ++i)
		constraints.add_entry (dofs_on_children[row],
				       dofs_on_mother[i],
				       fe.constraints()(row,i));
	    };
	};
};

#endif



template <int dim, typename Number>
void DoFTools::distribute_cell_to_dof_vector (const DoFHandler<dim> &dof_handler,
					      const Vector<Number>  &cell_data,
					      Vector<double>        &dof_data,
					      const unsigned int     component)
{
  const Triangulation<dim> &tria = dof_handler.get_tria();
  const FiniteElement<dim> &fe   = dof_handler.get_fe();
  
  Assert (cell_data.size()==tria.n_active_cells(),
	  ExcWrongSize (cell_data.size(), tria.n_active_cells()));
  Assert (dof_data.size()==dof_handler.n_dofs(),
	  ExcWrongSize (dof_data.size(), dof_handler.n_dofs()));
  Assert (component < fe.n_components(),
	  ExcInvalidComponent(component, fe.n_components()));

				   // store a flag whether we should care
				   // about different components. this is
				   // just a simplification, we could ask
				   // for this at every single place
				   // equally well
  const bool consider_components = (fe.n_components() != 1);
  
				   // count how often we have added a value
				   // in the sum for each dof
  vector<unsigned char> touch_count (dof_handler.n_dofs(), 0);

  DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
					endc = dof_handler.end();
  unsigned int present_cell = 0;
  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  vector<unsigned int> dof_indices (dofs_per_cell);

  for (; cell!=endc; ++cell, ++present_cell) 
    {
      cell->get_dof_indices (dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
					 // consider this dof only if it
					 // is the right component. if there
					 // is only one component, short cut
					 // the test
	if (!consider_components ||
	    (fe.system_to_component_index(i).first == component))
	  {
					     // sum up contribution of the
					     // present_cell to this dof
	    dof_data(dof_indices[i]) += cell_data(present_cell);
					     // note that we added another
					     // summand
	    ++touch_count[dof_indices[i]];
	  };
    };
  
				   // compute the mean value on all the
				   // dofs by dividing with the number
				   // of summands.
  for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
    {
				       // assert that each dof was used
				       // at least once. this needs not be
				       // the case if the vector has more than
				       // one component
      Assert (consider_components || (touch_count[i]!=0),
	      ExcInternalError());
      if (touch_count[i] != 0)
	dof_data(i) /=  touch_count[i];
    };
};


template<int dim>
void
DoFTools::extract_dofs (const DoFHandler<dim> &dof,
			const vector<bool>    &component_select,
			vector<bool>          &selected_dofs)
{
  const FiniteElement<dim> &fe = dof.get_fe();
  Assert(component_select.size() == fe.n_components(),
	 ExcDimensionMismatch(component_select.size(), fe.n_components()));
  Assert(selected_dofs.size() == dof.n_dofs(),
	 ExcDimensionMismatch(selected_dofs.size(), dof.n_dofs()));

				   // preset all values by false
  fill_n (selected_dofs.begin(), dof.n_dofs(), false);
  
  vector<unsigned int> indices(fe.dofs_per_cell);
  
  DoFHandler<dim>::active_cell_iterator c;
  for (c = dof.begin_active() ; c != dof.end() ; ++ c)
    {
      c->get_dof_indices(indices);
      for (unsigned int i=0;i<fe.dofs_per_cell;++i)
	{
	  const unsigned int component = fe.system_to_component_index(i).first;

	  if (component_select[component] == true)
	    selected_dofs[indices[i]] = true;
	}
    }
}


template<int dim>
void
DoFTools::extract_level_dofs(const unsigned int       level,
			     const MGDoFHandler<dim> &dof,
			     const vector<bool>      &component_select,
			     vector<bool>            &selected_dofs)
{
  const FiniteElement<dim>& fe = dof.get_fe();
  Assert(component_select.size() == fe.n_components(),
	 ExcDimensionMismatch(component_select.size(), fe.n_components()));
  Assert(selected_dofs.size() == dof.n_dofs(level),
	 ExcDimensionMismatch(selected_dofs.size(), dof.n_dofs(level)));

				   // preset all values by false
  fill_n (selected_dofs.begin(), dof.n_dofs(level), false);

  vector<unsigned int> indices(fe.dofs_per_cell);
  
  MGDoFHandler<dim>::cell_iterator c;
  for (c = dof.begin(level) ; c != dof.end(level) ; ++ c)
    {
      c->get_mg_dof_indices(indices);
      for (unsigned int i=0;i<fe.dofs_per_cell;++i)
	{
	  const unsigned int component = fe.system_to_component_index(i).first;
	  if (component_select[component]  == true)
	    selected_dofs[indices[i]] = true;
	}
    }
}



#if deal_II_dimension != 1

template <int dim>
void
DoFTools::extract_boundary_dofs (const DoFHandler<dim> &dof_handler,
				 const vector<bool>    &component_select,
				 vector<bool>          &selected_dofs,
				 const set<unsigned char> &boundary_indicators)
{
  Assert (component_select.size() == dof_handler.get_fe().n_components(),
	  ExcWrongSize (component_select.size(),
			dof_handler.get_fe().n_components()));
  Assert (boundary_indicators.find (255) == boundary_indicators.end(),
	  ExcInvalidBoundaryIndicator());

				   // let's see whether we have to
				   // check for certain boundary
				   // indicators or whether we can
				   // accept all
  const bool check_boundary_indicator = (boundary_indicators.size() != 0);
  
				   // clear and reset array by default
				   // values
  selected_dofs.clear ();
  selected_dofs.resize (dof_handler.n_dofs(), false);
  vector<unsigned int> face_dof_indices (dof_handler.get_fe().dofs_per_face);
  for (DoFHandler<dim>::active_cell_iterator cell=dof_handler.begin_active();
       cell!=dof_handler.end(); ++cell)
    for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
      if (cell->at_boundary(face))
	if (! check_boundary_indicator ||
	    (boundary_indicators.find (cell->face(face)->boundary_indicator())
	     != boundary_indicators.end()))
	{
	  cell->face(face)->get_dof_indices (face_dof_indices);
	  for (unsigned int i=0; i<dof_handler.get_fe().dofs_per_face; ++i)
	    if (component_select[dof_handler.get_fe().
				face_system_to_component_index(i).first] == true)
	      selected_dofs[face_dof_indices[i]] = true;
	};
};


#else

template <>
void
DoFTools::extract_boundary_dofs (const DoFHandler<1> &dof_handler,
				 const vector<bool>  &component_select,
				 vector<bool>        &selected_dofs)
{
  Assert (component_select.size() == dof_handler.get_fe().n_components(),
	  ExcWrongSize (component_select.size(),
			dof_handler.get_fe().n_components()));
	  
				   // clear and reset array by default
				   // values
  selected_dofs.clear ();
  selected_dofs.resize (dof_handler.n_dofs(), false);

  Assert (dof_handler.get_fe().dofs_per_face == dof_handler.get_fe().dofs_per_vertex,
	  ExcInternalError());
  
				   // loop over coarse grid cells
  for (DoFHandler<1>::cell_iterator cell=dof_handler.begin(0);
       cell!=dof_handler.end(0); ++cell)
    {
				       // check left-most vertex
      if (cell->neighbor(0) == dof_handler.end())
	for (unsigned int i=0; i<dof_handler.get_fe().dofs_per_face; ++i)
	  if (component_select[dof_handler.get_fe().
			      face_system_to_component_index(i).first] == true)
	    selected_dofs[cell->vertex_dof_index(0,i)] = true;
				       // check right-most vertex
      if (cell->neighbor(1) == dof_handler.end())
	for (unsigned int i=0; i<dof_handler.get_fe().dofs_per_face; ++i)
	  if (component_select[dof_handler.get_fe().
			      face_system_to_component_index(i).first] == true)
	    selected_dofs[cell->vertex_dof_index(1,i)] = true;
    };
};


#endif



#if deal_II_dimension == 1

template <>
void
DoFTools::extract_hanging_node_dofs (const DoFHandler<1> &dof_handler,
				     vector<bool>        &selected_dofs)
{
  Assert(selected_dofs.size() == dof_handler.n_dofs(),
	 ExcDimensionMismatch(selected_dofs.size(), dof_handler.n_dofs()));
				   // preset all values by false
  fill_n (selected_dofs.begin(), dof_handler.n_dofs(), false);

				   // there are no hanging nodes in 1d
};

#endif



#if deal_II_dimension == 2

template <>
void
DoFTools::extract_hanging_node_dofs (const DoFHandler<2> &dof_handler,
				     vector<bool>        &selected_dofs)
{
  const unsigned int dim = 2;
  
  Assert(selected_dofs.size() == dof_handler.n_dofs(),
	 ExcDimensionMismatch(selected_dofs.size(), dof_handler.n_dofs()));
				   // preset all values by false
  fill_n (selected_dofs.begin(), dof_handler.n_dofs(), false);

  const FiniteElement<dim> &fe   = dof_handler.get_fe();

				   // this function is similar to the
				   // make_sparsity_pattern function,
				   // see there for more information
  DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
					endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
      if (cell->face(face)->has_children()) 
	{
	  const DoFHandler<dim>::line_iterator line = cell->face(face);

	  for (unsigned int dof=0; dof!=fe.dofs_per_vertex; ++dof)
	    selected_dofs[line->child(0)->vertex_dof_index(1,dof)] = true;
	  
	  for (unsigned int child=0; child<2; ++child)
	    for (unsigned int dof=0; dof!=fe.dofs_per_line; ++dof)
	      selected_dofs[line->child(child)->dof_index(dof)] = true;
	};
};

#endif



#if deal_II_dimension == 3

template <>
void
DoFTools::extract_hanging_node_dofs (const DoFHandler<3> &dof_handler,
				     vector<bool>        &selected_dofs)
{
  const unsigned int dim = 3;

  Assert(selected_dofs.size() == dof_handler.n_dofs(),
	 ExcDimensionMismatch(selected_dofs.size(), dof_handler.n_dofs()));
				   // preset all values by false
  fill_n (selected_dofs.begin(), dof_handler.n_dofs(), false);

  const FiniteElement<dim> &fe   = dof_handler.get_fe();
  
				   // this function is similar to the
				   // make_sparsity_pattern function,
				   // see there for more information

  DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
					endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
      if (cell->face(f)->has_children()) 
	{
	  const DoFHandler<dim>::face_iterator face = cell->face(f);
	  
	  for (unsigned int dof=0; dof!=fe.dofs_per_vertex; ++dof)
	    selected_dofs[face->child(0)->vertex_dof_index(2,dof)] = true;
	  
					   // dof numbers on the centers of
					   // the lines bounding this face
	  for (unsigned int line=0; line<4; ++line)
	    for (unsigned int dof=0; dof!=fe.dofs_per_vertex; ++dof)
	      selected_dofs[face->line(line)->child(0)->vertex_dof_index(1,dof)] = true;
	  
					   // next the dofs on the lines interior
					   // to the face; the order of these
					   // lines is laid down in the
					   // FiniteElement class documentation
	  for (unsigned int dof=0; dof<fe.dofs_per_line; ++dof)
	    selected_dofs[face->child(0)->line(1)->dof_index(dof)] = true;
	  for (unsigned int dof=0; dof<fe.dofs_per_line; ++dof)
	    selected_dofs[face->child(1)->line(2)->dof_index(dof)] = true;
	  for (unsigned int dof=0; dof<fe.dofs_per_line; ++dof)
	    selected_dofs[face->child(2)->line(3)->dof_index(dof)] = true;
	  for (unsigned int dof=0; dof<fe.dofs_per_line; ++dof)
	    selected_dofs[face->child(3)->line(0)->dof_index(dof)] = true;
	  
					   // dofs on the bordering lines
	  for (unsigned int line=0; line<4; ++line)
	    for (unsigned int child=0; child<2; ++child)
	      for (unsigned int dof=0; dof!=fe.dofs_per_line; ++dof)
		selected_dofs[face->line(line)->child(child)->dof_index(dof)] = true;
	
					   // finally, for the dofs interior
					   // to the four child faces
	  for (unsigned int child=0; child<4; ++child)
	    for (unsigned int dof=0; dof!=fe.dofs_per_quad; ++dof)
	      selected_dofs[face->child(child)->dof_index(dof)] = true;
	};
};

#endif



template <int dim>
void
DoFTools::compute_intergrid_constraints (const DoFHandler<dim>              &coarse_grid,
					 const unsigned int                  coarse_component,
					 const DoFHandler<dim>              &fine_grid,
					 const unsigned int                  fine_component,
					 const InterGridMap<DoFHandler,dim> &coarse_to_fine_grid_map,
					 ConstraintMatrix                   &constraints)
{
				   // aliases to the finite elements
				   // used by the dof handlers:
  const FiniteElement<dim> &coarse_fe = coarse_grid.get_fe(),
			   &fine_fe   = fine_grid.get_fe();

				   // global numbers of dofs
  const unsigned int n_coarse_dofs = coarse_grid.n_dofs(),
		     n_fine_dofs   = fine_grid.n_dofs();

				   // local numbers of dofs
  const unsigned int fine_dofs_per_cell   = fine_fe.dofs_per_cell;

				   // alias the number of dofs per
				   // cell belonging to the
				   // coarse_component which is to be
				   // the restriction of the fine
				   // grid:
  const unsigned int coarse_dofs_per_cell_component
    = coarse_fe.base_element(coarse_fe.component_to_base(coarse_component)).dofs_per_cell;
  

				   // Try to find out whether the
				   // grids stem from the same coarse
				   // grid. This is a rather crude
				   // test, but better than nothing
  Assert (coarse_grid.get_tria().n_cells(0) == fine_grid.get_tria().n_cells(0),
	  ExcGridsDontMatch());

				   // check whether the map correlates
				   // the right objects
  Assert (&coarse_to_fine_grid_map.get_source_grid() == &coarse_grid,
	  ExcGridsDontMatch ());
  Assert (&coarse_to_fine_grid_map.get_destination_grid() == &fine_grid,
	  ExcGridsDontMatch ());
  
  
				   // check whether component numbers
				   // are valid
  Assert (coarse_component < coarse_fe.n_components(),
	  ExcInvalidComponent (coarse_component, coarse_fe.n_components()));
  Assert (fine_component < fine_fe.n_components(),
	  ExcInvalidComponent (fine_component, fine_fe.n_components()));
				   // check whether respective finite
				   // elements are equal
  Assert (coarse_fe.base_element (coarse_fe.component_to_base(coarse_component))
	  ==
	  fine_fe.base_element (fine_fe.component_to_base(fine_component)),
	  ExcFiniteElementsDontMatch());

#ifdef DEBUG
				   // if in debug mode, check whether
				   // the coarse grid is indeed
				   // coarser everywhere than the fine
				   // grid
  for (typename DoFHandler<dim>::active_cell_iterator cell=coarse_grid.begin_active();
       cell != coarse_grid.end(); ++cell)
    Assert (cell->level() <= coarse_to_fine_grid_map[cell]->level(),
	    ExcGridNotCoarser());
#endif

  

/*
 * From here on: the term `parameter' refers to the selected component
 * on the coarse grid and its analogon on the fine grid. The naming of
 * variables containing this term is due to the fact that
 * `selected_component' is longer, but also due to the fact that the
 * code of this function was initially written for a program where the
 * component which we wanted to match between grids was actually the
 * `parameter' variable.
 *
 * Likewise, the terms `parameter grid' and `state grid' refer to the
 * coarse and fine grids, respectively.
 *
 * Changing the names of variables would in principle be a good idea,
 * but would not make things simpler and would be another source of
 * errors. If anyone feels like doing so: patches would be welcome!
 */


  
				   // set up vectors of cell-local
				   // data; each vector represents one
				   // degree of freedom of the
				   // coarse-grid variable in the
				   // fine-grid element
  vector<Vector<double> > parameter_dofs (coarse_dofs_per_cell_component,
					  Vector<double>(fine_dofs_per_cell));
				   // for each coarse dof: find its
				   // position within the fine element
				   // and set this value to one in the
				   // respective vector (all other values
				   // are zero by construction)
  for (unsigned int local_coarse_dof=0;
       local_coarse_dof<coarse_dofs_per_cell_component;
       ++local_coarse_dof)
    parameter_dofs[local_coarse_dof]
      (fine_fe.component_to_system_index (fine_component,local_coarse_dof)) = 1.;


				   // find out how many DoFs there are
				   // on the grids belonging to the
				   // components we want to match
  unsigned int n_parameters_on_fine_grid=0;
  if (true)
    {
				       // have a flag for each dof on
				       // the fine grid and set it
				       // to true if this is an
				       // interesting dof. finally count
				       // how many true's there
      vector<bool> dof_is_interesting (fine_grid.n_dofs(), false);
      vector<unsigned int>  local_dof_indices (fine_fe.dofs_per_cell);
      
      for (DoFHandler<dim>::active_cell_iterator cell=fine_grid.begin_active();
	   cell!=fine_grid.end(); ++cell)
	{
	  cell->get_dof_indices (local_dof_indices);
	  for (unsigned int i=0; i<fine_fe.dofs_per_cell; ++i)
	    if (fine_fe.system_to_component_index(i).first == fine_component)
	      dof_is_interesting[local_dof_indices[i]] = true;
	};

      n_parameters_on_fine_grid = count (dof_is_interesting.begin(),
					 dof_is_interesting.end(),
					 true);
    };

				   // get an array in which we store
				   // which dof on the coarse grid is
				   // a parameter and which is not
  vector<bool> coarse_dof_is_parameter (coarse_grid.n_dofs());
  if (true)
    {
      vector<bool> mask (coarse_grid.get_fe().n_components(),
			 false);
      mask[coarse_component] = true;
      extract_dofs (coarse_grid, mask, coarse_dof_is_parameter);
    };
  
  

				   // store the weights with which a dof
				   // on the parameter grid contributes
				   // to a dof on the fine grid. see the
				   // long doc below for more info
				   //
				   // allocate as many rows as there are
				   // parameter dofs on the coarse grid
				   // and as many columns as there are
				   // parameter dofs on the fine grid.
				   //
				   // weight_mapping is used to map the
				   // global (fine grid) parameter dof
				   // indices to the columns
				   //
				   // note that the `weights' array
				   // can take up huge amounts of
				   // memory, and in particular is
				   // roughly quadratic in the memory
				   // consumption!
  FullMatrix<float> weights (n_coarse_dofs,
			     n_parameters_on_fine_grid);
				   // this is this mapping. there is one
				   // entry for each dof on the fine grid;
				   // if it is a parameter dof, then its
				   // value is the column in weights for
				   // that parameter dof, if it is any
				   // other dof, then its value is -1,
				   // indicating an error
  vector<int> weight_mapping (n_fine_dofs, -1);

				   // set up this mapping
  if (true)
    {
      vector<unsigned int> local_dof_indices(fine_fe.dofs_per_cell);
      unsigned int next_free_index=0;
      for (typename DoFHandler<dim>::active_cell_iterator cell=fine_grid.begin_active();
	   cell != fine_grid.end(); ++cell)
	{
	  cell->get_dof_indices (local_dof_indices);
	  for (unsigned int i=0; i<fine_fe.dofs_per_cell; ++i)
					     // if this DoF is a
					     // parameter dof and has
					     // not yet been numbered,
					     // then do so
	    if ((fine_fe.system_to_component_index(i).first == fine_component) &&
		(weight_mapping[local_dof_indices[i]] == -1))
	      {
		weight_mapping[local_dof_indices[i]] = next_free_index;
		++next_free_index;
	      };
	};

      Assert (next_free_index == n_parameters_on_fine_grid,
	      ExcInternalError());
    };

  
				   // for each cell on the parameter grid:
				   // find out which degrees of freedom on the
				   // fine grid correspond in which way to
				   // the degrees of freedom on the parameter
				   // grid
				   //
				   // do this in a separate function
				   // to allow for multithreading
				   // there. see this function also if
				   // you want to read more
				   // information on the algorithm
				   // used.
  compute_intergrid_weights (coarse_grid, coarse_component,
			     coarse_to_fine_grid_map, parameter_dofs,
			     weight_mapping, weights);

				   // ok, now we have all weights for each
				   // dof on the fine grid. if in debug
				   // mode lets see if everything went smooth,
				   // i.e. each dof has sum of weights one
				   //
				   // in other words this means that
				   // if the sum of all shape
				   // functions on the parameter grid
				   // is one (which is always the
				   // case), then the representation
				   // on the state grid should be as
				   // well (division of unity)
				   //
				   // if the parameter grid has more
				   // than one component, then the
				   // respective dofs of the other
				   // components have sum of weights
				   // zero, of course. we do not
				   // explicitely ask which component
				   // a dof belongs to, but this at
				   // least tests some errors
#ifdef DEBUG
  for (unsigned int col=0; col<weights.n(); ++col)
    {
      double sum=0;
      for (unsigned int row=0; row<weights.m(); ++row)
	sum += weights(row,col);
      Assert ((sum==1) ||
	      ((coarse_fe.n_components()>1) && (sum==0)), ExcInternalError());
    };
#endif


				   // now we know that the weights in
				   // each row constitute a
				   // constraint. enter this into the
				   // constraints object
				   //
				   // first task: for each parameter
				   // dof on the parameter grid, find
				   // a representant on the fine,
				   // global grid. this is possible
				   // since we use conforming finite
				   // element. we take this
				   // representant to be the first
				   // element in this row with weight
				   // identical to one. the
				   // representant will become an
				   // unconstrained degree of freedom,
				   // while all others will be
				   // constrained to this dof (and
				   // possibly others)
  vector<int> representants(weights.m(), -1);
  for (unsigned int parameter_dof=0; parameter_dof<weights.m(); ++parameter_dof)
    if (coarse_dof_is_parameter[parameter_dof] == true)
      {
	unsigned int column=0;
	for (; column<weights.n(); ++column)
	  if (weights(parameter_dof,column) == 1)
	    break;
	Assert (column < weights.n(), ExcInternalError());
	
					 // now we know in which column of
					 // weights the representant is, but
					 // we don't know its global index. get
					 // it using the inverse operation of
					 // the weight_mapping
	unsigned int global_dof=0;
	for (; global_dof<weight_mapping.size(); ++global_dof)
	  if (weight_mapping[global_dof] == static_cast<int>(column))
	    break;
	Assert (global_dof < weight_mapping.size(), ExcInternalError());
	
					 // now enter the representants global
					 // index into our list
	representants[parameter_dof] = global_dof;
      }
    else
      {
					 // consistency check: if this
					 // is no parameter dof on the
					 // coarse grid, then the
					 // respective row must be
					 // empty!
	for (unsigned int col=0; col<weights.n(); ++col)
	  Assert (weights(parameter_dof,col) == 0, ExcInternalError());
      };
  


				   // note for people that want to
				   // optimize this function: the
				   // largest part of the computing
				   // time is spent in the following,
				   // rather innocent block of
				   // code. basically, it must be the
				   // ConstraintMatrix::add_entry call
				   // which takes the bulk of the
				   // time, but it is not known to the
				   // author how to make it faster...
  vector<pair<unsigned int,double> > constraint_line;
  for (unsigned int global_dof=0; global_dof<n_fine_dofs; ++global_dof)
    if (weight_mapping[global_dof] != -1)
				       // this global dof is a parameter
				       // dof, so it may carry a constraint
				       // note that for each global dof,
				       // the sum of weights shall be one,
				       // so we can find out whether this
				       // dof is constrained in the following
				       // way: if the only weight in this row
				       // is a one, and the representant for
				       // the parameter dof of the line in
				       // which this one is is the present
				       // dof, then we consider this dof
				       // to be unconstrained. otherwise,
				       // all other dofs are constrained
      {
	unsigned int first_used_row=0;
	for (; first_used_row<weights.m(); ++first_used_row)
	  if (weights(first_used_row,weight_mapping[global_dof]) != 0)
	    break;

	if ((weights(first_used_row,weight_mapping[global_dof]) == 1) &&
	    (representants[first_used_row] == static_cast<int>(global_dof)))
					   // dof unconstrained or
					   // constrained to itself
					   // (in case this cell is
					   // mapped to itself, rather
					   // than to children of
					   // itself)
	  continue;

					 // otherwise enter all constraints
	constraints.add_line (global_dof);

	constraint_line.clear ();
	for (unsigned int row=first_used_row; row<weights.m(); ++row)
	  if (weights(row,weight_mapping[global_dof]) != 0)
	    constraint_line.push_back (make_pair(representants[row],
						 weights(row,weight_mapping[global_dof])));

	constraints.add_entries (global_dof, constraint_line);
      };
};



template <int dim>
void
DoFTools::compute_intergrid_weights (const DoFHandler<dim>              &coarse_grid,
				     const unsigned int                  coarse_component,
				     const InterGridMap<DoFHandler,dim> &coarse_to_fine_grid_map,
				     const vector<Vector<double> >      &parameter_dofs,
				     const vector<int>                  &weight_mapping,
				     FullMatrix<float>                  &weights)
{
				   // simply distribute the range of
				   // cells to different threads
  typedef typename DoFHandler<dim>::active_cell_iterator active_cell_iterator;
  vector<pair<active_cell_iterator,active_cell_iterator> >
    cell_intervals = Threads::split_range<active_cell_iterator> (coarse_grid.begin_active(),
								 coarse_grid.end(),
								 multithread_info.n_default_threads);

  Threads::ThreadManager thread_manager;
  for (unsigned int i=0; i<multithread_info.n_default_threads; ++i)
    Threads::spawn (thread_manager,
		    Threads::encapsulate (&DoFTools::template compute_intergrid_weights_1<dim>)
		    .collect_args (coarse_grid, coarse_component,
				   coarse_to_fine_grid_map, parameter_dofs,
				   weight_mapping, weights,
				   cell_intervals[i].first,
				   cell_intervals[i].second));

				   // wait for the threads to finish
  thread_manager.wait ();
};



template <int dim>
void
DoFTools::compute_intergrid_weights_1 (const DoFHandler<dim>              &coarse_grid,
				       const unsigned int                  coarse_component,
				       const InterGridMap<DoFHandler,dim> &coarse_to_fine_grid_map,
				       const vector<Vector<double> >      &parameter_dofs,
				       const vector<int>                  &weight_mapping,
				       FullMatrix<float>                  &weights,
				       const typename DoFHandler<dim>::active_cell_iterator &begin,
				       const typename DoFHandler<dim>::active_cell_iterator &end)
{
				   // aliases to the finite elements
				   // used by the dof handlers:
  const FiniteElement<dim> &coarse_fe = coarse_grid.get_fe();    

  const unsigned int coarse_dofs_per_cell_component
    = coarse_fe.base_element(coarse_fe.component_to_base(coarse_component)).dofs_per_cell;
  
				   // for each cell on the parameter grid:
				   // find out which degrees of freedom on the
				   // fine grid correspond in which way to
				   // the degrees of freedom on the parameter
				   // grid
				   //
				   // since for continuous FEs some
				   // dofs exist on more than one
				   // cell, we have to track which
				   // ones were already visited. the
				   // problem is that if we visit a
				   // dof first on one cell and
				   // compute its weight with respect
				   // to some global dofs to be
				   // non-zero, and later visit the
				   // dof again on another cell and
				   // (since we are on another cell)
				   // recompute the weights with
				   // respect to the same dofs as
				   // above to be zero now, we have to
				   // preserve them. we therefore
				   // overwrite all weights if they
				   // are nonzero and do not enforce
				   // zero weights since that might be
				   // only due to the fact that we are
				   // on another cell.
				   //
				   // example:
				   // coarse grid
				   //  |     |     |
				   //  *-----*-----*
				   //  | cell|cell |
				   //  |  1  |  2  |
				   //  |     |     |
				   //  0-----1-----*
				   //
				   // fine grid
				   //  |  |  |  |  |
				   //  *--*--*--*--*
				   //  |  |  |  |  |
				   //  *--*--*--*--*
				   //  |  |  |  |  |
				   //  *--x--y--*--*
				   //
				   // when on cell 1, we compute the
				   // weights of dof 'x' to be 1/2
				   // from parameter dofs 0 and 1,
				   // respectively. however, when
				   // later we are on cell 2, we again
				   // compute the prolongation of
				   // shape function 1 restricted to
				   // cell 2 to the globla grid and
				   // find that the weight of global
				   // dof 'x' now is zero. however, we
				   // should not overwrite the old
				   // value.
				   //
				   // we therefore always only set
				   // nonzero values. why adding up is
				   // not useful: dof 'y' would get
				   // weight 1 from parameter dof 1 on
				   // both cells 1 and 2, but the
				   // correct weight is nevertheless
				   // only 1.

				   // vector to hold the representation of
				   // a single degree of freedom on the
				   // coarse grid (for the selected fe)
				   // on the fine grid
  const unsigned int n_fine_dofs = weight_mapping.size();
  Vector<double> global_parameter_representation (n_fine_dofs);
  
  typename DoFHandler<dim>::active_cell_iterator cell;
  vector<unsigned int> parameter_dof_indices (coarse_fe.dofs_per_cell);
  
  for (cell=begin; cell!=end; ++cell)
    {
				       // get the global indices of the
				       // parameter dofs on this parameter
				       // grid cell
      cell->get_dof_indices (parameter_dof_indices);

				       // loop over all dofs on this
				       // cell and check whether they
				       // are interesting for us
      for (unsigned int local_dof=0;
	   local_dof<coarse_fe.dofs_per_cell;
	   ++local_dof)
	if (coarse_fe.system_to_component_index(local_dof).first
	    ==
	    coarse_component)
	  {
					     // the how-many-th
					     // parameter is this on
					     // this cell?
	    const unsigned int local_parameter_dof
	      = coarse_fe.system_to_component_index(local_dof).second;
	    
	    global_parameter_representation.clear ();
	    
					     // distribute the representation of
					     // @p{local_parameter_dof} on the
					     // parameter grid cell @p{cell} to
					     // the global data space
	    coarse_to_fine_grid_map[cell]->
	      set_dof_values_by_interpolation (parameter_dofs[local_parameter_dof],
					       global_parameter_representation);
					     // now that we've got the global
					     // representation of each parameter
					     // dof, we've only got to clobber the
					     // non-zero entries in that vector and
					     // store the result
					     //
					     // what we have learned: if entry @p{i}
					     // of the global vector holds the value
					     // @p{v[i]}, then this is the weight with
					     // which the present dof contributes
					     // to @p{i}. there may be several such
					     // @p{i}s and their weights' sum should
					     // be one. Then, @p{v[i]} should be
					     // equal to @p{\sum_j w_{ij} p[j]} with
					     // @p{p[j]} be the values of the degrees
					     // of freedom on the coarse grid. we
					     // can thus compute constraints which
					     // link the degrees of freedom @p{v[i]}
					     // on the fine grid to those on the
					     // coarse grid, @p{p[j]}. Now to use
					     // these as real constraints, rather
					     // than as additional equations, we
					     // have to identify representants
					     // among the @p{i} for each @p{j}. this will
					     // be done by simply taking the first
					     // @p{i} for which @p{w_{ij}==1}.
	    for (unsigned int i=0; i<global_parameter_representation.size(); ++i)
					       // set this weight if it belongs
					       // to a parameter dof.
	      if (weight_mapping[i] != -1)
		{
						   // only overwrite old
						   // value if not by
						   // zero
		  if (global_parameter_representation(i) != 0)
		    {
		      const unsigned int wi = parameter_dof_indices[local_dof],
					 wj = weight_mapping[i];
		      weights(wi,wj) = global_parameter_representation(i);
		    };
		}
	      else
		Assert (global_parameter_representation(i) == 0,
			ExcInternalError());
	  };
    };
};




#if deal_II_dimension == 1

template <>
void DoFTools::map_dof_to_boundary_indices (const DoFHandler<1> &dof_handler,
					    vector<unsigned int> &)
{
  Assert (&dof_handler.get_fe() != 0, ExcNoFESelected());
  Assert (false, ExcNotImplemented());
};



template <>
void DoFTools::map_dof_to_boundary_indices (const DoFHandler<1> &dof_handler,
					    const set<unsigned char> &,
					    vector<unsigned int> &)
{
  Assert (&dof_handler.get_fe() != 0, ExcNoFESelected());
  Assert (false, ExcNotImplemented());
};


#else


template <int dim>
void DoFTools::map_dof_to_boundary_indices (const DoFHandler<dim> &dof_handler,
					    vector<unsigned int>  &mapping)
{
  Assert (&dof_handler.get_fe() != 0, ExcNoFESelected());

  mapping.clear ();
  mapping.insert (mapping.end(), dof_handler.n_dofs(),
		  DoFHandler<dim>::invalid_dof_index);
  
  const unsigned int dofs_per_face = dof_handler.get_fe().dofs_per_face;
  vector<unsigned int> dofs_on_face(dofs_per_face);
  unsigned int next_boundary_index = 0;
  
  typename DoFHandler<dim>::active_face_iterator face = dof_handler.begin_active_face(),
						 endf = dof_handler.end_face();
  for (; face!=endf; ++face)
    if (face->at_boundary()) 
      {
	face->get_dof_indices (dofs_on_face);
	for (unsigned int i=0; i<dofs_per_face; ++i)
	  if (mapping[dofs_on_face[i]] == DoFHandler<dim>::invalid_dof_index)
	    mapping[dofs_on_face[i]] = next_boundary_index++;
      };

  Assert (next_boundary_index == dof_handler.n_boundary_dofs(),
	  ExcInternalError());
};



template <int dim>
void DoFTools::map_dof_to_boundary_indices (const DoFHandler<dim>    &dof_handler,
					    const set<unsigned char> &boundary_indicators,
					    vector<unsigned int>     &mapping)
{
  Assert (&dof_handler.get_fe() != 0, ExcNoFESelected());
  Assert (boundary_indicators.find (255) == boundary_indicators.end(),
	  ExcInvalidBoundaryIndicator());

  mapping.clear ();
  mapping.insert (mapping.end(), dof_handler.n_dofs(),
		  DoFHandler<dim>::invalid_dof_index);

				   // return if there is nothing to do
  if (boundary_indicators.size() == 0)
    return;
  
  const unsigned int dofs_per_face = dof_handler.get_fe().dofs_per_face;
  vector<unsigned int> dofs_on_face(dofs_per_face);
  unsigned int next_boundary_index = 0;
  
  typename DoFHandler<dim>::active_face_iterator face = dof_handler.begin_active_face(),
						 endf = dof_handler.end_face();
  for (; face!=endf; ++face)
    if (boundary_indicators.find (face->boundary_indicator()) !=
	boundary_indicators.end())
      {
	face->get_dof_indices (dofs_on_face);
	for (unsigned int i=0; i<dofs_per_face; ++i)
	  if (mapping[dofs_on_face[i]] == DoFHandler<dim>::invalid_dof_index)
	    mapping[dofs_on_face[i]] = next_boundary_index++;
      };

  Assert (next_boundary_index == dof_handler.n_boundary_dofs(boundary_indicators),
	  ExcInternalError());
};

#endif



// explicit instantiations
#if deal_II_dimension > 1
template void
DoFTools::make_flux_sparsity_pattern (const DoFHandler<deal_II_dimension>& dof,
				      SparsityPattern    &sparsity);
template void
DoFTools::make_flux_sparsity_pattern (const DoFHandler<deal_II_dimension>& dof,
				      SparsityPattern    &,
				      const FullMatrix<double>&,
				      const FullMatrix<double>&);
template void
DoFTools::make_boundary_sparsity_pattern (const DoFHandler<deal_II_dimension>& dof,
					  const vector<unsigned int>  &,
					  SparsityPattern    &);
template void
DoFTools::make_boundary_sparsity_pattern (const DoFHandler<deal_II_dimension>& dof,
					  const DoFHandler<deal_II_dimension>::FunctionMap  &boundary_indicators,
					  const vector<unsigned int>  &dof_to_boundary_mapping,
					  SparsityPattern    &sparsity);
#endif

template void
DoFTools::make_sparsity_pattern (const DoFHandler<deal_II_dimension>& dof,
				 SparsityPattern    &sparsity);

template void
DoFTools::make_sparsity_pattern (const DoFHandler<deal_II_dimension> &dof,
				 BlockSparsityPattern                &sparsity);

template void 
DoFTools::make_sparsity_pattern (const DoFHandler<deal_II_dimension>& dof,
				 const vector<vector<bool> > &mask,
				 SparsityPattern             &sparsity);

template void 
DoFTools::make_sparsity_pattern (const DoFHandler<deal_II_dimension>& dof,
				 const vector<vector<bool> > &mask,
				 BlockSparsityPattern        &sparsity);

template
void
DoFTools::distribute_cell_to_dof_vector (const DoFHandler<deal_II_dimension> &dof_handler,
					 const Vector<float>  &cell_data,
					 Vector<double>       &dof_data);

template
void
DoFTools::distribute_cell_to_dof_vector (const DoFHandler<deal_II_dimension> &dof_handler,
					 const Vector<double> &cell_data,
					 Vector<double>       &dof_data);


template void DoFTools::extract_dofs(const DoFHandler<deal_II_dimension>& dof,
				     const vector<bool>& component_select,
				     vector<bool>& selected_dofs);
template void DoFTools::extract_level_dofs(unsigned int level,
					   const MGDoFHandler<deal_II_dimension>& dof,
					   const vector<bool>& component_select,
					   vector<bool>& selected_dofs);


#if deal_II_dimension != 1
template
void
DoFTools::extract_boundary_dofs (const DoFHandler<deal_II_dimension> &,
				 const vector<bool>                  &,
				 vector<bool>                        &);
#endif 

template
void
DoFTools::compute_intergrid_constraints (const DoFHandler<deal_II_dimension> &,
					 const unsigned int                   ,
					 const DoFHandler<deal_II_dimension> &,
					 const unsigned int                   ,
					 const InterGridMap<DoFHandler,deal_II_dimension> &,
					 ConstraintMatrix                    &);



#if deal_II_dimension != 1

template
void
DoFTools::map_dof_to_boundary_indices (const DoFHandler<deal_II_dimension> &,
				       vector<unsigned int> &);

template
void
DoFTools::map_dof_to_boundary_indices (const DoFHandler<deal_II_dimension> &,
				       const set<unsigned char> &,
				       vector<unsigned int> &);

#endif 
