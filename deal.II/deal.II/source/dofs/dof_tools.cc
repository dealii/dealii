// $Id$


#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_handler.h>
#include <dofs/mg_dof_handler.h>
#include <dofs/dof_accessor.h>
#include <dofs/mg_dof_accessor.h>
#include <dofs/dof_constraints.h>
#include <fe/fe.h>
#include <fe/fe_system.h>
#include <dofs/dof_tools.h>
#include <lac/sparsity_pattern.h>
#include <lac/vector.h>

#include <algorithm>


template <int dim>
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
  vector<int> dofs_on_this_cell(dofs_per_cell);
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


template <int dim>
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
  
  
  vector<int> dofs_on_this_cell(dofs_per_cell);
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
					       const vector<int>  &,
					       SparsityPattern    &)
{
    Assert (false, ExcInternalError());
};



template <>
void
DoFTools::make_boundary_sparsity_pattern (const DoFHandler<1>&,
					  const DoFHandler<1>::FunctionMap  &,
					  const vector<int>  &,
					  SparsityPattern    &)
{
  Assert (false, ExcInternalError());
}

#endif



template <int dim>
void
DoFTools::make_boundary_sparsity_pattern (const DoFHandler<dim>& dof,
					  const vector<int>  &dof_to_boundary_mapping,
					  SparsityPattern    &sparsity)
{
  const unsigned int n_dofs = dof.n_dofs();

  Assert (dof_to_boundary_mapping.size() == n_dofs, ExcInternalError());
  Assert (sparsity.n_rows() == dof.n_boundary_dofs(),
	  ExcDimensionMismatch (sparsity.n_rows(), dof.n_boundary_dofs()));
  Assert (sparsity.n_cols() == dof.n_boundary_dofs(),
	  ExcDimensionMismatch (sparsity.n_cols(), dof.n_boundary_dofs()));
  Assert (*max_element(dof_to_boundary_mapping.begin(),
		       dof_to_boundary_mapping.end()) == (signed int)sparsity.n_rows()-1,
	  ExcInternalError());

  const unsigned int dofs_per_face = dof.get_fe().dofs_per_face;
  vector<int> dofs_on_this_face(dofs_per_face);
  DoFHandler<dim>::active_face_iterator face = dof.begin_active_face(),
		       endf = dof.end_face();
  for (; face!=endf; ++face)
    if (face->at_boundary())
      {
	face->get_dof_indices (dofs_on_this_face);

					 // make sure all dof indices have a
					 // boundary index
	Assert (*min_element(dofs_on_this_face.begin(),
			     dofs_on_this_face.end()) >=0,
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
					       const vector<int>  &dof_to_boundary_mapping,
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
  Assert (*max_element(dof_to_boundary_mapping.begin(),
		       dof_to_boundary_mapping.end()) == (signed int)sparsity.n_rows()-1,
	  ExcInternalError());

  const unsigned int dofs_per_face = dof.get_fe().dofs_per_face;
  vector<int> dofs_on_this_face(dofs_per_face);
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
			     dofs_on_this_face.end()) >=0,
		ExcInternalError());
					 // make sparsity pattern for this cell
	for (unsigned int i=0; i<dofs_per_face; ++i)
	  for (unsigned int j=0; j<dofs_per_face; ++j)
	    sparsity.add (dof_to_boundary_mapping[dofs_on_this_face[i]],
			  dof_to_boundary_mapping[dofs_on_this_face[j]]);
      };
};




//TODO: Check this function for potential of optimization. (G)

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
	  if (! cell->at_boundary(face) )
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
  const Triangulation<dim> &tria = dof_handler.get_tria();
  const FiniteElement<dim> &fe   = dof_handler.get_fe();
  
				   // first mark all faces which are subject
				   // to constraints. We do so by looping
				   // over all active cells and checking
				   // whether any of the faces are refined
				   // which can only be from the neighboring
				   // cell because this one is active. In that
				   // case, the face is subject to constraints
  DoFHandler<dim>::line_iterator line = dof_handler.begin_line(),
				 endl = dof_handler.end_line();
  for (; line!=endl; ++line)
    line->clear_user_flag ();
  
  Triangulation<dim>::active_cell_iterator cell = tria.begin_active(),
					   endc = tria.end();
  for (; cell!=endc; ++cell)
    for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
      if (cell->face(face)->has_children()) 
	cell->face(face)->set_user_flag();
	  
  

  
				   // loop over all lines; only on lines
				   // there can be constraints.
  for (line = dof_handler.begin_line(); line != endl; ++line)
				     // if dofs on this line are subject
				     // to constraints
    if (line->user_flag_set() == true)
      {
					 // reserve space to gather
					 // the dof numbers. We could
					 // get them when we need them,
					 // but it seems easier to gather
					 // them only once.
	vector<int> dofs_on_mother;
	vector<int> dofs_on_children;
	dofs_on_mother.reserve (2*fe.dofs_per_vertex+
				fe.dofs_per_line);
	dofs_on_children.reserve (fe.dofs_per_vertex+
				  2*fe.dofs_per_line);

	Assert(2*fe.dofs_per_vertex+fe.dofs_per_line ==
	       fe.constraints().n(),
	       ExcDimensionMismatch(2*fe.dofs_per_vertex+
				      fe.dofs_per_line,
				      fe.constraints().n()));
	Assert(fe.dofs_per_vertex+2*fe.dofs_per_line ==
	       fe.constraints().m(),
	       ExcDimensionMismatch(3*fe.dofs_per_vertex+
				      2*fe.dofs_per_line,
				      fe.constraints().m()));
	
					 // fill the dofs indices. Use same
					 // enumeration scheme as in
					 // #FiniteElement::constraints()#
	for (unsigned int vertex=0; vertex<2; ++vertex)
	  for (unsigned int dof=0; dof!=fe.dofs_per_vertex; ++dof)
	    dofs_on_mother.push_back (line->vertex_dof_index(vertex,dof));
	for (unsigned int dof=0; dof!=fe.dofs_per_line; ++dof)
	  dofs_on_mother.push_back (line->dof_index(dof));

	for (unsigned int dof=0; dof!=fe.dofs_per_vertex; ++dof)
	  dofs_on_children.push_back (line->child(0)->vertex_dof_index(1,dof));
	for (unsigned int child=0; child<2; ++child)
	  for (unsigned int dof=0; dof!=fe.dofs_per_line; ++dof)
	    dofs_on_children.push_back (line->child(child)->dof_index(dof));

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
  const Triangulation<dim> &tria = dof_handler.get_tria();
  const FiniteElement<dim> &fe   = dof_handler.get_fe();
  
				   // first mark all faces which are subject
				   // to constraints. We do so by looping
				   // over all active cells and checking
				   // whether any of the faces are refined
				   // which can only be from the neighboring
				   // cell because this one is active. In that
				   // case, the face is subject to constraints
  DoFHandler<dim>::face_iterator face = dof_handler.begin_face(),
				 endf = dof_handler.end_face();
  for (; face!=endf; ++face)
    face->clear_user_flag ();

  Triangulation<dim>::active_cell_iterator cell = tria.begin_active(),
					   endc = tria.end();
  for (; cell!=endc; ++cell)
    for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
      if (cell->face(face)->has_children()) 
	cell->face(face)->set_user_flag();
	  
  
				   // loop over all faces; only on faces
				   // there can be constraints.
  for (face=dof_handler.begin_face(); face != endf; ++face)
				     // if dofs on this line are subject
				     // to constraints
    if (face->user_flag_set() == true)
      {
					 // reserve space to gather
					 // the dof numbers. We could
					 // get them when we need them,
					 // but it seems easier to gather
					 // them only once.
	vector<int> dofs_on_mother;
	vector<int> dofs_on_children;
	dofs_on_mother.reserve (4*fe.dofs_per_vertex+
				4*fe.dofs_per_line+
				fe.dofs_per_quad);
	dofs_on_children.reserve (5*fe.dofs_per_vertex+
				  12*fe.dofs_per_line+
				  4*fe.dofs_per_quad);

	Assert(4*fe.dofs_per_vertex+
	       4*fe.dofs_per_line+
	       fe.dofs_per_quad
	       ==
	       fe.constraints().n(),
	       ExcDimensionMismatch(4*fe.dofs_per_vertex+
				      4*fe.dofs_per_line+
				      fe.dofs_per_quad,
				      fe.constraints().n()));
	Assert(5*fe.dofs_per_vertex+
	       12*fe.dofs_per_line+
	       4*fe.dofs_per_quad
	       ==
	       fe.constraints().m(),
	       ExcDimensionMismatch(5*fe.dofs_per_vertex+
				      12*fe.dofs_per_line+
				      4*fe.dofs_per_quad,
				      fe.constraints().m()));
	
					 // fill the dofs indices. Use same
					 // enumeration scheme as in
					 // #FiniteElement::constraints()#
	for (unsigned int vertex=0; vertex<4; ++vertex)
	  for (unsigned int dof=0; dof!=fe.dofs_per_vertex; ++dof)
	    dofs_on_mother.push_back (face->vertex_dof_index(vertex,dof));
	for (unsigned int line=0; line<4; ++line)
	  for (unsigned int dof=0; dof!=fe.dofs_per_line; ++dof)
	    dofs_on_mother.push_back (face->line(line)->dof_index(dof));
	for (unsigned int dof=0; dof!=fe.dofs_per_quad; ++dof)
	  dofs_on_mother.push_back (face->dof_index(dof));

					 // dof numbers on vertex at the center
					 // of the face, which is vertex 2 of
					 // child zero, or vertex 3 of child 1
					 // or vertex 0 of child 2 or vertex 1
					 // of child 3. We're a bit cautious and
					 // check this (also an additional safety
					 // check for the internal states of the
					 // library)
	Assert ((face->child(0)->vertex_dof_index(2,0) ==
		 face->child(1)->vertex_dof_index(3,0)) &&
		(face->child(0)->vertex_dof_index(2,0) ==
		 face->child(2)->vertex_dof_index(0,0)) &&
		(face->child(0)->vertex_dof_index(2,0) ==
		 face->child(3)->vertex_dof_index(1,0)),
		ExcInternalError());
	for (unsigned int dof=0; dof!=fe.dofs_per_vertex; ++dof)
	  dofs_on_children.push_back (face->child(0)->vertex_dof_index(2,dof));
	
					 // dof numbers on the centers of
					 // the lines bounding this face
	for (unsigned int line=0; line<4; ++line)
	  for (unsigned int dof=0; dof!=fe.dofs_per_vertex; ++dof)
	    dofs_on_children.push_back (face->line(line)->child(0)->vertex_dof_index(1,dof));

					 // next the dofs on the lines interior
					 // to the face; the order of these
					 // lines is laid down in the
					 // FiniteElement class documentation
	for (unsigned int dof=0; dof<fe.dofs_per_line; ++dof)
	  dofs_on_children.push_back (face->child(0)->line(1)->dof_index(dof));
	for (unsigned int dof=0; dof<fe.dofs_per_line; ++dof)
	  dofs_on_children.push_back (face->child(1)->line(2)->dof_index(dof));
	for (unsigned int dof=0; dof<fe.dofs_per_line; ++dof)
	  dofs_on_children.push_back (face->child(2)->line(3)->dof_index(dof));
	for (unsigned int dof=0; dof<fe.dofs_per_line; ++dof)
	  dofs_on_children.push_back (face->child(3)->line(0)->dof_index(dof));

					 // dofs on the bordering lines
	for (unsigned int line=0; line<4; ++line)
	  for (unsigned int child=0; child<2; ++child)
	    for (unsigned int dof=0; dof!=fe.dofs_per_line; ++dof)
	      dofs_on_children.push_back (face->line(line)->child(child)->dof_index(dof));
	
					 // finally, for the dofs interior
					 // to the four child faces
	for (unsigned int child=0; child<4; ++child)
	  for (unsigned int dof=0; dof!=fe.dofs_per_quad; ++dof)
	    dofs_on_children.push_back (face->child(child)->dof_index(dof));

	Assert (dofs_on_children.size() ==
	       fe.constraints().m(),
	       ExcDimensionMismatch(dofs_on_children.size(),
				      fe.constraints().m()));
	Assert (dofs_on_mother.size() ==
	       fe.constraints().n(),
	       ExcDimensionMismatch(dofs_on_mother.size(),
				      fe.constraints().n()));

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
  vector<int> dof_indices (dofs_per_cell);

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
template void
DoFTools::make_flux_sparsity_pattern (const DoFHandler<deal_II_dimension>& dof,
				      SparsityPattern    &sparsity);
template void
DoFTools::make_boundary_sparsity_pattern (const DoFHandler<deal_II_dimension>& dof,
					  const vector<int>  &,
					  SparsityPattern    &);
template void
DoFTools::make_boundary_sparsity_pattern (const DoFHandler<deal_II_dimension>& dof,
					  const DoFHandler<deal_II_dimension>::FunctionMap  &boundary_indicators,
					  const vector<int>  &dof_to_boundary_mapping,
					  SparsityPattern    &sparsity);
#endif

template void
DoFTools::make_sparsity_pattern (const DoFHandler<deal_II_dimension>& dof,
					       SparsityPattern    &sparsity);

template void 
DoFTools::make_sparsity_pattern (const DoFHandler<deal_II_dimension>& dof,
				 const vector<vector<bool> > &mask,
				 SparsityPattern             &sparsity);

template
void
DoFTools::distribute_cell_to_dof_vector (const DoFHandler<deal_II_dimension> &dof_handler,
					 const Vector<float>  &cell_data,
					 Vector<double>       &dof_data) const;

template
void
DoFTools::distribute_cell_to_dof_vector (const DoFHandler<deal_II_dimension> &dof_handler,
					 const Vector<double> &cell_data,
					 Vector<double>       &dof_data) const;



template void DoFTools::extract_dofs(const DoFHandler<deal_II_dimension>& dof,
				     const vector<bool>& local_select,
				     vector<bool>& selected_dofs);
template void DoFTools::extract_level_dofs(unsigned int level,
					   const MGDoFHandler<deal_II_dimension>& dof,
					   const vector<bool>& local_select,
					   vector<bool>& selected_dofs);

