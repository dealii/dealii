/* $Id$ */
/* Copyright W. Bangerth, G. Kanschat University of Heidelberg, 1999 */


#include <fe/fe_system.h>





template <int dim>
FESystem<dim>::~FESystem ()
{
  for (unsigned i=0;i<base_elements.size();++i)
    {
      base_elements[i].first -> unsubscribe ();
      delete base_elements[i].first;
    }
};


template <int dim>
void
FESystem<dim>::build_cell_table()
{
  unsigned total_index = 0;
  for (unsigned base=0 ; base < n_base_elements() ; ++base)
    for (unsigned m = 0; m < element_multiplicity(base); ++m)
      component_to_base_table[total_index++] = base;

				   // Initialize index table
				   // Multi-component base elements have
				   // to be thought of.
  
				   // 1. Vertices
  total_index = 0;
  for (unsigned vertex_number= 0 ; vertex_number < GeometryInfo<dim>::vertices_per_cell ;
       ++vertex_number)
    {
      unsigned comp_start = 0;
      for(unsigned base = 0; base < n_base_elements() ;
	  ++base)
	{
	  for (unsigned m = 0; m < element_multiplicity(base); ++m)
	    {
	      for (unsigned local_index = 0 ;
		   local_index < base_element(base).dofs_per_vertex ;
		   ++local_index)
		{
		  system_to_component_table[total_index++]
		    = pair<unsigned,unsigned>
		    (comp_start+m,
		     vertex_number*base_element(base).dofs_per_vertex+local_index);
		}
	    }
	  comp_start += element_multiplicity(base);
	}
    }
  
				   // 2. Lines
  for (unsigned line_number= 0 ; ((line_number != GeometryInfo<dim>::lines_per_cell) &&
				  (GeometryInfo<dim>::lines_per_cell > 0));
       ++line_number)
    {
      unsigned comp_start = 0;
      for(unsigned base = 0; base < n_base_elements() ;
	  ++base)
	{
	  for (unsigned m = 0; m < element_multiplicity(base); ++m)
	    {
	      for (unsigned local_index = 0 ;
		   local_index < base_element(base).dofs_per_line ;
		   ++local_index)
		{
		  system_to_component_table[total_index++]
		    = pair<unsigned,unsigned>
		    (comp_start+m,
		     line_number*base_element(base).dofs_per_line
		     +local_index+base_element(base).first_line_index);
		}
	    }
	  comp_start += element_multiplicity(base);
	}
    }
  
				   // 3. Quads
  for (unsigned quad_number= 0 ; ((quad_number != GeometryInfo<dim>::quads_per_cell) &&
				  (GeometryInfo<dim>::quads_per_cell > 0));
       ++quad_number)
    {
      unsigned comp_start = 0;
      for(unsigned base = 0; base < n_base_elements() ;
	  ++base)
	{
	  for (unsigned m = 0; m < element_multiplicity(base); ++m)
	    {
	      for (unsigned local_index = 0 ;
		   local_index < base_element(base).dofs_per_quad ;
		   ++local_index)
		{
		  system_to_component_table[total_index++]
		    = pair<unsigned,unsigned>
		    (comp_start+m,
		     quad_number*base_element(base).dofs_per_quad
		     +local_index+base_element(base).first_quad_index);
		}
	    }
	  comp_start += element_multiplicity(base);
	}
    }
  
				   // 4. Hex
  for (unsigned hex_number= 0 ; ((hex_number != GeometryInfo<dim>::hexes_per_cell) &&
				 (GeometryInfo<dim>::hexes_per_cell > 0));
       ++hex_number)
    {
      unsigned comp_start = 0;
      for(unsigned base = 0; base < n_base_elements() ;
	  ++base)
	{
	  for (unsigned m = 0; m < element_multiplicity(base); ++m)
	    {
	      for (unsigned local_index = 0 ;
		   local_index < base_element(base).dofs_per_hex ;
		   ++local_index)
		{
		  system_to_component_table[total_index++]
		    = pair<unsigned,unsigned>
		    (comp_start+m,
		     hex_number*base_element(base).dofs_per_hex
		     +local_index+base_element(base).first_hex_index);
		}
	    }
	  comp_start += element_multiplicity(base);
	  
	}
    }
				   // Inintialize mapping from component
				   // to base element
				   // Initialize mapping from components to
				   // linear index. Fortunately, this is
				   // the inverse of what we just did.
  for (unsigned comp=0 ; comp<n_components ; ++comp)
    component_to_system_table[comp]
      .resize(base_element(component_to_base_table[comp]).total_dofs);

  for (unsigned sys=0 ; sys < total_dofs ; ++sys)
    component_to_system_table[system_to_component_table[sys].first]
      [system_to_component_table[sys].second] = sys;
}


template <int dim>
void
FESystem<dim>::build_face_table()
{
  unsigned total_index = 0;
  for (unsigned base=0 ; base < n_base_elements() ; ++base)
    for (unsigned m = 0; m < element_multiplicity(base); ++m)
      component_to_base_table[total_index++] = base;

				   // Initialize index table
				   // Multi-component base elements have
				   // to be thought of.
  
				   // 1. Vertices
  total_index = 0;
  for (unsigned vertex_number= 0 ; vertex_number < GeometryInfo<dim>::vertices_per_face ;
       ++vertex_number)
    {
      unsigned comp_start = 0;
      for(unsigned base = 0; base < n_base_elements() ;
	  ++base)
	{
	  for (unsigned m = 0; m < element_multiplicity(base); ++m)
	    {
	      for (unsigned local_index = 0 ;
		   local_index < base_element(base).dofs_per_vertex ;
		   ++local_index)
		{
		  face_system_to_component_table[total_index++]
		    = pair<unsigned,unsigned>
		    (comp_start+m,
		     vertex_number*base_element(base).dofs_per_vertex+local_index);
		}
	    }
	  comp_start += element_multiplicity(base);
	}
    }
  Assert (total_index <= face_system_to_component_table.size(),
	  ExcInternalError());
  
				   // 2. Lines
  for (unsigned line_number= 0 ; ((line_number != GeometryInfo<dim>::lines_per_face) &&
				  (GeometryInfo<dim>::lines_per_cell > 0));
       ++line_number)
    {
      unsigned comp_start = 0;
      for(unsigned base = 0; base < n_base_elements() ;
	  ++base)
	{
	  for (unsigned m = 0; m < element_multiplicity(base); ++m)
	    {
	      for (unsigned local_index = 0 ;
		   local_index < base_element(base).dofs_per_line ;
		   ++local_index)
		{
		  face_system_to_component_table[total_index++]
		    = pair<unsigned,unsigned>
		    (comp_start+m,
		     line_number*base_element(base).dofs_per_line
		     +local_index+base_element(base).first_face_line_index);
		}
	    }
	  comp_start += element_multiplicity(base);
	}
    }
  Assert (total_index <= face_system_to_component_table.size(),
	  ExcInternalError());
  
				   // 3. Quads
  for (unsigned quad_number= 0 ; ((quad_number != GeometryInfo<dim>::quads_per_face) &&
				  (GeometryInfo<dim>::quads_per_cell > 0));
       ++quad_number)
    {
      unsigned comp_start = 0;
      for(unsigned base = 0; base < n_base_elements() ;
	  ++base)
	{
	  for (unsigned m = 0; m < element_multiplicity(base); ++m)
	    {
	      for (unsigned local_index = 0 ;
		   local_index < base_element(base).dofs_per_quad ;
		   ++local_index)
		{
		  face_system_to_component_table[total_index++]
		    = pair<unsigned,unsigned>
		    (comp_start+m,
		     quad_number*base_element(base).dofs_per_quad
		     +local_index+base_element(base).first_face_quad_index);
		}
	    }
	  comp_start += element_multiplicity(base);
	}
    }
  Assert (total_index <= face_system_to_component_table.size(),
	  ExcInternalError());
  
				   // Inintialize mapping from component
				   // to base element
				   // Initialize mapping from components to
				   // linear index. Fortunately, this is
				   // the inverse of what we just did.
  for (unsigned comp=0 ; comp<n_components ; ++comp)
    face_component_to_system_table[comp]
      .resize(base_element(component_to_base_table[comp]).total_dofs);

  for (unsigned sys=0 ; sys < dofs_per_face ; ++sys)
    face_component_to_system_table[face_system_to_component_table[sys].first]
      [face_system_to_component_table[sys].second] = sys;
}



template <int dim>
void FESystem<dim>::build_interface_constraints () 
{
				   // the layout of the constraints matrix is
				   // described in the FiniteElement class. you
				   // may want to look there first before trying
				   // to understand the following, especially
				   // the mapping of the #n# index.
				   //
				   // in order to map it to the fe-system class,
				   // we have to know which base element a
				   // degree of freedom within a vertex, line,
				   // etc belongs to. this can be accomplished
				   // by the system_to_component_index
				   // function in conjunction with the
				   // numbers first_{line,quad,...}_index
  for (unsigned int n=0; n<interface_constraints.n(); ++n)
    for (unsigned int m=0; m<interface_constraints.m(); ++m)
      {
					 // for the pair (n,m) find out which
					 // component they belong to and
					 // the number therein
					 //
					 // first value in pair is component,
					 // second is index
	const pair<unsigned int, unsigned int> n_index
	  = face_system_to_component_index (n);

	pair<unsigned int, unsigned int> m_index;
	switch (dim)
	  {
	    case 1:
	    {
					       // we should never get here!
					       // (in 1d, the constraints matrix
					       // should be of size zero)
	      Assert (false, ExcInternalError());
	      break;
	    };

	    case 2:
	    {
					       // the indices m=0..d_v-1 are
					       // from the center vertex.
					       // their order is the same
					       // as for the first vertex
					       // of the whole cell, so we
					       // can use the
					       // system_to_component_index
					       // function (using the
					       // face_s_t_c_i function would
					       // yield the same)
	      if (m < dofs_per_vertex)
		m_index = system_to_component_index (m);
	      else
						 // then come the two sets of
						 // line indices
		{
		  const unsigned int index_in_line
		    = (m-dofs_per_vertex) % dofs_per_line;
		  const unsigned int sub_line
		    = (m-dofs_per_vertex) / dofs_per_line;
		  Assert (sub_line < 2, ExcInternalError());
		  
						   // get the component by
						   // asking s_t_c_index and
						   // tweaking the index a bit
		  m_index.first = system_to_component_index
				  (GeometryInfo<2>::vertices_per_cell * dofs_per_vertex
				   + index_in_line).first;
						   // first find out the how-many'th
						   // line index of that component
						   // this was
		  m_index.second = (system_to_component_index
				    (GeometryInfo<2>::vertices_per_cell * dofs_per_vertex
				     + index_in_line).second
				    - base_element (component_to_base_table[m_index.first]).first_line_index)
						   // then add the number of dofs
						   // per vertex to get the index
						   // on the first line
				   + base_element(component_to_base_table[m_index.first]).dofs_per_vertex
						   // if on the second line: add
						   // some more
				   + base_element(component_to_base_table[m_index.first]).dofs_per_line * sub_line;
		};
	      break;
	    };

	    case 3:
	    {
					       // the indices
					       // m=0..5*d_v-1 are
					       // from the center and
					       // the four subline
					       // vertices.  their
					       // order is the same as
					       // for the first vertex
					       // of the whole cell,
					       // so we can use the
					       // simple arithmetic
	      if (m < 5*dofs_per_vertex)
		{
		  m_index.first  = system_to_component_index(m % dofs_per_vertex).first;
		  m_index.second = m / dofs_per_vertex;
		}
	      else
						 // then come the 12 sets of
						 // line indices
		if (m < 5*dofs_per_vertex + 12*dofs_per_line)
		  {   
		    const unsigned int index_in_line
		      = (m-5*dofs_per_vertex) % dofs_per_line;
		    const unsigned int sub_line
		      = (m-5*dofs_per_vertex) / dofs_per_line;
		    Assert (sub_line < 12, ExcInternalError());
		  
						     // get the component by
						     // asking s_t_c_index and
						     // tweaking the index a bit
		    m_index.first = system_to_component_index
				    (GeometryInfo<3>::vertices_per_cell * dofs_per_vertex
				     + index_in_line).first;
		    
						     // first find out the how-many'th
						     // line index of that component
						     // this was
		    m_index.second = (system_to_component_index
				      (GeometryInfo<3>::vertices_per_cell * dofs_per_vertex
				       + index_in_line).second
				      - base_element (component_to_base_table[m_index.first]).first_line_index)
								      // then add the number of dofs
								      // for the five vertices to get
								      // the index on the first line
				     + 5*base_element(component_to_base_table[m_index.first]).dofs_per_vertex
								      // and correct for the
								      // how-many'th line
				     + base_element(component_to_base_table[m_index.first]).dofs_per_line * sub_line;
		  }
		else
						   // on one of the four sub-quads
		  {   
		    const unsigned int index_in_quad
		      = (m-5*dofs_per_vertex-12*dofs_per_line) % dofs_per_line;
		    const unsigned int sub_quad
		      = (m-5*dofs_per_vertex-12*dofs_per_line) / dofs_per_line;
		    Assert (sub_quad < 4, ExcInternalError());
		  
						     // get the component by
						     // asking s_t_c_index and
						     // tweaking the index a bit
		    m_index.first = system_to_component_index
				    (GeometryInfo<3>::vertices_per_cell * dofs_per_vertex
				     + GeometryInfo<3>::lines_per_cell * dofs_per_line
				     + index_in_quad).first;
		    
						     // first find out the how-many'th
						     // quad index of that component
						     // this was
		    m_index.second = (system_to_component_index
				      (GeometryInfo<3>::vertices_per_cell * dofs_per_vertex
				       + GeometryInfo<3>::lines_per_cell * dofs_per_line
				       + index_in_quad).second
				      - base_element (component_to_base_table[m_index.first]).first_quad_index)
								      // then add the number of dofs
								      // for the five vertices and 12 lines
								      // to get the index on the first quad
				     + 5*base_element(component_to_base_table[m_index.first]).dofs_per_vertex
				     + 12*base_element(component_to_base_table[m_index.first]).dofs_per_line
								      // and correct for the
								      // how-many'th line
				     + base_element(component_to_base_table[m_index.first]).dofs_per_quad * sub_quad;
		  };
	      
	      break;
	    };
		  
	    default:
		  Assert (false, ExcNotImplemented());
	  };

					 // now that we gathered all information:
					 // use it to build the matrix. note
					 // that if n and m belong to different
					 // components, there definitely will be
					 // no coupling
	if (n_index.first == m_index.first)
	  interface_constraints(m,n)
	    = (base_element(component_to_base_table[n_index.first])
	       .constraints()(m_index.second,
			      n_index.second));
      };
};





template <int dim>
void FESystem<dim>::initialize ()
{
  build_cell_table();
  build_face_table();
				   // distribute the matrices of the base
				   // finite elements to the matrices of
				   // this object
  for (unsigned int component=0; component<n_components; ++component)
				     // transform restriction and
				     // prolongation matrices
    for (unsigned int i=0; i<base_element(component_to_base_table[component]).total_dofs; ++i)
      for (unsigned int j=0; j<base_element(component_to_base_table[component]).total_dofs; ++j)
					 // only fill block diagonals, no
					 // intermixing of subelements
	for (unsigned int child=0; child<GeometryInfo<dim>::children_per_cell; ++child)
	  {
	    restriction[child] (component_to_system_index (component,i),
				component_to_system_index (component, j))
	      = base_element(component_to_base_table[component]).restrict(child)(i,j);
	    prolongation[child] (component_to_system_index (component,i),
				 component_to_system_index (component, j))
	      = base_element(component_to_base_table[component]).prolongate(child)(i,j);
	  };


				   // now set up the interface constraints.
				   // this is kind'o hairy, so don't try
				   // to do it dimension independent
  build_interface_constraints ();
};




#if deal_II_dimension == 1

template <>
FiniteElementData<1>
FESystem<1>::multiply_dof_numbers (const FiniteElementData<1> &fe_data,
				   const unsigned int          N)
{
  return FiniteElementData<1> (fe_data.dofs_per_vertex * N,
			       fe_data.dofs_per_line * N,
			       fe_data.n_transform_functions,
			       fe_data.n_components * N);
};


template <>
FiniteElementData<1>
FESystem<1>::multiply_dof_numbers (const FiniteElementData<1> &fe1,
				   const unsigned int          N1,
				   const FiniteElementData<1> &fe2,
				   const unsigned int          N2)
{
  return FiniteElementData<1> (fe1.dofs_per_vertex * N1 + fe2.dofs_per_vertex * N2 ,
			       fe1.dofs_per_line * N1 + fe2.dofs_per_line * N2 ,
			       fe1.n_transform_functions,
			       fe1.n_components * N1 + fe2.n_components * N2 );
};


template <>
FiniteElementData<1>
FESystem<1>::multiply_dof_numbers (const FiniteElementData<1> &fe1,
				   const unsigned int          N1,
				   const FiniteElementData<1> &fe2,
				   const unsigned int          N2,
				   const FiniteElementData<1> &fe3,
				   const unsigned int          N3)
{
  return FiniteElementData<1> (fe1.dofs_per_vertex * N1
			       + fe2.dofs_per_vertex * N2
			       + fe3.dofs_per_vertex * N3,
			       fe1.dofs_per_line * N1
			       + fe2.dofs_per_line * N2
			       + fe3.dofs_per_line * N3,
			       fe1.n_transform_functions,
			       fe1.n_components * N1
			       + fe2.n_components * N2
			       + fe3.n_components * N3);
};

#endif


#if deal_II_dimension == 2

template <>
FiniteElementData<2>
FESystem<2>::multiply_dof_numbers (const FiniteElementData<2> &fe_data,
				   const unsigned int          N)
{
  return FiniteElementData<2> (fe_data.dofs_per_vertex * N,
			       fe_data.dofs_per_line * N,
			       fe_data.dofs_per_quad * N,
			       fe_data.n_transform_functions,
			       fe_data.n_components * N);
};

template <>
FiniteElementData<2>
FESystem<2>::multiply_dof_numbers (const FiniteElementData<2> &fe1,
				   const unsigned int          N1,
				   const FiniteElementData<2> &fe2,
				   const unsigned int          N2)
{
  return FiniteElementData<2> (fe1.dofs_per_vertex * N1 + fe2.dofs_per_vertex * N2 ,
			       fe1.dofs_per_line * N1 + fe2.dofs_per_line * N2 ,
			       fe1.dofs_per_quad * N1 + fe2.dofs_per_quad * N2 ,
			       fe1.n_transform_functions,
			       fe1.n_components * N1 + fe2.n_components * N2 );
};


template <>
FiniteElementData<2>
FESystem<2>::multiply_dof_numbers (const FiniteElementData<2> &fe1,
				   const unsigned int          N1,
				   const FiniteElementData<2> &fe2,
				   const unsigned int          N2,
				   const FiniteElementData<2> &fe3,
				   const unsigned int          N3)
{
  return FiniteElementData<2> (fe1.dofs_per_vertex * N1
			       + fe2.dofs_per_vertex * N2
			       + fe3.dofs_per_vertex * N3 ,
			       fe1.dofs_per_line * N1
			       + fe2.dofs_per_line * N2
			       + fe3.dofs_per_line * N3 ,
			       fe1.dofs_per_quad * N1
			       + fe2.dofs_per_quad * N2
			       + fe3.dofs_per_quad * N3 ,
			       fe1.n_transform_functions,
			       fe1.n_components * N1
			       + fe2.n_components * N2
			       + fe3.n_components * N3 );
};

#endif


#if deal_II_dimension == 3

template <>
FiniteElementData<3>
FESystem<3>::multiply_dof_numbers (const FiniteElementData<3> &fe_data,
				   const unsigned int          N)
{
  return FiniteElementData<3> (fe_data.dofs_per_vertex * N,
			       fe_data.dofs_per_line * N,
			       fe_data.dofs_per_quad * N,
			       fe_data.dofs_per_hex * N,
			       fe_data.n_transform_functions,
			       fe_data.n_components * N);
};

template <>
FiniteElementData<3>
FESystem<3>::multiply_dof_numbers (const FiniteElementData<3> &fe1,
				   const unsigned int          N1,
				   const FiniteElementData<3> &fe2,
				   const unsigned int          N2)
{
  return FiniteElementData<3> (fe1.dofs_per_vertex * N1 + fe2.dofs_per_vertex * N2 ,
			       fe1.dofs_per_line * N1 + fe2.dofs_per_line * N2 ,
			       fe1.dofs_per_quad * N1 + fe2.dofs_per_quad * N2 ,
			       fe1.dofs_per_hex * N1 + fe2.dofs_per_hex * N2 ,
			       fe1.n_transform_functions,
			       fe1.n_components * N1 + fe2.n_components * N2 );
};


template <>
FiniteElementData<3>
FESystem<3>::multiply_dof_numbers (const FiniteElementData<3> &fe1,
				   const unsigned int          N1,
				   const FiniteElementData<3> &fe2,
				   const unsigned int          N2,
				   const FiniteElementData<3> &fe3,
				   const unsigned int          N3)
{
  return FiniteElementData<3> (fe1.dofs_per_vertex * N1
			       + fe2.dofs_per_vertex * N2
			       + fe3.dofs_per_vertex * N3 ,
			       fe1.dofs_per_line * N1
			       + fe2.dofs_per_line * N2
			       + fe3.dofs_per_line * N3 ,
			       fe1.dofs_per_quad * N1
			       + fe2.dofs_per_quad * N2
			       + fe3.dofs_per_quad * N3 ,
			       fe1.dofs_per_hex * N1
			       + fe2.dofs_per_hex * N2
			       + fe3.dofs_per_hex * N3 ,
			       fe1.n_transform_functions,
			       fe1.n_components * N1
			       + fe2.n_components * N2
			       + fe3.n_components * N3 );
};

#endif


template <int dim>
vector<bool>
FESystem<dim>::compute_restriction_is_additive_flags (const FiniteElement<dim> &fe,
						      const unsigned int n_elements) 
{
  vector<bool> tmp;
  for (unsigned int i=0; i<n_elements; ++i)
    for (unsigned int component=0; component<fe.n_components; ++component)
      tmp.push_back (fe.restriction_is_additive (component));
  return tmp;
};


template <int dim>
vector<bool>
FESystem<dim>::compute_restriction_is_additive_flags (const FiniteElement<dim> &fe1,
						      const unsigned int        N1,
						      const FiniteElement<dim> &fe2,
						      const unsigned int        N2) 
{
  vector<bool> tmp;
  for (unsigned int i=0; i<N1; ++i)
    for (unsigned int component=0; component<fe1.n_components; ++component)
      tmp.push_back (fe1.restriction_is_additive (component));
  for (unsigned int i=0; i<N2; ++i)
    for (unsigned int component=0; component<fe2.n_components; ++component)
      tmp.push_back (fe2.restriction_is_additive (component));
  return tmp;
};


template <int dim>
vector<bool>
FESystem<dim>::compute_restriction_is_additive_flags (const FiniteElement<dim> &fe1,
						      const unsigned int        N1,
						      const FiniteElement<dim> &fe2,
						      const unsigned int        N2,
						      const FiniteElement<dim> &fe3,
						      const unsigned int        N3) 
{
  vector<bool> tmp;
  for (unsigned int i=0; i<N1; ++i)
    for (unsigned int component=0; component<fe1.n_components; ++component)
      tmp.push_back (fe1.restriction_is_additive (component));
  for (unsigned int i=0; i<N2; ++i)
    for (unsigned int component=0; component<fe2.n_components; ++component)
      tmp.push_back (fe2.restriction_is_additive (component));
  for (unsigned int i=0; i<N3; ++i)
    for (unsigned int component=0; component<fe3.n_components; ++component)
      tmp.push_back (fe3.restriction_is_additive (component));
  return tmp;
};

		


template <int dim>
double
FESystem<dim>::shape_value (const unsigned int i,
			    const Point<dim>  &p) const
{
  Assert((i<total_dofs), ExcIndexRange(i, 0, total_dofs));

  pair<unsigned,unsigned> comp = system_to_component_index(i);
  
  return base_element(component_to_base_table[comp.first])
    .shape_value(comp.second, p);
};



template <int dim>
Tensor<1,dim>
FESystem<dim>::shape_grad (const unsigned int  i,
			   const Point<dim>   &p) const
{
  Assert((i<total_dofs), ExcIndexRange(i, 0, total_dofs));

  pair<unsigned,unsigned> comp = system_to_component_index(i);
  
  return base_element(component_to_base_table[comp.first])
    .shape_grad(comp.second, p);
};



template <int dim>
Tensor<2,dim>
FESystem<dim>::shape_grad_grad (const unsigned int  i,
				const Point<dim>   &p) const
 {
  Assert((i<total_dofs), ExcIndexRange(i, 0, total_dofs));


  pair<unsigned,unsigned> comp = system_to_component_index(i);
  
  return base_element(component_to_base_table[comp.first])
    .shape_grad_grad(comp.second, p);
};



template <int dim>
void FESystem<dim>::get_unit_support_points (
  vector<Point<dim> > &unit_support_points) const
{
  Assert(unit_support_points.size() == total_dofs,
	 ExcWrongFieldDimension (unit_support_points.size(), total_dofs));

  vector<Point<dim> > base_unit_support_points (base_element(0).total_dofs);
  unsigned int component = 0;
  for (unsigned int base_el=0 ; base_el<n_base_elements(); ++base_el)
    {
      const unsigned int base_element_total_dofs
	=base_element(base_el).total_dofs;
 
      base_unit_support_points.resize(base_element_total_dofs);
      base_element(base_el).get_unit_support_points (base_unit_support_points);
      for (unsigned int n = 0 ; n < element_multiplicity(base_el); ++n)
	{
	  for (unsigned int i=0; i<base_element_total_dofs; ++i)
	    {
	      unit_support_points[component_to_system_index(component,i)]
		= base_unit_support_points[i];
	    }
	  ++component;
	}
    }

				   // An alternative version

				   // base unit support points
//   vector<vector<Point<dim> > > base_us_points(n_base_elements());
//   for (unsigned int base_el=0 ; base_el<n_base_elements(); ++base_el)
//     {
//       const unsigned int base_element_total_dofs
// 	=base_element(base_el).total_dofs;

//       base_us_points[base_el].resize(base_element_total_dofs);
//       base_element(base_el).get_unit_support_points (base_us_points[base_el]);
//     }

//   for (unsigned int i=0; i<total_dofs; ++i)
//     {
//       const unsigned int comp=system_to_component_index(i).first,
// 		     base_dof=system_to_component_index(i).second,
// 		      base_el=component_to_base_table[comp];
      
//       unit_support_points[i]=base_us_points[base_el][base_dof];
//    }
};



template <int dim>
void FESystem<dim>::get_support_points (const DoFHandler<dim>::cell_iterator &cell,
					vector<Point<dim> > &support_points) const
{
  Assert(support_points.size() == total_dofs,
	 ExcWrongFieldDimension (support_points.size(), total_dofs));

  vector<Point<dim> > base_support_points (base_element(0).total_dofs);
  unsigned int component = 0;
  for (unsigned int base_el=0 ; base_el<n_base_elements(); ++base_el)
    {
      const unsigned int base_element_total_dofs
	=base_element(base_el).total_dofs;
      
      base_support_points.resize(base_element_total_dofs);
      base_element(base_el).get_support_points (cell, base_support_points);
      for (unsigned int n = 0 ; n < element_multiplicity(base_el); ++n)
	{
	  for (unsigned int i=0; i<base_element_total_dofs; ++i)
	    {
	      support_points[component_to_system_index(component,i)]
		= base_support_points[i];
	    }
	  ++component;
	}
    }
};



template <int dim>
void FESystem<dim>::get_face_support_points (const DoFHandler<dim>::face_iterator & face,
					     vector<Point<dim> > & support_points) const
{
  Assert (support_points.size() == dofs_per_face,
	  ExcWrongFieldDimension (support_points.size(), dofs_per_face));

  vector<Point<dim> > base_support_points (base_element(0).dofs_per_face);
  unsigned int comp = 0;
  for (unsigned int base=0 ; base<n_base_elements(); ++base)
    {
      base_support_points.resize(base_element(base).dofs_per_face);
      base_element(base).get_face_support_points (face, base_support_points);
      for (unsigned int inbase = 0 ; inbase < element_multiplicity(base); ++inbase)
	{
	  for (unsigned int i=0; i<base_element(base).dofs_per_face; ++i)
	    {
	      support_points[face_component_to_system_index(comp,i)]
		= base_support_points[i];
	    }
	  
	      ++comp;
	}
    }
}



template <int dim>
void FESystem<dim>::get_local_mass_matrix (const DoFHandler<dim>::cell_iterator &cell,
					   FullMatrix<double>  &local_mass_matrix) const
{
  Assert (local_mass_matrix.n() == total_dofs,
	  ExcWrongFieldDimension(local_mass_matrix.n(),total_dofs));
  Assert (local_mass_matrix.m() == total_dofs,
	  ExcWrongFieldDimension(local_mass_matrix.m(),total_dofs));

				   // track which component we are
				   // presently working with, since we
				   // only have the number of the base
				   // element and the number within
				   // its multiplicity
  unsigned int component = 0;  
  for (unsigned int base_el=0; base_el<n_base_elements(); ++base_el)
    {
				       // first get the local mass matrix for
				       // the base object
      const unsigned int base_element_total_dofs=base_element(base_el).total_dofs;
      FullMatrix<double> base_mass_matrix (base_element_total_dofs,
					   base_element_total_dofs);
      base_element(base_el).get_local_mass_matrix (cell, base_mass_matrix);
      
				       // now distribute it to the mass matrix
				       // of this object
      const unsigned int el_multiplicity=element_multiplicity(base_el);
      for (unsigned int n=0; n<el_multiplicity; ++n)
	{
	  for (unsigned int i=0; i<base_element_total_dofs; ++i)
	    for (unsigned int j=0; j<base_element_total_dofs; ++j)
					     // only fill diagonals of the blocks
	      local_mass_matrix (component_to_system_index(component,i),
				 component_to_system_index(component,j))
		= base_mass_matrix (i,j);
	  ++component;
	};
    };
  Assert (component == n_components, ExcInternalError());
};



template <int dim>
double FESystem<dim>::shape_value_transform (const unsigned int i,
					     const Point<dim>  &p) const
{
  return base_elements[0].first->shape_value_transform(i,p);
};


template <int dim>
Tensor<1,dim> FESystem<dim>::shape_grad_transform (const unsigned int i,
						   const Point<dim>  &p) const
{
  return base_elements[0].first->shape_grad_transform (i, p);
};



template <int dim>
void FESystem<dim>::get_face_jacobians (const DoFHandler<dim>::face_iterator &face,
					const vector<Point<dim-1> > &unit_points,
					vector<double>      &face_jacobi_determinants) const
{
  base_elements[0].first->get_face_jacobians (face, unit_points,
					      face_jacobi_determinants);
};



template <int dim>
void FESystem<dim>::get_subface_jacobians (const DoFHandler<dim>::face_iterator &face,
					   const unsigned int           subface_no,
					   const vector<Point<dim-1> > &unit_points,
					   vector<double>      &face_jacobi_determinants) const
{
  base_elements[0].first->get_subface_jacobians (face, subface_no, unit_points,
					   face_jacobi_determinants);
};




template <int dim>
void FESystem<dim>::get_normal_vectors (const DoFHandler<dim>::cell_iterator &cell,
					const unsigned int          face_no,
					const vector<Point<dim-1> > &unit_points,
					vector<Point<dim> >         &normal_vectors) const
{
  base_elements[0].first->get_normal_vectors (cell, face_no, unit_points,
					normal_vectors);
};



template <int dim>
void FESystem<dim>::get_normal_vectors (const DoFHandler<dim>::cell_iterator &cell,
					const unsigned int          face_no,
					const unsigned int          subface_no,
					const vector<Point<dim-1> > &unit_points,
					vector<Point<dim> >         &normal_vectors) const
{
  base_elements[0].first->get_normal_vectors (cell, face_no, subface_no, unit_points,
				    normal_vectors);
};


template <int dim>
void
FESystem<dim>::fill_fe_values (const DoFHandler<dim>::cell_iterator &cell,
			       const vector<Point<dim> >            &unit_points,
			       vector<Tensor<2,dim> >               &jacobians,
			       const bool              compute_jacobians,
			       vector<Tensor<3,dim> > &jacobians_grad,
			       const bool              compute_jacobians_grad,
			       vector<Point<dim> > &support_points,
			       const bool           compute_support_points,
			       vector<Point<dim> > &q_points,
			       const bool           compute_q_points,
			       const FullMatrix<double>  &shape_values_transform,
			       const vector<vector<Tensor<1,dim> > > &shape_grad_transform) const
{
  vector<Point<dim> > supp(base_elements[0].first->total_dofs);

  base_elements[0].first->fill_fe_values (cell, unit_points, jacobians, compute_jacobians,
					  jacobians_grad, compute_jacobians_grad,
					  support_points, compute_support_points,
					  q_points, compute_q_points,
					  shape_values_transform, shape_grad_transform);
  
  if (compute_support_points)
    {
      unsigned component = 0;
      
      for (unsigned m=0 ; m < element_multiplicity(0) ; ++ m)
	{
	  for (unsigned i=0 ; i < base_element(0).total_dofs ; ++i)
	    support_points[component_to_system_index(component,i)] = supp[i];
	  ++component;
	}
      for (unsigned base=1 ; base < n_base_elements() ; ++base)
	{
	  supp.resize(base_elements[base].first->total_dofs);
	  base_elements[base].first->fill_fe_values (cell, unit_points, jacobians, false,
						     jacobians_grad, false,
						     supp, true,
						     q_points, false,
						     shape_values_transform, shape_grad_transform);
	  
	  for (unsigned m=0 ; m < element_multiplicity(base) ; ++ m)
	    {
	      for (unsigned i=0 ; i < base_element(base).total_dofs ; ++i)
		support_points[component_to_system_index(component,i)] = supp[i];
	      ++component;
	    }
	}    
    }
}

    

// explicit instantiations
template class FESystem<deal_II_dimension>;



