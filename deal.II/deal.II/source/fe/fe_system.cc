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
void FESystem<dim>::initialize ()
{
  build_cell_table();
  build_face_table();
				   // distribute the matrices of the base
				   // finite element to the matrices of
				   // this object
/*  for (unsigned int i=0; i<base_element->total_dofs; ++i)
    for (unsigned int j=0; j<base_element->total_dofs; ++j)
      for (unsigned int n=0; n<n_sub_elements; ++n)
					 // only fill diagonals of the blocks
	{
	  for (unsigned int child=0; child<GeometryInfo<dim>::children_per_cell; ++child)
	    {
	      restriction[child] (i*n_sub_elements + n,
				  j*n_sub_elements + n)
		= base_element->restrict(child)(i,j);
	      prolongation[child] (i*n_sub_elements + n,
				   j*n_sub_elements + n)
		= base_element->prolongate(child)(i,j);
	    };

	  interface_constraints (i*n_sub_elements + n,
				 j*n_sub_elements + n)
	    = base_element->constraints()(i,j);
	};
*/};




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


template <int dim>
double
FESystem<dim>::shape_value (const unsigned int i,
			    const Point<dim>  &p) const
{
  Assert((i<total_dofs), ExcInvalidIndex(i));

  pair<unsigned,unsigned> comp = system_to_component_index(i);
  
  return base_element(component_to_base_table[comp.first])
    .shape_value(comp.second, p);
};



template <int dim>
Tensor<1,dim>
FESystem<dim>::shape_grad (const unsigned int  i,
			   const Point<dim>   &p) const
{
  Assert((i<total_dofs), ExcInvalidIndex(i));

  pair<unsigned,unsigned> comp = system_to_component_index(i);
  
  return base_element(component_to_base_table[comp.first])
    .shape_grad(comp.second, p);
};



template <int dim>
Tensor<2,dim>
FESystem<dim>::shape_grad_grad (const unsigned int  i,
				const Point<dim>   &p) const
 {
  Assert((i<total_dofs), ExcInvalidIndex(i));


  pair<unsigned,unsigned> comp = system_to_component_index(i);
  
  return base_element(component_to_base_table[comp.first])
    .shape_grad_grad(comp.second, p);
};



template <int dim>
void FESystem<dim>::get_unit_support_points (vector<Point<dim> > &/*support_points*/) const
{
  Assert(false, ExcNotImplemented());
/*
  Assert (support_points.size() == total_dofs,
	  ExcWrongFieldDimension (support_points.size(), total_dofs));

  vector<Point<dim> > base_support_points (base_element->total_dofs);
  base_element->get_unit_support_points (base_support_points);
  
  for (unsigned int i=0; i<base_element->total_dofs; ++i) 
    for (unsigned int n=0; n<n_sub_elements; ++n)
      support_points[i*n_sub_elements+n] = base_support_points[i];
*/
};



template <int dim>
void FESystem<dim>::get_support_points (const DoFHandler<dim>::cell_iterator &/*cell*/,
					const Boundary<dim> &/*boundary*/,
					vector<Point<dim> > &/*support_points*/) const
{
  Assert(false, ExcNotImplemented());
/*
  Assert (support_points.size() == total_dofs,
	  ExcWrongFieldDimension (support_points.size(), total_dofs));

  vector<Point<dim> > base_support_points (base_element->total_dofs);
  base_element->get_support_points (cell, boundary, base_support_points);
  
  for (unsigned int i=0; i<base_element->total_dofs; ++i) 
    for (unsigned int n=0; n<n_sub_elements; ++n)
      support_points[i*n_sub_elements+n] = base_support_points[i];
*/
};



template <int dim>
void FESystem<dim>::get_face_support_points (const DoFHandler<dim>::face_iterator & face,
					     const Boundary<dim> & boundary,
					     vector<Point<dim> > & support_points) const
{
  Assert (support_points.size() == dofs_per_face,
	  ExcWrongFieldDimension (support_points.size(), dofs_per_face));

  vector<Point<dim> > base_support_points (base_element(0).dofs_per_face);
  unsigned int comp = 0;
  for (unsigned int base=0 ; base<n_base_elements(); ++base)
    {
      base_support_points.resize(base_element(base).dofs_per_face);
      base_element(base).get_face_support_points (face, boundary, base_support_points);
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
void FESystem<dim>::get_local_mass_matrix (const DoFHandler<dim>::cell_iterator &/*cell*/,
					   const Boundary<dim> &/*boundary*/,
					   FullMatrix<double>  &/*local_mass_matrix*/) const
{
  Assert(false, ExcNotImplemented());
/*
  Assert (local_mass_matrix.n() == total_dofs,
	  ExcWrongFieldDimension(local_mass_matrix.n(),total_dofs));
  Assert (local_mass_matrix.m() == total_dofs,
	  ExcWrongFieldDimension(local_mass_matrix.m(),total_dofs));

				   // first get the local mass matrix for
				   // the base object
  FullMatrix<double> base_mass_matrix (base_element.total_dofs,
			     base_element.total_dofs);
  base_element.get_local_mass_matrix (cell, boundary, base_mass_matrix);


				   // now distribute it to the mass matrix
				   // of this object
  for (unsigned int i=0; i<base_element.total_dofs; ++i)
    for (unsigned int j=0; j<base_element.total_dofs; ++j)
      for (unsigned int n=0; n<n_sub_elements; ++n)
					 // only fill diagonals of the blocks
	local_mass_matrix (i*n_sub_elements + n,
			   j*n_sub_elements + n) = base_mass_matrix (i,j);
*/
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
					const Boundary<dim>         &boundary,
					const vector<Point<dim-1> > &unit_points,
					vector<double>      &face_jacobi_determinants) const
{
  base_elements[0].first->get_face_jacobians (face, boundary, unit_points,
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
					const Boundary<dim>         &boundary,
					const vector<Point<dim-1> > &unit_points,
					vector<Point<dim> >         &normal_vectors) const
{
  base_elements[0].first->get_normal_vectors (cell, face_no, boundary, unit_points,
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
			       const vector<vector<Tensor<1,dim> > > &shape_grad_transform,
			       const Boundary<dim> &boundary) const
{
  vector<Point<dim> > supp(base_elements[0].first->total_dofs);

  base_elements[0].first->fill_fe_values (cell, unit_points, jacobians, compute_jacobians,
					  jacobians_grad, compute_jacobians_grad,
					  support_points, compute_support_points,
					  q_points, compute_q_points,
					  shape_values_transform, shape_grad_transform, boundary);
  
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
						     shape_values_transform, shape_grad_transform, boundary);
	  
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



