/* $Id$ */
/* Copyright W. Bangerth, University of Heidelberg, 1990 */


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
void FESystem<dim>::initialize ()
{
				   // Initialize index table
				   // Multi-component base elements have to be thought of.
				   // 1. Vertices
  unsigned total_index = 0;
  for (unsigned vertex_number= 0 ; vertex_number < GeometryInfo<2>::vertices_per_cell ;
       ++vertex_number)
    {
      unsigned comp_start = 0;
      for(unsigned comp = 0; comp < n_component_elements() ;
	  ++comp)
	{
	  for (unsigned m = 0; m < element_multiplicity(comp); ++m)
	    {
	      for (unsigned local_index = 0 ;
		   local_index < base_element(comp).dofs_per_vertex ;
		   ++local_index)
		{
		  system_to_component_table[total_index++]
		    = pair<unsigned,unsigned>
		    (comp_start+m,
		     vertex_number*base_element(comp).dofs_per_vertex+local_index);
		}
	    }
	  comp_start += element_multiplicity(comp);
	}
    }
  
				   // 2. Lines
  for (unsigned line_number= 0 ; line_number < GeometryInfo<2>::lines_per_cell ;
       ++line_number)
    {
      unsigned comp_start = 0;
      for(unsigned comp = 0; comp < n_component_elements() ;
	  ++comp)
	{
	  for (unsigned m = 0; m < element_multiplicity(comp); ++m)
	    {
	      for (unsigned local_index = 0 ;
		   local_index < base_element(comp).dofs_per_line ;
		   ++local_index)
		{
		  system_to_component_table[total_index++]
		    = pair<unsigned,unsigned>
		    (comp_start+m,
		     line_number*base_element(comp).dofs_per_line
		     +local_index+base_element(comp).first_line_index);
		}
	    }
	  comp_start += element_multiplicity(comp);
	}
    }
  
				   // 3. Quads
  for (unsigned quad_number= 0 ; quad_number < GeometryInfo<2>::quads_per_cell ;
       ++quad_number)
    {
      unsigned comp_start = 0;
      for(unsigned comp = 0; comp < n_component_elements() ;
	  ++comp)
	{
	  for (unsigned m = 0; m < element_multiplicity(comp); ++m)
	    {
	      for (unsigned local_index = 0 ;
		   local_index < base_element(comp).dofs_per_quad ;
		   ++local_index)
		{
		  system_to_component_table[total_index++]
		    = pair<unsigned,unsigned>
		    (comp_start+m,
		     quad_number*base_element(comp).dofs_per_quad
		     +local_index+base_element(comp).first_quad_index);
		}
	    }
	  comp_start += element_multiplicity(comp);
	}
    }
  
				   // 4. Hex
  for (unsigned hex_number= 0 ; hex_number < GeometryInfo<2>::hexes_per_cell ;
       ++hex_number)
    {
      unsigned comp_start = 0;
      for(unsigned comp = 0; comp < n_component_elements() ;
	  ++comp, comp_start += element_multiplicity(comp))
	{
	  for (unsigned m = 0; m < element_multiplicity(comp); ++m)
	    {
	      for (unsigned local_index = 0 ;
		   local_index < base_element(comp).dofs_per_hex ;
		   ++local_index)
		{
		  system_to_component_table[total_index++]
		    = pair<unsigned,unsigned>
		    (comp_start+m,
		     hex_number*base_element(comp).dofs_per_hex
		     +local_index+base_element(comp).first_hex_index);
		}
	    }
	}
    }
  
  

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
			       fe1.n_transform_functions + fe2.n_transform_functions ,
			       fe1.n_components * N1 + fe2.n_components * N2 );
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
			       fe1.n_transform_functions + fe2.n_transform_functions ,
			       fe1.n_components * N1 + fe2.n_components * N2 );
};

#endif


template <int dim>
double FESystem<dim>::shape_value (const unsigned int i,
				   const Point<dim>  &p) const
{
  Assert(false, ExcNotImplemented());
  return 0.;
  
  Assert((i<total_dofs), ExcInvalidIndex(i));


//  return base_element->shape_value (i / n_sub_elements, p);
};



template <int dim>
Tensor<1,dim>
FESystem<dim>::shape_grad (const unsigned int  i,
			   const Point<dim>   &p) const
{
  Assert(false, ExcNotImplemented());
  return Tensor<1,dim>();

  Assert((i<total_dofs), ExcInvalidIndex(i));

//  return base_element->shape_grad (i / n_sub_elements, p);
};



template <int dim>
Tensor<2,dim>
FESystem<dim>::shape_grad_grad (const unsigned int  i,
				const Point<dim>   &p) const
 {
  Assert(false, ExcNotImplemented());
  return Tensor<2,dim>();

  Assert((i<total_dofs), ExcInvalidIndex(i));

//  return base_element->shape_grad_grad (i / n_sub_elements, p);
};



template <int dim>
void FESystem<dim>::get_unit_support_points (vector<Point<dim> > &support_points) const
{
/*  Assert (support_points.size() == total_dofs,
	  ExcWrongFieldDimension (support_points.size(), total_dofs));

  vector<Point<dim> > base_support_points (base_element->total_dofs);
  base_element->get_unit_support_points (base_support_points);
  
  for (unsigned int i=0; i<base_element->total_dofs; ++i) 
    for (unsigned int n=0; n<n_sub_elements; ++n)
      support_points[i*n_sub_elements+n] = base_support_points[i];
*/
};



template <int dim>
void FESystem<dim>::get_support_points (const DoFHandler<dim>::cell_iterator &cell,
					const Boundary<dim> &boundary,
					vector<Point<dim> > &support_points) const
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
void FESystem<dim>::get_face_support_points (const DoFHandler<dim>::face_iterator &face,
					     const Boundary<dim> &boundary,
					     vector<Point<dim> > &support_points) const
{
  Assert(false, ExcNotImplemented());
/*
  Assert (support_points.size() == dofs_per_face,
	  ExcWrongFieldDimension (support_points.size(), dofs_per_face));

  vector<Point<dim> > base_support_points (base_element->dofs_per_face);
  base_element->get_face_support_points (face, boundary, base_support_points);
  
  for (unsigned int i=0; i<base_element->dofs_per_face; ++i) 
    for (unsigned int n=0; n<n_sub_elements; ++n)
      support_points[i*n_sub_elements+n] = base_support_points[i];
*/
};



template <int dim>
void FESystem<dim>::get_local_mass_matrix (const DoFHandler<dim>::cell_iterator &cell,
					   const Boundary<dim> &boundary,
					   dFMatrix            &local_mass_matrix) const
{
  Assert(false, ExcNotImplemented());
/*
  Assert (local_mass_matrix.n() == total_dofs,
	  ExcWrongFieldDimension(local_mass_matrix.n(),total_dofs));
  Assert (local_mass_matrix.m() == total_dofs,
	  ExcWrongFieldDimension(local_mass_matrix.m(),total_dofs));

				   // first get the local mass matrix for
				   // the base object
  dFMatrix base_mass_matrix (base_element->total_dofs,
			     base_element->total_dofs);
  base_element->get_local_mass_matrix (cell, boundary, base_mass_matrix);


				   // now distribute it to the mass matrix
				   // of this object
  for (unsigned int i=0; i<base_element->total_dofs; ++i)
    for (unsigned int j=0; j<base_element->total_dofs; ++j)
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
  Assert(false, ExcNotImplemented());
  return 0.;
  
//  return base_element->shape_value_transform (i, p);
};


template <int dim>
Tensor<1,dim> FESystem<dim>::shape_grad_transform (const unsigned int i,
						   const Point<dim>  &p) const
{
  Assert(false, ExcNotImplemented());
  return Tensor<1,dim>();

//  return base_element->shape_grad_transform (i, p);
};



template <int dim>
void FESystem<dim>::get_face_jacobians (const DoFHandler<dim>::face_iterator &face,
					const Boundary<dim>         &boundary,
					const vector<Point<dim-1> > &unit_points,
					vector<double>      &face_jacobi_determinants) const
{
  Assert(false, ExcNotImplemented());

//  base_element->get_face_jacobians (face, boundary, unit_points, face_jacobi_determinants);
};



template <int dim>
void FESystem<dim>::get_subface_jacobians (const DoFHandler<dim>::face_iterator &face,
					   const unsigned int           subface_no,
					   const vector<Point<dim-1> > &unit_points,
					   vector<double>      &face_jacobi_determinants) const
{
  Assert(false, ExcNotImplemented());

//  base_element->get_subface_jacobians (face, subface_no, unit_points, face_jacobi_determinants);
};




template <int dim>
void FESystem<dim>::get_normal_vectors (const DoFHandler<dim>::cell_iterator &cell,
					const unsigned int          face_no,
					const Boundary<dim>         &boundary,
					const vector<Point<dim-1> > &unit_points,
					vector<Point<dim> >         &normal_vectors) const
{
  Assert(false, ExcNotImplemented());

//  base_element->get_normal_vectors (cell, face_no, boundary, unit_points, normal_vectors);
};



template <int dim>
void FESystem<dim>::get_normal_vectors (const DoFHandler<dim>::cell_iterator &cell,
					const unsigned int          face_no,
					const unsigned int          subface_no,
					const vector<Point<dim-1> > &unit_points,
					vector<Point<dim> >         &normal_vectors) const
{
  Assert(false, ExcNotImplemented());

//  base_element->get_normal_vectors (cell, face_no, subface_no, unit_points, normal_vectors);
};





// explicit instantiations
template class FESystem<deal_II_dimension>;



