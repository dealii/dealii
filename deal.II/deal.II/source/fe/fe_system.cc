/* $Id$ */
/* Copyright W. Bangerth, University of Heidelberg, 1990 */


#include <fe/fe_system.h>





template <int dim>
FESystem<dim>::~FESystem ()
{
  base_element->unsubscribe ();
  delete base_element;
};



template <int dim>
void FESystem<dim>::initialize_matrices ()
{
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
				   const unsigned int          N) {
  return FiniteElementData<1> (fe_data.dofs_per_vertex * N,
			       fe_data.dofs_per_line * N,
			       fe_data.n_transform_functions,
			       fe_data.n_components * N);
};

#endif


#if deal_II_dimension == 2

template <>
FiniteElementData<2>
FESystem<2>::multiply_dof_numbers (const FiniteElementData<2> &fe_data,
				   const unsigned int          N) {
  return FiniteElementData<2> (fe_data.dofs_per_vertex * N,
			       fe_data.dofs_per_line * N,
			       fe_data.dofs_per_quad * N,
			       fe_data.n_transform_functions,
			       fe_data.n_components * N);
};

#endif



template <int dim>
double FESystem<dim>::shape_value (const unsigned int i,
				   const Point<dim>  &p) const {
  Assert((i<total_dofs), ExcInvalidIndex(i));

  return base_element->shape_value (i / n_sub_elements, p);
};



template <int dim>
Tensor<1,dim>
FESystem<dim>::shape_grad (const unsigned int  i,
			   const Point<dim>   &p) const {
  Assert((i<total_dofs), ExcInvalidIndex(i));

  return base_element->shape_grad (i / n_sub_elements, p);
};



template <int dim>
Tensor<2,dim>
FESystem<dim>::shape_grad_grad (const unsigned int  i,
				const Point<dim>   &p) const {
  Assert((i<total_dofs), ExcInvalidIndex(i));

  return base_element->shape_grad_grad (i / n_sub_elements, p);
};



template <int dim>
void FESystem<dim>::get_unit_support_points (vector<Point<dim> > &support_points) const {
  Assert (support_points.size() == total_dofs,
	  ExcWrongFieldDimension (support_points.size(), total_dofs));

  vector<Point<dim> > base_support_points (base_element->total_dofs);
  base_element->get_unit_support_points (base_support_points);
  
  for (unsigned int i=0; i<base_element->total_dofs; ++i) 
    for (unsigned int n=0; n<n_sub_elements; ++n)
      support_points[i*n_sub_elements+n] = base_support_points[i];
};



template <int dim>
void FESystem<dim>::get_support_points (const DoFHandler<dim>::cell_iterator &cell,
					const Boundary<dim> &boundary,
					vector<Point<dim> > &support_points) const {
  Assert (support_points.size() == total_dofs,
	  ExcWrongFieldDimension (support_points.size(), total_dofs));

  vector<Point<dim> > base_support_points (base_element->total_dofs);
  base_element->get_support_points (cell, boundary, base_support_points);
  
  for (unsigned int i=0; i<base_element->total_dofs; ++i) 
    for (unsigned int n=0; n<n_sub_elements; ++n)
      support_points[i*n_sub_elements+n] = base_support_points[i];
};



template <int dim>
void FESystem<dim>::get_face_support_points (const DoFHandler<dim>::face_iterator &face,
					     const Boundary<dim> &boundary,
					     vector<Point<dim> > &support_points) const {
  Assert (support_points.size() == dofs_per_face,
	  ExcWrongFieldDimension (support_points.size(), dofs_per_face));

  vector<Point<dim> > base_support_points (base_element->dofs_per_face);
  base_element->get_face_support_points (face, boundary, base_support_points);
  
  for (unsigned int i=0; i<base_element->dofs_per_face; ++i) 
    for (unsigned int n=0; n<n_sub_elements; ++n)
      support_points[i*n_sub_elements+n] = base_support_points[i];
};



template <int dim>
void FESystem<dim>::get_local_mass_matrix (const DoFHandler<dim>::cell_iterator &cell,
					   const Boundary<dim> &boundary,
					   dFMatrix            &local_mass_matrix) const {
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
};



template <int dim>
double FESystem<dim>::shape_value_transform (const unsigned int i,
					     const Point<dim>  &p) const {
  return base_element->shape_value_transform (i, p);
};



template <int dim>
Tensor<1,dim> FESystem<dim>::shape_grad_transform (const unsigned int i,
						   const Point<dim>  &p) const {
  return base_element->shape_grad_transform (i, p);
};



template <int dim>
void FESystem<dim>::get_face_jacobians (const DoFHandler<dim>::face_iterator &face,
					const Boundary<dim>         &boundary,
					const vector<Point<dim-1> > &unit_points,
					vector<double>      &face_jacobi_determinants) const {
  base_element->get_face_jacobians (face, boundary, unit_points, face_jacobi_determinants);
};



template <int dim>
void FESystem<dim>::get_subface_jacobians (const DoFHandler<dim>::face_iterator &face,
					   const unsigned int           subface_no,
					   const vector<Point<dim-1> > &unit_points,
					   vector<double>      &face_jacobi_determinants) const {
  base_element->get_subface_jacobians (face, subface_no, unit_points, face_jacobi_determinants);
};




template <int dim>
void FESystem<dim>::get_normal_vectors (const DoFHandler<dim>::cell_iterator &cell,
					const unsigned int          face_no,
					const Boundary<dim>         &boundary,
					const vector<Point<dim-1> > &unit_points,
					vector<Point<dim> >         &normal_vectors) const {
  base_element->get_normal_vectors (cell, face_no, boundary, unit_points, normal_vectors);
};



template <int dim>
void FESystem<dim>::get_normal_vectors (const DoFHandler<dim>::cell_iterator &cell,
					const unsigned int          face_no,
					const unsigned int          subface_no,
					const vector<Point<dim-1> > &unit_points,
					vector<Point<dim> >         &normal_vectors) const {
  base_element->get_normal_vectors (cell, face_no, subface_no, unit_points, normal_vectors);
};





// explicit instantiations
template class FESystem<deal_II_dimension>;



