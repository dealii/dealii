/* $Id$ */
/* Copyright W. Bangerth, University of Heidelberg, 1998 */


#include <fe/fe.h>
#include <fe/fe_values.h>
#include <base/quadrature.h>
#include <grid/tria_iterator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_boundary.h>
#include <grid/dof_accessor.h>
#include <lac/dvector.h>


/*------------------------------- FEValuesBase ---------------------------*/


template <int dim>
FEValuesBase<dim>::FEValuesBase (const unsigned int n_q_points,
				 const unsigned int n_support_points,
				 const unsigned int n_dofs,
				 const unsigned int n_transform_functions,
				 const unsigned int n_values_arrays,
				 const UpdateFlags update_flags,
				 const FiniteElement<dim> &fe) :
		n_quadrature_points (n_q_points),
		total_dofs (n_dofs),
		n_transform_functions (n_transform_functions),
		shape_values (n_values_arrays, dFMatrix(n_dofs, n_q_points)),
		shape_gradients (n_dofs, vector<Tensor<1,dim> >(n_q_points)),
		shape_2nd_derivatives (n_dofs, vector<Tensor<2,dim> >(n_q_points)),
		weights (n_q_points, 0),
		JxW_values (n_q_points, 0),
		quadrature_points (n_q_points, Point<dim>()),
		support_points (n_support_points, Point<dim>()),
		jacobi_matrices (n_q_points, Tensor<2,dim>()),
		jacobi_matrices_grad (n_q_points, Tensor<3,dim>()),
		shape_values_transform (n_values_arrays,
					dFMatrix(n_transform_functions,
						 n_quadrature_points)),
		selected_dataset (0),
		update_flags (update_flags),
		fe(&fe)        {};




template <int dim>
double FEValuesBase<dim>::shape_value (const unsigned int i,
				       const unsigned int j) const {
  Assert (selected_dataset<shape_values.size(),
	  ExcInvalidIndex (selected_dataset, shape_values.size()));
  Assert (i<shape_values[selected_dataset].m(),
	  ExcInvalidIndex (i, shape_values[selected_dataset].m()));
  Assert (j<shape_values[selected_dataset].n(),
	  ExcInvalidIndex (j, shape_values[selected_dataset].n()));

  return shape_values[selected_dataset](i,j);
};



template <int dim>
void FEValuesBase<dim>::get_function_values (const dVector  &fe_function,
					     vector<double> &values) const
{
  Assert (fe.n_components == 1,
	  ExcWrongNoOfComponents());
  Assert (selected_dataset<shape_values.size(),
	  ExcInvalidIndex (selected_dataset, shape_values.size()));
  Assert (values.size() == n_quadrature_points,
	  ExcWrongVectorSize(values.size(), n_quadrature_points));

				   // get function values of dofs
				   // on this cell
  dVector dof_values (total_dofs);
  present_cell->get_dof_values (fe_function, dof_values);

				   // initialize with zero
  fill_n (values.begin(), n_quadrature_points, 0);

				   // add up contributions of trial
				   // functions
  for (unsigned int point=0; point<n_quadrature_points; ++point)
    for (unsigned int shape_func=0; shape_func<total_dofs; ++shape_func)
      values[point] += (dof_values(shape_func) *
			shape_values[selected_dataset](shape_func, point));
};



template <int dim>
const Tensor<1,dim> &
FEValuesBase<dim>::shape_grad (const unsigned int i,
			       const unsigned int j) const {
  Assert (i<shape_gradients.size(),
	  ExcInvalidIndex (i, shape_gradients.size()));
  Assert (j<shape_gradients[i].size(),
	  ExcInvalidIndex (j, shape_gradients[i].size()));
  Assert (update_flags & update_gradients, ExcAccessToUninitializedField());

  return shape_gradients[i][j];
};



template <int dim>
void FEValuesBase<dim>::get_function_grads (const dVector          &fe_function,
					    vector<Tensor<1,dim> > &gradients) const {
  Assert (gradients.size() == n_quadrature_points,
	  ExcWrongVectorSize(gradients.size(), n_quadrature_points));

				   // get function values of dofs
				   // on this cell
  dVector dof_values (total_dofs);
  present_cell->get_dof_values (fe_function, dof_values);

				   // initialize with zero
  fill_n (gradients.begin(), n_quadrature_points, Tensor<1,dim>());

				   // add up contributions of trial
				   // functions
  for (unsigned int point=0; point<n_quadrature_points; ++point)
    for (unsigned int shape_func=0; shape_func<total_dofs; ++shape_func)
      {
	Tensor<1,dim> tmp(shape_gradients[shape_func][point]);
	tmp *= dof_values(shape_func);
	gradients[point] += tmp;
      };
};



template <int dim>
const Tensor<2,dim> &
FEValuesBase<dim>::shape_2nd_derivative (const unsigned int i,
					 const unsigned int j) const {
  Assert (i<shape_2nd_derivatives.size(),
	  ExcInvalidIndex (i, shape_2nd_derivatives.size()));
  Assert (j<shape_2nd_derivatives[i].size(),
	  ExcInvalidIndex (j, shape_2nd_derivatives[i].size()));
  Assert (update_flags & update_second_derivatives, ExcAccessToUninitializedField());

  return shape_2nd_derivatives[i][j];
};



template <int dim>
void FEValuesBase<dim>::get_function_2nd_derivatives (const dVector &fe_function,
						      vector<Tensor<2,dim> > &second_derivatives) const {
  Assert (second_derivatives.size() == n_quadrature_points,
	  ExcWrongVectorSize(second_derivatives.size(), n_quadrature_points));

				   // get function values of dofs
				   // on this cell
  dVector dof_values (total_dofs);
  present_cell->get_dof_values (fe_function, dof_values);

				   // initialize with zero
  fill_n (second_derivatives.begin(), n_quadrature_points, Tensor<2,dim>());

				   // add up contributions of trial
				   // functions
  for (unsigned int point=0; point<n_quadrature_points; ++point)
    for (unsigned int shape_func=0; shape_func<total_dofs; ++shape_func)
      {
	Tensor<2,dim> tmp(shape_2nd_derivatives[shape_func][point]);
	tmp *= dof_values(shape_func);
	second_derivatives[point] += tmp;
      };
};



template <int dim>
const Point<dim> &
FEValuesBase<dim>::quadrature_point (const unsigned int i) const {
  Assert (i<quadrature_points.size(), ExcInvalidIndex(i,quadrature_points.size()));
  Assert (update_flags & update_q_points, ExcAccessToUninitializedField());
  
  return quadrature_points[i];
};



template <int dim>
const Point<dim> &
FEValuesBase<dim>::support_point (const unsigned int i) const {
  Assert (i<support_points.size(), ExcInvalidIndex(i, support_points.size()));
  Assert (update_flags & update_support_points, ExcAccessToUninitializedField());
  
  return support_points[i];
};



template <int dim>
double FEValuesBase<dim>::JxW (const unsigned int i) const {
  Assert (i<JxW_values.size(), ExcInvalidIndex(i, JxW_values.size()));
  Assert (update_flags & update_JxW_values, ExcAccessToUninitializedField());
  
  return JxW_values[i];
};




/*------------------------------- FEValues -------------------------------*/

template <int dim>
FEValues<dim>::FEValues (const FiniteElement<dim> &fe,
			 const Quadrature<dim>    &quadrature,
			 const UpdateFlags         update_flags) :
		FEValuesBase<dim> (quadrature.n_quadrature_points,
				   fe.total_dofs,
				   fe.total_dofs,
				   fe.n_transform_functions,
				   1,
				   update_flags,
				   fe),
		unit_shape_gradients(fe.total_dofs,
				     vector<Tensor<1,dim> >(quadrature.n_quadrature_points)),
		unit_shape_2nd_derivatives(fe.total_dofs,
					   vector<Tensor<2,dim> >(quadrature.n_quadrature_points)),
		unit_shape_gradients_transform(fe.n_transform_functions,
					       vector<Tensor<1,dim> >(quadrature.n_quadrature_points)),
		unit_quadrature_points(quadrature.get_quad_points())
{
  Assert ((update_flags & update_normal_vectors) == false,
	  ExcInvalidUpdateFlag());

  for (unsigned int i=0; i<fe.total_dofs; ++i)
    for (unsigned int j=0; j<n_quadrature_points; ++j) 
      {
	shape_values[0](i,j) = fe.shape_value(i, unit_quadrature_points[j]);
	unit_shape_gradients[i][j]
	  = fe.shape_grad(i, unit_quadrature_points[j]);
	unit_shape_2nd_derivatives[i][j]
	  = fe.shape_grad_grad(i, unit_quadrature_points[j]);
      };

  for (unsigned int i=0; i<n_transform_functions; ++i)
    for (unsigned int j=0; j<n_quadrature_points; ++j)
      {
	shape_values_transform[0] (i,j)
	  = fe.shape_value_transform (i, unit_quadrature_points[j]);
	unit_shape_gradients_transform[i][j]
	  = fe.shape_grad_transform(i, unit_quadrature_points[j]);
      };
  
  weights = quadrature.get_weights ();
};







template <int dim>
void FEValues<dim>::reinit (const typename DoFHandler<dim>::cell_iterator &cell,
			    const Boundary<dim>                           &boundary) {
  present_cell = cell;
				   // fill jacobi matrices and real
				   // quadrature points
  if ((update_flags & update_jacobians)          ||
      (update_flags & update_JxW_values)         ||
      (update_flags & update_q_points)           ||
      (update_flags & update_gradients)          ||
      (update_flags & update_second_derivatives) ||
      (update_flags & update_support_points))
    fe->fill_fe_values (cell,
			unit_quadrature_points,
			jacobi_matrices,
			update_flags & (update_jacobians  |
					update_JxW_values |
					update_gradients  |
					update_second_derivatives),
			jacobi_matrices_grad,
			update_flags & update_second_derivatives,
			support_points,
			update_flags & update_support_points,
			quadrature_points,
			update_flags & update_q_points,
			shape_values_transform[0], unit_shape_gradients_transform,
			boundary);
  
				   // compute gradients on real element if
				   // requested
  if (update_flags & update_gradients) 
    for (unsigned int i=0; i<fe->total_dofs; ++i)
      for (unsigned int j=0; j<n_quadrature_points; ++j)
	for (unsigned int s=0; s<dim; ++s)
	  {
	    shape_gradients[i][j][s] = 0;
	    
					     // (grad psi)_s =
					     // (grad_{\xi\eta})_b J_{bs}
					     // with J_{bs}=(d\xi_b)/(dx_s)
	    for (unsigned int b=0; b<dim; ++b)
	      shape_gradients[i][j][s]
		+=
		unit_shape_gradients[i][j][b] * jacobi_matrices[j][b][s];
	  };

  
  if (update_flags & update_second_derivatives)
    for (unsigned int i=0; i<fe->total_dofs; ++i)
      for (unsigned int j=0; j<n_quadrature_points; ++j)
	{
	  static Tensor<2,dim> tmp1, tmp2;

					   // tmp1 := (d_k d_l phi) J_lj
	  contract (tmp1, unit_shape_2nd_derivatives[i][j], jacobi_matrices[j]);
					   // tmp2 := tmp1_kj J_ki
	  contract (tmp2, tmp1, 1, jacobi_matrices[j], 1);


					   // second part:
					   // tmp1 := (d_k J_lj) (d_l phi)
	  contract (tmp1, jacobi_matrices_grad[j], 2, unit_shape_gradients[i][j]);
					   // tmp1_kj J_ki
	  contract (shape_2nd_derivatives[i][j],
		    jacobi_matrices[j], 1,
		    tmp1, 1);

					   // add up first contribution
	  shape_2nd_derivatives[i][j] += tmp2;
	};
  
  
				   // compute Jacobi determinants in
				   // quadrature points.
				   // refer to the general doc for
				   // why we take the inverse of the
				   // determinant
  if (update_flags & update_JxW_values) 
    for (unsigned int i=0; i<n_quadrature_points; ++i)
      JxW_values[i] = weights[i] / determinant(jacobi_matrices[i]);
};





/*------------------------------- FEFaceValuesBase --------------------------*/


template <int dim>
FEFaceValuesBase<dim>::FEFaceValuesBase (const unsigned int n_q_points,
					 const unsigned int n_support_points,
					 const unsigned int n_dofs,
					 const unsigned int n_transform_functions,
					 const unsigned int n_faces_or_subfaces,
					 const UpdateFlags         update_flags,
					 const FiniteElement<dim> &fe) :
		FEValuesBase<dim> (n_q_points,
				   n_support_points,
				   n_dofs,
				   n_transform_functions,
				   n_faces_or_subfaces,
				   update_flags,
				   fe),
		unit_shape_gradients (n_faces_or_subfaces,
				      vector<vector<Tensor<1,dim> > >(n_dofs,
								   vector<Tensor<1,dim> >(n_q_points))),
		unit_shape_2nd_derivatives(n_faces_or_subfaces,
					   vector<vector<Tensor<2,dim> > >(n_dofs,
									   vector<Tensor<2,dim> >(n_q_points))),
		unit_shape_gradients_transform (n_faces_or_subfaces,
						vector<vector<Tensor<1,dim> > >(n_transform_functions,
										vector<Tensor<1,dim> >(n_q_points))),
		unit_face_quadrature_points (n_q_points, Point<dim-1>()),
		unit_quadrature_points (n_faces_or_subfaces,
					vector<Point<dim> >(n_q_points, Point<dim>())),
		face_jacobi_determinants (n_q_points, 0),
		normal_vectors (n_q_points)
{};



template <int dim>
const Point<dim> &
FEFaceValuesBase<dim>::normal_vector (const unsigned int i) const {
  Assert (i<normal_vectors.size(), ExcInvalidIndex(i, normal_vectors.size()));
  Assert (update_flags & update_normal_vectors,
	  ExcAccessToUninitializedField());
  
  return normal_vectors[i];
};






/*------------------------------- FEFaceValues -------------------------------*/


template <int dim>
FEFaceValues<dim>::FEFaceValues (const FiniteElement<dim> &fe,
				 const Quadrature<dim-1>  &quadrature,
				 const UpdateFlags         update_flags) :
		FEFaceValuesBase<dim> (quadrature.n_quadrature_points,
				       fe.dofs_per_face,
				       fe.total_dofs,
				       fe.n_transform_functions,
				       GeometryInfo<dim>::faces_per_cell,
				       update_flags,
				       fe)
{
  unit_face_quadrature_points = quadrature.get_quad_points();
  weights = quadrature.get_weights ();  

  				   // set up an array of the unit points
				   // on the given face, but in coordinates
				   // of the space with #dim# dimensions.
				   // the points are still on the unit
				   // cell, not on the real cell.
  for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
    QProjector<dim>::project_to_face (quadrature, face, unit_quadrature_points[face]);

  for (unsigned int i=0; i<fe.total_dofs; ++i)
    for (unsigned int j=0; j<n_quadrature_points; ++j) 
      for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
	{
	  shape_values[face](i,j)
	    = fe.shape_value(i, unit_quadrature_points[face][j]);
	  unit_shape_gradients[face][i][j]
	    = fe.shape_grad(i, unit_quadrature_points[face][j]);
	  unit_shape_2nd_derivatives[face][i][j]
	    = fe.shape_grad_grad(i, unit_quadrature_points[face][j]);
	};

  for (unsigned int i=0; i<n_transform_functions; ++i)
    for (unsigned int j=0; j<n_quadrature_points; ++j)
      for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
	{
	  shape_values_transform[face] (i,j)
	    = fe.shape_value_transform (i, unit_quadrature_points[face][j]);
	  unit_shape_gradients_transform[face][i][j]
	    = fe.shape_grad_transform(i, unit_quadrature_points[face][j]);
	};
};



template <int dim>
void FEFaceValues<dim>::reinit (const typename DoFHandler<dim>::cell_iterator &cell,
				const unsigned int                             face_no,
				const Boundary<dim>                           &boundary) {
  present_cell  = cell;
  selected_dataset = face_no;
				   // fill jacobi matrices and real
				   // quadrature points
  if ((update_flags & update_jacobians)          ||
      (update_flags & update_JxW_values)         ||
      (update_flags & update_q_points)           ||
      (update_flags & update_gradients)          ||
      (update_flags & update_second_derivatives) ||
      (update_flags & update_support_points)     ||
      (update_flags & update_normal_vectors)     ||
      (update_flags & update_JxW_values))
    fe->fill_fe_face_values (cell,
			     face_no,
			     unit_face_quadrature_points,
			     unit_quadrature_points[face_no],
			     jacobi_matrices,
			     update_flags & (update_jacobians |
					     update_gradients |
					     update_JxW_values |
					     update_second_derivatives),
			     jacobi_matrices_grad,
			     update_flags & update_second_derivatives,
			     support_points,
			     update_flags & update_support_points,
			     quadrature_points,
			     update_flags & update_q_points,
			     face_jacobi_determinants,
			     update_flags & update_JxW_values,
			     normal_vectors,
			     update_flags & update_normal_vectors,
			     shape_values_transform[face_no],
			     unit_shape_gradients_transform[face_no],
			     boundary);

				   // compute gradients on real element if
				   // requested
  if (update_flags & update_gradients) 
    for (unsigned int i=0; i<fe->total_dofs; ++i)
      {
	fill_n (shape_gradients[i].begin(),
		n_quadrature_points,
		Tensor<1,dim>());
	for (unsigned int j=0; j<n_quadrature_points; ++j) 
	  for (unsigned int s=0; s<dim; ++s)
					     // (grad psi)_s =
					     // (grad_{\xi\eta})_b J_{bs}
					     // with J_{bs}=(d\xi_b)/(dx_s)
	    for (unsigned int b=0; b<dim; ++b)
	      shape_gradients[i][j][s]
		+= (unit_shape_gradients[face_no][i][j][b] *
		    jacobi_matrices[j][b][s]);
      };


  if (update_flags & update_second_derivatives)
    for (unsigned int i=0; i<fe->total_dofs; ++i)
      for (unsigned int j=0; j<n_quadrature_points; ++j)
	{
	  static Tensor<2,dim> tmp1, tmp2;

					   // tmp1 := (d_k d_l phi) J_lj
	  contract (tmp1, unit_shape_2nd_derivatives[face_no][i][j], jacobi_matrices[j]);
					   // tmp2 := tmp1_kj J_ki
	  contract (tmp2, tmp1, 1, jacobi_matrices[j], 1);


					   // second part:
					   // tmp1 := (d_k J_lj) (d_l phi)
	  contract (tmp1,
		    jacobi_matrices_grad[j], 2,
		    unit_shape_gradients[face_no][i][j]);
					   // tmp1_kj J_ki
	  contract (shape_2nd_derivatives[i][j],
		    jacobi_matrices[j], 1,
		    tmp1, 1);

					   // add up first contribution
	  shape_2nd_derivatives[i][j] += tmp2;
	};
  
  
				   // compute Jacobi determinants in
				   // quadrature points.
				   // refer to the general doc for
				   // why we take the inverse of the
				   // determinant
  if (update_flags & update_JxW_values) 
    for (unsigned int i=0; i<n_quadrature_points; ++i)
      JxW_values[i] = weights[i] * face_jacobi_determinants[i];
};






/*------------------------------- FESubFaceValues -------------------------------*/


template <int dim>
FESubfaceValues<dim>::FESubfaceValues (const FiniteElement<dim> &fe,
				       const Quadrature<dim-1>  &quadrature,
				       const UpdateFlags         update_flags) :
		FEFaceValuesBase<dim> (quadrature.n_quadrature_points,
				       0,
				       fe.total_dofs,
				       fe.n_transform_functions,
				       GeometryInfo<dim>::faces_per_cell * GeometryInfo<dim>::subfaces_per_face,
				       update_flags,
				       fe)
{
  Assert ((update_flags & update_support_points) == false,
	  ExcInvalidUpdateFlag());
  
  unit_face_quadrature_points = quadrature.get_quad_points();
  weights = quadrature.get_weights ();  

  				   // set up an array of the unit points
				   // on the given face, but in coordinates
				   // of the space with #dim# dimensions.
				   // the points are still on the unit
				   // cell, not on the real cell.
  for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
    for (unsigned int subface=0; subface<GeometryInfo<dim>::subfaces_per_face; ++subface)
      QProjector<dim>::project_to_subface (quadrature,
					   face, subface,
					   unit_quadrature_points[face*(1<<(dim-1))+subface]);

  for (unsigned int i=0; i<fe.total_dofs; ++i)
    for (unsigned int j=0; j<n_quadrature_points; ++j) 
      for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
	for (unsigned int subface=0; subface<GeometryInfo<dim>::subfaces_per_face; ++subface)
	  {
	    shape_values[face*GeometryInfo<dim>::subfaces_per_face+subface](i,j)
	      = fe.shape_value(i, unit_quadrature_points[face *
							GeometryInfo<dim>::
							subfaces_per_face+subface][j]);
	    unit_shape_gradients[face*GeometryInfo<dim>::subfaces_per_face+subface][i][j]
	      = fe.shape_grad(i, unit_quadrature_points[face *
						       GeometryInfo<dim>::
						       subfaces_per_face+subface][j]);
	    unit_shape_2nd_derivatives[face*GeometryInfo<dim>::subfaces_per_face+subface][i][j]
	      = fe.shape_grad_grad(i, unit_quadrature_points[face *
							    GeometryInfo<dim>::
							    subfaces_per_face+subface][j]);
	  };
  for (unsigned int i=0; i<n_transform_functions; ++i)
    for (unsigned int j=0; j<n_quadrature_points; ++j)
      for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
	for (unsigned int subface=0; subface<GeometryInfo<dim>::subfaces_per_face; ++subface)
	  {
	    shape_values_transform[face*GeometryInfo<dim>::subfaces_per_face+subface] (i,j)
	      = fe.shape_value_transform (i, unit_quadrature_points[face *
								   GeometryInfo<dim>::
								   subfaces_per_face +
								   subface][j]);
	    unit_shape_gradients_transform[face *
					  GeometryInfo<dim>::subfaces_per_face +
					  subface][i][j]
	      = fe.shape_grad_transform(i, unit_quadrature_points[face *
								 GeometryInfo<dim>::
								 subfaces_per_face +
								 subface][j]);
	  };
};



template <int dim>
void FESubfaceValues<dim>::reinit (const typename DoFHandler<dim>::cell_iterator &cell,
				   const unsigned int         face_no,
				   const unsigned int         subface_no,
				   const Boundary<dim>       &boundary) {
  Assert (cell->face(face_no)->at_boundary() == false,
	  ExcReinitCalledWithBoundaryFace());
  
  present_cell  = cell;
  selected_dataset = face_no*(1<<(dim-1)) + subface_no;
				   // fill jacobi matrices and real
				   // quadrature points
  if ((update_flags & update_jacobians)          ||
      (update_flags & update_JxW_values)         ||
      (update_flags & update_q_points)           ||
      (update_flags & update_gradients)          ||
      (update_flags & update_second_derivatives) ||
      (update_flags & update_normal_vectors)     ||
      (update_flags & update_JxW_values))
    fe->fill_fe_subface_values (cell,
				face_no,
				subface_no,
				unit_face_quadrature_points,
				unit_quadrature_points[selected_dataset],
				jacobi_matrices,
				update_flags & (update_jacobians |
						update_gradients |
						update_JxW_values|
						update_second_derivatives),
				jacobi_matrices_grad,
				update_flags & update_second_derivatives,
				quadrature_points,
				update_flags & update_q_points,
				face_jacobi_determinants,
				update_flags & update_JxW_values,
				normal_vectors,
				update_flags & update_normal_vectors,
				shape_values_transform[selected_dataset],
				unit_shape_gradients_transform[selected_dataset],
				boundary);

				   // compute gradients on real element if
				   // requested
  if (update_flags & update_gradients) 
    for (unsigned int i=0; i<fe->total_dofs; ++i) 
      {
	fill_n (shape_gradients[i].begin(),
		n_quadrature_points,
		Tensor<1,dim>());
	for (unsigned int j=0; j<n_quadrature_points; ++j) 
	  for (unsigned int s=0; s<dim; ++s)
					     // (grad psi)_s =
					     // (grad_{\xi\eta})_b J_{bs}
					     // with J_{bs}=(d\xi_b)/(dx_s)
	    for (unsigned int b=0; b<dim; ++b)
	      shape_gradients[i][j][s]
		+= (unit_shape_gradients[selected_dataset][i][j][b] *
		    jacobi_matrices[j][b][s]);
      };

  if (update_flags & update_second_derivatives)
    for (unsigned int i=0; i<fe->total_dofs; ++i)
      for (unsigned int j=0; j<n_quadrature_points; ++j)
	{
	  static Tensor<2,dim> tmp1, tmp2;

					   // tmp1 := (d_k d_l phi) J_lj
	  contract (tmp1,
		    unit_shape_2nd_derivatives[selected_dataset][i][j],
		    jacobi_matrices[j]);
					   // tmp2 := tmp1_kj J_ki
	  contract (tmp2, tmp1, 1, jacobi_matrices[j], 1);


					   // second part:
					   // tmp1 := (d_k J_lj) (d_l phi)
	  contract (tmp1,
		    jacobi_matrices_grad[j], 2,
		    unit_shape_gradients[selected_dataset][i][j]);
					   // tmp1_kj J_ki
	  contract (shape_2nd_derivatives[i][j],
		    jacobi_matrices[j], 1,
		    tmp1, 1);

					   // add up first contribution
	  shape_2nd_derivatives[i][j] += tmp2;
	};
  
  
				   // compute Jacobi determinants in
				   // quadrature points.
				   // refer to the general doc for
				   // why we take the inverse of the
				   // determinant
  if (update_flags & update_JxW_values) 
    for (unsigned int i=0; i<n_quadrature_points; ++i)
      JxW_values[i] = weights[i] * face_jacobi_determinants[i];
};





/*------------------------------- Explicit Instantiations -------------*/

template class FEValuesBase<deal_II_dimension>;
template class FEValues<deal_II_dimension>;

#if deal_II_dimension >= 2
template class FEFaceValuesBase<deal_II_dimension>;
template class FEFaceValues<deal_II_dimension>;
template class FESubfaceValues<deal_II_dimension>;
#endif
