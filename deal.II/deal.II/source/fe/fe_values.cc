//----------------------------  fe_values.cc  ---------------------------
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
//----------------------------  fe_values.cc  ---------------------------


#include <fe/fe.h>
#include <fe/fe_values.h>
#include <base/quadrature.h>
#include <grid/tria_iterator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_boundary.h>
#include <dofs/dof_accessor.h>
#include <lac/vector.h>
#include <lac/block_vector.h>


/*------------------------------- FEValuesBase ---------------------------*/


template <int dim>
FEValuesBase<dim>::FEValuesBase (const unsigned int n_q_points,
				 const unsigned int n_support_points,
				 const unsigned int dofs_per_cell,
				 const unsigned int n_transform_functions,
				 const unsigned int n_values_arrays,
				 const UpdateFlags update_flags,
				 const FiniteElement<dim> &fe)
		:
		n_quadrature_points (n_q_points),
		dofs_per_cell (dofs_per_cell),
		n_transform_functions (n_transform_functions),
		shape_values (n_values_arrays, FullMatrix<double>(dofs_per_cell, n_q_points)),
		shape_gradients (dofs_per_cell, vector<Tensor<1,dim> >(n_q_points)),
		shape_2nd_derivatives (dofs_per_cell, vector<Tensor<2,dim> >(n_q_points)),
		weights (n_q_points, 0),
		JxW_values (n_q_points, 0),
		quadrature_points (n_q_points, Point<dim>()),
		support_points (n_support_points, Point<dim>()),
		jacobi_matrices (n_q_points, Tensor<2,dim>()),
		jacobi_matrices_grad (n_q_points, Tensor<3,dim>()),
		shape_values_transform (n_values_arrays,
					FullMatrix<double>(n_transform_functions,
							   n_quadrature_points)),
		selected_dataset (0),
		update_flags (update_flags),
		fe(&fe)
{};



template <int dim>
double FEValuesBase<dim>::shape_value (const unsigned int i,
				       const unsigned int j) const
{
  Assert (update_flags & update_values, ExcAccessToUninitializedField());
  Assert (selected_dataset<shape_values.size(),
	  ExcIndexRange (selected_dataset, 0, shape_values.size()));
  Assert (i<shape_values[selected_dataset].m(),
	  ExcIndexRange (i, 0, shape_values[selected_dataset].m()));
  Assert (j<shape_values[selected_dataset].n(),
	  ExcIndexRange (j, 0, shape_values[selected_dataset].n()));

  return shape_values[selected_dataset](i,j);
};



template <int dim>
template <class InputVector, typename number>
void FEValuesBase<dim>::get_function_values (const InputVector &fe_function,
					     vector<number>    &values) const
{
  Assert (fe->n_components() == 1,
	  ExcWrongNoOfComponents());
  Assert (selected_dataset<shape_values.size(),
	  ExcIndexRange (selected_dataset, 0, shape_values.size()));
  Assert (values.size() == n_quadrature_points,
	  ExcWrongVectorSize(values.size(), n_quadrature_points));

				   // get function values of dofs
				   // on this cell
  Vector<typename InputVector::value_type> dof_values (dofs_per_cell);
  if (present_cell->active())
    present_cell->get_dof_values (fe_function, dof_values);
  else
    present_cell->get_interpolated_dof_values(fe_function, dof_values);

				   // initialize with zero
  fill_n (values.begin(), n_quadrature_points, 0);

				   // add up contributions of trial
				   // functions
  for (unsigned int point=0; point<n_quadrature_points; ++point)
    for (unsigned int shape_func=0; shape_func<dofs_per_cell; ++shape_func)
      values[point] += (dof_values(shape_func) *
			shape_values[selected_dataset](shape_func, point));
};



template <int dim>
template <class InputVector, typename number>
void FEValuesBase<dim>::get_function_values (const InputVector       &fe_function,
					     vector<Vector<number> > &values) const
{
  Assert (n_quadrature_points == values.size(),
	  ExcWrongVectorSize(values.size(), n_quadrature_points));
  Assert (selected_dataset<shape_values.size(),
	  ExcIndexRange (selected_dataset, 0, shape_values.size()));
  for (unsigned i=0;i<values.size();++i)
    Assert (values[i].size() == fe->n_components(),
	    ExcWrongNoOfComponents());

				   // get function values of dofs
				   // on this cell
  Vector<typename InputVector::value_type> dof_values (dofs_per_cell);
  if (present_cell->active())
    present_cell->get_dof_values (fe_function, dof_values);
  else
    present_cell->get_interpolated_dof_values(fe_function, dof_values);
  
				   // initialize with zero
  for (unsigned i=0;i<values.size();++i)
    fill_n (values[i].begin(), values[i].size(), 0);

				   // add up contributions of trial
				   // functions
  for (unsigned int point=0; point<n_quadrature_points; ++point)
    for (unsigned int shape_func=0; shape_func<dofs_per_cell; ++shape_func)
      values[point](fe->system_to_component_index(shape_func).first)
	+= (dof_values(shape_func) * shape_values[selected_dataset](shape_func, point));
};



template <int dim>
const Tensor<1,dim> &
FEValuesBase<dim>::shape_grad (const unsigned int i,
			       const unsigned int j) const
{
  Assert (i<shape_gradients.size(),
	  ExcIndexRange (i, 0, shape_gradients.size()));
  Assert (j<shape_gradients[i].size(),
	  ExcIndexRange (j, 0, shape_gradients[i].size()));
  Assert (update_flags & update_gradients, ExcAccessToUninitializedField());

  return shape_gradients[i][j];
};



template <int dim>
template <class InputVector>
void FEValuesBase<dim>::get_function_grads (const InputVector      &fe_function,
					    vector<Tensor<1,dim> > &gradients) const
{
  Assert (fe->n_components() == 1,
	  ExcWrongNoOfComponents());
  Assert (gradients.size() == n_quadrature_points,
	  ExcWrongVectorSize(gradients.size(), n_quadrature_points));

				   // get function values of dofs
				   // on this cell
  Vector<typename InputVector::value_type> dof_values (dofs_per_cell);
  if (present_cell->active())
    present_cell->get_dof_values (fe_function, dof_values);
  else
    present_cell->get_interpolated_dof_values(fe_function, dof_values);

				   // initialize with zero
  fill_n (gradients.begin(), n_quadrature_points, Tensor<1,dim>());

				   // add up contributions of trial
				   // functions
  for (unsigned int point=0; point<n_quadrature_points; ++point)
    for (unsigned int shape_func=0; shape_func<dofs_per_cell; ++shape_func)
      {
	Tensor<1,dim> tmp(shape_gradients[shape_func][point]);
	tmp *= dof_values(shape_func);
	gradients[point] += tmp;
      };
};



template <int dim>
template <class InputVector>
void FEValuesBase<dim>::get_function_grads (const InputVector               &fe_function,
					    vector<vector<Tensor<1,dim> > > &gradients) const
{
  Assert (n_quadrature_points == gradients.size(),
	  ExcWrongNoOfComponents());
  Assert (selected_dataset<shape_values.size(),
	  ExcIndexRange (selected_dataset, 0, shape_values.size()));
  for (unsigned i=0;i<gradients.size();++i)
    Assert (gradients[i].size() == fe->n_components(),
	    ExcWrongVectorSize(gradients[i].size(), fe->n_components()));

				   // get function values of dofs
				   // on this cell
  Vector<typename InputVector::value_type> dof_values (dofs_per_cell);
  if (present_cell->active())
    present_cell->get_dof_values (fe_function, dof_values);
  else
    present_cell->get_interpolated_dof_values(fe_function, dof_values);

				   // initialize with zero
  for (unsigned i=0;i<gradients.size();++i)
    fill_n (gradients[i].begin(), gradients[i].size(), Tensor<1,dim>());

				   // add up contributions of trial
				   // functions
  for (unsigned int point=0; point<n_quadrature_points; ++point)
    for (unsigned int shape_func=0; shape_func<dofs_per_cell; ++shape_func)
      {
	Tensor<1,dim> tmp(shape_gradients[shape_func][point]);
	tmp *= dof_values(shape_func);
	gradients[point][fe->system_to_component_index(shape_func).first]
	  += tmp;
      };
};



template <int dim>
const Tensor<2,dim> &
FEValuesBase<dim>::shape_2nd_derivative (const unsigned int i,
					 const unsigned int j) const
{
  Assert (i<shape_2nd_derivatives.size(),
	  ExcIndexRange (i, 0, shape_2nd_derivatives.size()));
  Assert (j<shape_2nd_derivatives[i].size(),
	  ExcIndexRange (j, 0, shape_2nd_derivatives[i].size()));
  Assert (update_flags & update_second_derivatives, ExcAccessToUninitializedField());

  return shape_2nd_derivatives[i][j];
};



template <int dim>
template <class InputVector>
void FEValuesBase<dim>::get_function_2nd_derivatives (const InputVector      &fe_function,
						      vector<Tensor<2,dim> > &second_derivatives) const
{
  Assert (fe->n_components() == 1,
	  ExcWrongNoOfComponents());
  Assert (second_derivatives.size() == n_quadrature_points,
	  ExcWrongVectorSize(second_derivatives.size(), n_quadrature_points));

				   // get function values of dofs
				   // on this cell
  Vector<typename InputVector::value_type> dof_values (dofs_per_cell);
  if (present_cell->active())
    present_cell->get_dof_values (fe_function, dof_values);
  else
    present_cell->get_interpolated_dof_values(fe_function, dof_values);

				   // initialize with zero
  fill_n (second_derivatives.begin(), n_quadrature_points, Tensor<2,dim>());

				   // add up contributions of trial
				   // functions
  for (unsigned int point=0; point<n_quadrature_points; ++point)
    for (unsigned int shape_func=0; shape_func<dofs_per_cell; ++shape_func)
      {
	Tensor<2,dim> tmp(shape_2nd_derivatives[shape_func][point]);
	tmp *= dof_values(shape_func);
	second_derivatives[point] += tmp;
      };
};



template <int dim>
const Point<dim> &
FEValuesBase<dim>::quadrature_point (const unsigned int i) const
{
  Assert (i<quadrature_points.size(), ExcIndexRange(i, 0, quadrature_points.size()));
  Assert (update_flags & update_q_points, ExcAccessToUninitializedField());
  
  return quadrature_points[i];
};



template <int dim>
const Point<dim> &
FEValuesBase<dim>::support_point (const unsigned int i) const
{
  Assert (i<support_points.size(), ExcIndexRange(i, 0, support_points.size()));
  Assert (update_flags & update_support_points, ExcAccessToUninitializedField());
  
  return support_points[i];
};



template <int dim>
double FEValuesBase<dim>::JxW (const unsigned int i) const
{
  Assert (i<JxW_values.size(), ExcIndexRange(i, 0, JxW_values.size()));
  Assert (update_flags & update_JxW_values, ExcAccessToUninitializedField());
  
  return JxW_values[i];
};


/*------------------------------- FEValues -------------------------------*/

template <int dim>
FEValues<dim>::FEValues (const FiniteElement<dim> &fe,
			 const Quadrature<dim>    &quadrature,
			 const UpdateFlags         update_flags)
		:
		FEValuesBase<dim> (quadrature.n_quadrature_points,
				   fe.dofs_per_cell,
				   fe.dofs_per_cell,
				   fe.transform_functions,
				   1,
				   update_flags,
				   fe),
		unit_shape_gradients(fe.dofs_per_cell,
				     vector<Tensor<1,dim> >(quadrature.n_quadrature_points)),
		unit_shape_2nd_derivatives(fe.dofs_per_cell,
					   vector<Tensor<2,dim> >(quadrature.n_quadrature_points)),
		unit_shape_gradients_transform(fe.n_transform_functions(),
					       vector<Tensor<1,dim> >(quadrature.n_quadrature_points)),
		unit_quadrature_points(quadrature.get_points())
{
  Assert ((update_flags & update_normal_vectors) == false,
	  ExcInvalidUpdateFlag());

  for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
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
void FEValues<dim>::reinit (const typename DoFHandler<dim>::cell_iterator &cell)
{
  present_cell = cell;

				   // assert that the finite elements
				   // passed to the constructor and
				   // used by the DoFHandler used by
				   // this cell, are the same
  Assert (static_cast<const FiniteElementData<dim>&>(*fe)
	  ==
	  static_cast<const FiniteElementData<dim>&>(cell->get_dof_handler().get_fe()),
	  ExcFEDontMatch());
  
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
			shape_values_transform[0], unit_shape_gradients_transform);
  
				   // compute gradients on real element if
				   // requested
  if (update_flags & update_gradients) 
    for (unsigned int i=0; i<fe->dofs_per_cell; ++i)
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
  
  Tensor<2,dim> tmp1, tmp2;
  if (update_flags & update_second_derivatives)
    for (unsigned int i=0; i<fe->dofs_per_cell; ++i)
      for (unsigned int j=0; j<n_quadrature_points; ++j)
	{
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
					 const unsigned int dofs_per_cell,
					 const unsigned int n_transform_functions,
					 const unsigned int n_faces_or_subfaces,
					 const UpdateFlags         update_flags,
					 const FiniteElement<dim> &fe)
		:
		FEValuesBase<dim> (n_q_points,
				   n_support_points,
				   dofs_per_cell,
				   n_transform_functions,
				   n_faces_or_subfaces,
				   update_flags,
				   fe),
		unit_shape_gradients (n_faces_or_subfaces,
				      vector<vector<Tensor<1,dim> > >(dofs_per_cell,
								   vector<Tensor<1,dim> >(n_q_points))),
		unit_shape_2nd_derivatives(n_faces_or_subfaces,
					   vector<vector<Tensor<2,dim> > >(dofs_per_cell,
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
FEFaceValuesBase<dim>::normal_vector (const unsigned int i) const
{
  Assert (i<normal_vectors.size(), ExcIndexRange(i, 0, normal_vectors.size()));
  Assert (update_flags & update_normal_vectors,
	  ExcAccessToUninitializedField());
  
  return normal_vectors[i];
};


/*------------------------------- FEFaceValues -------------------------------*/


template <int dim>
FEFaceValues<dim>::FEFaceValues (const FiniteElement<dim> &fe,
				 const Quadrature<dim-1>  &quadrature,
				 const UpdateFlags         update_flags)
		:
		FEFaceValuesBase<dim> (quadrature.n_quadrature_points,
				       fe.dofs_per_face,
				       fe.dofs_per_cell,
				       fe.n_transform_functions(),
				       GeometryInfo<dim>::faces_per_cell,
				       update_flags,
				       fe)
{
  unit_face_quadrature_points = quadrature.get_points();
  weights = quadrature.get_weights ();  

  				   // set up an array of the unit points
				   // on the given face, but in coordinates
				   // of the space with #dim# dimensions.
				   // the points are still on the unit
				   // cell, not on the real cell.
  for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
    QProjector<dim>::project_to_face (quadrature, face, unit_quadrature_points[face]);

  for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
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
				const unsigned int                             face_no)
{
  present_cell  = cell;
  selected_dataset = face_no;

				   // assert that the finite elements
				   // passed to the constructor and
				   // used by the DoFHandler used by
				   // this cell, are the same
  Assert (static_cast<const FiniteElementData<dim>&>(*fe)
	  ==
	  static_cast<const FiniteElementData<dim>&>(cell->get_dof_handler().get_fe()),
	  ExcFEDontMatch());
  Assert (face_no < GeometryInfo<dim>::faces_per_cell,
	  ExcIndexRange (face_no, 0, GeometryInfo<dim>::faces_per_cell));
  
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
			     unit_shape_gradients_transform[face_no]);

				   // compute gradients on real element if
				   // requested
  if (update_flags & update_gradients) 
    for (unsigned int i=0; i<fe->dofs_per_cell; ++i)
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


  Tensor<2,dim> tmp1, tmp2;
  if (update_flags & update_second_derivatives)
    for (unsigned int i=0; i<fe->dofs_per_cell; ++i)
      for (unsigned int j=0; j<n_quadrature_points; ++j)
	{
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
				       const UpdateFlags         update_flags)
		:
		FEFaceValuesBase<dim> (quadrature.n_quadrature_points,
				       0,
				       fe.dofs_per_cell,
				       fe.n_transform_functions(),
				       GeometryInfo<dim>::faces_per_cell * GeometryInfo<dim>::subfaces_per_face,
				       update_flags,
				       fe)
{
  Assert ((update_flags & update_support_points) == false,
	  ExcInvalidUpdateFlag());
  
  unit_face_quadrature_points = quadrature.get_points();
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

  for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
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
				   const unsigned int         subface_no)
{
  Assert (cell->face(face_no)->at_boundary() == false,
	  ExcReinitCalledWithBoundaryFace());
  
  present_cell  = cell;
  selected_dataset = face_no*(1<<(dim-1)) + subface_no;

				   // assert that the finite elements
				   // passed to the constructor and
				   // used by the DoFHandler used by
				   // this cell, are the same
  Assert (static_cast<const FiniteElementData<dim>&>(*fe)
	  ==
	  static_cast<const FiniteElementData<dim>&>(cell->get_dof_handler().get_fe()),
	  ExcFEDontMatch());
  Assert (face_no < GeometryInfo<dim>::faces_per_cell,
	  ExcIndexRange (face_no, 0, GeometryInfo<dim>::faces_per_cell));
  Assert (subface_no < GeometryInfo<dim>::subfaces_per_face,
	  ExcIndexRange (subface_no, 0, GeometryInfo<dim>::subfaces_per_face));
    
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
				unit_shape_gradients_transform[selected_dataset]);

				   // compute gradients on real element if
				   // requested
  if (update_flags & update_gradients) 
    for (unsigned int i=0; i<fe->dofs_per_cell; ++i) 
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

  Tensor<2,dim> tmp1, tmp2;
  if (update_flags & update_second_derivatives)
    for (unsigned int i=0; i<fe->dofs_per_cell; ++i)
      for (unsigned int j=0; j<n_quadrature_points; ++j)
	{
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


//-----------------------------------------------------------------------------

template
void FEValuesBase<deal_II_dimension>::get_function_values (const Vector<double> &,
					     vector<double>       &) const;
template
void FEValuesBase<deal_II_dimension>::get_function_values (const Vector<float> &,
					     vector<float>      &) const;
template
void FEValuesBase<deal_II_dimension>::get_function_values (const BlockVector<2,double> &,
					     vector<double>      &) const;
template
void FEValuesBase<deal_II_dimension>::get_function_values (const BlockVector<3,double> &,
					     vector<double>      &) const;

//-----------------------------------------------------------------------------

template
void FEValuesBase<deal_II_dimension>::get_function_values (const Vector<double> &,
					     vector<Vector<double> > &) const;
template
void FEValuesBase<deal_II_dimension>::get_function_values (const Vector<float> &,
					     vector<Vector<float> > &) const;
template
void FEValuesBase<deal_II_dimension>::get_function_values (const BlockVector<2,double> &,
					     vector<Vector<double> >     &) const;
template
void FEValuesBase<deal_II_dimension>::get_function_values (const BlockVector<3,double> &,
					     vector<Vector<double> >     &) const;

//-----------------------------------------------------------------------------

template
void FEValuesBase<deal_II_dimension>::get_function_grads (const Vector<double> &,
					    vector<Tensor<1,deal_II_dimension> > &) const;
template
void FEValuesBase<deal_II_dimension>::get_function_grads (const Vector<float> &,
					     vector<Tensor<1,deal_II_dimension> > &) const;
template
void FEValuesBase<deal_II_dimension>::get_function_grads (const BlockVector<2,double> &,
					     vector<Tensor<1,deal_II_dimension> > &) const;
template
void FEValuesBase<deal_II_dimension>::get_function_grads (const BlockVector<3,double> &,
					     vector<Tensor<1,deal_II_dimension> > &) const;

//-----------------------------------------------------------------------------

template
void FEValuesBase<deal_II_dimension>::get_function_grads (const Vector<double> &,
					     vector<vector<Tensor<1,deal_II_dimension> > > &) const;
template
void FEValuesBase<deal_II_dimension>::get_function_grads (const Vector<float> &,
					     vector<vector<Tensor<1,deal_II_dimension> > > &) const;
template
void FEValuesBase<deal_II_dimension>::get_function_grads (const BlockVector<2,double> &,
					     vector<vector<Tensor<1,deal_II_dimension> > > &) const;
template
void FEValuesBase<deal_II_dimension>::get_function_grads (const BlockVector<3,double> &,
					     vector<vector<Tensor<1,deal_II_dimension> > > &) const;

//-----------------------------------------------------------------------------

template
void FEValuesBase<deal_II_dimension>::get_function_2nd_derivatives (const Vector<double> &,
					    vector<Tensor<2,deal_II_dimension> > &) const;
template
void FEValuesBase<deal_II_dimension>::get_function_2nd_derivatives (const Vector<float> &,
					     vector<Tensor<2,deal_II_dimension> > &) const;
template
void FEValuesBase<deal_II_dimension>::get_function_2nd_derivatives (const BlockVector<2,double> &,
					     vector<Tensor<2,deal_II_dimension> > &) const;
template
void FEValuesBase<deal_II_dimension>::get_function_2nd_derivatives (const BlockVector<3,double> &,
					     vector<Tensor<2,deal_II_dimension> > &) const;
