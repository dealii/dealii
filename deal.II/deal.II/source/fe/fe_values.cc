/* $Id$ */

#include <fe/fe.h>
#include <fe/fe_values.h>
#include <fe/quadrature.h>
#include <grid/tria_iterator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_boundary.h>
#include <grid/dof_accessor.h>



/*------------------------------- FEValuesBase ---------------------------*/


template <int dim>
FEValuesBase<dim>::FEValuesBase (const unsigned int n_q_points,
				 const unsigned int n_ansatz_points,
				 const unsigned int n_dofs,
				 const unsigned int n_values_arrays,
				 const UpdateFlags update_flags) :
		n_quadrature_points (n_q_points),
		total_dofs (n_dofs),
		shape_values (n_values_arrays, dFMatrix(n_dofs, n_q_points)),
		shape_gradients (n_dofs, vector<Point<dim> >(n_q_points)),
		weights (n_q_points, 0),
		JxW_values (n_q_points, 0),
		quadrature_points (n_q_points, Point<dim>()),
		ansatz_points (n_ansatz_points, Point<dim>()),
		jacobi_matrices (n_q_points, dFMatrix(dim,dim)),
		selected_dataset (0),
		update_flags (update_flags) {};




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
					     vector<double> &values) const {
  Assert (selected_dataset<shape_values.size(),
	  ExcInvalidIndex (selected_dataset, shape_values.size()));
  Assert (values.size() == n_quadrature_points,
	  ExcWrongVectorSize(values.size(), n_quadrature_points));

				   // get function values of dofs
				   // on this cell
  vector<double> dof_values (total_dofs, 0);
  present_cell->get_dof_values (fe_function, dof_values);

				   // initialize with zero
  fill_n (values.begin(), n_quadrature_points, 0);

				   // add up contributions of ansatz
				   // functions
  for (unsigned int point=0; point<n_quadrature_points; ++point)
    for (unsigned int shape_func=0; shape_func<total_dofs; ++shape_func)
      values[point] += (dof_values[shape_func] *
			shape_values[selected_dataset](shape_func, point));
};



template <int dim>
const Point<dim> &
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
void FEValuesBase<dim>::get_function_grads (const dVector       &fe_function,
					    vector<Point<dim> > &gradients) const {
  Assert (gradients.size() == n_quadrature_points,
	  ExcWrongVectorSize(gradients.size(), n_quadrature_points));

				   // get function values of dofs
				   // on this cell
  vector<double> dof_values (total_dofs, 0);
  present_cell->get_dof_values (fe_function, dof_values);

				   // initialize with zero
  fill_n (gradients.begin(), n_quadrature_points, Point<dim>());

				   // add up contributions of ansatz
				   // functions
  for (unsigned int point=0; point<n_quadrature_points; ++point)
    for (unsigned int shape_func=0; shape_func<total_dofs; ++shape_func)
      gradients[point] += (dof_values[shape_func] *
			   shape_gradients[shape_func][point]);
};



template <int dim>
const Point<dim> & FEValuesBase<dim>::quadrature_point (const unsigned int i) const {
  Assert (i<n_quadrature_points, ExcInvalidIndex(i, n_quadrature_points));
  Assert (update_flags & update_q_points, ExcAccessToUninitializedField());
  
  return quadrature_points[i];
};



template <int dim>
const Point<dim> & FEValuesBase<dim>::ansatz_point (const unsigned int i) const {
  Assert (i<ansatz_points.size(), ExcInvalidIndex(i, ansatz_points.size()));
  Assert (update_flags & update_ansatz_points, ExcAccessToUninitializedField());
  
  return ansatz_points[i];
};



template <int dim>
double FEValuesBase<dim>::JxW (const unsigned int i) const {
  Assert (i<n_quadrature_points, ExcInvalidIndex(i, n_quadrature_points));
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
				   1,
				   update_flags),
		unit_shape_gradients(fe.total_dofs,
				     vector<Point<dim> >(quadrature.n_quadrature_points)),
		unit_quadrature_points(quadrature.get_quad_points())
{
  for (unsigned int i=0; i<fe.total_dofs; ++i)
    for (unsigned int j=0; j<n_quadrature_points; ++j) 
      {
	shape_values[0](i,j) = fe.shape_value(i, quadrature.quad_point(j));
	unit_shape_gradients[i][j]
	  = fe.shape_grad(i, quadrature.quad_point(j));
      };

  for (unsigned int i=0; i<n_quadrature_points; ++i) 
    {
      weights[i] = quadrature.weight(i);
    };
};







template <int dim>
void FEValues<dim>::reinit (const typename DoFHandler<dim>::cell_iterator &cell,
			    const FiniteElement<dim>                      &fe,
			    const Boundary<dim>                           &boundary) {
  present_cell = cell;
				   // fill jacobi matrices and real
				   // quadrature points
  if ((update_flags & update_jacobians) ||
      (update_flags & update_q_points)  ||
      (update_flags & update_ansatz_points))
    fe.fill_fe_values (cell,
		       unit_quadrature_points,
		       jacobi_matrices,
		       update_flags & update_jacobians,
		       ansatz_points,
		       update_flags & update_ansatz_points,
		       quadrature_points,
		       update_flags & update_q_points,
		       boundary);

				   // compute gradients on real element if
				   // requested
  if (update_flags & update_gradients) 
    {
      Assert (update_flags & update_jacobians, ExcCannotInitializeField());
      
      for (unsigned int i=0; i<fe.total_dofs; ++i)
	for (unsigned int j=0; j<n_quadrature_points; ++j)
	  for (unsigned int s=0; s<dim; ++s)
	    {
	      shape_gradients[i][j](s) = 0;
	      
					       // (grad psi)_s =
					       // (grad_{\xi\eta})_b J_{bs}
					       // with J_{bs}=(d\xi_b)/(dx_s)
	      for (unsigned int b=0; b<dim; ++b)
		shape_gradients[i][j](s)
		  +=
		  unit_shape_gradients[i][j](b) * jacobi_matrices[j](b,s);
	    };
    };
  
  
				   // compute Jacobi determinants in
				   // quadrature points.
				   // refer to the general doc for
				   // why we take the inverse of the
				   // determinant
  if (update_flags & update_JxW_values) 
    {
      Assert (update_flags & update_jacobians, ExcCannotInitializeField());
      for (unsigned int i=0; i<n_quadrature_points; ++i)
	JxW_values[i] = weights[i] / jacobi_matrices[i].determinant();
    };
};





/*------------------------------- FEFaceValues -------------------------------*/


template <int dim>
FEFaceValues<dim>::FEFaceValues (const FiniteElement<dim> &fe,
				 const Quadrature<dim-1>  &quadrature,
				 const UpdateFlags         update_flags) :
		FEValuesBase<dim> (quadrature.n_quadrature_points,
				   fe.dofs_per_face,
				   fe.total_dofs,
				   2*dim,
				   update_flags),
		unit_quadrature_points(quadrature.get_quad_points()),
		face_jacobi_determinants (quadrature.n_quadrature_points,0),
		normal_vectors (quadrature.n_quadrature_points,Point<dim>())
{
  for (unsigned int face=0; face<2*dim; ++face)
    {
      unit_shape_gradients[face].resize (fe.total_dofs,
					 vector<Point<dim> >(quadrature.
							     n_quadrature_points));
      global_unit_quadrature_points[face].resize (quadrature.n_quadrature_points,
						  Point<dim>());
    };

  				   // set up an array of the unit points
				   // on the given face, but in coordinates
				   // of the space with #dim# dimensions.
				   // the points are still on the unit
				   // cell.
  for (unsigned int face=0; face<2*dim; ++face)
    for (unsigned int p=0; p<n_quadrature_points; ++p)
      switch (dim) 
	{
	  case 2:
	  {
	    
	    switch (face)
	      {
		case 0:
		      global_unit_quadrature_points[face][p]
			= Point<dim>(unit_quadrature_points[p](0),0);
		      break;	   
		case 1:
		      global_unit_quadrature_points[face][p]
			= Point<dim>(1,unit_quadrature_points[p](0));
		      break;	   
		case 2:
		      global_unit_quadrature_points[face][p]
			= Point<dim>(unit_quadrature_points[p](0),1);
		      break;	   
		case 3:
		      global_unit_quadrature_points[face][p]
			= Point<dim>(0,unit_quadrature_points[p](0));
		      break;
		default:
		      Assert (false, ExcInternalError());
	      };
	    
	    break;
	  };
	  default:
		Assert (false, ExcNotImplemented());
	};

  for (unsigned int i=0; i<n_quadrature_points; ++i) 
    weights[i] = quadrature.weight(i);

  for (unsigned int face=0; face<2*dim; ++face)
    for (unsigned int i=0; i<fe.total_dofs; ++i)
      for (unsigned int j=0; j<n_quadrature_points; ++j) 
	{
	  shape_values[face](i,j)
	    = fe.shape_value(i, global_unit_quadrature_points[face][j]);
	  unit_shape_gradients[face][i][j]
	    = fe.shape_grad(i, global_unit_quadrature_points[face][j]);
	};
};



template <int dim>
const Point<dim> & FEFaceValues<dim>::normal_vector (const unsigned int i) const {
  Assert (i<normal_vectors.size(), ExcInvalidIndex(i, normal_vectors.size()));
  Assert (update_flags & update_normal_vectors,
	  ExcAccessToUninitializedField());
  
  return normal_vectors[i];
};



template <int dim>
void FEFaceValues<dim>::reinit (const typename DoFHandler<dim>::cell_iterator &cell,
				const unsigned int                             face_no,
				const FiniteElement<dim>                      &fe,
				const Boundary<dim>                           &boundary) {
  present_cell  = cell;
  selected_dataset = face_no;
				   // fill jacobi matrices and real
				   // quadrature points
  if ((update_flags & update_jacobians) ||
      (update_flags & update_q_points)  ||
      (update_flags & update_ansatz_points) ||
      (update_flags & update_JxW_values))
    fe.fill_fe_face_values (cell,
			    face_no,
			    unit_quadrature_points,
			    global_unit_quadrature_points[face_no],
			    jacobi_matrices,
			    update_flags & update_jacobians,
			    ansatz_points,
			    update_flags & update_ansatz_points,
			    quadrature_points,
			    update_flags & update_q_points,
			    face_jacobi_determinants,
			    update_flags & update_JxW_values,
			    normal_vectors,
			    update_flags & update_normal_vectors,
			    boundary);

				   // compute gradients on real element if
				   // requested
  if (update_flags & update_gradients) 
    {
      Assert (update_flags & update_jacobians, ExcCannotInitializeField());
      
      for (unsigned int i=0; i<fe.total_dofs; ++i)
	for (unsigned int j=0; j<n_quadrature_points; ++j) 
	  {
	    shape_gradients[i][j] = Point<dim>();
	    
	    for (unsigned int s=0; s<dim; ++s)
					       // (grad psi)_s =
					       // (grad_{\xi\eta})_b J_{bs}
					       // with J_{bs}=(d\xi_b)/(dx_s)
	      for (unsigned int b=0; b<dim; ++b)
		shape_gradients[i][j](s)
		  += (unit_shape_gradients[face_no][i][j](b) *
		      jacobi_matrices[j](b,s));
	  };
    };
  
  
				   // compute Jacobi determinants in
				   // quadrature points.
				   // refer to the general doc for
				   // why we take the inverse of the
				   // determinant
  if (update_flags & update_JxW_values) 
    {
      Assert (update_flags & update_jacobians,
	      ExcCannotInitializeField());
      for (unsigned int i=0; i<n_quadrature_points; ++i)
	JxW_values[i] = weights[i] * face_jacobi_determinants[i];
    };
};





/*------------------------------- Explicit Instantiations -------------*/

template class FEValuesBase<1>;
template class FEValuesBase<2>;

template class FEValues<1>;
template class FEValues<2>;

template class FEFaceValues<2>;


