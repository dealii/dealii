/* $Id$ */

#include <fe/fe.h>
#include <fe/fe_values.h>
#include <fe/quadrature.h>
#include <grid/tria_iterator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_boundary.h>



/*------------------------------- FEValues -------------------------------*/


template <int dim>
FEValues<dim>::FEValues (const FiniteElement<dim> &fe,
			 const Quadrature<dim>    &quadrature,
			 const UpdateFields        update_flags) :
		n_quadrature_points(quadrature.n_quadrature_points),
		total_dofs(fe.total_dofs),
		shape_values(fe.total_dofs, quadrature.n_quadrature_points),
		shape_gradients(fe.total_dofs,
				vector<Point<dim> >(quadrature.n_quadrature_points)),
		unit_shape_gradients(fe.total_dofs,
				     vector<Point<dim> >(quadrature.n_quadrature_points)),
		weights(quadrature.n_quadrature_points, 0),
		JxW_values(quadrature.n_quadrature_points, 0),
		quadrature_points(quadrature.n_quadrature_points, Point<dim>()),
		unit_quadrature_points(quadrature.n_quadrature_points, Point<dim>()),
		ansatz_points (fe.total_dofs, Point<dim>()),
		jacobi_matrices (quadrature.n_quadrature_points,
				 dFMatrix(dim,dim)),
		update_flags (update_flags)
{
  for (unsigned int i=0; i<fe.total_dofs; ++i)
    for (unsigned int j=0; j<n_quadrature_points; ++j) 
      {
	shape_values(i,j) = fe.shape_value(i, quadrature.quad_point(j));
	unit_shape_gradients[i][j]
	  = fe.shape_grad(i, quadrature.quad_point(j));
      };

  for (unsigned int i=0; i<n_quadrature_points; ++i) 
    {
      weights[i] = quadrature.weight(i);
      unit_quadrature_points[i] = quadrature.quad_point(i);
    };
};



template <int dim>
double FEValues<dim>::shape_value (const unsigned int i,
				   const unsigned int j) const {
  Assert (i<shape_values.m(), ExcInvalidIndex (i, shape_values.m()));
  Assert (j<shape_values.n(), ExcInvalidIndex (j, shape_values.n()));

  return shape_values(i,j);
};



template <int dim>
const Point<dim> &
FEValues<dim>::shape_grad (const unsigned int i,
			   const unsigned int j) const {
  Assert (i<shape_values.m(), ExcInvalidIndex (i, shape_values.m()));
  Assert (j<shape_values.n(), ExcInvalidIndex (j, shape_values.n()));
  Assert (update_flags | update_gradients, ExcAccessToUninitializedField());

  return shape_gradients[i][j];
};



template <int dim>
const Point<dim> & FEValues<dim>::quadrature_point (const unsigned int i) const {
  Assert (i<n_quadrature_points, ExcInvalidIndex(i, n_quadrature_points));
  Assert (update_flags | update_q_points, ExcAccessToUninitializedField());
  
  return quadrature_points[i];
};



template <int dim>
const Point<dim> & FEValues<dim>::ansatz_point (const unsigned int i) const {
  Assert (i<ansatz_points.size(), ExcInvalidIndex(i, ansatz_points.size()));
  Assert (update_flags | update_ansatz_points, ExcAccessToUninitializedField());
  
  return ansatz_points[i];
};



template <int dim>
double FEValues<dim>::JxW (const unsigned int i) const {
  Assert (i<n_quadrature_points, ExcInvalidIndex(i, n_quadrature_points));
  Assert (update_flags | update_JxW_values, ExcAccessToUninitializedField());
  
  return JxW_values[i];
};



template <int dim>
void FEValues<dim>::reinit (const typename Triangulation<dim>::cell_iterator &cell,
			    const FiniteElement<dim>                         &fe,
			    const Boundary<dim>                              &boundary) {
				   // fill jacobi matrices and real
				   // quadrature points
  if ((update_flags | update_jacobians) ||
      (update_flags | update_q_points))
    fe.fill_fe_values (cell,
		       unit_quadrature_points,
		       jacobi_matrices,
		       update_flags | update_jacobians,
		       ansatz_points,
		       update_flags | update_ansatz_points,
		       quadrature_points,
		       update_flags | update_q_points,
		       boundary);

				   // compute gradients on real element if
				   // requested
  if (update_flags | update_gradients) 
    {
      Assert (update_flags | update_jacobians, ExcCannotInitializeField());
      
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
  if (update_flags | update_JxW_values) 
    {
      Assert (update_flags | update_jacobians, ExcCannotInitializeField());
      for (unsigned int i=0; i<n_quadrature_points; ++i)
	JxW_values[i] = weights[i] / jacobi_matrices[i].determinant();
    };
};





/*------------------------------- Explicit Instantiations -------------*/

template class FEValues<1>;
template class FEValues<2>;

