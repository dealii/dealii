/* $Id$ */
/* Copyright W. Bangerth, University of Heidelberg, 1998 */


#include "poisson.h"




#if deal_II_dimension == 1

void PoissonEquation<1>::assemble (dFMatrix            &cell_matrix,
				   dVector             &rhs,
				   const FEValues<1>   &fe_values,
				   const Triangulation<1>::cell_iterator &) const {
  for (unsigned int point=0; point<fe_values.n_quadrature_points; ++point)
    for (unsigned int i=0; i<fe_values.total_dofs; ++i) 
      {
	for (unsigned int j=0; j<fe_values.total_dofs; ++j)
	  cell_matrix(i,j) += (fe_values.shape_grad(i,point) *
			       fe_values.shape_grad(j,point)) *
			      fe_values.JxW(point);
	rhs(i) += fe_values.shape_value(i,point) *
		  right_hand_side(fe_values.quadrature_point(point)) *
		  fe_values.JxW(point);
      };
};

#endif



#if deal_II_dimension == 2

void PoissonEquation<2>::assemble (dFMatrix            &cell_matrix,
				   dVector             &rhs,
				   const FEValues<2>   &fe_values,
				   const Triangulation<2>::cell_iterator &) const {
  const vector<vector<Point<2> > >&gradients = fe_values.get_shape_grads ();
  const dFMatrix       &values    = fe_values.get_shape_values ();
  vector<double>        rhs_values (fe_values.n_quadrature_points);
  const vector<double> &weights   = fe_values.get_JxW_values ();

  right_hand_side.value_list (fe_values.get_quadrature_points(), rhs_values);
   
  for (unsigned int point=0; point<fe_values.n_quadrature_points; ++point)
    for (unsigned int i=0; i<fe_values.total_dofs; ++i) 
      {
	for (unsigned int j=0; j<fe_values.total_dofs; ++j)
	  cell_matrix(i,j) += (gradients[i][point] *
			       gradients[j][point]) *
			      weights[point];
	rhs(i) += values(i,point) *
		  rhs_values[point] *
		  weights[point];
      };
};

#endif



template <int dim>
void PoissonEquation<dim>::assemble (dFMatrix            &,
				     const FEValues<dim> &,
				     const Triangulation<dim>::cell_iterator &) const {
  Assert (false, ExcPureVirtualFunctionCalled());
};



template <int dim>
void PoissonEquation<dim>::assemble (dVector             &,
				     const FEValues<dim> &,
				     const Triangulation<dim>::cell_iterator &) const {
  Assert (false, ExcPureVirtualFunctionCalled());
};






template class PoissonEquation<2>;
