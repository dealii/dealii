/* $Id$ */


#include <basic/function.h>
#include <grid/dof.h>
#include <grid/dof_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/dof_constraints.h>
#include <fe/fe.h>
#include <fe/fe_values.h>
#include <fe/quadrature.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <lac/dvector.h>
#include <lac/dsmatrix.h>

#include "../../../mia/control.h"
#include "../../../mia/vectormemory.h"
#include "../../../mia/cg.h"

#include <numeric>
#include <algorithm>
#include <cmath>




inline double sqr (const double x) {
  return x*x;
};


template <int dim>
inline double sqr_point (const Point<dim> &p) {
  return p.square();
};






template <int dim>
void VectorTools<dim>::interpolate (const DoFHandler<dim>    &dof,
				      const FiniteElement<dim> &fe,
				      const Boundary<dim>      &boundary,
				      const Function<dim>      &function,
				      dVector                  &vec) {
  DoFHandler<dim>::active_cell_iterator cell = dof.begin_active(),
					endc = dof.end();
  vector<int>         dofs_on_cell (fe.total_dofs);
  vector<double>      dof_values_on_cell (fe.total_dofs);
  vector<Point<dim> > ansatz_points (fe.total_dofs);
  for (; cell!=endc; ++cell) 
    {
				       // for each cell:
				       // get location of finite element
				       // off-points
      fe.get_ansatz_points (cell, boundary, ansatz_points);
				       // get function values at these points
      function.value_list (ansatz_points, dof_values_on_cell);
				       // get indices of the dofs on this cell
      cell->get_dof_indices (dofs_on_cell);
				       // distribute function values to the
				       // whole vector
      for (unsigned int i=0; i<fe.total_dofs; ++i)
	vec(dofs_on_cell[i]) = dof_values_on_cell[i];
    };
};




void VectorTools<1>::project (const DoFHandler<1>    &,
			      const ConstraintMatrix &,
			      const FiniteElement<1> &,
			      const Quadrature<1>    &,
			      const Quadrature<0>    &,
			      const Boundary<1>      &,
			      const Function<1>      &,
			      const bool              ,
			      dVector                &) {
				   // this function should easily be implemented
				   // using the template below. However some
				   // changes have to be made since faces don't
				   // exist in 1D. Maybe integrate the creation of
				   // zero boundary values into the
				   // project_boundary_values function?
  Assert (false, ExcNotImplemented());
};




template <int dim>
void VectorTools<dim>::project (const DoFHandler<dim>    &dof,
				const ConstraintMatrix   &constraints,
				const FiniteElement<dim> &fe,
				const Quadrature<dim>    &q,
				const Quadrature<dim-1>  &q_boundary,
				const Boundary<dim>      &boundary,
				const Function<dim>      &function,
				const bool                has_zero_boundary,
				dVector                  &vec) {
  				   // make up boundary values
  map<int,double> boundary_values;

  if (has_zero_boundary == false) 
    {
				       // set up a list of boundary functions for
				       // the different boundary parts. We want the
				       // #function# to hold on all parts of the
				       // boundary
      FunctionMap boundary_functions;
      for (unsigned char c=0; c<255; ++c)
	boundary_functions[c] = &function;
      project_boundary_values (dof, boundary_functions, fe, q_boundary,
			       boundary, boundary_values);
    }
  else
				     // no need to project boundary values
    {
      DoFHandler<dim>::active_face_iterator face = dof.begin_active_face(),
					    endf = dof.end_face();
      vector<int> face_dof_indices (fe.dofs_per_face);
      for (; face!=endf; ++face)
	{
	  face->get_dof_indices (face_dof_indices);
	  for (unsigned int i=0; i<fe.dofs_per_face; ++i)
					     // enter zero boundary values
					     // for all boundary nodes
	    boundary_values[face_dof_indices[i]] = 0.;
	};
    };
  
      
				   // set up mass matrix and right hand side
  vec.reinit (dof.n_dofs());
  dSMatrixStruct sparsity(dof.n_dofs(),
			  dof.n_dofs(),
			  dof.max_couplings_between_dofs());
  dof.make_sparsity_pattern (sparsity);
  constraints.condense (sparsity);
  
  dSMatrix mass_matrix (sparsity);
  dVector tmp (mass_matrix.n());
  MatrixCreator<dim>::create_mass_matrix (dof, fe, q, boundary,
					  mass_matrix, function, tmp);

  constraints.condense (mass_matrix);
  MatrixTools<dim>::apply_boundary_values (boundary_values,
					   mass_matrix, vec, tmp);

  int    max_iter  = 1000;
  double tolerance = 1.e-16;
  Control                          control1(max_iter,tolerance);
  PrimitiveVectorMemory<dVector>   memory(tmp.size());
  CG<dSMatrix,dVector>             cg(control1,memory);

				   // solve
  cg (mass_matrix, vec, tmp);
				   // distribute solution
  constraints.distribute (vec);
};



void
VectorTools<1>::interpolate_boundary_values (const DoFHandler<1> &,
					     const FunctionMap &,
					     const FiniteElement<1> &,
					     const Boundary<1> &,
					     map<int,double>   &) {
  Assert (false, ExcNotImplemented());
};




template <int dim>
void
VectorTools<dim>::interpolate_boundary_values (const DoFHandler<dim> &dof,
					       const FunctionMap     &dirichlet_bc,
					       const FiniteElement<dim> &fe,
					       const Boundary<dim>      &boundary,
					       map<int,double>   &boundary_values) {
  Assert (dirichlet_bc.find(255) == dirichlet_bc.end(),
	  ExcInvalidBoundaryIndicator());
				   // use two face iterators, since we need
				   // a DoF-iterator for the dof indices, but
				   // a Tria-iterator for the fe object
  DoFHandler<dim>::active_face_iterator face = dof.begin_active_face(),
					endf = dof.end_face();
  
  FunctionMap::const_iterator function_ptr;

				   // field to store the indices of dofs
  vector<int>         face_dofs (fe.dofs_per_face);
  vector<Point<dim> > dof_locations (face_dofs.size(), Point<dim>());
  vector<double>      dof_values (fe.dofs_per_face);
	
  for (; face!=endf; ++face)
    if ((function_ptr = dirichlet_bc.find(face->boundary_indicator())) !=
	dirichlet_bc.end()) 
				       // face is subject to one of the
				       // bc listed in #dirichlet_bc#
      {
					 // get indices, physical location and
					 // boundary values of dofs on this
					 // face
	face->get_dof_indices (face_dofs);
	fe.get_face_ansatz_points (face, boundary, dof_locations);
	function_ptr->second->value_list (dof_locations, dof_values);

					 // enter into list
	for (unsigned int i=0; i<face_dofs.size(); ++i)
	  boundary_values[face_dofs[i]] = dof_values[i];
      };
};



template <int dim>
void
VectorTools<dim>::project_boundary_values (const DoFHandler<dim>    &dof,
					   const FunctionMap        &boundary_functions,
					   const FiniteElement<dim> &fe,
					   const Quadrature<dim-1>  &q,
					   const Boundary<dim>      &boundary,
					   map<int,double>   &boundary_values) {
  vector<int>    dof_to_boundary_mapping;
  dof.map_dof_to_boundary_indices (dof_to_boundary_mapping);
  
				   // set up sparsity structure
  dSMatrixStruct sparsity(dof.n_boundary_dofs(boundary_functions),
			  dof.max_couplings_between_boundary_dofs());
  dof.make_boundary_sparsity_pattern (boundary_functions, dof_to_boundary_mapping,
				      sparsity);

				   // note: for three or more dimensions, there
				   // may be constrained nodes on the boundary
				   // in this case the boundary mass matrix has
				   // to be condensed and the solution is to
				   // be distributed afterwards, which is not
				   // yet implemented
  if (dim<3)
    sparsity.compress();
  else
    Assert (false, ExcNotImplemented());
  

				   // make mass matrix and right hand side
  dSMatrix       mass_matrix(sparsity);
  dVector        rhs(sparsity.n_rows());
  

  MatrixTools<dim>::create_boundary_mass_matrix (dof, fe, q, boundary,
						 mass_matrix, boundary_functions,
						 rhs, dof_to_boundary_mapping);
  
};




template <int dim>
void VectorTools<dim>::integrate_difference (const DoFHandler<dim>    &dof,
					     const dVector            &fe_function,
					     const Function<dim>      &exact_solution,
					     dVector                  &difference,
					     const Quadrature<dim>    &q,
					     const FiniteElement<dim> &fe,
					     const NormType           &norm,
					     const Boundary<dim>      &boundary) {
  Assert (fe == dof.get_selected_fe(), ExcInvalidFE());
  
  difference.reinit (dof.get_tria().n_active_cells());
  
  UpdateFlags update_flags = UpdateFlags (update_q_points  |
					  update_JxW_values);
  if ((norm==H1_seminorm) || (norm==H1_norm))
    update_flags = UpdateFlags (update_flags | update_gradients);
  FEValues<dim> fe_values(fe, q, update_flags);
  
				   // loop over all cells
  DoFHandler<dim>::active_cell_iterator cell = dof.begin_active(),
					endc = dof.end();
  for (unsigned int index=0; cell != endc; ++cell, ++index)
    {
      double diff=0;
				       // initialize for this cell
      fe_values.reinit (cell, fe, boundary);

      switch (norm) 
	{
	  case mean:
	  case L1_norm:
	  case L2_norm:
	  case Linfty_norm:
	  case H1_norm:
	  {
					     // we need the finite element
					     // function \psi at the different
					     // integration points. Compute
					     // it like this:
					     // \psi(x_j)=\sum_i v_i \phi_i(x_j)
					     // with v_i the nodal values of the
					     // fe_function and \phi_i(x_j) the
					     // matrix of the ansatz function
					     // values at the integration point
					     // x_j. Then the vector
					     // of the \psi(x_j) is v*Phi with
					     // v being the vector of nodal
					     // values on this cell and Phi
					     // the matrix.
					     //
					     // we then need the difference:
					     // reference_function(x_j)-\psi_j
					     // and assign that to the vector
					     // \psi.
	    const unsigned int n_q_points = q.n_quadrature_points;
	    vector<double>   psi (n_q_points);

					     // in praxi: first compute
					     // exact fe_function vector
	    exact_solution.value_list (fe_values.get_quadrature_points(),
				       psi);
					     // then subtract finite element
					     // fe_function
	    if (true) 
	      {
		vector<double> function_values (n_q_points, 0);
		fe_values.get_function_values (fe_function, function_values);

		transform (psi.begin(), psi.end(),
			   function_values.begin(),
			   psi.begin(),
			   minus<double>());
	      };	    

					     // for L1_norm and Linfty_norm:
					     // take absolute
					     // value, for the L2_norm take
					     // square of psi
	    switch (norm) 
	      {
		case mean:
		      break;
		case L1_norm:
		case Linfty_norm:
		      transform (psi.begin(), psi.end(),
				 psi.begin(), ptr_fun(fabs));
		      break;
		case L2_norm:
		case H1_norm:
		      transform (psi.begin(), psi.end(),
				 psi.begin(), ptr_fun(sqr));
		      break;
		default:
		      Assert (false, ExcNotImplemented());
	      };

					     // ok, now we have the integrand,
					     // let's compute the integral,
					     // which is
					     // sum_j psi_j JxW_j
					     // (or |psi_j| or |psi_j|^2
	    switch (norm) 
	      {
		case mean:
		case L1_norm:
		      diff = inner_product (psi.begin(), psi.end(),
					    fe_values.get_JxW_values().begin(),
					    0.0);
		      break;
		case L2_norm:
		case H1_norm:
		      diff = sqrt(inner_product (psi.begin(), psi.end(),
						 fe_values.get_JxW_values().begin(),
						 0.0));
		      break;
		case Linfty_norm:
		      diff = *max_element (psi.begin(), psi.end());
		      break;
		default:
		      Assert (false, ExcNotImplemented());
	      };

					     // note: the H1_norm uses the result
					     // of the L2_norm and control goes
					     // over to the next case statement!
	    if (norm != H1_norm)
	      break;
	  };

	  case H1_seminorm:
	  {
					     // note: the computation of the
					     // H1_norm starts at the previous
					     // case statement, but continues
					     // here!

					     // for H1_norm: re-square L2_norm.
	    diff = sqr(diff);

					     // same procedure as above, but now
					     // psi is a vector of gradients
	    const unsigned int n_q_points = q.n_quadrature_points;
	    vector<Point<dim> >   psi (n_q_points);

					     // in praxi: first compute
					     // exact fe_function vector
	    exact_solution.gradient_list (fe_values.get_quadrature_points(),
					  psi);
	    
					     // then subtract finite element
					     // fe_function
	    if (true) 
	      {
		vector<Point<dim> > function_grads (n_q_points, Point<dim>());
		fe_values.get_function_grads (fe_function, function_grads);

		transform (psi.begin(), psi.end(),
			   function_grads.begin(),
			   psi.begin(),
			   minus<Point<dim> >());
	      };
					     // take square of integrand
	    vector<double> psi_square (psi.size(), 0.0);
	    for (unsigned int i=0; i<n_q_points; ++i)
	      psi_square[i] = sqr_point(psi[i]);

					     // add seminorm to L_2 norm or
					     // to zero
	    diff += inner_product (psi_square.begin(), psi_square.end(),
				   fe_values.get_JxW_values().begin(),
				   0.0);
	    diff = sqrt(diff);

	    break;
	  };
					     
	  default:
		Assert (false, ExcNotImplemented());
	};

      
				       // append result of this cell
				       // to the end of the vector
      difference(index) = diff;
    };
};


template VectorTools<1>;
template VectorTools<2>;
