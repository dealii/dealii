/* $Id$ */

#include <numerics/assembler.h>
#include <numerics/base.h>
#include <grid/dof_constraints.h>
#include <grid/tria_iterator.h>
#include <basic/data_io.h>
#include <basic/function.h>
#include <fe/fe.h>
#include <fe/quadrature.h>

#include "../../../mia/control.h"
#include "../../../mia/vectormemory.h"
#include "../../../mia/cg.h"

#include <map>
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
ProblemBase<dim>::ProblemBase () :
		tria(0),
		dof_handler(0),
		system_sparsity(1,1,1),   // dummy initialisation, is later reinit'd
		system_matrix()           // dummy initialisation, is later reinit'd
{};



template <int dim>
void ProblemBase<dim>::set_tria_and_dof (Triangulation<dim> *t,
					 DoFHandler<dim>    *d) {
  tria        = t;
  dof_handler = d;

  Assert ((tria!=0) && (dof_handler!=0), ExcNoTriaSelected());
  Assert (tria == &dof_handler->get_tria(), ExcDofAndTriaDontMatch());
};



template <int dim>
void ProblemBase<dim>::clear () {
  tria        = 0;
  dof_handler = 0;
  system_sparsity.reinit (1,1,1);
  system_matrix.clear ();
  right_hand_side.reinit (1);
  solution.reinit (1);
  constraints.clear ();
};



template <int dim>
ProblemBase<dim>::~ProblemBase () {};



template <int dim>
void ProblemBase<dim>::assemble (const Equation<dim>      &equation,
				 const Quadrature<dim>    &quadrature,
				 const FiniteElement<dim> &fe,
				 const UpdateFlags         update_flags,
				 const DirichletBC        &dirichlet_bc,
				 const Boundary<dim>      &boundary) {
  Assert ((tria!=0) && (dof_handler!=0), ExcNoTriaSelected());
  
  system_sparsity.reinit (dof_handler->n_dofs(),
			  dof_handler->n_dofs(),
			  dof_handler->max_couplings_between_dofs());
  right_hand_side.reinit (dof_handler->n_dofs());
  
				   // make up sparsity pattern and
				   // compress with constraints
  constraints.clear ();
  dof_handler->make_constraint_matrix (constraints);
  dof_handler->make_sparsity_pattern (system_sparsity);
  constraints.condense (system_sparsity);

				   // reinite system matrix
  system_matrix.reinit (system_sparsity);
				   // reinit solution vector, preset
				   // with zeroes.
  solution.reinit (dof_handler->n_dofs());
  
				   // create assembler
  AssemblerData<dim> data (*dof_handler,
			   true, true, //assemble matrix and rhs
			   system_matrix,
			   right_hand_side,
			   quadrature,
			   fe,
			   update_flags,
			   boundary);
  active_assemble_iterator assembler (tria,
				      tria->begin_active()->level(),
				      tria->begin_active()->index(),
				      &data);
				   // loop over all cells, fill matrix and rhs
  do 
    {
      assembler->assemble (equation);
    }
  while ((++assembler).state() == valid);

				   // condense system matrix in-place
  constraints.condense (system_matrix);

				   // condense right hand side in-place
  constraints.condense (right_hand_side);

				   // apply Dirichlet bc as described
				   // in the docs
  apply_dirichlet_bc (system_matrix, solution,
		      right_hand_side,
		      dirichlet_bc, fe, boundary);
  
};




template <int dim>
void ProblemBase<dim>::solve () {
  Assert ((tria!=0) && (dof_handler!=0), ExcNoTriaSelected());
  
  int    max_iter  = 4000;
  double tolerance = 1.e-16;
  
  Control                          control1(max_iter,tolerance);
  PrimitiveVectorMemory<dVector>   memory(right_hand_side.n());
  CG<dSMatrix,dVector>             cg(control1,memory);

				   // solve
  cg (system_matrix, solution, right_hand_side);
				   // distribute solution
  constraints.distribute (solution);
};



template <int dim>
void ProblemBase<dim>::integrate_difference (const Function<dim>      &exact_solution,
					     dVector                  &difference,
					     const Quadrature<dim>    &q,
					     const FiniteElement<dim> &fe,
					     const NormType           &norm,
					     const Boundary<dim>      &boundary) const {
  Assert ((tria!=0) && (dof_handler!=0), ExcNoTriaSelected());  
  Assert (fe == dof_handler->get_selected_fe(), ExcInvalidFE());

  difference.reinit (tria->n_active_cells());
  
  UpdateFlags update_flags = UpdateFlags (update_q_points  |
					  update_jacobians |
					  update_JxW_values);
  if ((norm==H1_seminorm) || (norm==H1_norm))
    update_flags = UpdateFlags (update_flags | update_gradients);
  FEValues<dim> fe_values(fe, q, update_flags);
  
				   // loop over all cells
  DoFHandler<dim>::active_cell_iterator cell = dof_handler->begin_active(),
					endc = dof_handler->end();
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
					     // solution and \phi_i(x_j) the
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
	    const unsigned int n_dofs = fe.total_dofs;
	    const unsigned int n_q_points = q.n_quadrature_points;
	    const dFMatrix & shape_values = fe_values.get_shape_values();
	    vector<double>   dof_values;
	    cell->get_dof_values (solution, dof_values);
	    
	    vector<double>   psi;

					     // in praxi: first compute
					     // exact solution vector
	    exact_solution.value_list (fe_values.get_quadrature_points(),
				       psi);
					     // then subtract finite element
					     // solution
	    for (unsigned int j=0; j<n_q_points; ++j) 
	      for (unsigned int i=0; i<n_dofs; ++i)
		psi[j] -= dof_values[i]*shape_values(i,j);

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
	    const unsigned int n_dofs = fe.total_dofs;
	    const unsigned int n_q_points = q.n_quadrature_points;
	    const vector<vector<Point<dim> > > & shape_grads = fe_values.get_shape_grads();
	    vector<double>   dof_values;
	    cell->get_dof_values (solution, dof_values);
	    
	    vector<Point<dim> >   psi;

					     // in praxi: first compute
					     // exact solution vector
	    exact_solution.gradient_list (fe_values.get_quadrature_points(),
					  psi);
	    
					     // then subtract finite element
					     // solution
	    for (unsigned int j=0; j<n_q_points; ++j) 
	      for (unsigned int i=0; i<n_dofs; ++i)
		psi[j] -= dof_values[i]*shape_grads[i][j];

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



template <int dim>
void ProblemBase<dim>::fill_data (DataOut<dim> &out) const {
  Assert ((tria!=0) && (dof_handler!=0), ExcNoTriaSelected());
  
  out.clear_data_vectors ();
  out.attach_dof_handler (*dof_handler);

  pair<char*,char*> solution_name = get_solution_name ();
  out.add_data_vector (solution,
		       solution_name.first, solution_name.second);
};



template <int dim>
pair<char*,char*> ProblemBase<dim>::get_solution_name () const {
  return pair<char*,char*>("solution", "");
};




template <int dim>
void ProblemBase<dim>::apply_dirichlet_bc (dSMatrix &matrix,
					   dVector  &solution,
					   dVector  &right_hand_side,
					   const DirichletBC &dirichlet_bc,
					   const FiniteElement<dim> &fe,
					   const Boundary<dim>      &boundary) {
  Assert ((tria!=0) && (dof_handler!=0), ExcNoTriaSelected());
  Assert (dirichlet_bc.find(255) == dirichlet_bc.end(),
	  ExcInvalidBoundaryIndicator());

				   // first make up a list of dofs subject
				   // to any boundary condition and which
				   // value they take; if a node occurs
				   // with two bc (e.g. a corner node, with
				   // the lines in 2D being subject to
				   // different bc's), the last value is taken
  map<int,double> boundary_values;
  make_boundary_value_list (dirichlet_bc, fe, boundary, boundary_values);

  map<int,double>::const_iterator dof, endd;
  const unsigned int n_dofs   = matrix.m();
  const dSMatrixStruct &sparsity = matrix.get_sparsity_pattern();
  const unsigned int *sparsity_rowstart = sparsity.get_rowstart_indices();
  const int          *sparsity_colnums  = sparsity.get_column_numbers();

  for (dof=boundary_values.begin(), endd=boundary_values.end(); dof != endd; ++dof)
    {
				       // for each boundary dof:

				       // set entries of this line
				       // to zero
      for (unsigned int j=sparsity_rowstart[dof->first];
	   j<sparsity_rowstart[dof->first+1]; ++j)
	if (sparsity_colnums[j] != dof->first)
					   // if not main diagonal entry
	  matrix.global_entry(j) = 0.;
      
				       // set right hand side to
				       // wanted value: if main diagonal
				       // entry nonzero, don't touch it
				       // and scale rhs accordingly. If
				       // zero, take the first main
				       // diagonal entry we can find, or
				       // one if no nonzero main diagonal
				       // element exists. Normally, however,
				       // the main diagonal entry should
				       // not be zero.
				       //
				       // store the new rhs entry to make
				       // the gauss step more efficient
      double new_rhs;
      if (matrix.diag_element(dof->first) != 0.0)
	new_rhs = right_hand_side(dof->first)
		= dof->second * matrix.diag_element(dof->first);
      else
	{
	  double first_diagonal_entry = 1;
	  for (unsigned int i=0; i<n_dofs; ++i)
	    if (matrix.diag_element(i) != 0)
	      {
		first_diagonal_entry = matrix.diag_element(i);
		break;
	      };
	  
	  matrix.set(dof->first, dof->first,
		     first_diagonal_entry);
	  new_rhs = right_hand_side(dof->first)
		  = dof->second * first_diagonal_entry;
	};
      
				       // store the only nonzero entry
				       // of this line for the Gauss
				       // elimination step
      const double diagonal_entry = matrix.diag_element(dof->first);

				       // do the Gauss step
      for (unsigned int row=0; row<n_dofs; ++row) 
	for (unsigned int j=sparsity_rowstart[row];
	     j<sparsity_rowstart[row+1]; ++j)
	  if ((sparsity_colnums[j] == (signed int)dof->first) &&
	      ((signed int)row != dof->first))
					     // this line has an entry
					     // in the regarding column
					     // but this is not the main
					     // diagonal entry
	    {
					       // correct right hand side
	      right_hand_side(row) -= matrix.global_entry(j)/diagonal_entry *
				      new_rhs;
	      
					       // set matrix entry to zero
	      matrix.global_entry(j) = 0.;
	    };
      
      
				       // preset solution vector
      solution(dof->first) = dof->second;
    };
};





void
ProblemBase<1>::make_boundary_value_list (const DirichletBC &,
					  const FiniteElement<1> &,
					  const Boundary<1> &,
					  map<int,double>   &) const {
  Assert ((tria!=0) && (dof_handler!=0), ExcNoTriaSelected());  
  Assert (false, ExcNotImplemented());
};




template <int dim>
void
ProblemBase<dim>::make_boundary_value_list (const DirichletBC        &dirichlet_bc,
					    const FiniteElement<dim> &fe,
					    const Boundary<dim>      &boundary,
					    map<int,double>   &boundary_values) const {
  Assert ((tria!=0) && (dof_handler!=0), ExcNoTriaSelected());

				   // use two face iterators, since we need
				   // a DoF-iterator for the dof indices, but
				   // a Tria-iterator for the fe object
  DoFHandler<dim>::active_face_iterator face = dof_handler->begin_active_face(),
					endf = dof_handler->end_face();
  
  DirichletBC::const_iterator function_ptr;

				   // field to store the indices of dofs
				   // initialize once to get the size right
				   // for the following fields.
  vector<int>         face_dofs;
  face->get_dof_indices (face_dofs);
  vector<Point<dim> > dof_locations (face_dofs.size(), Point<dim>());
  vector<double>      dof_values;
	
  for (; face!=endf; ++face)
    if ((function_ptr = dirichlet_bc.find(face->boundary_indicator())) !=
	dirichlet_bc.end()) 
				       // face is subject to one of the
				       // bc listed in #dirichlet_bc#
      {
					 // get indices, physical location and
					 // boundary values of dofs on this
					 // face
	face_dofs.erase (face_dofs.begin(), face_dofs.end());
	dof_values.erase (dof_values.begin(), dof_values.end());
	face->get_dof_indices (face_dofs);
	fe.get_face_ansatz_points (face, boundary, dof_locations);
	function_ptr->second->value_list (dof_locations, dof_values);

					 // enter into list
	for (unsigned int i=0; i<face_dofs.size(); ++i)
	  boundary_values[face_dofs[i]] = dof_values[i];
      };
};

  


// explicit instantiations
template class ProblemBase<1>;
template class ProblemBase<2>;
