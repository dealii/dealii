/* $Id$ */

#include <numerics/assembler.h>
#include <numerics/base.h>
#include <grid/dof_constraints.h>
#include <grid/tria_iterator.h>
#include <basic/data_io.h>
#include <basic/function.h>

#include "../../../mia/control.h"
#include "../../../mia/vectormemory.h"
#include "../../../mia/cg.h"

#include <algo.h>




inline double sqr (double x) {
  return x*x;
};






template <int dim>
ProblemBase<dim>::ProblemBase (Triangulation<dim> *tria,
			       DoFHandler<dim>    *dof) :
		tria(tria),
		dof_handler(dof),
		system_sparsity(1,1,1),   // dummy initialisation, is later reinit'd
		system_matrix()           // dummy initialisation, is later reinit'd
{
  Assert (tria == &dof->get_tria(), ExcDofAndTriaDontMatch());
};



template <int dim>
ProblemBase<dim>::~ProblemBase () {};



template <int dim>
void ProblemBase<dim>::assemble (const Equation<dim>      &equation,
				 const Quadrature<dim>    &quadrature,
				 const FiniteElement<dim> &fe,
				 const DirichletBC        &dirichlet_bc) {
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
			   *this,
			   quadrature,
			   fe);
  TriaActiveIterator<dim, Assembler<dim> > assembler (tria,
						      tria->begin_active()->level(),
						      tria->begin_active()->index(),
						      &data);

				   // loop over all cells, fill matrix and rhs
  do 
    {
      assembler->assemble (equation);
    }
  while ((++assembler).state() == valid);

				   // apply Dirichlet bc as described
				   // in the docs
  apply_dirichlet_bc (system_matrix, solution,
		      right_hand_side, dirichlet_bc);
  
				   // condense system matrix in-place
  constraints.condense (system_matrix);

				   // condense right hand side in-place
  constraints.condense (right_hand_side);
};




template <int dim>
void ProblemBase<dim>::solve () {
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
					     vector<double>           &difference,
					     const Quadrature<dim>    &q,
					     const FiniteElement<dim> &fe,
					     const NormType           &norm) const {
  Assert (fe == dof_handler->get_selected_fe(), ExcInvalidFE());

  difference.erase (difference.begin(), difference.end());
  difference.reserve (tria->n_cells());

  FEValues<dim> fe_values(fe, q);
  
				   // loop over all cells
				   // (we need an iterator on the triangulation for
				   // fe_values.reinit, but we also need an iterator
				   // on the dofhandler; we could generate the first
				   // out of the second or vice versa each time
				   // we need them, but conversion seems to be slow.
				   // therefore, we keep two iterators, hoping that
				   // incrementing is faster than conversion...)
  DoFHandler<dim>::active_cell_iterator cell = dof_handler->begin_active(),
					endc = dof_handler->end();
  Triangulation<dim>::active_cell_iterator tria_cell = dof_handler->get_tria().begin();
  
  for (; cell != endc; ++cell, ++tria_cell) 
    {
      double diff=0;
				       // initialize for this cell
      fe_values.reinit (tria_cell, fe);

      switch (norm) 
	{
	  case mean:
	  case L1_norm:
	  case L2_norm:
	  case Linfty_norm:
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
	    for (unsigned int j=0; j<n_dofs; ++j)
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
				 psi.begin(), ptr_fun(abs));
		      break;
		case L2_norm:
		      transform (psi.begin(), psi.end(),
				 psi.begin(), ptr_fun(sqr));
		      break;
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
		case L2_norm:
		      diff = inner_product (psi.begin(), psi.end(),
					    fe_values.get_JxW_values().begin(),
					    0.0);
		      break;
		case Linfty_norm:
		      diff = *max_element (psi.begin(), psi.end());
		      break;
	      };

	    if (norm==L2_norm)
	      diff = sqrt(diff);
	    
	    break;
	  };

	  default:
		Assert (false, ExcNotImplemented());
	};

      
				       // append result of this cell
				       // to the end of the vector
      difference.push_back (diff);
    };
};



template <int dim>
void ProblemBase<dim>::fill_data (DataOut<dim> &out) const {
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
					   const DirichletBC &dirichlet_bc) {
//  l;
};

  


// explicit instantiations
template class ProblemBase<1>;
template class ProblemBase<2>;
