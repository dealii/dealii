/*      $Id$                 */

#include <basic/function.h>
#include <grid/dof.h>
#include <grid/dof_accessor.h>
#include <grid/tria_iterator.h>
#include <fe/quadrature.h>
#include <fe/fe.h>
#include <fe/fe_values.h>
#include <numerics/matrices.h>
#include <numerics/assembler.h>
#include <lac/dsmatrix.h>



template <int dim>
void MatrixCreator<dim>::create_mass_matrix (const DoFHandler<dim>    &dof,
					     const FiniteElement<dim> &fe,
					     const Quadrature<dim>    &q,
					     const Boundary<dim>      &boundary,
					     dSMatrix                 &matrix,
					     const Function<dim> * const a) {
  dVector dummy;    // no entries, should give an error if accessed
  UpdateFlags update_flags = UpdateFlags(update_jacobians |
					 update_JxW_values);
  if (a != 0)
    update_flags = UpdateFlags (update_flags | update_q_points);
  const AssemblerData<dim> data (dof,
				 true, false,  // assemble matrix but not rhs
				 matrix, dummy,
				 q, fe, update_flags, boundary);
  TriaActiveIterator<dim, Assembler<dim> >
    assembler (const_cast<Triangulation<dim>*>(&dof.get_tria()),
	       dof.get_tria().begin_active()->level(),
	       dof.get_tria().begin_active()->index(),
	       &data);
  MassMatrix<dim> equation(0,a);
  do 
    {
      assembler->assemble (equation);
    }
  while ((++assembler).state() == valid);
};




template <int dim>
void MatrixCreator<dim>::create_mass_matrix (const DoFHandler<dim>    &dof,
					     const FiniteElement<dim> &fe,
					     const Quadrature<dim>    &q,
					     const Boundary<dim>      &boundary,
					     dSMatrix                 &matrix,
					     const Function<dim>      &rhs,
					     dVector                  &rhs_vector,
					     const Function<dim> * const a) {
  UpdateFlags update_flags = UpdateFlags(update_q_points |
					 update_jacobians |
					 update_JxW_values);
  const AssemblerData<dim> data (dof,
				 true, true,
				 matrix, rhs_vector,
				 q, fe,	 update_flags, boundary);
  TriaActiveIterator<dim, Assembler<dim> >
    assembler (const_cast<Triangulation<dim>*>(&dof.get_tria()),
	       dof.get_tria().begin_active()->level(),
	       dof.get_tria().begin_active()->index(),
	       &data);
  MassMatrix<dim> equation(&rhs,a);
  do 
    {
      assembler->assemble (equation);
    }
  while ((++assembler).state() == valid);
};




template <int dim>
void MatrixCreator<dim>::create_laplace_matrix (const DoFHandler<dim>    &dof,
						const FiniteElement<dim> &fe,
						const Quadrature<dim>    &q,
						const Boundary<dim>      &boundary,
						dSMatrix                 &matrix,
						const Function<dim> * const a) {
  dVector dummy;   // no entries, should give an error if accessed
  UpdateFlags update_flags = UpdateFlags(update_gradients |
					 update_jacobians |
					 update_JxW_values);
  if (a != 0)
    update_flags = UpdateFlags(update_flags | update_q_points);
  const AssemblerData<dim> data (dof,
				 true, false,  // assemble matrix but not rhs
				 matrix, dummy,
				 q, fe,	 update_flags, boundary);
  TriaActiveIterator<dim, Assembler<dim> >
    assembler (const_cast<Triangulation<dim>*>(&dof.get_tria()),
	       dof.get_tria().begin_active()->level(),
	       dof.get_tria().begin_active()->index(),
	       &data);
  LaplaceMatrix<dim> equation (0, a);
  do 
    {
      assembler->assemble (equation);
    }
  while ((++assembler).state() == valid);
};



template <int dim>
void MatrixCreator<dim>::create_laplace_matrix (const DoFHandler<dim>    &dof,
						const FiniteElement<dim> &fe,
						const Quadrature<dim>    &q,
						const Boundary<dim>      &boundary,
						dSMatrix                 &matrix,
						const Function<dim>      &rhs,
						dVector                  &rhs_vector,
						const Function<dim> * const a) {
  UpdateFlags update_flags = UpdateFlags(update_q_points  |
					 update_gradients |
					 update_jacobians |
					 update_JxW_values);
  const AssemblerData<dim> data (dof,
				 true, true,
				 matrix, rhs_vector,
				 q, fe,
				 update_flags,
				 boundary);
  TriaActiveIterator<dim, Assembler<dim> >
    assembler (const_cast<Triangulation<dim>*>(&dof.get_tria()),
	       dof.get_tria().begin_active()->level(),
	       dof.get_tria().begin_active()->index(),
	       &data);
  LaplaceMatrix<dim> equation (&rhs, a);
  do 
    {
      assembler->assemble (equation);
    }
  while ((++assembler).state() == valid);
};






template <int dim>
void MatrixTools<dim>::apply_boundary_values (const map<int,double> &boundary_values,
					      dSMatrix  &matrix,
					      dVector   &solution,
					      dVector   &right_hand_side) {

  map<int,double>::const_iterator dof  = boundary_values.begin(),
				  endd = boundary_values.end();
  const unsigned int n_dofs             = matrix.m();
  const dSMatrixStruct &sparsity        = matrix.get_sparsity_pattern();
  const unsigned int *sparsity_rowstart = sparsity.get_rowstart_indices();
  const int          *sparsity_colnums  = sparsity.get_column_numbers();

  for (; dof != endd; ++dof)
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
MatrixTools<1>::interpolate_boundary_values (const DoFHandler<1> &,
					     const FunctionMap &,
					     const FiniteElement<1> &,
					     const Boundary<1> &,
					     map<int,double>   &) {
  Assert (false, ExcNotImplemented());
};




template <int dim>
void
MatrixTools<dim>::interpolate_boundary_values (const DoFHandler<dim> &dof,
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
MassMatrix<dim>::MassMatrix (const Function<dim> * const rhs,
			     const Function<dim> * const a) :
		Equation<dim> (1),
		right_hand_side (rhs),
		coefficient (a)   {};



template <int dim>
void MassMatrix<dim>::assemble (dFMatrix            &cell_matrix,
				const FEValues<dim> &fe_values,
				const typename Triangulation<dim>::cell_iterator &) const {
  const dFMatrix       &values    = fe_values.get_shape_values ();
  const vector<double> &weights   = fe_values.get_JxW_values ();

  if (coefficient != 0)
    {
      vector<double> coefficient_values (fe_values.n_quadrature_points);
      coefficient->value_list (fe_values.get_quadrature_points(),
			       coefficient_values);
      for (unsigned int point=0; point<fe_values.n_quadrature_points; ++point)
	for (unsigned int i=0; i<fe_values.total_dofs; ++i) 
	  for (unsigned int j=0; j<fe_values.total_dofs; ++j)
	    cell_matrix(i,j) += (values(i,point) *
				 values(j,point) *
				 weights[point] *
				 coefficient_values[point]);
    }
  else
    for (unsigned int point=0; point<fe_values.n_quadrature_points; ++point)
      for (unsigned int i=0; i<fe_values.total_dofs; ++i) 
	for (unsigned int j=0; j<fe_values.total_dofs; ++j)
	  cell_matrix(i,j) += (values(i,point) *
			       values(j,point) *
			       weights[point]);
};



template <int dim>
void MassMatrix<dim>::assemble (dFMatrix            &cell_matrix,
				dVector             &rhs,
				const FEValues<dim> &fe_values,
				const Triangulation<dim>::cell_iterator &) const {
  Assert (right_hand_side != 0, ExcNoRHSSelected());
  
  const dFMatrix       &values    = fe_values.get_shape_values ();
  const vector<double> &weights   = fe_values.get_JxW_values ();
  vector<double>        rhs_values (fe_values.n_quadrature_points);
  right_hand_side->value_list (fe_values.get_quadrature_points(), rhs_values);

  if (coefficient != 0)
    {
      vector<double> coefficient_values (fe_values.n_quadrature_points);
      coefficient->value_list (fe_values.get_quadrature_points(),
			       coefficient_values);
      for (unsigned int point=0; point<fe_values.n_quadrature_points; ++point)
	for (unsigned int i=0; i<fe_values.total_dofs; ++i) 
	  {
	    for (unsigned int j=0; j<fe_values.total_dofs; ++j)
	      cell_matrix(i,j) += (values(i,point) *
				   values(j,point) *
				   weights[point] *
				   coefficient_values[point]);
	    rhs(i) += values(i,point) *
		      rhs_values[point] *
		      weights[point];
	  };
    }
  else
    for (unsigned int point=0; point<fe_values.n_quadrature_points; ++point)
      for (unsigned int i=0; i<fe_values.total_dofs; ++i) 
	{
	  for (unsigned int j=0; j<fe_values.total_dofs; ++j)
	    cell_matrix(i,j) += (values(i,point) *
				 values(j,point) *
				 weights[point]);
	  rhs(i) += values(i,point) *
		    rhs_values[point] *
		    weights[point];
	};
};



template <int dim>
void MassMatrix<dim>::assemble (dVector             &rhs,
				const FEValues<dim> &fe_values,
				const Triangulation<dim>::cell_iterator &) const {
  Assert (right_hand_side != 0, ExcNoRHSSelected());
  
  const dFMatrix       &values    = fe_values.get_shape_values ();
  const vector<double> &weights   = fe_values.get_JxW_values ();
  vector<double>        rhs_values(fe_values.n_quadrature_points);
  right_hand_side->value_list (fe_values.get_quadrature_points(), rhs_values);

  for (unsigned int point=0; point<fe_values.n_quadrature_points; ++point)
    for (unsigned int i=0; i<fe_values.total_dofs; ++i) 
      rhs(i) += values(i,point) *
		rhs_values[point] *
		weights[point];
};





template <int dim>
LaplaceMatrix<dim>::LaplaceMatrix (const Function<dim> * const rhs,
				   const Function<dim> * const a) :
		Equation<dim> (1),
		right_hand_side (rhs),
		coefficient (a) {};


template <int dim>
void LaplaceMatrix<dim>::assemble (dFMatrix            &cell_matrix,
				   dVector             &rhs,
				   const FEValues<dim> &fe_values,
				   const Triangulation<dim>::cell_iterator &) const {
  Assert (right_hand_side != 0, ExcNoRHSSelected());
  
  const vector<vector<Point<dim> > >&gradients = fe_values.get_shape_grads ();
  const dFMatrix       &values    = fe_values.get_shape_values ();
  vector<double>        rhs_values(fe_values.n_quadrature_points);
  const vector<double> &weights   = fe_values.get_JxW_values ();
  right_hand_side->value_list (fe_values.get_quadrature_points(), rhs_values);

  if (coefficient != 0)
    {
      vector<double> coefficient_values(fe_values.n_quadrature_points);
      coefficient->value_list (fe_values.get_quadrature_points(),
			       coefficient_values);
      for (unsigned int point=0; point<fe_values.n_quadrature_points; ++point)
	for (unsigned int i=0; i<fe_values.total_dofs; ++i) 
	  {
	    for (unsigned int j=0; j<fe_values.total_dofs; ++j)
	      cell_matrix(i,j) += (gradients[i][point] *
				   gradients[j][point]) *
				  weights[point] *
				  coefficient_values[point];
	    rhs(i) += values(i,point) *
		      rhs_values[point] *
		      weights[point];
	  };
    }
  else
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



template <int dim>
void LaplaceMatrix<dim>::assemble (dFMatrix            &cell_matrix,
				   const FEValues<dim> &fe_values,
				   const Triangulation<dim>::cell_iterator &) const {
  const vector<vector<Point<dim> > >&gradients = fe_values.get_shape_grads ();
  const vector<double> &weights   = fe_values.get_JxW_values ();
   
  if (coefficient != 0)
    {
      vector<double> coefficient_values(fe_values.n_quadrature_points);
      coefficient->value_list (fe_values.get_quadrature_points(),
			       coefficient_values);
      for (unsigned int point=0; point<fe_values.n_quadrature_points; ++point)
	for (unsigned int i=0; i<fe_values.total_dofs; ++i) 
	  for (unsigned int j=0; j<fe_values.total_dofs; ++j)
	    cell_matrix(i,j) += (gradients[i][point] *
				 gradients[j][point]) *
				weights[point] *
				coefficient_values[point];
    }
  else
    for (unsigned int point=0; point<fe_values.n_quadrature_points; ++point)
      for (unsigned int i=0; i<fe_values.total_dofs; ++i) 
	for (unsigned int j=0; j<fe_values.total_dofs; ++j)
	  cell_matrix(i,j) += (gradients[i][point] *
			       gradients[j][point]) *
			      weights[point];
};



template <int dim>
void LaplaceMatrix<dim>::assemble (dVector             &rhs,
				   const FEValues<dim> &fe_values,
				   const Triangulation<dim>::cell_iterator &) const {
  Assert (right_hand_side != 0, ExcNoRHSSelected());
  
  const dFMatrix       &values    = fe_values.get_shape_values ();
  const vector<double> &weights   = fe_values.get_JxW_values ();
  vector<double>        rhs_values(fe_values.n_quadrature_points);
  right_hand_side->value_list (fe_values.get_quadrature_points(), rhs_values);
   
  for (unsigned int point=0; point<fe_values.n_quadrature_points; ++point)
    for (unsigned int i=0; i<fe_values.total_dofs; ++i) 
      rhs(i) += values(i,point) *
		rhs_values[point] *
		weights[point];
};








template class MatrixCreator<1>;
template class MatrixCreator<2>;
template class MatrixTools<1>;
template class MatrixTools<2>;
template class MassMatrix<1>;
template class MassMatrix<2>;
template class LaplaceMatrix<1>;
template class LaplaceMatrix<2>;

