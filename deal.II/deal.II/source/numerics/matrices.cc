/*      $Id$                 */

#include <basic/function.h>
#include <grid/dof.h>
#include <grid/dof_accessor.h>
#include <grid/tria_iterator.h>
#include <fe/quadrature.h>
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
template class MassMatrix<1>;
template class MassMatrix<2>;
template class LaplaceMatrix<1>;
template class LaplaceMatrix<2>;

