/*      $Id$                 */

#include <base/function.h>
#include <grid/dof.h>
#include <grid/dof_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/geometry_info.h>
#include <base/quadrature.h>
#include <fe/fe.h>
#include <fe/fe_values.h>
#include <numerics/matrices.h>
#include <numerics/assembler.h>
#include <lac/dvector.h>
#include <lac/dsmatrix.h>

#include <algorithm>
#include <set>
#include <cmath>



template <int dim>
void MatrixCreator<dim>::create_mass_matrix (const DoFHandler<dim>    &dof,
					     const FiniteElement<dim> &fe,
					     const Quadrature<dim>    &q,
					     const Boundary<dim>      &boundary,
					     dSMatrix                 &matrix,
					     const Function<dim> * const a) {
  dVector dummy;    // no entries, should give an error if accessed
  UpdateFlags update_flags = update_JxW_values;
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
void MatrixCreator<dim>::create_mass_matrix (const DoFHandler<dim>    &dof,
					     const FiniteElement<dim> &fe,
					     const Boundary<dim>      &boundary,
					     dSMatrix                 &matrix) {
  const unsigned int total_dofs = fe.total_dofs;
  
  dFMatrix    local_mass_matrix (total_dofs, total_dofs);
  vector<int> dofs_on_this_cell (total_dofs);
  
  DoFHandler<dim>::active_cell_iterator cell = dof.begin_active(),
					endc = dof.end();
  for (; cell!=endc; ++cell) 
    {
      cell->get_dof_indices (dofs_on_this_cell);
      fe.get_local_mass_matrix (cell, boundary, local_mass_matrix);
      
      for (unsigned int i=0; i<total_dofs; ++i)
	for (unsigned int j=0; j<total_dofs; ++j)
	  matrix.add (dofs_on_this_cell[i], dofs_on_this_cell[j],
		      local_mass_matrix(i,j));
    };
};



#if deal_II_dimension == 1

template <>
void MatrixCreator<1>::create_boundary_mass_matrix (const DoFHandler<1>    &,
						    const FiniteElement<1> &,
						    const Quadrature<0>    &,
						    const Boundary<1>      &,
						    dSMatrix               &,
						    const FunctionMap      &,
						    dVector                &,
						    vector<int>            &,
						    const Function<1>      *) {
  Assert (false, ExcNotImplemented());
};

#endif



template <int dim>
void MatrixCreator<dim>::create_boundary_mass_matrix (const DoFHandler<dim>    &dof,
						      const FiniteElement<dim> &fe,
						      const Quadrature<dim-1>  &q,
						      const Boundary<dim>      &boundary,
						      dSMatrix                 &matrix,
						      const FunctionMap        &rhs,
						      dVector                  &rhs_vector,
						      vector<int>              &dof_to_boundary_mapping,
						      const Function<dim>      *a) {
  Assert (matrix.n() == dof.n_boundary_dofs(rhs), ExcInternalError());
  Assert (matrix.n() == matrix.m(), ExcInternalError());
  Assert (matrix.n() == rhs_vector.size(), ExcInternalError());
  Assert (rhs.size() != 0, ExcInternalError());
  Assert (dof.get_fe() == fe, ExcInternalError());
  Assert (dof_to_boundary_mapping.size() == dof.n_dofs(), ExcInternalError());
  Assert (*max_element(dof_to_boundary_mapping.begin(),dof_to_boundary_mapping.end()) ==
	  (signed int)matrix.n()-1,
	  ExcInternalError());
  
  const unsigned int dofs_per_cell = fe.total_dofs,
		     dofs_per_face = fe.dofs_per_face;
  dFMatrix cell_matrix(dofs_per_cell, dofs_per_cell);
  dVector  cell_vector(dofs_per_cell);
  
  
  UpdateFlags update_flags = UpdateFlags (update_JxW_values | update_q_points);
  FEFaceValues<dim> fe_values (fe, q, update_flags);
  
  DoFHandler<dim>::active_cell_iterator cell = dof.begin_active (),
					endc = dof.end ();
  for (; cell!=endc; ++cell)
    for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
				       // check if this face is on that part of
				       // the boundary we are interested in
      if (rhs.find(cell->face(face)->boundary_indicator()) != rhs.end())
	{
	  cell_matrix.clear ();
	  cell_vector.clear ();
	  
	  fe_values.reinit (cell, face, boundary);

	  const dFMatrix       &values    = fe_values.get_shape_values ();
	  const vector<double> &weights   = fe_values.get_JxW_values ();
	  vector<double>        rhs_values (fe_values.n_quadrature_points);
	  rhs.find(cell->face(face)->boundary_indicator())
	    ->second->value_list (fe_values.get_quadrature_points(), rhs_values);
	  
	  if (a != 0)
	    {
	      vector<double> coefficient_values (fe_values.n_quadrature_points);
	      a->value_list (fe_values.get_quadrature_points(), coefficient_values);
	      for (unsigned int point=0; point<fe_values.n_quadrature_points; ++point)
		for (unsigned int i=0; i<fe_values.total_dofs; ++i) 
		  {
		    for (unsigned int j=0; j<fe_values.total_dofs; ++j)
		      cell_matrix(i,j) += (values(i,point) *
					   values(j,point) *
					   weights[point] *
					   coefficient_values[point]);
		    cell_vector(i) += values(i,point) *
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
		  cell_vector(i) += values(i,point) *
				    rhs_values[point] *
				    weights[point];
		};


					   // now transfer cell matrix and vector
					   // to the whole boundary matrix
					   //
					   // in the following: dof[i] holds the
					   // global index of the i-th degree of
					   // freedom on the present cell. If it
					   // is also a dof on the boundary, it
					   // must be a nonzero entry in the
					   // dof_to_boundary_mapping and then
					   // the boundary index of this dof is
					   // dof_to_boundary_mapping[dof[i]].
					   //
					   // if dof[i] is not on the boundary,
					   // it should be zero on the boundary
					   // therefore on all quadrature
					   // points and finally all of its
					   // entries in the cell matrix and
					   // vector should be zero. If not, we
					   // throw an error (note: because of
					   // the evaluation of the shape
					   // functions only up to machine
					   // precision, the term "must be zero"
					   // really should mean: "should be
					   // very small". since this is only an
					   // assertion and not part of the
					   // code, we may choose "very small"
					   // quite arbitrarily)
					   //
					   // the main problem here is that the
					   // matrix or vector entry should also
					   // be zero if the degree of freedom
					   // dof[i] is on the boundary, but not
					   // on the present face, i.e. on
					   // another face of the same cell also
					   // on the boundary. We can therefore
					   // not rely on the
					   // dof_to_boundary_mapping[dof[i]]
					   // being !=-1, we really have to
					   // determine whether dof[i] is a
					   // dof on the present face. We do so
					   // by getting the dofs on the
					   // face into #dofs_on_face_vector#,
					   // a vector as always. Usually,
					   // searching in a vector is
					   // inefficient, so we copy the dofs
					   // into a set, which enables binary
					   // searches.
	  vector<int> dofs (dofs_per_cell);
	  cell->get_dof_indices (dofs);

	  vector<int> dofs_on_face_vector (dofs_per_face);
	  cell->face(face)->get_dof_indices (dofs_on_face_vector);
	  set<int> dofs_on_face (dofs_on_face_vector.begin(),
				 dofs_on_face_vector.end());
#ifdef DEBUG
					   // in debug mode: compute an element
					   // in the matrix which is
					   // guaranteed to belong to a boundary
					   // dof. We do this to check that the
					   // entries in the cell matrix are
					   // guaranteed to be zero if the
					   // respective dof is not on the
					   // boundary. Since because of
					   // round-off, the actual
					   // value of the matrix entry may be
					   // only close to zero, we assert that
					   // it is small relative to an element
					   // which is guaranteed to be nonzero.
					   // (absolute smallness does not
					   // suffice since the size of the
					   // domain scales in here)
					   //
					   // for this purpose we seek the
					   // diagonal of the matrix, where there
					   // must be an element belonging to
					   // the boundary. we take the maximum
					   // diagonal entry.
	  double max_diag_entry = 0;
	  for (unsigned int i=0; i<dofs_per_cell; ++i)
	    if (fabs(cell_matrix(i,i)) > max_diag_entry)
	      max_diag_entry = fabs(cell_matrix(i,i));
#endif  
	  
	  for (unsigned int i=0; i<dofs_per_cell; ++i)
	    for (unsigned int j=0; j<dofs_per_cell; ++j)
	      if ((dofs_on_face.find(dofs[i]) != dofs_on_face.end()) &&
		  (dofs_on_face.find(dofs[j]) != dofs_on_face.end()))
		matrix.add(dof_to_boundary_mapping[dofs[i]],
			   dof_to_boundary_mapping[dofs[j]],
			   cell_matrix(i,j));
	      else
		{
						   // compare here for relative
						   // smallness
		  Assert (fabs(cell_matrix(i,j)) <= 1e-10 * max_diag_entry,
			  ExcInternalError ());
		};
	  
	  for (unsigned int j=0; j<dofs_per_cell; ++j)
	    if (dofs_on_face.find(dofs[j]) != dofs_on_face.end())
	      rhs_vector(dof_to_boundary_mapping[dofs[j]]) += cell_vector(j);
	    else
	      {
						   // compare here for relative
						   // smallness
		Assert (fabs(cell_vector(j)) <= 1e-10 * max_diag_entry,
			ExcInternalError());
	      };
	};
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
  Assert (matrix.n() == matrix.m(),
	  ExcDimensionsDontMatch(matrix.n(), matrix.m()));
  Assert (matrix.n() == right_hand_side.size(),
	  ExcDimensionsDontMatch(matrix.n(), right_hand_side.size()));
  Assert (matrix.n() == solution.size(),
	  ExcDimensionsDontMatch(matrix.n(), solution.size()));
				   // if no boundary values are to be applied
				   // simply return
  if (boundary_values.size() == 0)
    return;
  
  
  map<int,double>::const_iterator dof  = boundary_values.begin(),
				  endd = boundary_values.end();
  const unsigned int n_dofs             = matrix.m();
  const dSMatrixStruct &sparsity        = matrix.get_sparsity_pattern();
  const unsigned int *sparsity_rowstart = sparsity.get_rowstart_indices();
  const int          *sparsity_colnums  = sparsity.get_column_numbers();

  for (; dof != endd; ++dof)
    {
      const int dof_number = dof->first;
				       // for each boundary dof:
      
				       // set entries of this line
				       // to zero except for the diagonal
				       // entry. Note that the diagonal
				       // entry is always the first one
				       // for square matrices, i.e.
				       // we shall not set
				       // matrix.global_entry(
				       //     sparsity_rowstart[dof.first])
      const unsigned int last = sparsity_rowstart[dof_number+1];
      for (unsigned int j=sparsity_rowstart[dof_number]+1; j<last; ++j)
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
      if (matrix.diag_element(dof_number) != 0.0)
	new_rhs = right_hand_side(dof_number)
		= dof->second * matrix.diag_element(dof_number);
      else
	{
	  double first_diagonal_entry = 1;
	  for (unsigned int i=0; i<n_dofs; ++i)
	    if (matrix.diag_element(i) != 0)
	      {
		first_diagonal_entry = matrix.diag_element(i);
		break;
	      };

					   // use the dSMatrix:: as a workaround
					   // for a bug in egcs
	  matrix.dSMatrix::set(dof_number, dof_number,
			       first_diagonal_entry);
	  new_rhs = right_hand_side(dof_number)
		  = dof->second * first_diagonal_entry;
	};
      
				       // store the only nonzero entry
				       // of this line for the Gauss
				       // elimination step
      const double diagonal_entry = matrix.diag_element(dof_number);

				       // do the Gauss step
      for (unsigned int row=0; row<n_dofs; ++row) 
	{
					   // we need not handle the
					   // row we have already cleared
	  if ((signed int)row == dof_number)
	    continue;

					   // check whether the line has
					   // an entry in the row corresponding
					   // to the dof presently worked with.
					   // note again: the first entry is
					   // the diagonal entry which we
					   // cannot be interested in; following
					   // are the other entries in sorted
					   // order, so we can use a binary
					   // search
					   //
					   // if this row contains an element
					   // for this dof, *p==dof_number
	  const int * p = lower_bound(&sparsity_colnums[sparsity_rowstart[row]+1],
				      &sparsity_colnums[sparsity_rowstart[row+1]],
				      dof_number);
					   // check whether this line has
					   // an entry in the regarding column
					   // (check for ==dof_number and
					   // != next_row, since if
					   // row==dof_number-1, *p is a
					   // past-the-end pointer but points
					   // to dof_number anyway...)
	  if ((*p == dof_number) &&
	      (p != &sparsity_colnums[sparsity_rowstart[row+1]]))
					     // this line has an entry
					     // in the regarding column
	    {
	      const unsigned int global_entry
		= (p - &sparsity_colnums[sparsity_rowstart[0]]);
	      
					       // correct right hand side
	      right_hand_side(row) -= matrix.global_entry(global_entry) /
				      diagonal_entry * new_rhs;
	      
					       // set matrix entry to zero
	      matrix.global_entry(global_entry) = 0.;
	    };
	};
      
      
      
				       // preset solution vector
      solution(dof_number) = dof->second;
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
  const unsigned int total_dofs = fe_values.total_dofs,
		     n_q_points = fe_values.n_quadrature_points;

  Assert (cell_matrix.n() == total_dofs,
	  ExcWrongSize(cell_matrix.n(), total_dofs));
  Assert (cell_matrix.m() == total_dofs,
	  ExcWrongSize(cell_matrix.m(), total_dofs));
  Assert (cell_matrix.all_zero(), ExcObjectNotEmpty());
  
  const dFMatrix       &values    = fe_values.get_shape_values ();
  const vector<double> &weights   = fe_values.get_JxW_values ();

  
  if (coefficient != 0)
    {
      vector<double> coefficient_values (fe_values.n_quadrature_points);
      coefficient->value_list (fe_values.get_quadrature_points(),
			       coefficient_values);
      for (unsigned int i=0; i<total_dofs; ++i) 
	for (unsigned int j=0; j<total_dofs; ++j)
	  for (unsigned int point=0; point<n_q_points; ++point)
	    cell_matrix(i,j) += (values(i,point) *
				 values(j,point) *
				 weights[point] *
				 coefficient_values[point]);
    }
  else
    for (unsigned int i=0; i<total_dofs; ++i) 
      for (unsigned int j=0; j<total_dofs; ++j)
	for (unsigned int point=0; point<n_q_points; ++point)
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

  const unsigned int total_dofs = fe_values.total_dofs,
		     n_q_points = fe_values.n_quadrature_points;

  Assert (cell_matrix.n() == total_dofs,
	  ExcWrongSize(cell_matrix.n(), total_dofs));
  Assert (cell_matrix.m() == total_dofs,
	  ExcWrongSize(cell_matrix.m(), total_dofs));
  Assert (rhs.size() == total_dofs,
	  ExcWrongSize(rhs.size(), total_dofs));
  Assert (cell_matrix.all_zero(), ExcObjectNotEmpty());
  Assert (rhs.all_zero(), ExcObjectNotEmpty());

  const dFMatrix       &values    = fe_values.get_shape_values ();
  const vector<double> &weights   = fe_values.get_JxW_values ();
  vector<double>        rhs_values (fe_values.n_quadrature_points);
  right_hand_side->value_list (fe_values.get_quadrature_points(), rhs_values);

  if (coefficient != 0)
    {
      vector<double> coefficient_values (n_q_points);
      coefficient->value_list (fe_values.get_quadrature_points(),
			       coefficient_values);
      for (unsigned int point=0; point<n_q_points; ++point)
	for (unsigned int i=0; i<total_dofs; ++i) 
	  {
	    for (unsigned int j=0; j<total_dofs; ++j)
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
    for (unsigned int point=0; point<n_q_points; ++point)
      for (unsigned int i=0; i<total_dofs; ++i) 
	{
	  for (unsigned int j=0; j<total_dofs; ++j)
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

  const unsigned int total_dofs = fe_values.total_dofs,
		     n_q_points = fe_values.n_quadrature_points;

  Assert (rhs.size() == total_dofs, ExcWrongSize(rhs.size(), total_dofs));
  Assert (rhs.all_zero(), ExcObjectNotEmpty());

  const dFMatrix       &values    = fe_values.get_shape_values ();
  const vector<double> &weights   = fe_values.get_JxW_values ();
  vector<double>        rhs_values(fe_values.n_quadrature_points);
  right_hand_side->value_list (fe_values.get_quadrature_points(), rhs_values);

  for (unsigned int point=0; point<n_q_points; ++point)
    for (unsigned int i=0; i<total_dofs; ++i) 
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
  
  const unsigned int total_dofs = fe_values.total_dofs,
		     n_q_points = fe_values.n_quadrature_points;

  Assert (cell_matrix.n() == total_dofs,
	  ExcWrongSize(cell_matrix.n(), total_dofs));
  Assert (cell_matrix.m() == total_dofs,
	  ExcWrongSize(cell_matrix.m(), total_dofs));
  Assert (rhs.size() == total_dofs,
	  ExcWrongSize(rhs.size(), total_dofs));
  Assert (cell_matrix.all_zero(), ExcObjectNotEmpty());
  Assert (rhs.all_zero(), ExcObjectNotEmpty());

  const vector<vector<Point<dim> > >&gradients = fe_values.get_shape_grads ();
  const dFMatrix       &values    = fe_values.get_shape_values ();
  vector<double>        rhs_values(fe_values.n_quadrature_points);
  const vector<double> &weights   = fe_values.get_JxW_values ();
  right_hand_side->value_list (fe_values.get_quadrature_points(), rhs_values);

  if (coefficient != 0)
    {
      vector<double> coefficient_values(n_q_points);
      coefficient->value_list (fe_values.get_quadrature_points(),
			       coefficient_values);
      for (unsigned int point=0; point<n_q_points; ++point)
	for (unsigned int i=0; i<total_dofs; ++i) 
	  {
	    for (unsigned int j=0; j<total_dofs; ++j)
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
    for (unsigned int point=0; point<n_q_points; ++point)
      for (unsigned int i=0; i<total_dofs; ++i) 
	{
	  for (unsigned int j=0; j<total_dofs; ++j)
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
  const unsigned int total_dofs = fe_values.total_dofs,
		     n_q_points = fe_values.n_quadrature_points;

  Assert (cell_matrix.n() == total_dofs,
	  ExcWrongSize(cell_matrix.n(), total_dofs));
  Assert (cell_matrix.m() == total_dofs,
	  ExcWrongSize(cell_matrix.m(), total_dofs));
  Assert (cell_matrix.all_zero(), ExcObjectNotEmpty());
  
  const vector<vector<Point<dim> > >&gradients = fe_values.get_shape_grads ();
  const vector<double> &weights   = fe_values.get_JxW_values ();
   
  if (coefficient != 0)
    {
      vector<double> coefficient_values(n_q_points);
      coefficient->value_list (fe_values.get_quadrature_points(),
			       coefficient_values);
      for (unsigned int point=0; point<n_q_points; ++point)
	for (unsigned int i=0; i<total_dofs; ++i) 
	  for (unsigned int j=0; j<total_dofs; ++j)
	    cell_matrix(i,j) += (gradients[i][point] *
				 gradients[j][point]) *
				weights[point] *
				coefficient_values[point];
    }
  else
    for (unsigned int point=0; point<n_q_points; ++point)
      for (unsigned int i=0; i<total_dofs; ++i) 
	for (unsigned int j=0; j<total_dofs; ++j)
	  cell_matrix(i,j) += (gradients[i][point] *
			       gradients[j][point]) *
			      weights[point];
};



template <int dim>
void LaplaceMatrix<dim>::assemble (dVector             &rhs,
				   const FEValues<dim> &fe_values,
				   const Triangulation<dim>::cell_iterator &) const {
  Assert (right_hand_side != 0, ExcNoRHSSelected());

  const unsigned int total_dofs = fe_values.total_dofs,
		     n_q_points = fe_values.n_quadrature_points;

  Assert (rhs.size() == total_dofs, ExcWrongSize(rhs.size(), total_dofs));
  Assert (rhs.all_zero(), ExcObjectNotEmpty());

  const dFMatrix       &values    = fe_values.get_shape_values ();
  const vector<double> &weights   = fe_values.get_JxW_values ();
  vector<double>        rhs_values(fe_values.n_quadrature_points);
  right_hand_side->value_list (fe_values.get_quadrature_points(), rhs_values);
   
  for (unsigned int point=0; point<n_q_points; ++point)
    for (unsigned int i=0; i<total_dofs; ++i) 
      rhs(i) += values(i,point) *
		rhs_values[point] *
		weights[point];
}

template<int dim> void
MatrixCreator<dim>::create_interpolation_matrix(const FiniteElement<dim> &high,
						const FiniteElement<dim> &low,
						dFMatrix& result)
{
  result.reinit (low.total_dofs, high.total_dofs);

  vector<Point<dim> > unit_support_points (high.total_dofs);
  high.get_unit_support_points (unit_support_points);
  
  for (unsigned int i=0; i<low.total_dofs; ++i)
    for (unsigned int j=0; j<high.total_dofs; ++j)
      result(i,j) = high.shape_value (j, unit_support_points[i]);
}




template class MatrixCreator<deal_II_dimension>;
template class MatrixTools<deal_II_dimension>;
template class MassMatrix<deal_II_dimension>;
template class LaplaceMatrix<deal_II_dimension>;

