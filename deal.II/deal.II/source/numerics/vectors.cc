/* $Id$ */
/* Copyright W. Bangerth, University of Heidelberg, 1998 */



#include <base/function.h>
#include <base/tensorfunction.h>
#include <grid/dof.h>
#include <grid/dof_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/dof_constraints.h>
#include <fe/fe.h>
#include <fe/fe_values.h>
#include <base/quadrature.h>
#include <numerics/assembler.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <lac/vector.h>
#include <lac/sparsematrix.h>
#include <lac/precondition.h>
#include <lac/solver_cg.h>
#include <lac/vector_memory.h>

#include <numeric>
#include <algorithm>
#include <cmath>




inline double sqr (const double x) {
  return x*x;
};


template <int dim>
inline double sqr_point (const Tensor<1,dim> &p) {
  return p * p;
};






template <int dim>
void VectorTools<dim>::interpolate (const DoFHandler<dim>    &dof,
				    const Function<dim>      &function,
				    Vector<double>           &vec)
{
  const FiniteElement<dim> &fe = dof.get_fe();
  
  DoFHandler<dim>::active_cell_iterator cell = dof.begin_active(),
					endc = dof.end();
  vector<int>         dofs_on_cell (fe.total_dofs);
  vector<double>      dof_values_on_cell (fe.total_dofs);
  vector<Point<dim> > support_points (fe.total_dofs);
  for (; cell!=endc; ++cell) 
    {
				       // for each cell:
				       // get location of finite element
				       // off-points
      fe.get_support_points (cell, support_points);
				       // get function values at these points
      function.value_list (support_points, dof_values_on_cell);
				       // get indices of the dofs on this cell
      cell->get_dof_indices (dofs_on_cell);
				       // distribute function values to the
				       // whole vector
      for (unsigned int i=0; i<fe.total_dofs; ++i)
	vec(dofs_on_cell[i]) = dof_values_on_cell[i];
    };
};



template <int dim>
void VectorTools<dim>::interpolate (const DoFHandler<dim>    &dof,
				    const VectorFunction<dim>&vectorfunction,
				    Vector<double>           &vec)
{
  const FiniteElement<dim> &fe = dof.get_fe();
  
  DoFHandler<dim>::active_cell_iterator cell = dof.begin_active(),
					endc = dof.end();

				   // For FESystems many of the
				   // unit_support_points will
				   // appear multiply, as a point
				   // may be unit_support_point
				   // for several of the components
				   // of the system.
				   // The following is rather
				   // complicated as it is
				   // avoided to evaluate
				   // the vectorfunction multiply at
				   // the same point on a cell.
  vector<Point<dim> > unit_support_points (fe.total_dofs);
  fe.get_unit_support_points(unit_support_points);

				   // The following works well
				   // if #dofs_per_cell<=1# as then
				   // the multiple support_points
				   // are placed one after another.

				   // find the support points 
				   // on a cell that
				   // are multiply mentioned in 
				   // #unit_support_points#.
				   // Mark the first representative
				   // of each multiply mentioned
				   // support point by setting
				   // #true# in the boolean vector 
				   // #is_representative_point#.
//   vector<bool>  is_representative_point(fe.total_dofs, false);
//   is_representative_point[0]=true;
//   unsigned int n_rep_points=1;
//   for (unsigned int last_rep_point=0, i=1; i<fe.total_dofs; ++i)
//     {
//       if (unit_support_points[i] != unit_support_points[last_rep_point])
// 	{
// 	  is_representative_point[i] = true;
// 	  last_rep_point=i;
// 	  ++n_rep_points;
// 	}
//    };

//   vector<int>         dofs_on_cell (fe.total_dofs);
//   vector<Point<dim> > support_points (fe.total_dofs);

//   vector<Point<dim> > rep_points (n_rep_points);
//   vector<Vector<double> > function_values_at_rep_points (
//     n_rep_points, Vector<double>(fe.n_components));

//   for (; cell!=endc; ++cell)
//     {
// 				       // for each cell:
// 				       // get location of finite element
// 				       // off-points (support_points)
//       fe.get_support_points (cell, support_points);

// 				       // pick out the representative
// 				       // support points
//       unsigned int j=0;
//       for (unsigned int i=0; i<fe.total_dofs; ++i)
// 	if (is_representative_point[i])
// 	  rep_points[j++]=support_points[i];
//       Assert(j == n_rep_points, ExcInternalError());

// 				       // get function values at these points
//       vectorfunction.value_list (rep_points, function_values_at_rep_points);
  
// 					     // get indices of the dofs on this cell
//       cell->get_dof_indices (dofs_on_cell);

// 				       // distribute function values to the
// 				       // whole vector
//       int last_rep_point = -1;
// 				       // it holds `is_representative_point[0]=true'
// 				       // therefore the first #last_rep_point# is 0
// 				       // and we need to start with
// 				       // `last_rep_point = -1'
//       for (unsigned int i=0; i<fe.total_dofs; ++i)
// 	{
// 	  if (is_representative_point[i])
// 	    ++last_rep_point;

// 	  const unsigned int component
// 	    = fe.system_to_component_index(i).first;
// 	  vec(dofs_on_cell[i])
// 	    = function_values_at_rep_points[last_rep_point](component);
// 	} 
//     }

				   // The following is more general.
				   // It also works if #dofs_per_cell>1#,
				   // i.e. it is usable also for systems
				   // including
				   // FEQ3, FEQ4, FEDG_Qx.

				   // Find the support points 
				   // on a cell that
				   // are multiply mentioned in 
				   // #unit_support_points#.
				   // Mark the first representative
				   // of each multiply mentioned
				   // support point by appending its
				   // dof index to #dofs_of_rep_points#.
				   // Each multiple point gets to know
				   // the dof index of its representative
				   // point by the #dof_to_rep_dof_table#.

				   // the following vector collects all dofs i,
				   // 0<=i<fe.total_dofs, for that
				   // unit_support_points[i] 
				   // is a representative one.
  vector<unsigned int> dofs_of_rep_points;
				   // the following table converts a dof i
				   // to the dof of the representative
				   // point.
  vector<unsigned int> dof_to_rep_dof_table;
  unsigned int n_rep_points=0;
  for (unsigned int i=0; i<fe.total_dofs; ++i)
    {
      bool representative=true;
				       // the following loop is looped
				       // the other way round to get
				       // the minimal effort of
				       // O(fe.total_dofs) for multiple
				       // support points that are placed
				       // one after the other.
      for (unsigned int j=dofs_of_rep_points.size(); j>0; --j)
	if (unit_support_points[i] 
	    == unit_support_points[dofs_of_rep_points[j-1]])
	  {
	    dof_to_rep_dof_table.push_back(j-1);
	    representative=false;
	    break;
	  }
      
      if (representative)
	{
	  dofs_of_rep_points.push_back(i);
	  dof_to_rep_dof_table.push_back(i);
	  ++n_rep_points;
	}
    }
  Assert(dofs_of_rep_points.size()==n_rep_points, ExcInternalError());
  Assert(dof_to_rep_dof_table.size()==fe.total_dofs, ExcInternalError());

  cout << "n_rep_points=" << n_rep_points << endl;

  vector<int>         dofs_on_cell (fe.total_dofs);
  vector<Point<dim> > support_points (fe.total_dofs);

  vector<Point<dim> > rep_points (n_rep_points);
  vector<Vector<double> > function_values_at_rep_points (
    n_rep_points, Vector<double>(fe.n_components));

  for (; cell!=endc; ++cell)
    {
				       // for each cell:
				       // get location of finite element
				       // off-points (support_points)
      fe.get_support_points (cell, support_points);
      
				       // pick out the representative
				       // support points
      for (unsigned int j=0; j<dofs_of_rep_points.size(); ++j)
	rep_points[j]=support_points[dofs_of_rep_points[j]];

				       // get function values at these points
      vectorfunction.value_list (rep_points, function_values_at_rep_points);

				       // get indices of the dofs on this cell
      cell->get_dof_indices (dofs_on_cell);

				       // distribute the function values to
				       // the global vector
      for (unsigned int i=0; i<fe.total_dofs; ++i)
	{
	  const unsigned int component
	    = fe.system_to_component_index(i).first;
	  const unsigned int rep_dof=dof_to_rep_dof_table[i];
	  vec(dofs_on_cell[i])
	    = function_values_at_rep_points[rep_dof](component);
	}
    }
}


template <int dim> void
VectorTools<dim>::interpolate(const DoFHandler<dim>    &high_dof,
			      const DoFHandler<dim>    &low_dof,
			      const FullMatrix<double>        &transfer,
			      const Vector<double>            &high,
			      Vector<double>                  &low)
{
  Vector<double> cell_high(high_dof.get_fe().total_dofs);
  Vector<double> cell_low(low_dof.get_fe().total_dofs);
  
  DoFHandler<dim>::active_cell_iterator h = high_dof.begin_active();
  DoFHandler<dim>::active_cell_iterator l = low_dof.begin_active();
  
  for(; h != high_dof.end(); ++h, ++l)
  {
    h->get_dof_values(high, cell_high);
    transfer.vmult(cell_low, cell_high);
    l->distribute_local_to_global(cell_low, low);
  }
}



#if deal_II_dimension == 1

template <>
void VectorTools<1>::project (const DoFHandler<1>    &,
			      const ConstraintMatrix &,
			      const Quadrature<1>    &,
			      const Function<1>      &,
			      Vector<double>         &,
			      const bool              ,
			      const Quadrature<0>    &,
			      const bool              ) {
				   // this function should easily be implemented
				   // using the template below. However some
				   // changes have to be made since faces don't
				   // exist in 1D. Maybe integrate the creation of
				   // zero boundary values into the
				   // project_boundary_values function?
  Assert (false, ExcNotImplemented());
};

#endif




template <int dim>
void VectorTools<dim>::project (const DoFHandler<dim>    &dof,
				const ConstraintMatrix   &constraints,
				const Quadrature<dim>    &quadrature,
				const Function<dim>      &function,
				Vector<double>           &vec,
				const bool                enforce_zero_boundary,
				const Quadrature<dim-1>  &q_boundary,
				const bool                project_to_boundary_first) {
  const FiniteElement<dim> &fe = dof.get_fe();

				   // make up boundary values
  map<int,double> boundary_values;

  if (enforce_zero_boundary == true) 
				     // no need to project boundary values, but
				     // enforce homogeneous boundary values
				     // anyway
    {
      DoFHandler<dim>::active_face_iterator face = dof.begin_active_face(),
					    endf = dof.end_face();
      vector<int> face_dof_indices (fe.dofs_per_face);
      for (; face!=endf; ++face)
	if (face->at_boundary())
	  {
	    face->get_dof_indices (face_dof_indices);
	    for (unsigned int i=0; i<fe.dofs_per_face; ++i)
					       // enter zero boundary values
					       // for all boundary nodes
	      boundary_values[face_dof_indices[i]] = 0.;
	  };
    }
  else
				     // no homogeneous boundary values
    if (project_to_boundary_first == true)
				       // boundary projection required
      {
					 // set up a list of boundary functions for
					 // the different boundary parts. We want the
					 // #function# to hold on all parts of the
					 // boundary
	FunctionMap boundary_functions;
	for (unsigned char c=0; c<255; ++c)
	  boundary_functions[c] = &function;
	project_boundary_values (dof, boundary_functions, q_boundary,
				 boundary_values);
      };

  
      
				   // set up mass matrix and right hand side
  vec.reinit (dof.n_dofs());
  SparseMatrixStruct sparsity(dof.n_dofs(),
			      dof.n_dofs(),
			      dof.max_couplings_between_dofs());
  dof.make_sparsity_pattern (sparsity);
  constraints.condense (sparsity);
  
  SparseMatrix<double> mass_matrix (sparsity);
  Vector<double> tmp (mass_matrix.n());

				   // try to assemble the mass matrix by exact
				   // integration. if this is not supported,
				   // then use quadrature
  try 
    {
      MatrixCreator<dim>::create_mass_matrix (dof, mass_matrix);
    }
  catch (FiniteElement<dim>::ExcComputationNotUseful)
    {
      MatrixCreator<dim>::create_mass_matrix (dof, quadrature, mass_matrix);
    };
  
  VectorTools<dim>::create_right_hand_side (dof, quadrature, function, tmp);

  constraints.condense (mass_matrix);
  constraints.condense (tmp);
  MatrixTools<dim>::apply_boundary_values (boundary_values,
					   mass_matrix, vec, tmp);

  SolverControl                    control(1000,1e-16);
  PrimitiveVectorMemory<Vector<double> >   memory;
  SolverCG<SparseMatrix<double>,Vector<double> >       cg(control,memory);

  PreconditionRelaxation<SparseMatrix<double>, Vector<double> >
    prec(mass_matrix,
	 &SparseMatrix<double>::template precondition_SSOR<double>,
	 1.2);
				   // solve
  cg.solve (mass_matrix, vec, tmp, prec);
  
				   // distribute solution
  constraints.distribute (vec);
};



template <int dim>
void VectorTools<dim>::create_right_hand_side (const DoFHandler<dim>    &dof,
					       const Quadrature<dim>    &quadrature,
					       const Function<dim>      &rhs,
					       Vector<double>           &rhs_vector) {
  UpdateFlags update_flags = UpdateFlags(update_q_points |
					 update_JxW_values);
  SparseMatrix<double> dummy;
  const Assembler<dim>::AssemblerData data (dof,
					    false, true,
					    dummy, rhs_vector,
					    quadrature, update_flags);
  TriaActiveIterator<dim, Assembler<dim> >
    assembler (const_cast<Triangulation<dim>*>(&dof.get_tria()),
	       dof.get_tria().begin_active()->level(),
	       dof.get_tria().begin_active()->index(),
	       &data);
  MassMatrix<dim> equation(&rhs,0);
  do 
    {
      assembler->assemble (equation);
    }
  while ((++assembler).state() == valid);
};



#if deal_II_dimension == 1

template <>
void
VectorTools<1>::interpolate_boundary_values (const DoFHandler<1> &dof,
					     const FunctionMap &dirichlet_bc,
					     map<int,double>   &boundary_values)
{
  Assert (dirichlet_bc.find(255) == dirichlet_bc.end(),
	  ExcInvalidBoundaryIndicator());

  const FiniteElement<1> &fe = dof.get_fe();
  Assert (fe.dofs_per_vertex == 1, ExcInvalidFE());

				   // check whether boundary values at the
				   // left boundary of the line are requested
  if (dirichlet_bc.find(0) != dirichlet_bc.end())
    {
				       // first find the leftmost active
				       // cell by first traversing the coarse
				       // grid to its left end and then going
				       // to the children
      DoFHandler<1>::cell_iterator leftmost_cell = dof.begin(0);
      while (leftmost_cell->neighbor(0).state() == valid)
	leftmost_cell = leftmost_cell->neighbor(0);

      while (leftmost_cell->has_children())
	leftmost_cell = leftmost_cell->child(0);

				       // now set the value of the leftmost
				       // degree of freedom
      boundary_values[leftmost_cell->vertex_dof_index(0,0)]
	= dirichlet_bc.find(0)->second->operator()(leftmost_cell->vertex(0));
    };

				   // same for the right boundary of
				   // the line are requested
  if (dirichlet_bc.find(1) != dirichlet_bc.end())
    {
				       // first find the leftmost active
				       // cell by first traversing the coarse
				       // grid to its left end and then going
				       // to the children
      DoFHandler<1>::cell_iterator rightmost_cell = dof.last(0);
      while (rightmost_cell->neighbor(1).state() == valid)
	rightmost_cell = rightmost_cell->neighbor(1);

      while (rightmost_cell->has_children())
	rightmost_cell = rightmost_cell->child(1);

				       // now set the value of the rightmost
				       // degree of freedom
      boundary_values[rightmost_cell->vertex_dof_index(1,0)]
	= dirichlet_bc.find(1)->second->operator()(rightmost_cell->vertex(1));
    };
  
};



template <>
void VectorTools<1>::interpolate_boundary_values (const DoFHandler<1> &,
						  const VectorFunctionMap&,
						  map<int,double>&)
{
  Assert (false, ExcNotImplemented());
};

#endif



template <int dim>
void
VectorTools<dim>::interpolate_boundary_values (const DoFHandler<dim> &dof,
					       const FunctionMap     &dirichlet_bc,
					       map<int,double>   &boundary_values) {
  Assert (dirichlet_bc.find(255) == dirichlet_bc.end(),
	  ExcInvalidBoundaryIndicator());

  const FiniteElement<dim> &fe = dof.get_fe();

				   // use two face iterators, since we need
				   // a DoF-iterator for the dof indices, but
				   // a Tria-iterator for the fe object
  DoFHandler<dim>::active_face_iterator face = dof.begin_active_face(),
					endf = dof.end_face();
  
  typename FunctionMap::const_iterator function_ptr;

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
	fe.get_face_support_points (face, dof_locations);
	function_ptr->second->value_list (dof_locations, dof_values);

					 // enter into list
	for (unsigned int i=0; i<face_dofs.size(); ++i)
	  boundary_values[face_dofs[i]] = dof_values[i];
      };
};



template <int dim>
void
VectorTools<dim>::interpolate_boundary_values (const DoFHandler<dim> &dof,
					       const VectorFunctionMap     &dirichlet_bc,
					       map<int,double>   &boundary_values)
{
  Assert (dirichlet_bc.find(255) == dirichlet_bc.end(),
	  ExcInvalidBoundaryIndicator());

  const FiniteElement<dim> &fe = dof.get_fe();

				   // use two face iterators, since we need
				   // a DoF-iterator for the dof indices, but
				   // a Tria-iterator for the fe object
  DoFHandler<dim>::active_face_iterator face = dof.begin_active_face(),
					endf = dof.end_face();
  
  typename VectorFunctionMap::const_iterator function_ptr;

				   // field to store the indices of dofs
  vector<int>         face_dofs (fe.dofs_per_face);
  vector<Point<dim> > dof_locations (face_dofs.size(), Point<dim>());
  vector< Vector<double> > dof_values (fe.dofs_per_face, Vector<double>(fe.n_components));
	
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
	fe.get_face_support_points (face, dof_locations);
	function_ptr->second->value_list (dof_locations, dof_values);

					 // enter into list

	for (unsigned int i=0; i<face_dofs.size(); ++i)
	  {
	    pair<unsigned int, unsigned int>
	      index = fe.face_system_to_component_index(i);
	    double s = dof_values[i](index.first);
	    if (s != HUGE_VAL)
	      boundary_values[face_dofs[i]] = s;
	  }
      }
}



template <int dim>
void
VectorTools<dim>::project_boundary_values (const DoFHandler<dim>    &dof,
					   const FunctionMap        &boundary_functions,
					   const Quadrature<dim-1>  &q,
					   map<int,double>          &boundary_values) {
  vector<int>    dof_to_boundary_mapping;
  dof.map_dof_to_boundary_indices (boundary_functions, dof_to_boundary_mapping);
  
				   // set up sparsity structure
  SparseMatrixStruct sparsity(dof.n_boundary_dofs(boundary_functions),
			      dof.max_couplings_between_boundary_dofs());
  dof.make_boundary_sparsity_pattern (boundary_functions, dof_to_boundary_mapping,
				      sparsity);

				   // note: for three or more dimensions, there
				   // may be constrained nodes on the boundary
				   // in this case the boundary mass matrix has
				   // to be condensed and the solution is to
				   // be distributed afterwards, which is not
				   // yet implemented. The reason for this is
				   // that we cannot simply use the #condense#
				   // family of functions, since the matrices
				   // and vectors do not use the global
				   // numbering but rather the boundary
				   // numbering, i.e. the condense function
				   // needs to use another indirection. There
				   // should be not many technical problems,
				   // but it needs to be implemented
  if (dim<3)
    sparsity.compress();
  else
    Assert (false, ExcNotImplemented());
  

				   // make mass matrix and right hand side
  SparseMatrix<double> mass_matrix(sparsity);
  Vector<double>       rhs(sparsity.n_rows());
  

  MatrixTools<dim>::create_boundary_mass_matrix (dof, q, 
						 mass_matrix, boundary_functions,
						 rhs, dof_to_boundary_mapping);

				   // same thing as above: if dim>=3 we need
				   // to consider constraints
  Assert (dim<3, ExcNotImplemented());

  
  Vector<double> boundary_projection (rhs.size());

  SolverControl                    control(1000, 1e-16);
  PrimitiveVectorMemory<Vector<double> >   memory;
  SolverCG<SparseMatrix<double>,Vector<double> >       cg(control,memory);

  PreconditionRelaxation<SparseMatrix<double>, Vector<double> >
    prec(mass_matrix,
	 &SparseMatrix<double>::template precondition_SSOR<double>,
	 1.2);
				   // solve
  cg.solve (mass_matrix, boundary_projection, rhs, prec);

				   // fill in boundary values
  for (unsigned int i=0; i<dof_to_boundary_mapping.size(); ++i)
    if (dof_to_boundary_mapping[i] != -1)
				       // this dof is on one of the
				       // interesting boundary parts
				       //
				       // remember: #i# is the global dof
				       // number, #dof_to_boundary_mapping[i]#
				       // is the number on the boundary and
				       // thus in the solution vector
      boundary_values[i] = boundary_projection(dof_to_boundary_mapping[i]);
};




template <int dim>
void VectorTools<dim>::integrate_difference (const DoFHandler<dim>    &dof,
					     const Vector<double>     &fe_function,
					     const Function<dim>      &exact_solution,
					     Vector<float>            &difference,
					     const Quadrature<dim>    &q,
					     const NormType           &norm) {
  const FiniteElement<dim> &fe = dof.get_fe();
    
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
      fe_values.reinit (cell);

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
					     // matrix of the trial function
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
	    vector<Tensor<1,dim> >   psi (n_q_points);

					     // in praxi: first compute
					     // exact fe_function vector
	    exact_solution.gradient_list (fe_values.get_quadrature_points(),
					  psi);
	    
					     // then subtract finite element
					     // fe_function
	    if (true) 
	      {
		vector<Tensor<1,dim> > function_grads (n_q_points, Tensor<1,dim>());
		fe_values.get_function_grads (fe_function, function_grads);

		transform (psi.begin(), psi.end(),
			   function_grads.begin(),
			   psi.begin(),
			   minus<Tensor<1,dim> >());
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



template <int dim>
void
VectorTools<dim>::integrate_difference (const DoFHandler<dim>    &dof,
					const Vector<double>     &fe_function,
					const VectorFunction<dim>&exact_solution,
					Vector<float>            &difference,
					const Quadrature<dim>    &q,
					const NormType           &norm,
					const Function<dim>      *weight)
{
  Assert(norm != mean , ExcNotUseful());

  const FiniteElement<dim> &fe = dof.get_fe();
  
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
      fe_values.reinit (cell);
      
      switch (norm)
 	{
 	  case mean:
		break;
 	  case L1_norm:
 	  case L2_norm:
 	  case Linfty_norm:
 	  case H1_norm:
 	  {
 	    const unsigned int n_q_points = q.n_quadrature_points;
 	    vector<Vector<double> >  psi (n_q_points);

 					     // first compute the exact solution
					     // (vectors) at the quadrature points
 	    exact_solution.value_list (fe_values.get_quadrature_points(), psi);
 					     // then subtract finite element
 					     // fe_function
 	    if (true) 
 	      {
 		vector< Vector<double> > function_values (
		  n_q_points, Vector<double>(fe.n_components));

 		fe_values.get_function_values (fe_function, function_values);

		for (unsigned int q=0; q<n_q_points; ++q)
		  psi[q] -= function_values[q];
 	      };	    

					     // for L1_norm, Linfty_norm, L2_norm
					     // and H1_norm take square of the
					     // vectors psi[q]. Afterwards
 					     // for L1_norm and Linfty_norm:
 					     // take square root to get finally
					     // the (euclidean) vector norm.
					     // Use psi_scalar to store the squares
					     // of the vectors or the vector norms
					     // respectively.
 	    vector<double>  psi_scalar (n_q_points);
 	    switch (norm)
	      {
		case mean:
		      break;
		case L1_norm:
 		case Linfty_norm:
 		case L2_norm:
 		case H1_norm:
		      for (unsigned int q=0; q<n_q_points; ++q)
			psi_scalar[q]=psi[q].norm_sqr();
		      
		      if (norm == L1_norm || norm == Linfty_norm)
			transform (psi_scalar.begin(), psi_scalar.end(),
				   psi_scalar.begin(), ptr_fun(sqrt));
 		      break;
 		default:
 		      Assert (false, ExcNotImplemented());
 	      };

					     // now weight the values
					     // at the quadrature points due
					     // to the weighting function
	    if (weight)
	      {
		vector<double> w(n_q_points);
		weight->value_list(fe_values.get_quadrature_points(),w);
		for (unsigned int q=0; q<n_q_points; ++q)
		  psi_scalar[q]*=w[q];
	      }

 					     // ok, now we have the integrand,
 					     // let's compute the integral,
 					     // which is
 					     // sum_j psi_j weight_j JxW_j
 					     // (or |psi_j| or |psi_j|^2
 	    switch (norm)
 	      {
 		case mean:
		      break;      
 		case L1_norm:
 		case L2_norm:
		case H1_norm:
		      diff = inner_product (psi_scalar.begin(), psi_scalar.end(),
					    fe_values.get_JxW_values().begin(),
					    0.0);
		      if (norm == L2_norm)
			diff=sqrt(diff);

 		      break;
 		case Linfty_norm:
 		      diff = *max_element (psi_scalar.begin(), psi_scalar.end());
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
					     // Until now, #diff# includes the
					     // square of the L2_norm.

 					     // same procedure as above, but now
 					     // psi is a vector of Jacobians
					     // i.e. psi is a vector of vectors of
					     // gradients.
 	    const unsigned int n_q_points = q.n_quadrature_points;
 	    vector<vector<Tensor<1,dim> > >   psi (
	      n_q_points, vector<Tensor<1,dim> >(fe.n_components, Tensor<1,dim>()));
	    
 					     // in praxi: first compute
 					     // exact fe_function vector
 	    exact_solution.gradient_list (fe_values.get_quadrature_points(), psi);

 					     // then subtract finite element
 					     // function_grads
 	    if (true) 
 	      {
 		vector<vector<Tensor<1,dim> > > function_grads (
		  n_q_points, vector<Tensor<1,dim> >(fe.n_components, Tensor<1,dim>()));
 		fe_values.get_function_grads (fe_function, function_grads);

		for (unsigned int q=0; q<n_q_points; ++q)
		  for (unsigned int k=0; k<fe.n_components; ++k)
		    psi[q][k] -= function_grads[q][k];
 	      };
 					     // take square of integrand
 	    vector<double> psi_square (psi.size(), 0.0);
 	    for (unsigned int q=0; q<n_q_points; ++q)
	      for (unsigned int k=0; k<fe.n_components; ++k)
		psi_square[q] += sqr_point(psi[q][k]);

					     // now weight the values
					     // at the quadrature points due
					     // to the weighting function
	    if (weight)
	      {
		vector<double> w(n_q_points);
		weight->value_list(fe_values.get_quadrature_points(),w);
		for (unsigned int q=0; q<n_q_points; ++q)
		  psi_square[q]*=w[q];
	      }

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



// explicit instantiations
template VectorTools<deal_II_dimension>;
