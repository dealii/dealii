//----------------------------  matrices.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  matrices.cc  ---------------------------


#include <base/function.h>
#include <base/thread_management.h>
#include <base/multithread_info.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/geometry_info.h>
#include <base/quadrature.h>
#include <fe/fe.h>
#include <fe/fe_values.h>
#include <numerics/matrices.h>
#include <lac/vector.h>
#include <lac/block_vector.h>
#include <lac/sparse_matrix.h>
#include <lac/block_sparse_matrix.h>
#include <fe/mapping_q1.h>

#include <algorithm>
#include <set>
#include <cmath>


// if necessary try to work around a bug in the IBM xlC compiler
#ifdef XLC_WORK_AROUND_STD_BUG
using namespace std;
#endif




template <int dim>
inline
MatrixCreator::IteratorRange<dim>::
IteratorRange (const active_cell_iterator &first,
						  const active_cell_iterator &second)
		:
		first (first),
		second (second)
{};



template <int dim>
inline
MatrixCreator::IteratorRange<dim>::IteratorRange (const iterator_pair &ip)
		:
		first (ip.first),
		second (ip.second)
{};




template <int dim>
void MatrixCreator::create_mass_matrix (const Mapping<dim>       &mapping,
					const DoFHandler<dim>    &dof,
					const Quadrature<dim>    &q,
					SparseMatrix<double>     &matrix,
					const Function<dim> * const coefficient)
{
  const unsigned int n_threads = multithread_info.n_default_threads;
  Threads::ThreadManager thread_manager;

				   // define starting and end point
				   // for each thread
  typedef typename DoFHandler<dim>::active_cell_iterator active_cell_iterator;  
  const std::vector<std::pair<active_cell_iterator,active_cell_iterator> >
    thread_ranges = Threads::split_range<active_cell_iterator> (dof.begin_active(),
								 dof.end(), n_threads);

				   // mutex to synchronise access to
				   // the matrix
  Threads::ThreadMutex mutex;
  
				   // then assemble in parallel
  for (unsigned int thread=0; thread<n_threads; ++thread)
    Threads::spawn (thread_manager,
		    Threads::encapsulate(&MatrixCreator::template create_mass_matrix_1<dim>)
		    .collect_args (mapping, dof, q, matrix, coefficient,
				   thread_ranges[thread], mutex));
  thread_manager.wait ();  
};



// TODO:[RH] extend this function to use vector valued coefficient functions for system elements.
template <int dim>
void MatrixCreator::create_mass_matrix_1 (const Mapping<dim>       &mapping,
					  const DoFHandler<dim>    &dof,
					  const Quadrature<dim>    &q,
					  SparseMatrix<double>     &matrix,
					  const Function<dim> * const coefficient,
					  const IteratorRange<dim>  range,
					  Threads::ThreadMutex     &mutex)
{
  UpdateFlags update_flags = UpdateFlags(update_values | update_JxW_values);
  if (coefficient != 0)
    update_flags = UpdateFlags (update_flags | update_q_points);

  FEValues<dim> fe_values (mapping, dof.get_fe(), q, update_flags);
    
  const unsigned int dofs_per_cell = fe_values.dofs_per_cell,
		     n_q_points    = fe_values.n_quadrature_points;
  const FiniteElement<dim>    &fe  = fe_values.get_fe();
  const unsigned int n_components  = fe.n_components();

  FullMatrix<double>  cell_matrix (dofs_per_cell, dofs_per_cell);
  std::vector<double> coefficient_values (n_q_points);
  
  std::vector<unsigned int> dof_indices (dofs_per_cell);
  
  typename DoFHandler<dim>::active_cell_iterator cell = range.first;
  for (; cell!=range.second; ++cell)
    {
      fe_values.reinit (cell);
      
      cell_matrix.clear ();
      cell->get_dof_indices (dof_indices);
      
      const FullMatrix<double>  &values    = fe_values.get_shape_values ();
      const std::vector<double> &weights   = fe_values.get_JxW_values ();
      
      if (coefficient != 0)
	{
	  coefficient->value_list (fe_values.get_quadrature_points(),
				   coefficient_values);
	  for (unsigned int point=0; point<n_q_points; ++point)
	    for (unsigned int i=0; i<dofs_per_cell; ++i) 
	      for (unsigned int j=0; j<dofs_per_cell; ++j)
		if ((n_components==1) ||
		    (fe.system_to_component_index(i).first ==
		     fe.system_to_component_index(j).first))
		  cell_matrix(i,j) += (values(i,point) *
				       values(j,point) *
				       weights[point] *
				       coefficient_values[point]);
	}
      else
	for (unsigned int point=0; point<n_q_points; ++point)
	  for (unsigned int i=0; i<dofs_per_cell; ++i) 
	    for (unsigned int j=0; j<dofs_per_cell; ++j)
	      if ((n_components==1) ||
		  (fe.system_to_component_index(i).first ==
		   fe.system_to_component_index(j).first))
		cell_matrix(i,j) += (values(i,point) *
				     values(j,point) *
				     weights[point]);

				       // transfer everything into the
				       // global object
      mutex.acquire ();
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	for (unsigned int j=0; j<dofs_per_cell; ++j)
	  if ((n_components==1) ||
	      (fe.system_to_component_index(i).first ==
	       fe.system_to_component_index(j).first))
	    matrix.add (dof_indices[i], dof_indices[j],
			cell_matrix(i,j));
      mutex.release ();
    };
};


template <int dim>
void MatrixCreator::create_mass_matrix (const DoFHandler<dim>    &dof,
					const Quadrature<dim>    &q,
					SparseMatrix<double>     &matrix,
					const Function<dim> * const coefficient)
{
  Assert (DEAL_II_COMPAT_MAPPING, ExcCompatibility("mapping"));
  static const MappingQ1<dim> mapping;
  create_mass_matrix(mapping, dof, q, matrix, coefficient);
}


template <int dim>
void MatrixCreator::create_mass_matrix (const Mapping<dim>       &mapping,
					const DoFHandler<dim>    &dof,
					const Quadrature<dim>    &q,
					SparseMatrix<double>     &matrix,
					const Function<dim>      &rhs,
					Vector<double>           &rhs_vector,
					const Function<dim> * const coefficient)
{
  const unsigned int n_threads = multithread_info.n_default_threads;
  Threads::ThreadManager thread_manager;

				   // define starting and end point
				   // for each thread
  typedef typename DoFHandler<dim>::active_cell_iterator active_cell_iterator;  
  std::vector<std::pair<active_cell_iterator,active_cell_iterator> > thread_ranges
    = Threads::split_range<active_cell_iterator> (dof.begin_active(),
						  dof.end(), n_threads);

				   // mutex to synchronise access to
				   // the matrix
  Threads::ThreadMutex mutex;
  
				   // then assemble in parallel
  for (unsigned int thread=0; thread<n_threads; ++thread)
    Threads::spawn (thread_manager,
		    Threads::encapsulate(&MatrixCreator::template
					 create_mass_matrix_2<dim>)
		    .collect_args (mapping, dof, q, matrix, rhs,
				   rhs_vector, coefficient,
				   thread_ranges[thread], mutex));
  thread_manager.wait ();  
};



// TODO:[RH] extend this function to use vector valued coefficient functions for system elements.
template <int dim>
void
MatrixCreator::create_mass_matrix_2 (const Mapping<dim>       &mapping,
				     const DoFHandler<dim>    &dof,
				     const Quadrature<dim>    &q,
				     SparseMatrix<double>     &matrix,
				     const Function<dim>      &rhs,
				     Vector<double>           &rhs_vector,
				     const Function<dim> * const coefficient,
				     const IteratorRange<dim>  range,
				     Threads::ThreadMutex     &mutex)
{
  UpdateFlags update_flags = UpdateFlags(update_values |
					 update_q_points |
					 update_JxW_values);
  if (coefficient != 0)
    update_flags = UpdateFlags (update_flags | update_q_points);

  FEValues<dim> fe_values (mapping, dof.get_fe(), q, update_flags);
    
  const unsigned int dofs_per_cell = fe_values.dofs_per_cell,
		     n_q_points    = fe_values.n_quadrature_points;
  const FiniteElement<dim>    &fe  = fe_values.get_fe();
  const unsigned int n_components  = fe.n_components();

  FullMatrix<double>  cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>      local_rhs (dofs_per_cell);
  std::vector<double> rhs_values (fe_values.n_quadrature_points);
  std::vector<double> coefficient_values (n_q_points);
  
  std::vector<unsigned int> dof_indices (dofs_per_cell);
  
  typename DoFHandler<dim>::active_cell_iterator cell = range.first;
  for (; cell!=range.second; ++cell)
    {
      fe_values.reinit (cell);
      
      cell_matrix.clear ();
      local_rhs.clear ();
      cell->get_dof_indices (dof_indices);
      
      const FullMatrix<double>  &values    = fe_values.get_shape_values ();
      const std::vector<double> &weights   = fe_values.get_JxW_values ();
      rhs.value_list (fe_values.get_quadrature_points(), rhs_values);
      
      if (coefficient != 0)
	{
	  coefficient->value_list (fe_values.get_quadrature_points(),
				   coefficient_values);
	  for (unsigned int point=0; point<n_q_points; ++point)
	    for (unsigned int i=0; i<dofs_per_cell; ++i) 
	      {
		for (unsigned int j=0; j<dofs_per_cell; ++j)
		  if ((n_components==1) ||
		      (fe.system_to_component_index(i).first ==
		       fe.system_to_component_index(j).first))
		    cell_matrix(i,j) += (values(i,point) *
					 values(j,point) *
					 weights[point] *
					 coefficient_values[point]);
		local_rhs(i) += values(i,point) *
				rhs_values[point] *
				weights[point];
	      };
	}
      else
	for (unsigned int point=0; point<n_q_points; ++point)
	  for (unsigned int i=0; i<dofs_per_cell; ++i) 
	    {
	      for (unsigned int j=0; j<dofs_per_cell; ++j)
		if ((n_components==1) ||
		    (fe.system_to_component_index(i).first ==
		     fe.system_to_component_index(j).first))
		  cell_matrix(i,j) += (values(i,point) *
				       values(j,point) *
				       weights[point]);
	      local_rhs(i) += values(i,point) *
			      rhs_values[point] *
			      weights[point];
	    };

				       // transfer everything into the
				       // global object
      mutex.acquire ();
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	for (unsigned int j=0; j<dofs_per_cell; ++j)
	  if ((n_components==1) ||
	      (fe.system_to_component_index(i).first ==
	       fe.system_to_component_index(j).first))
	    matrix.add (dof_indices[i], dof_indices[j],
			cell_matrix(i,j));
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	rhs_vector(dof_indices[i]) += local_rhs(i);
      mutex.release ();
    };
};


template <int dim>
void MatrixCreator::create_mass_matrix (const DoFHandler<dim>    &dof,
					const Quadrature<dim>    &q,
					SparseMatrix<double>     &matrix,
					const Function<dim>      &rhs,
					Vector<double>           &rhs_vector,
					const Function<dim> * const coefficient)
{
  Assert (DEAL_II_COMPAT_MAPPING, ExcCompatibility("mapping"));
  static const MappingQ1<dim> mapping;
  create_mass_matrix(mapping, dof, q, matrix, rhs, rhs_vector, coefficient);
}


#if deal_II_dimension == 1

void MatrixCreator::create_boundary_mass_matrix (const Mapping<1>          &,
						 const DoFHandler<1>       &,
						 const Quadrature<0>       &,
						 SparseMatrix<double>      &,
						 const FunctionMap<1>::type&,
						 Vector<double>            &,
						 std::vector<unsigned int> &,
						 const Function<1>         * const)
{
				   // what would that be in 1d? the
				   // identity matrix on the boundary
				   // dofs?
  Assert (false, ExcNotImplemented());
};


#endif


// TODO:[RH] extend this function to use vector valued coefficient functions for system elements.
template <int dim>
void
MatrixCreator::create_boundary_mass_matrix (const Mapping<dim>        &mapping,
					    const DoFHandler<dim>     &dof,
					    const Quadrature<dim-1>   &q,
					    SparseMatrix<double>      &matrix,
					    const typename FunctionMap<dim>::type         &boundary_functions,
					    Vector<double>            &rhs_vector,
					    std::vector<unsigned int> &dof_to_boundary_mapping,
					    const Function<dim> * const a)
{
  const unsigned int n_threads = multithread_info.n_default_threads;
  Threads::ThreadManager thread_manager;

				   // define starting and end point
				   // for each thread
  typedef typename DoFHandler<dim>::active_cell_iterator active_cell_iterator;  
  std::vector<std::pair<active_cell_iterator,active_cell_iterator> > thread_ranges
    = Threads::split_range<active_cell_iterator> (dof.begin_active(),
						  dof.end(), n_threads);

				   // mutex to synchronise access to
				   // the matrix
  Threads::ThreadMutex mutex;
  
				   // then assemble in parallel
  for (unsigned int thread=0; thread<n_threads; ++thread)
    Threads::spawn (thread_manager,
		    Threads::encapsulate(&MatrixCreator::template 
					 create_boundary_mass_matrix_1<dim>)
		    .collect_args (mapping, dof, q, matrix,
				   boundary_functions, rhs_vector,
				   dof_to_boundary_mapping, a,
				   thread_ranges[thread], mutex));
  thread_manager.wait ();  
};



template <int dim>
void
MatrixCreator::
create_boundary_mass_matrix_1 (const Mapping<dim>        &mapping,
			       const DoFHandler<dim>     &dof,
			       const Quadrature<dim-1>   &q,
			       SparseMatrix<double>      &matrix,
			       const typename FunctionMap<dim>::type &boundary_functions,
			       Vector<double>            &rhs_vector,
			       std::vector<unsigned int> &dof_to_boundary_mapping,
			       const Function<dim> * const a,
			       const IteratorRange<dim>   range,
			       Threads::ThreadMutex      &mutex)
{
  const FiniteElement<dim> &fe = dof.get_fe();
  const unsigned int n_components  = fe.n_components();
  const bool         fe_is_system  = (n_components != 1);
  
  Assert (matrix.n() == dof.n_boundary_dofs(boundary_functions), ExcInternalError());
  Assert (matrix.n() == matrix.m(), ExcInternalError());
  Assert (matrix.n() == rhs_vector.size(), ExcInternalError());
  Assert (boundary_functions.size() != 0, ExcInternalError());
  Assert (dof.get_fe() == fe, ExcInternalError());
  Assert (dof_to_boundary_mapping.size() == dof.n_dofs(),
	  ExcInternalError());
  Assert (n_components == boundary_functions.begin()->second->n_components,
	  ExcComponentMismatch());
#ifdef DEBUG
  if (true)
    {
      unsigned int max_element = 0;
      for (std::vector<unsigned int>::const_iterator i=dof_to_boundary_mapping.begin();
	   i!=dof_to_boundary_mapping.end(); ++i)
	if ((*i != DoFHandler<dim>::invalid_dof_index) &&
	    (*i > max_element))
	  max_element = *i;
      Assert (max_element  == matrix.n()-1, ExcInternalError());
    };
#endif
  
  const unsigned int dofs_per_cell = fe.dofs_per_cell,
		     dofs_per_face = fe.dofs_per_face;
  
  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_vector(dofs_per_cell);


  UpdateFlags update_flags = UpdateFlags (update_values     |
					  update_JxW_values |
					  update_q_points);
  FEFaceValues<dim> fe_values (mapping, fe, q, update_flags);

				   // two variables for the coefficient,
				   // one for the two cases indicated in
				   // the name
  std::vector<double>          coefficient_values_scalar (fe_values.n_quadrature_points);
  std::vector<Vector<double> > coefficient_values_system (fe_values.n_quadrature_points,
							  Vector<double>(n_components));

  std::vector<double>          rhs_values_scalar (fe_values.n_quadrature_points);
  std::vector<Vector<double> > rhs_values_system (fe_values.n_quadrature_points,
						  Vector<double>(n_components));

  std::vector<unsigned int> dofs (dofs_per_cell);
  std::vector<unsigned int> dofs_on_face_vector (dofs_per_face);

				   // for each dof on the cell, have a
				   // flag whether it is on the face
  std::vector<bool>         dof_is_on_face(dofs_per_cell);

  typename DoFHandler<dim>::active_cell_iterator cell = range.first;
  for (; cell!=range.second; ++cell)
    for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
				       // check if this face is on that part of
				       // the boundary we are interested in
      if (boundary_functions.find(cell->face(face)->boundary_indicator()) !=
	  boundary_functions.end())
	{
	  cell_matrix.clear ();
	  cell_vector.clear ();
	  
	  fe_values.reinit (cell, face);

	  const FullMatrix<double>  &values  = fe_values.get_shape_values ();
	  const std::vector<double> &weights = fe_values.get_JxW_values ();

	  if (fe_is_system)
					     // FE has several components
	    {
	      boundary_functions.find(cell->face(face)->boundary_indicator())
		->second->vector_value_list (fe_values.get_quadrature_points(),
					     rhs_values_system);

	      if (a != 0)
		{
		  a->vector_value_list (fe_values.get_quadrature_points(),
					coefficient_values_system);
		  for (unsigned int point=0; point<fe_values.n_quadrature_points; ++point)
		    for (unsigned int i=0; i<fe_values.dofs_per_cell; ++i) 
		      {
			for (unsigned int j=0; j<fe_values.dofs_per_cell; ++j)
			  if (fe.system_to_component_index(i).first ==
			      fe.system_to_component_index(j).first)
			    {
			      cell_matrix(i,j)
				+= (values(i,point) *
				    values(j,point) *
				    weights[point] *
				    coefficient_values_system[point](
				      fe.system_to_component_index(i).first));
			    };
			
			cell_vector(i) += values(i,point) *
					  rhs_values_system[point](
					    fe.system_to_component_index(i).first) *
					  weights[point];
		      };
		}
	      else
		for (unsigned int point=0; point<fe_values.n_quadrature_points; ++point)
		  for (unsigned int i=0; i<fe_values.dofs_per_cell; ++i) 
		    {
		      for (unsigned int j=0; j<fe_values.dofs_per_cell; ++j)
			if (fe.system_to_component_index(i).first ==
			    fe.system_to_component_index(j).first)
			  {
			    cell_matrix(i,j) += (values(i,point) *
						 values(j,point) *
						 weights[point]);
			  };
		      
		      cell_vector(i) += values(i,point) *
					rhs_values_system[point](
					  fe.system_to_component_index(i).first) *
					weights[point];
		    };
	    }
	  else
					     // FE is a scalar one
	    {
	      boundary_functions.find(cell->face(face)->boundary_indicator())
		->second->value_list (fe_values.get_quadrature_points(), rhs_values_scalar);

	      if (a != 0)
		{
		  a->value_list (fe_values.get_quadrature_points(),
				 coefficient_values_scalar);
		  for (unsigned int point=0; point<fe_values.n_quadrature_points; ++point)
		    for (unsigned int i=0; i<fe_values.dofs_per_cell; ++i) 
		      {
			for (unsigned int j=0; j<fe_values.dofs_per_cell; ++j)
			  cell_matrix(i,j) += (values(i,point) *
					       values(j,point) *
					       weights[point] *
					       coefficient_values_scalar[point]);
			cell_vector(i) += values(i,point) *
					  rhs_values_scalar[point] *
					  weights[point];
		      };
		}
	      else
		for (unsigned int point=0; point<fe_values.n_quadrature_points; ++point)
		  for (unsigned int i=0; i<fe_values.dofs_per_cell; ++i) 
		    {
		      for (unsigned int j=0; j<fe_values.dofs_per_cell; ++j)
			cell_matrix(i,j) += (values(i,point) *
					     values(j,point) *
					     weights[point]);
		      cell_vector(i) += values(i,point) *
					rhs_values_scalar[point] *
					weights[point];
		    };
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
					   // face into @p{dofs_on_face_vector},
					   // a vector as always. Usually,
					   // searching in a vector is
					   // inefficient, so we copy the dofs
					   // into a set, which enables binary
					   // searches.
	  cell->get_dof_indices (dofs);
	  cell->face(face)->get_dof_indices (dofs_on_face_vector);

					   // check for each of the
					   // dofs on this cell
					   // whether it is on the
					   // face
	  for (unsigned int i=0; i<dofs_per_cell; ++i)
	    dof_is_on_face[i] = (std::find(dofs_on_face_vector.begin(),
					   dofs_on_face_vector.end(),
					   dofs[i])
				 !=
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
	    if (std::fabs(cell_matrix(i,i)) > max_diag_entry)
	      max_diag_entry = std::fabs(cell_matrix(i,i));
#endif  

	  mutex.acquire ();
	  for (unsigned int i=0; i<dofs_per_cell; ++i)
	    for (unsigned int j=0; j<dofs_per_cell; ++j)
	      if (dof_is_on_face[i] && dof_is_on_face[j])
		matrix.add(dof_to_boundary_mapping[dofs[i]],
			   dof_to_boundary_mapping[dofs[j]],
			   cell_matrix(i,j));
	      else
		{
//TODO:[?] We assume here that shape functions that have support also on the boundary also have their support point on the boundary (by having their indices in dof_is_on_face set true). This is not true, however, e.g. for discontinuous elements.
						   // compare here for relative
						   // smallness
		  Assert (std::fabs(cell_matrix(i,j)) <= 1e-10 * max_diag_entry,
			  ExcInternalError ());
		};
	  
	  for (unsigned int j=0; j<dofs_per_cell; ++j)
	    if (dof_is_on_face[j])
	      rhs_vector(dof_to_boundary_mapping[dofs[j]]) += cell_vector(j);
	    else
	      {
						 // compare here for relative
						 // smallness
		Assert (std::fabs(cell_vector(j)) <= 1e-10 * max_diag_entry,
			ExcInternalError());
	      };
	  mutex.release ();
	};
};


template <int dim>
void MatrixCreator::create_boundary_mass_matrix (const DoFHandler<dim>     &dof,
						 const Quadrature<dim-1>   &q,
						 SparseMatrix<double>      &matrix,
						 const typename FunctionMap<dim>::type &rhs,
						 Vector<double>            &rhs_vector,
						 std::vector<unsigned int> &dof_to_boundary_mapping,
						 const Function<dim> * const a)
{
  Assert (DEAL_II_COMPAT_MAPPING, ExcCompatibility("mapping"));
  static const MappingQ1<dim> mapping;
  create_boundary_mass_matrix(mapping, dof, q, matrix, rhs, rhs_vector,
			      dof_to_boundary_mapping, a);
}



template <int dim>
void MatrixCreator::create_laplace_matrix (const Mapping<dim>       &mapping,
					   const DoFHandler<dim>    &dof,
					   const Quadrature<dim>    &q,
					   SparseMatrix<double>     &matrix,
					   const Function<dim> * const coefficient)
{
  const unsigned int n_threads = multithread_info.n_default_threads;
  Threads::ThreadManager thread_manager;

				   // define starting and end point
				   // for each thread
  typedef typename DoFHandler<dim>::active_cell_iterator active_cell_iterator;  
  std::vector<std::pair<active_cell_iterator,active_cell_iterator> > thread_ranges
    = Threads::split_range<active_cell_iterator> (dof.begin_active(),
						  dof.end(), n_threads);

				   // mutex to synchronise access to
				   // the matrix
  Threads::ThreadMutex mutex;
  
				   // then assemble in parallel
  for (unsigned int thread=0; thread<n_threads; ++thread)
    Threads::spawn (thread_manager,
		    Threads::encapsulate(&MatrixCreator::template
					 create_laplace_matrix_1<dim>)
		    .collect_args (mapping, dof, q, matrix, coefficient,
				   thread_ranges[thread], mutex));
  thread_manager.wait ();  
};



// TODO:[RH] extend this function to use vector valued coefficient functions for system elements.
template <int dim>
void MatrixCreator::create_laplace_matrix_1 (const Mapping<dim>       &mapping,
					     const DoFHandler<dim>    &dof,
					     const Quadrature<dim>    &q,
					     SparseMatrix<double>     &matrix,
					     const Function<dim> * const coefficient,
					     const IteratorRange<dim>  range,
					     Threads::ThreadMutex     &mutex)
{
  UpdateFlags update_flags = UpdateFlags(update_JxW_values |
					 update_gradients);
  if (coefficient != 0)
    update_flags = UpdateFlags (update_flags | update_q_points);

  FEValues<dim> fe_values (mapping, dof.get_fe(), q, update_flags);
    
  const unsigned int dofs_per_cell = fe_values.dofs_per_cell,
		     n_q_points    = fe_values.n_quadrature_points;
  const FiniteElement<dim>    &fe  = fe_values.get_fe();
  const unsigned int n_components  = fe.n_components();

  FullMatrix<double>  cell_matrix (dofs_per_cell, dofs_per_cell);
  std::vector<double> coefficient_values (n_q_points);
  
  std::vector<unsigned int> dof_indices (dofs_per_cell);
  
  typename DoFHandler<dim>::active_cell_iterator cell = range.first;
  for (; cell!=range.second; ++cell)
    {
      fe_values.reinit (cell);
      
      cell_matrix.clear ();
      cell->get_dof_indices (dof_indices);
      
      const std::vector<std::vector<Tensor<1,dim> > >
	&grads   = fe_values.get_shape_grads ();
      const std::vector<double> &weights = fe_values.get_JxW_values ();
      
      if (coefficient != 0)
	{
	  coefficient->value_list (fe_values.get_quadrature_points(),
				   coefficient_values);
	  for (unsigned int point=0; point<n_q_points; ++point)
	    for (unsigned int i=0; i<dofs_per_cell; ++i) 
	      for (unsigned int j=0; j<dofs_per_cell; ++j)
		if ((n_components==1) ||
		    (fe.system_to_component_index(i).first ==
		     fe.system_to_component_index(j).first))
		  cell_matrix(i,j) += (grads[i][point] *
				       grads[j][point] *
				       weights[point] *
				       coefficient_values[point]);
	}
      else
	for (unsigned int point=0; point<n_q_points; ++point)
	  for (unsigned int i=0; i<dofs_per_cell; ++i) 
	    for (unsigned int j=0; j<dofs_per_cell; ++j)
	      if ((n_components==1) ||
		  (fe.system_to_component_index(i).first ==
		   fe.system_to_component_index(j).first))
		cell_matrix(i,j) += (grads[i][point] *
				     grads[j][point] *
				     weights[point]);

				       // transfer everything into the
				       // global object
      mutex.acquire ();
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	for (unsigned int j=0; j<dofs_per_cell; ++j)
	  if ((n_components==1) ||
	      (fe.system_to_component_index(i).first ==
	       fe.system_to_component_index(j).first))
	    matrix.add (dof_indices[i], dof_indices[j],
			cell_matrix(i,j));
      mutex.release ();
    };
};



template <int dim>
void MatrixCreator::create_laplace_matrix (const DoFHandler<dim>    &dof,
					   const Quadrature<dim>    &q,
					   SparseMatrix<double>     &matrix,
					   const Function<dim> * const coefficient)
{
  Assert (DEAL_II_COMPAT_MAPPING, ExcCompatibility("mapping"));
  static const MappingQ1<dim> mapping;
  create_laplace_matrix(mapping, dof, q, matrix, coefficient);
}


//TODO:[GK,RH] maybe recreate this function
/*

  template <int dim>
  void MatrixCreator::create_level_laplace_matrix (unsigned int level,
  const MGDoFHandler<dim>    &dof,
  const Quadrature<dim>    &q,
  SparseMatrix<float>     &matrix,
  const Function<dim> * const a)
  {
  Vector<double> dummy;   // no entries, should give an error if accessed
  UpdateFlags update_flags = UpdateFlags(update_gradients |
  update_JxW_values);
  if (a != 0)
  update_flags = UpdateFlags(update_flags | update_q_points);
  const Assembler<dim>::AssemblerData data (dof,
  true, false,  // assemble matrix but not rhs
  matrix, dummy,
  q, update_flags);
  TriaIterator<dim, Assembler<dim> >
  assembler (const_cast<Triangulation<dim>*>(&dof.get_tria()),
  dof.get_tria().begin(level)->level(),
  dof.get_tria().begin(level)->index(),
  &data);
  LaplaceMatrix<dim> equation (0, a);
  do 
  {
  assembler->assemble (equation);
  ++assembler
  }
  while ( (assembler.state()==valid) && (assembler->level() == level) );
  };

*/



template <int dim>
void MatrixCreator::create_laplace_matrix (const Mapping<dim>       &mapping,
					   const DoFHandler<dim>    &dof,
					   const Quadrature<dim>    &q,
					   SparseMatrix<double>     &matrix,
					   const Function<dim>      &rhs,
					   Vector<double>           &rhs_vector,
					   const Function<dim> * const coefficient)
{
  const unsigned int n_threads = multithread_info.n_default_threads;
  Threads::ThreadManager thread_manager;

				   // define starting and end point
				   // for each thread
  typedef typename DoFHandler<dim>::active_cell_iterator active_cell_iterator;  
  std::vector<std::pair<active_cell_iterator,active_cell_iterator> > thread_ranges
    = Threads::split_range<active_cell_iterator> (dof.begin_active(),
						  dof.end(), n_threads);

				   // mutex to synchronise access to
				   // the matrix
  Threads::ThreadMutex mutex;
  
				   // then assemble in parallel
  for (unsigned int thread=0; thread<n_threads; ++thread)
    Threads::spawn (thread_manager,
		    Threads::encapsulate(&MatrixCreator::template
					 create_laplace_matrix_2<dim>)
		    .collect_args (mapping, dof, q, matrix, rhs,
				   rhs_vector, coefficient,
				   thread_ranges[thread], mutex));
  thread_manager.wait ();  
};



// TODO:[RH] extend this function to use vector valued coefficient functions for system elements.
template <int dim>
void
MatrixCreator::create_laplace_matrix_2 (const Mapping<dim>       &mapping,
					const DoFHandler<dim>    &dof,
					const Quadrature<dim>    &q,
					SparseMatrix<double>     &matrix,
					const Function<dim>      &rhs,
					Vector<double>           &rhs_vector,
					const Function<dim> * const coefficient,
					const IteratorRange<dim>  range,
					Threads::ThreadMutex     &mutex)
{
  UpdateFlags update_flags = UpdateFlags(update_values    |
					 update_gradients |
					 update_q_points  |
					 update_JxW_values);
  if (coefficient != 0)
    update_flags = UpdateFlags (update_flags | update_q_points);

  FEValues<dim> fe_values (mapping, dof.get_fe(), q, update_flags);
    
  const unsigned int dofs_per_cell = fe_values.dofs_per_cell,
		     n_q_points    = fe_values.n_quadrature_points;
  const FiniteElement<dim>    &fe  = fe_values.get_fe();
  const unsigned int n_components  = fe.n_components();

  FullMatrix<double>  cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>      local_rhs (dofs_per_cell);
  std::vector<double> rhs_values (fe_values.n_quadrature_points);
  std::vector<double> coefficient_values (n_q_points);
  
  std::vector<unsigned int> dof_indices (dofs_per_cell);
  
  typename DoFHandler<dim>::active_cell_iterator cell = range.first;
  for (; cell!=range.second; ++cell)
    {
      fe_values.reinit (cell);
      
      cell_matrix.clear ();
      local_rhs.clear ();
      cell->get_dof_indices (dof_indices);
      
      const FullMatrix<double>  &values    = fe_values.get_shape_values ();
      const std::vector<std::vector<Tensor<1,dim> > >
	&grads   = fe_values.get_shape_grads ();
      const std::vector<double> &weights   = fe_values.get_JxW_values ();
      rhs.value_list (fe_values.get_quadrature_points(), rhs_values);
      
      if (coefficient != 0)
	{
	  coefficient->value_list (fe_values.get_quadrature_points(),
				   coefficient_values);
	  for (unsigned int point=0; point<n_q_points; ++point)
	    for (unsigned int i=0; i<dofs_per_cell; ++i) 
	      {
		for (unsigned int j=0; j<dofs_per_cell; ++j)
		  if ((n_components==1) ||
		      (fe.system_to_component_index(i).first ==
		       fe.system_to_component_index(j).first))
		    cell_matrix(i,j) += (grads[i][point] *
					 grads[j][point] *
					 weights[point] *
					 coefficient_values[point]);
		local_rhs(i) += values(i,point) *
				rhs_values[point] *
				weights[point];
	      };
	}
      else
	for (unsigned int point=0; point<n_q_points; ++point)
	  for (unsigned int i=0; i<dofs_per_cell; ++i) 
	    {
	      for (unsigned int j=0; j<dofs_per_cell; ++j)
		if ((n_components==1) ||
		    (fe.system_to_component_index(i).first ==
		     fe.system_to_component_index(j).first))
		  cell_matrix(i,j) += (grads[i][point] *
				       grads[j][point] *
				       weights[point]);
	      local_rhs(i) += values(i,point) *
			      rhs_values[point] *
			      weights[point];
	    };

				       // transfer everything into the
				       // global object
      mutex.acquire ();
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	for (unsigned int j=0; j<dofs_per_cell; ++j)
	  if ((n_components==1) ||
	      (fe.system_to_component_index(i).first ==
	       fe.system_to_component_index(j).first))
	    matrix.add (dof_indices[i], dof_indices[j],
			cell_matrix(i,j));
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	rhs_vector(dof_indices[i]) += local_rhs(i);
      mutex.release ();
    };
};



template <int dim>
void MatrixCreator::create_laplace_matrix (const DoFHandler<dim>    &dof,
					   const Quadrature<dim>    &q,
					   SparseMatrix<double>     &matrix,
					   const Function<dim>      &rhs,
					   Vector<double>           &rhs_vector,
					   const Function<dim> * const coefficient)
{
  Assert (DEAL_II_COMPAT_MAPPING, ExcCompatibility("mapping"));
  static const MappingQ1<dim> mapping;  
  create_laplace_matrix(mapping, dof, q, matrix, rhs, rhs_vector, coefficient);
}






// explicit instantiations

template
void MatrixCreator::create_mass_matrix (const Mapping<deal_II_dimension>       &mapping,
					const DoFHandler<deal_II_dimension>    &dof,
					const Quadrature<deal_II_dimension>    &q,
					SparseMatrix<double>     &matrix,
					const Function<deal_II_dimension> * const coefficient);
template
void MatrixCreator::create_mass_matrix (const DoFHandler<deal_II_dimension>    &dof,
					const Quadrature<deal_II_dimension>    &q,
					SparseMatrix<double>     &matrix,
					const Function<deal_II_dimension> * const coefficient);
template
void MatrixCreator::create_mass_matrix (const Mapping<deal_II_dimension>       &mapping,
					const DoFHandler<deal_II_dimension>    &dof,
					const Quadrature<deal_II_dimension>    &q,
					SparseMatrix<double>     &matrix,
					const Function<deal_II_dimension>      &rhs,
					Vector<double>           &rhs_vector,
					const Function<deal_II_dimension> * const coefficient);
template
void MatrixCreator::create_mass_matrix (const DoFHandler<deal_II_dimension>    &dof,
					const Quadrature<deal_II_dimension>    &q,
					SparseMatrix<double>     &matrix,
					const Function<deal_II_dimension>      &rhs,
					Vector<double>           &rhs_vector,
					const Function<deal_II_dimension> * const coefficient);

#if deal_II_dimension != 1
template
void
MatrixCreator::create_boundary_mass_matrix (const Mapping<deal_II_dimension>        &mapping,
					    const DoFHandler<deal_II_dimension>     &dof,
					    const Quadrature<deal_II_dimension-1>   &q,
					    SparseMatrix<double>      &matrix,
					    const FunctionMap<deal_II_dimension>::type         &boundary_functions,
					    Vector<double>            &rhs_vector,
					    std::vector<unsigned int> &dof_to_boundary_mapping,
					    const Function<deal_II_dimension> * const a);
#endif

template
void MatrixCreator::create_boundary_mass_matrix (const DoFHandler<deal_II_dimension>     &dof,
						 const Quadrature<deal_II_dimension-1>   &q,
						 SparseMatrix<double>      &matrix,
						 const FunctionMap<deal_II_dimension>::type &rhs,
						 Vector<double>            &rhs_vector,
						 std::vector<unsigned int> &dof_to_boundary_mapping,
						 const Function<deal_II_dimension> * const a);
template
void MatrixCreator::create_laplace_matrix (const DoFHandler<deal_II_dimension>    &dof,
					   const Quadrature<deal_II_dimension>    &q,
					   SparseMatrix<double>     &matrix,
					   const Function<deal_II_dimension> * const coefficient);

template
void MatrixCreator::create_laplace_matrix (const Mapping<deal_II_dimension>       &mapping,
					   const DoFHandler<deal_II_dimension>    &dof,
					   const Quadrature<deal_II_dimension>    &q,
					   SparseMatrix<double>     &matrix,
					   const Function<deal_II_dimension> * const coefficient);
template
void MatrixCreator::create_laplace_matrix (const Mapping<deal_II_dimension>       &mapping,
					   const DoFHandler<deal_II_dimension>    &dof,
					   const Quadrature<deal_II_dimension>    &q,
					   SparseMatrix<double>     &matrix,
					   const Function<deal_II_dimension>      &rhs,
					   Vector<double>           &rhs_vector,
					   const Function<deal_II_dimension> * const coefficient);
template
void MatrixCreator::create_laplace_matrix (const DoFHandler<deal_II_dimension>    &dof,
					   const Quadrature<deal_II_dimension>    &q,
					   SparseMatrix<double>     &matrix,
					   const Function<deal_II_dimension>      &rhs,
					   Vector<double>           &rhs_vector,
					   const Function<deal_II_dimension> * const coefficient);




