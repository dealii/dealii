//----------------------------  error_estimator.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  error_estimator.cc  ---------------------------


#include <fe/fe.h>
#include <fe/fe_values.h>
#include <fe/fe_update_flags.h>
#include <base/thread_management.h>
#include <base/quadrature.h>
#include <base/quadrature_lib.h>
#include <numerics/error_estimator.h>
#include <dofs/dof_handler.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <grid/geometry_info.h>
#include <lac/vector.h>

#include <numeric>
#include <algorithm>
#include <cmath>
#include <vector>
#include <base/timer.h>



static
inline
double sqr (const double x)
{
  return x*x;
};




#if deal_II_dimension == 1

template <>
KellyErrorEstimator<1>::Data::Data(const DoFHandler<1>                 &,
				   const Quadrature<0>                 &,
				   const FunctionMap                   &,
				   const vector<const Vector<double>*> &,
				   const vector<bool>                  &,
				   const Function<1>                   *,
				   const unsigned int                   ,
				   FaceIntegrals                       &):
		dof_handler(* static_cast <const DoFHandler<1> *> (0)),
		quadrature(* static_cast <const Quadrature<0> *> (0)),
		neumann_bc(* static_cast <const FunctionMap *> (0)),
		solutions(* static_cast <const vector<const Vector<double>*> *> (0)),
		face_integrals (* static_cast<FaceIntegrals*> (0))
{
  Assert (false, ExcInternalError());
}

#else

template <int dim>
KellyErrorEstimator<dim>::Data::Data(const DoFHandler<dim>               &dof_handler,
				     const Quadrature<dim-1>             &quadrature,
				     const FunctionMap                   &neumann_bc,
				     const vector<const Vector<double>*> &solutions,
				     const vector<bool>                  &component_mask,
				     const Function<dim>                 *coefficients,
				     const unsigned int                   n_threads,
				     FaceIntegrals                       &face_integrals):
		dof_handler (dof_handler),
		quadrature (quadrature),
		neumann_bc (neumann_bc),
		solutions (solutions),
		component_mask (component_mask),
		coefficients (coefficients),
		n_threads (n_threads),
		n_solution_vectors (solutions.size()),
		face_integrals (face_integrals)
{
  const unsigned int n_components = dof_handler.get_fe().n_components();
  
  Assert (component_mask.size() == n_components, ExcInvalidComponentMask());
  Assert (count(component_mask.begin(), component_mask.end(), true) > 0,
	  ExcInvalidComponentMask());
  
  Assert ((coefficients == 0) ||
	  (coefficients->n_components == n_components) ||
	  (coefficients->n_components == 1),
	  ExcInvalidCoefficient());
  
  Assert (neumann_bc.find(255) == neumann_bc.end(),
	  ExcInvalidBoundaryIndicator());
  
  for (typename FunctionMap::const_iterator i=neumann_bc.begin(); i!=neumann_bc.end(); ++i)
    Assert (i->second->n_components == n_components, ExcInvalidBoundaryFunction());
  
  
				   // Init the size of a lot of vectors
				   // needed in the calculations once
				   // per thread.
  const unsigned int n_q_points = quadrature.n_quadrature_points;
  
  normal_vectors.resize(n_q_points);
  coefficient_values1.resize(n_q_points);
  coefficient_values.resize(n_q_points);
  JxW_values.resize(n_q_points);  

  phi.resize(n_solution_vectors);
  psi.resize(n_solution_vectors);
  neighbor_psi.resize(n_solution_vectors);

  for (unsigned int i=0; i<n_solution_vectors; ++i)
    {
      phi[i].resize(n_q_points);
      psi[i].resize(n_q_points);
      neighbor_psi[i].resize(n_q_points);

      for (unsigned int qp=0;qp<n_q_points;++qp)
	{
	  phi[i][qp].resize(n_components);
	  psi[i][qp].resize(n_components);
	  neighbor_psi[i][qp].resize(n_components);
	};
    };

  for (unsigned int qp=0;qp<n_q_points;++qp)
    coefficient_values[qp].reinit(n_components);
}

#endif



// the following function is still independent of dimension, but it
// calls dimension dependent functions
template <int dim>
void KellyErrorEstimator<dim>::estimate (const DoFHandler<dim>   &dof_handler,
					 const Quadrature<dim-1> &quadrature,
					 const FunctionMap       &neumann_bc,
					 const Vector<double>    &solution,
					 Vector<float>           &error,
					 const vector<bool>      &component_mask,
					 const Function<dim>     *coefficients,
					 unsigned int             n_threads)
{
				   // just pass on to the other function
  const vector<const Vector<double>*> solutions (1, &solution);
  vector<Vector<float>*>              errors (1, &error);
  estimate (dof_handler, quadrature, neumann_bc, solutions, errors,
	    component_mask, coefficients, n_threads);
};




#if deal_II_dimension == 1

template <>
void KellyErrorEstimator<1>::estimate_some (Data &, const unsigned int)
{
				   // in 1d, the @p{estimate} function
				   // does all the work
  Assert (false, ExcInternalError() );
}



template <>
void KellyErrorEstimator<1>::estimate (const DoFHandler<1>                 &dof_handler,
				       const Quadrature<0>                 &,
				       const FunctionMap                   &neumann_bc,
				       const vector<const Vector<double>*> &solutions,
				       vector<Vector<float>*>              &errors,
				       const vector<bool>                  &component_mask_,
				       const Function<1>                   *coefficient,
				       const unsigned int)
{
  const unsigned int n_components       = dof_handler.get_fe().n_components();
  const unsigned int n_solution_vectors = solutions.size();

				   // sanity checks
  Assert (solutions.size() > 0,
	  ExcNoSolutions());
  Assert (solutions.size() == errors.size(),
	  ExcIncompatibleNumberOfElements(solutions.size(), errors.size()));
  for (unsigned int n=0; n<solutions.size(); ++n)
    Assert (solutions[n]->size() == dof_handler.n_dofs(),
	    ExcInvalidSolutionVector());
  
				   // if no mask given: treat all components
  vector<bool> component_mask ((component_mask_.size() == 0)    ?
			       vector<bool>(n_components, true) :
			       component_mask_);
  Assert (component_mask.size() == n_components, ExcInvalidComponentMask());
  Assert (count(component_mask.begin(), component_mask.end(), true) > 0,
	  ExcInvalidComponentMask());
  
  Assert ((coefficient == 0) ||
	  (coefficient->n_components == n_components) ||
	  (coefficient->n_components == 1),
	  ExcInvalidCoefficient());

  for (FunctionMap::const_iterator i=neumann_bc.begin(); i!=neumann_bc.end(); ++i)
    Assert (i->second->n_components == n_components, ExcInvalidBoundaryFunction());


  const unsigned int dim=1;

				   // reserve one slot for each cell and set
				   // it to zero
  for (unsigned int n=0; n<n_solution_vectors; ++n)
    (*errors[n]).reinit (dof_handler.get_tria().n_active_cells());

				   // fields to get the gradients on
				   // the present and the neighbor cell.
				   //
				   // for the neighbor gradient, we
				   // need several auxiliary fields,
				   // depending on the way we get it
				   // (see below)
  vector<vector<vector<Tensor<1,1> > > >
    gradients_here (n_solution_vectors,
		    vector<vector<Tensor<1,1> > >(2, vector<Tensor<1,1> >(n_components)));
  vector<vector<vector<Tensor<1,1> > > >
    gradients_neighbor (gradients_here);
  vector<Vector<double> >
    grad_neighbor (n_solution_vectors, Vector<double>(n_components));

				   // reserve some space for
				   // coefficient values at one point.
				   // if there is no coefficient, then
				   // we fill it by unity once and for
				   // all and don't set it any more
  Vector<double> coefficient_values (n_components);
  if (coefficient == 0)
    for (unsigned int c=0; c<n_components; ++c)
      coefficient_values(c) = 1;
  
				   // loop over all
				   // cells. note that the error
				   // indicator is only a sum over the
				   // two contributions from the two
				   // vertices of each cell.
  QTrapez<1> quadrature;
  FEValues<1> fe_values (dof_handler.get_fe(), quadrature, update_gradients);
  active_cell_iterator cell = dof_handler.begin_active();
  for (unsigned int cell_index=0; cell != dof_handler.end(); ++cell, ++cell_index)
    {
      for (unsigned int n=0; n<n_solution_vectors; ++n)
	(*errors[n])(cell_index) = 0;
      
				       // loop over the two points bounding
				       // this line. n==0 is left point,
				       // n==1 is right point
      for (unsigned int n=0; n<2; ++n)
	{
					   // find right active neighbor
	  DoFHandler<dim>::cell_iterator neighbor = cell->neighbor(n);
	  if (neighbor.state() == valid)
	    while (neighbor->has_children())
	      neighbor = neighbor->child(n==0 ? 1 : 0);
      
					   // now get the gradients on the
					   // both sides of the point
	  fe_values.reinit (cell);

	  for (unsigned int s=0; s<n_solution_vectors; ++s)
	    fe_values.get_function_grads (*solutions[s], gradients_here[s]);

	  if (neighbor.state() == valid)
	    {
	      fe_values.reinit (neighbor);

	      for (unsigned int s=0; s<n_solution_vectors; ++s)
		fe_values.get_function_grads (*solutions[s],
					      gradients_neighbor[s]);

					       // extract the
					       // gradients of all the
					       // components. [0]
					       // means: x-derivative,
					       // which is the only
					       // one here
	      for (unsigned int s=0; s<n_solution_vectors; ++s)
		for (unsigned int c=0; c<n_components; ++c)
		  grad_neighbor[s](c) = gradients_neighbor[s][n==0 ? 1 : 0][c][0];
	    }
	  else
	    if (neumann_bc.find(n) != neumann_bc.end())
					       // if Neumann b.c., then fill
					       // the gradients field which
					       // will be used later on.
	      for (unsigned int s=0; s<n_solution_vectors; ++s)
		neumann_bc.find(n)->second->vector_value(cell->vertex(0),
							 grad_neighbor[s]);
	    else
					       // fill with zeroes.
	      for (unsigned int s=0; s<n_solution_vectors; ++s)
		grad_neighbor[s].clear ();

					   // if there is a
					   // coefficient, then
					   // evaluate it at the
					   // present position. if
					   // there is none, reuse the
					   // preset values.
	  if (coefficient != 0)
	    {
	      if (coefficient->n_components == 1)
		{
		  const double c_value = coefficient->value (cell->vertex(n));
		  for (unsigned int c=0; c<n_components; ++c)
		    coefficient_values(c) = c_value;
		}
	      else
		coefficient->vector_value(cell->vertex(n),
					  coefficient_values);
	    };


	  for (unsigned int s=0; s<n_solution_vectors; ++s)
	    for (unsigned int component=0; component<n_components; ++component)
	      if (component_mask[component] == true)
		{
						   // get gradient
						   // here. [0] means
						   // x-derivative
						   // (there is no
						   // other component
						   // in 1d)
		  const double grad_here = gradients_here[s][n][component][0];
		  
		  const double jump = ((grad_here - grad_neighbor[s](component)) *
				       coefficient_values(component));
		  (*errors[s])(cell_index) += jump*jump * cell->diameter();
		};
	};
      
      for (unsigned int s=0; s<n_solution_vectors; ++s)
	(*errors[s])(cell_index) = sqrt((*errors[s])(cell_index));
    };
};


#else // #if deal_II_dimension !=1


template <int dim>
void KellyErrorEstimator<dim>::estimate_some (Data               &data,
					      const unsigned int  this_thread) 
{
  const unsigned int n_solution_vectors = data.n_solution_vectors;
  
				   // make up a fe face values object for the
				   // restriction of the finite element function
				   // to a face, for the present cell and its
				   // neighbor. In principle we would only need
				   // one at a time, but this way we can
				   // have more fine grained access to what
				   // values really need to be computed (we
				   // need not compute all values on the
				   // neighbor cells, so using two objects
				   // gives us a performance gain).
  FEFaceValues<dim> fe_face_values_cell (data.dof_handler.get_fe(),
					 data.quadrature,
					 UpdateFlags(update_gradients      |
						     update_JxW_values     |
						     ((!data.neumann_bc.empty() ||
						       (data.coefficients != 0))  ?
						      update_q_points : 0) |
						     update_normal_vectors)); 
  FEFaceValues<dim> fe_face_values_neighbor (data.dof_handler.get_fe(),
					     data.quadrature,
					     update_gradients);
  FESubfaceValues<dim> fe_subface_values (data.dof_handler.get_fe(),
					  data.quadrature,
					  update_gradients);


  active_cell_iterator cell=data.dof_handler.begin_active();

				   // calculate the start cell for
				   // this thread. note that this way
				   // the threads work interleaved on
				   // successive cells, rather than on
				   // blocks of cells. the reason is
				   // that it takes vastly more time
				   // to work on cells with hanging
				   // nodes than on regular cells, but
				   // such cells are not evenly
				   // distributed across the range of
				   // cell iterators, so in order to
				   // have the different threads do
				   // approximately the same amount of
				   // work, we have to let them work
				   // interleaved to the effect of a
				   // pseudorandom distribution of the
				   // `hard' cells to the different
				   // threads.
  for (unsigned int t=0; (t<this_thread) && (cell!=data.dof_handler.end());
       ++t, ++cell);

  
				   // loop over all cells for this thread
				   // the iteration of cell is done at the end
  for (; cell!=data.dof_handler.end(); )
    {
      
				       // loop over all faces of this cell
      for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
	{
					   // if we already visited
					   // this face: do
					   // nothing. only check
					   // component for first
					   // solution vector, as we
					   // treat them all at the
					   // same time
	  if (data.face_integrals[cell->face(face_no)][0] >=0)
	    continue;


					   // if the neighboring cell is less
					   // refined than the present one, then
					   // do nothing since we integrate
					   // over the subfaces when we visit
					   // the coarse cells.
	  if (cell->at_boundary(face_no) == false)
	    if (cell->neighbor(face_no)->level() < cell->level())
	      continue;
	  
					   // if this face is part of the boundary
					   // but not of the neumann boundary
					   // -> nothing to do. However, to make
					   // things easier when summing up the
					   // contributions of the faces of cells,
					   // we enter this face into the list
					   // of faces with contribution zero.
	  const unsigned char boundary_indicator
	    = cell->face(face_no)->boundary_indicator();
	  if (cell->face(face_no)->at_boundary()
	      &&
	      data.neumann_bc.find(boundary_indicator)==data.neumann_bc.end()) 
	    {
	      for (unsigned int n=0; n<n_solution_vectors; ++n)
		data.face_integrals[cell->face(face_no)][n] = 0;
	      continue;
	    };


	  if (cell->face(face_no)->has_children() == false)
					     // if the face is a regular one, i.e.
					     // either on the other side there is
					     // nirvana (face is at boundary), or
					     // the other side's refinement level
					     // is the same as that of this side,
					     // then handle the integration of
					     // these both cases together
	    integrate_over_regular_face (data,
					 cell, face_no,
					 fe_face_values_cell,
					 fe_face_values_neighbor);
	  
	  else
					     // otherwise we need to do some
					     // special computations which do
					     // not fit into the framework of
					     // the above function
	    integrate_over_irregular_face (data,
					   cell, face_no,
					   fe_face_values_cell,
					   fe_subface_values);
	};

				       // go to next cell for this
				       // thread. note that the cells
				       // for each of the threads are
				       // interleaved.
      for (unsigned int t=0;
	   ((t<data.n_threads) && (cell!=data.dof_handler.end()));
	   ++t, ++cell);
    };
};



template <int dim>
void KellyErrorEstimator<dim>::estimate (const DoFHandler<dim>               &dof_handler,
					 const Quadrature<dim-1>             &quadrature,
					 const FunctionMap                   &neumann_bc,
					 const vector<const Vector<double>*> &solutions,
					 vector<Vector<float>*>              &errors,
					 const vector<bool>                  &component_mask,
					 const Function<dim>                 *coefficients,
					 unsigned int                         n_threads)
{
				   // sanity checks
  Assert (solutions.size() > 0,
	  ExcNoSolutions());
  Assert (solutions.size() == errors.size(),
	  ExcIncompatibleNumberOfElements(solutions.size(), errors.size()));
  for (unsigned int n=0; n<solutions.size(); ++n)
    Assert (solutions[n]->size() == dof_handler.n_dofs(),
	    ExcInvalidSolutionVector());
  
	  
				   // if NOT multithreaded, set n_threads to one
#ifndef DEAL_II_USE_MT
  n_threads = 1;
#endif
  Assert (n_threads > 0, ExcInternalError());
  
  const unsigned int n_solution_vectors = solutions.size();

  				   // Map of integrals indexed by
				   // the corresponding face. In this map
				   // we store the integrated jump of the
				   // gradient for each face.
				   // At the end of the function, we again
				   // loop over the cells and collect the
				   // contributions of the different faces
				   // of the cell.
				   // 
				   // The initial values for all faces
				   // are set to -10e20. It would cost
				   // a lot of time to synchronise the
				   // initialisation (i.e. the
				   // creation of new keys) of the map
				   // in multithreaded mode. Negative
				   // value indicates that the face
				   // has not yet been processed.
  vector<double> default_face_integrals (n_solution_vectors, -10e20);
  FaceIntegrals face_integrals;
  for (active_cell_iterator cell=dof_handler.begin_active(); cell!=dof_handler.end(); ++cell)
    for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
      face_integrals[cell->face(face_no)] = default_face_integrals;


				   // all the data needed in the error
				   // estimator by each of the threads
				   // is gathered in the following
				   // stuctures
				   //
				   // note that if no component mask
				   // was given, then treat all
				   // components
  vector<Data*> data_structures (n_threads);
  for (unsigned int i=0; i<n_threads; ++i)
    data_structures[i] = new Data (dof_handler,
				   quadrature,
				   neumann_bc,
				   solutions,
				   ((component_mask.size() == 0)    ?
				    vector<bool>(dof_handler.get_fe().n_components(),
						 true) :
				    component_mask),
				   coefficients,
				   n_threads,
				   face_integrals);
  
				   // split all cells into threads if
				   // multithreading is used and run
				   // the whole thing
  Threads::ThreadManager thread_manager;
  for (unsigned int i=0; i<n_threads; ++i)
    Threads::spawn (thread_manager,
		    Threads::encapsulate (&KellyErrorEstimator<dim>::estimate_some)
		    .collect_args (*data_structures[i], i));
  thread_manager.wait();

				   // delete the structures for the
				   // different threads again. the
				   // results are in the
				   // face_integrals object
  for (unsigned int i=0; i<n_threads; ++i)
    {
      delete data_structures[i];
      data_structures[i] = 0;
    };
  
  
				   // finally add up the contributions of the
				   // faces for each cell
  
				   // reserve one slot for each cell and set
				   // it to zero
  for (unsigned int n=0; n<n_solution_vectors; ++n)
    {
      (*errors[n]).reinit (dof_handler.get_tria().n_active_cells());
      for (unsigned int i=0;i<dof_handler.get_tria().n_active_cells();++i)
	(*errors[n])(i)=0;
    };

  unsigned int present_cell=0;
  
  for (active_cell_iterator cell=dof_handler.begin_active();
       cell!=dof_handler.end();
       ++cell, ++present_cell)
    {
				       // loop over all faces of this cell
      for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell;
	   ++face_no)
	{
	  Assert(face_integrals.find(cell->face(face_no)) != face_integrals.end(),
		 ExcInternalError());
	  const double factor = cell->diameter() / 24;
	  
	  for (unsigned int n=0; n<n_solution_vectors; ++n)
	    (*errors[n])(present_cell) += (face_integrals[cell->face(face_no)][n] *
					   factor);
	};

      for (unsigned int n=0; n<n_solution_vectors; ++n)
	(*errors[n])(present_cell) = sqrt((*errors[n])(present_cell));
    };
};

#endif


#if deal_II_dimension == 1

template <>
void KellyErrorEstimator<1>::integrate_over_regular_face (Data &,
							  const active_cell_iterator &,
							  const unsigned int      ,
							  FEFaceValues<1>        &,
							  FEFaceValues<1>        &)
{
  Assert (false, ExcInternalError());
};


template <>
void KellyErrorEstimator<1>::
integrate_over_irregular_face (Data &,
			       const active_cell_iterator &,
			       const unsigned int          ,
			       FEFaceValues<1>            &,
			       FESubfaceValues<1>         &)
{
  Assert (false, ExcInternalError());
};

#endif


template <int dim>
void KellyErrorEstimator<dim>::
integrate_over_regular_face (Data                       &data,
			     const active_cell_iterator &cell,
			     const unsigned int          face_no,
			     FEFaceValues<dim>          &fe_face_values_cell,
			     FEFaceValues<dim>          &fe_face_values_neighbor)
{
  const typename DoFHandler<dim>::face_iterator face = cell->face(face_no);
  const unsigned int n_q_points         = data.quadrature.n_quadrature_points,
		     n_components       = data.dof_handler.get_fe().n_components(),
		     n_solution_vectors = data.n_solution_vectors;
  
  
				   // initialize data of the restriction
				   // of this cell to the present face
  fe_face_values_cell.reinit (cell, face_no);
  
				   // get gradients of the finite element
				   // function on this cell
  for (unsigned int n=0; n<n_solution_vectors; ++n)
    fe_face_values_cell.get_function_grads (*data.solutions[n], data.psi[n]);
  
				   // now compute over the other side of
				   // the face
  if (face->at_boundary() == false)
				     // internal face; integrate jump
				     // of gradient across this face
    {
      Assert (cell->neighbor(face_no).state() == valid,
	      ExcInternalError());      
      
      const active_cell_iterator neighbor = cell->neighbor(face_no);
      
				       // find which number the current
				       // face has relative to the neighboring
				       // cell
      const unsigned int neighbor_neighbor = cell->neighbor_of_neighbor (face_no);
      Assert (neighbor_neighbor<GeometryInfo<dim>::faces_per_cell, ExcInternalError());
      
				       // get restriction of finite element
				       // function of @p{neighbor} to the
				       // common face.
      fe_face_values_neighbor.reinit (neighbor, neighbor_neighbor);
      
				       // get gradients on neighbor cell
      for (unsigned int n=0; n<n_solution_vectors; ++n)
	{
	  fe_face_values_neighbor.get_function_grads (*data.solutions[n],
						      data.neighbor_psi[n]);
      
					   // compute the jump in the gradients
	  for (unsigned int component=0; component<n_components; ++component)
	    for (unsigned int p=0; p<n_q_points; ++p)
	      data.psi[n][p][component] -= data.neighbor_psi[n][p][component];
	};
    };


				   // now psi contains the following:
				   // - for an internal face, psi=[grad u]
				   // - for a neumann boundary face,
				   //   psi=grad u
				   // each component being the
				   // mentioned value at one of the
				   // quadrature points
  
				   // next we have to multiply this with
				   // the normal vector. Since we have
				   // taken the difference of gradients
				   // for internal faces, we may chose
				   // the normal vector of one cell,
				   // taking that of the neighbor
				   // would only change the sign. We take
				   // the outward normal.
  
  data.normal_vectors=fe_face_values_cell.get_normal_vectors();
  
  for (unsigned int n=0; n<n_solution_vectors; ++n)
    for (unsigned int component=0; component<n_components; ++component)
      for (unsigned int point=0; point<n_q_points; ++point)
	data.phi[n][point][component] = data.psi[n][point][component]*
					data.normal_vectors[point];
  
				   // if a coefficient was given: use that
				   // to scale the jump in the gradient
  if (data.coefficients != 0)
    {
				       // scalar coefficient
      if (data.coefficients->n_components == 1)
	{
	  
	  data.coefficients->value_list (fe_face_values_cell.get_quadrature_points(),
					 data.coefficient_values1);
	  for (unsigned int n=0; n<n_solution_vectors; ++n)
	    for (unsigned int component=0; component<n_components; ++component)
	      for (unsigned int point=0; point<n_q_points; ++point)
		data.phi[n][point][component] *=
		  data.coefficient_values1[point];
	}
      else
					   // vector-valued coefficient
	{
	  data.coefficients->vector_value_list (fe_face_values_cell.get_quadrature_points(),
						data.coefficient_values);
	  for (unsigned int n=0; n<n_solution_vectors; ++n)
	    for (unsigned int component=0; component<n_components; ++component)
	      for (unsigned int point=0; point<n_q_points; ++point)
		data.phi[n][point][component] *=
		  data.coefficient_values[point](component);
	};
    };


  if (face->at_boundary() == true)
				     // neumann boundary face. compute
				     // difference between normal
				     // derivative and boundary function
    {
      const unsigned char boundary_indicator = face->boundary_indicator();
      
      Assert (data.neumann_bc.find(boundary_indicator) != data.neumann_bc.end(),
	      ExcInternalError ());
				       // get the values of the boundary
				       // function at the quadrature
				       // points
      
      vector<Vector<double> > g(n_q_points, Vector<double>(n_components));
      data.neumann_bc.find(boundary_indicator)->second
	->vector_value_list (fe_face_values_cell.get_quadrature_points(),
			     g);
      
      for (unsigned int n=0; n<n_solution_vectors; ++n)
	for (unsigned int component=0; component<n_components; ++component)
	  for (unsigned int point=0; point<n_q_points; ++point)
	    data.phi[n][point][component] -= g[point](component);
    };


				   // now phi contains the following:
				   // - for an internal face, phi=[a du/dn]
				   // - for a neumann boundary face,
				   //   phi=a du/dn-g
				   // each component being the
				   // mentioned value at one of the
				   // quadrature points

  data.JxW_values = fe_face_values_cell.get_JxW_values();
  
				   // take the square of the phi[i]
				   // for integration, and sum up
  vector<double> face_integral (n_solution_vectors, 0);
  for (unsigned int n=0; n<n_solution_vectors; ++n)
    for (unsigned int component=0; component<n_components; ++component)
      if (data.component_mask[component] == true)
	for (unsigned int p=0; p<n_q_points; ++p)
	  face_integral[n] += sqr(data.phi[n][p][component]) *
			      data.JxW_values[p];
  
  data.face_integrals[face] = face_integral;
};



template <int dim>
void KellyErrorEstimator<dim>::
integrate_over_irregular_face (Data                       &data,
			       const active_cell_iterator &cell,
			       const unsigned int          face_no,
			       FEFaceValues<dim>          &fe_face_values,
			       FESubfaceValues<dim>       &fe_subface_values)
{
  const typename DoFHandler<dim>::cell_iterator neighbor = cell->neighbor(face_no);
  const unsigned int n_q_points         = data.quadrature.n_quadrature_points,
		     n_components       = data.dof_handler.get_fe().n_components(),
		     n_solution_vectors = data.n_solution_vectors;

  Assert (neighbor.state() == valid, ExcInternalError());
  Assert (neighbor->has_children(), ExcInternalError());
				   // set up a vector of the gradients
				   // of the finite element function
				   // on this cell at the quadrature
				   // points
				   //
				   // let psi be a short name for
				   // [a grad u_h], where the second
				   // index be the component of the
				   // finite element, and the first
				   // index the number of the
				   // quadrature point
  
				   // store which number @p{cell} has in the
				   // list of neighbors of @p{neighbor}
  const unsigned int neighbor_neighbor = cell->neighbor_of_neighbor (face_no);
  Assert (neighbor_neighbor<GeometryInfo<dim>::faces_per_cell, ExcInternalError());
  
				   // loop over all subfaces
  for (unsigned int subface_no=0; subface_no<GeometryInfo<dim>::subfaces_per_face;
       ++subface_no)
    {
				       // get an iterator pointing to the
				       // cell behind the present subface
      const active_cell_iterator neighbor_child
	= neighbor->child(GeometryInfo<dim>::
			  child_cell_on_face(neighbor_neighbor,subface_no));
      Assert (neighbor_child->face(neighbor_neighbor) ==
	      cell->face(face_no)->child(subface_no),
	      ExcInternalError());
      Assert (!neighbor->child(GeometryInfo<dim>::
			       child_cell_on_face(neighbor_neighbor,subface_no))->has_children(),
	      ExcInternalError());
      
				       // restrict the finite element on the
				       // present cell to the subface and
				       // store the gradient of the solution
				       // in psi
      fe_subface_values.reinit (cell, face_no, subface_no);

      for (unsigned int n=0; n<n_solution_vectors; ++n)
	fe_subface_values.get_function_grads (*data.solutions[n], data.psi[n]);

				       // restrict the finite element on the
				       // neighbor cell to the common @p{subface}.
				       // store the gradient in @p{neighbor_psi}
      
      fe_face_values.reinit (neighbor_child, neighbor_neighbor);

      for (unsigned int n=0; n<n_solution_vectors; ++n)
	fe_face_values.get_function_grads (*data.solutions[n], data.neighbor_psi[n]);
      
				       // compute the jump in the gradients
      for (unsigned int n=0; n<n_solution_vectors; ++n)
	for (unsigned int component=0; component<n_components; ++component)
	  for (unsigned int p=0; p<n_q_points; ++p)
	    data.psi[n][p][component] -=
	      data.neighbor_psi[n][p][component];

				       // note that unlike for the
				       // case of regular faces
				       // (treated in the other
				       // function of this class), we
				       // have not to take care of
				       // boundary faces here, since
				       // they always are regular.
      
				       // next we have to multiply this with
				       // the normal vector. Since we have
				       // taken the difference of gradients
				       // for internal faces, we may chose
				       // the normal vector of one cell,
				       // taking that of the neighbor
				       // would only change the sign. We take
				       // the outward normal.
				       //
				       // let phi be the name of the integrand
     
      data.normal_vectors=fe_face_values.get_normal_vectors();


      for (unsigned int n=0; n<n_solution_vectors; ++n)
	for (unsigned int component=0; component<n_components; ++component)
	  for (unsigned int point=0; point<n_q_points; ++point)
	    data.phi[n][point][component] = (data.psi[n][point][component]*
					     data.normal_vectors[point]);
      
				       // if a coefficient was given: use that
				       // to scale the jump in the gradient
      if (data.coefficients != 0)
	{
					   // scalar coefficient
	  if (data.coefficients->n_components == 1)
	    {
	      data.coefficients->value_list (fe_face_values.get_quadrature_points(),
					     data.coefficient_values1);
	      for (unsigned int n=0; n<n_solution_vectors; ++n)
		for (unsigned int component=0; component<n_components; ++component)
		  for (unsigned int point=0; point<n_q_points; ++point)
		    data.phi[n][point][component] *=
		      data.coefficient_values1[point];
	    }
	  else
					     // vector-valued coefficient
	    {
	      data.coefficients->vector_value_list (fe_face_values.get_quadrature_points(),
						    data.coefficient_values);
	      for (unsigned int n=0; n<n_solution_vectors; ++n)
		for (unsigned int component=0; component<n_components; ++component)
		  for (unsigned int point=0; point<n_q_points; ++point)
		    data.phi[n][point][component] *=
		      data.coefficient_values[point](component);
	    };
	};

				       // get the weights for the
				       // integration. note that it
				       // does not matter whether we
				       // take the JxW values from the
				       // fe_face_values or the
				       // fe_subface_values, as the
				       // first is on the small
				       // neighbor cell, while the
				       // latter is on the refined
				       // face of the big cell here
      data.JxW_values = fe_face_values.get_JxW_values();
      
				       // take the square of the phi[i]
				       // for integration, and sum up
      vector<double> face_integral (n_solution_vectors, 0);
      for (unsigned int n=0; n<n_solution_vectors; ++n)
	for (unsigned int component=0; component<n_components; ++component)
	  if (data.component_mask[component] == true)
	    for (unsigned int p=0; p<n_q_points; ++p)
	      face_integral[n] += sqr(data.phi[n][p][component]) *
				  data.JxW_values[p];

      data.face_integrals[neighbor_child->face(neighbor_neighbor)] = face_integral;
    };


				   // finally loop over all subfaces to
				   // collect the contributions of the
				   // subfaces and store them with the
				   // mother face
  vector<double> sum (n_solution_vectors, 0);
  typename DoFHandler<dim>::face_iterator face = cell->face(face_no);
  for (unsigned int subface_no=0; subface_no<GeometryInfo<dim>::subfaces_per_face;
       ++subface_no) 
    {
      Assert (data.face_integrals.find(face->child(subface_no)) !=
	      data.face_integrals.end(),
	      ExcInternalError());
      Assert (data.face_integrals[face->child(subface_no)][0] >= 0,
	      ExcInternalError());

      for (unsigned int n=0; n<n_solution_vectors; ++n)
	sum[n] += data.face_integrals[face->child(subface_no)][n];
    };

  data.face_integrals[face] = sum;
};


// explicit instantiations

template class KellyErrorEstimator<deal_II_dimension>;
