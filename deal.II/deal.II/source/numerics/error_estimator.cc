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

				 // if multithreaded include
				 // ThreadManager
#ifdef DEAL_II_USE_MT
#  include <base/thread_management.h>
#endif



#if deal_II_dimension == 1

template <>
KellyErrorEstimator<1>::Data::Data(const DoFHandler<1>     &,
				   const Quadrature<0>     &,
				   const FunctionMap       &,
				   const Vector<double>    &,
				   vector<bool>            ,
				   const Function<1>       *,
				   unsigned int            ):
		dof(* static_cast <const DoFHandler<1> *> (0)),
		quadrature(* static_cast <const Quadrature<0> *> (0)),
		neumann_bc(* static_cast <const FunctionMap *> (0)),
		solution(* static_cast <const Vector<double> *> (0))
{
  Assert (false, ExcInternalError());
}

#else

template <int dim>
KellyErrorEstimator<dim>::Data::Data(const DoFHandler<dim>   &dof,
				     const Quadrature<dim-1> &quadrature,
				     const FunctionMap       &neumann_bc,
				     const Vector<double>    &solution,
				     vector<bool>            component_mask_,
				     const Function<dim>     *coefficients,
				     unsigned int            n_threads):
		dof(dof),
		quadrature(quadrature),
		neumann_bc(neumann_bc),
		solution(solution),
		coefficients(coefficients),
		n_threads(n_threads)
{
  n_components = dof.get_fe().n_components();
  
				   // if no mask given: treat all components
  component_mask = ((component_mask_.size() == 0)    ?
		    vector<bool>(n_components, true) :
		    component_mask_);
  
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
  
				   // the last cell, often needed
  endc=dof.end();
  
  n_q_points=quadrature.n_quadrature_points;
  
				   // Init the size of a lot of vectors
				   // needed in the calculations once
				   // per thread.
  phi.resize(n_threads);
  psi.resize(n_threads);
  neighbor_psi.resize(n_threads);
  normal_vectors.resize(n_threads);
  coefficient_values1.resize(n_threads);
  coefficient_values.resize(n_threads);
  JxW_values.resize(n_threads);
  
  for (unsigned int t=0;t<n_threads;++t)
    {
      phi[t].resize(n_q_points);
      psi[t].resize(n_q_points);
      neighbor_psi[t].resize(n_q_points);
      normal_vectors[t].resize(n_q_points);
      coefficient_values1[t].resize(n_q_points);
      coefficient_values[t].resize(n_q_points);
      JxW_values[t].resize(n_q_points);
      
      for (unsigned int qp=0;qp<n_q_points;++qp)
	{
	  phi[t][qp].resize(n_components);
	  psi[t][qp].resize(n_components);
	  neighbor_psi[t][qp].resize(n_components);
	  coefficient_values[t][qp].reinit(n_components);
	}
    }
}

#endif


inline static double sqr (const double x) {
  return x*x;
};


#if deal_II_dimension == 1

template <>
void KellyErrorEstimator<1>::estimate_some (Data &, const unsigned int)
{
  Assert (false, ExcInternalError() );
}


template <>
void KellyErrorEstimator<1>::estimate (const DoFHandler<1>  &dof,
				       const Quadrature<0>  &,
				       const FunctionMap    &neumann_bc,
				       const Vector<double> &solution,
				       Vector<float>        &error,
				       const vector<bool>   &component_mask_,
				       const Function<1>    *coefficient,
				       const unsigned int)
{
  const unsigned int n_components = dof.get_fe().n_components();

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
  error.reinit (dof.get_tria().n_active_cells());

				   // fields to get the gradients on
				   // the present and the neighbor cell.
				   //
				   // for the neighbor gradient, we
				   // need several auxiliary fields,
				   // depending on the way we get it
				   // (see below)
  vector<vector<Tensor<1,1> > > gradients_here (2, vector<Tensor<1,1> >(n_components));
  vector<vector<Tensor<1,1> > > gradients_neighbor (gradients_here);
  Vector<double>                grad_neighbor (n_components);

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
  FEValues<1> fe_values (dof.get_fe(), quadrature, update_gradients);
  DoFHandler<1>::active_cell_iterator cell = dof.begin_active();
  for (unsigned int cell_index=0; cell != dof.end(); ++cell, ++cell_index)
    {
      error(cell_index) = 0;
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
	  fe_values.get_function_grads (solution, gradients_here);

	  if (neighbor.state() == valid)
	    {
	      fe_values.reinit (neighbor);
	      fe_values.get_function_grads (solution, gradients_neighbor);

					       // extract the
					       // gradients of all the
					       // components. [0]
					       // means: x-derivative,
					       // which is the only
					       // one here
	      for (unsigned int c=0; c<n_components; ++c)
		grad_neighbor(c) = gradients_neighbor[n==0 ? 1 : 0][c][0];
	    }
	  else
	    if (neumann_bc.find(n) != neumann_bc.end())
					       // if Neumann b.c., then fill
					       // the gradients field which
					       // will be used later on.
	      neumann_bc.find(n)->second->vector_value(cell->vertex(0),
						       grad_neighbor);
	    else
					       // fill with zeroes.
	      grad_neighbor.clear ();

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


	  for (unsigned int component=0; component<n_components; ++component)
	    if (component_mask[component] == true)
	      {
						 // get gradient
						 // here. [0] means
						 // x-derivative
						 // (there is not
						 // other in 1d)
		const double grad_here = gradients_here[n][component][0];
	    
		const double jump = ((grad_here - grad_neighbor(component)) *
				     coefficient_values(component));
		error(cell_index) += jump*jump * cell->diameter();
	      };
	};
      
      error(cell_index) = sqrt(error(cell_index));
    };
};


#else // #if deal_II_dimension !=1


template <int dim>
void KellyErrorEstimator<dim>::estimate_some (Data               &data,
					      const unsigned int  this_thread) 
{
  
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
  FEFaceValues<dim> fe_face_values_cell (data.dof.get_fe(),
					 data.quadrature,
					 UpdateFlags(update_gradients      |
						     update_JxW_values     |
						     ((!data.neumann_bc.empty() ||
						       (data.coefficients != 0))  ?
						      update_q_points : 0) |
						     update_normal_vectors)); 
  FEFaceValues<dim> fe_face_values_neighbor (data.dof.get_fe(),
					     data.quadrature,
					     update_gradients);
  FESubfaceValues<dim> fe_subface_values (data.dof.get_fe(),
					  data.quadrature,
					  update_gradients);


  DoFHandler<dim>::active_cell_iterator cell=data.dof.begin_active();

				   // calculate the start cell for this
				   // thread. the enumeration is choosen
				   // in this strange way to generate a
				   // "random" distribution of the cells.
				   // if the sequence of the iterator would
				   // be used, the threads would take widely
				   // spread times to calculate their cells.
  for (unsigned int t=0;t<this_thread;++t,++cell);
				   // loop over all cells for this thread
				   // the iteration of cell is done at the end
  for (; cell!=data.endc; )
    {
      
				       // loop over all faces of this cell
      for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
	{
					   // if we already visited this
					   // face: do nothing
	  if (data.face_integrals[cell->face(face_no)] >=0)
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
	  if ((boundary_indicator != 255) &&
	      data.neumann_bc.find(boundary_indicator)==data.neumann_bc.end()) 
	    {
	      data.face_integrals[cell->face(face_no)] = 0;
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
					 this_thread,
					 cell, face_no,
					 fe_face_values_cell,
					 fe_face_values_neighbor);
	  
	  else
					     // otherwise we need to do some
					     // special computations which do
					     // not fit into the framework of
					     // the above function
	    integrate_over_irregular_face (data,
					   this_thread,cell, face_no,
					   fe_face_values_cell,
					   fe_subface_values);
	};

				       // next cell in this thread
      for (unsigned int t=0;((t<data.n_threads)&&(cell!=data.endc));++t,++cell)
	{};
    };
};



template <int dim>
void KellyErrorEstimator<dim>::estimate (const DoFHandler<dim>   &dof,
					 const Quadrature<dim-1> &quadrature,
					 const FunctionMap       &neumann_bc,
					 const Vector<double>    &solution,
					 Vector<float>           &error,
					 const vector<bool>      &component_mask,
					 const Function<dim>     *coefficients,
					 unsigned int            n_threads)
{
				   // if NOT multithreaded, set n_threads to one
#ifndef DEAL_II_USE_MT
  n_threads = 1;
#endif
  
				   // all the data needed in the error-
				   // estimator is gathered in this stuct.
  KellyErrorEstimator<dim>::Data data (dof,
				       quadrature,
				       neumann_bc,
				       solution,
				       component_mask,
				       coefficients,
				       n_threads);
  				   // map of integrals indexed by
				   // the corresponding face. In this map
				   // we store the integrated jump of the
				   // gradient for each face. By doing so,
				   // we can check whether we have already
				   // integrated along this face by testing
				   // whether the respective face is already
				   // a key in this map.
				   // At the end of the function, we again
				   // loop over the cells and collect the
				   // conrtibutions of the different faces
				   // of the cell.
				   // the values for all faces are set to
				   // -10e20. It would cost a lot of time
				   // to synchronise the initialisation
				   // of the map in multithreaded mode.
				   // negative value indicates that the
				   // face is not calculated.
  
  for (active_cell_iterator cell=data.dof.begin_active(); cell!=data.endc; ++cell)
    for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
      data.face_integrals[cell->face(face_no)]=-10e20;


				   // split all cells into threads
				   // if multithreading is used
#ifdef DEAL_II_USE_MT

  Threads::ThreadManager thread_manager;
  for (unsigned int i=0;i<data.n_threads; ++i)
    Threads::spawn (thread_manager,
		    Threads::encapsulate (&KellyErrorEstimator<dim>::estimate_some)
		    .collect_args (data, i));
  thread_manager.wait();
  
				   // ... ifdef DEAL_II_USE_MT
#else
				   // just one thread, calculate
				   // error on all cells
  KellyErrorEstimator::estimate_some(data,0);
  
#endif


				   // finally add up the contributions of the
				   // faces for each cell
  
				   // reserve one slot for each cell and set
				   // it to zero
  error.reinit (data.dof.get_tria().n_active_cells());
  for (unsigned int i=0;i<data.dof.get_tria().n_active_cells();++i)
    error(i)=0;

  unsigned int present_cell=0;
  
  for (active_cell_iterator cell=data.dof.begin_active();
       cell!=data.endc;
       ++cell, ++present_cell)
    {
				       // loop over all faces of this cell
      for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell;
	   ++face_no)
	{
	  Assert(data.face_integrals.find(cell->face(face_no)) != data.face_integrals.end(),
		 ExcInternalError());
	  error(present_cell) += (data.face_integrals[cell->face(face_no)] *
				  cell->diameter() / 24);
	};
      error(present_cell) = sqrt(error(present_cell));
    };
};

#endif


#if deal_II_dimension == 1

template <>
void KellyErrorEstimator<1>::integrate_over_regular_face (Data &,
							  const unsigned int,
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
			       const unsigned int,
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
			     const unsigned int          this_thread,
			     const active_cell_iterator &cell,
			     const unsigned int          face_no,
			     FEFaceValues<dim>          &fe_face_values_cell,
			     FEFaceValues<dim>          &fe_face_values_neighbor)
{
  const DoFHandler<dim>::face_iterator face = cell->face(face_no);
  
				   // initialize data of the restriction
				   // of this cell to the present face
  fe_face_values_cell.reinit (cell, face_no);
  
				   // get gradients of the finite element
				   // function on this cell
  fe_face_values_cell.get_function_grads (data.solution, data.psi[this_thread]);
  
				   // now compute over the other side of
				   // the face
  if (face->at_boundary() == false)
				     // internal face; integrate jump
				     // of gradient across this face
    {
      Assert (cell->neighbor(face_no).state() == valid,
	      ExcInternalError());      
      
      const DoFHandler<dim>::active_cell_iterator neighbor = cell->neighbor(face_no);
      
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
      fe_face_values_neighbor.get_function_grads (data.solution,
						  data.neighbor_psi[this_thread]);
      
				       // compute the jump in the gradients
      for (unsigned int component=0; component<data.n_components; ++component)
	for (unsigned int p=0; p<data.n_q_points; ++p)
	  data.psi[this_thread][p][component] -= data.neighbor_psi[this_thread][p][component];
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
  
  data.normal_vectors[this_thread]=fe_face_values_cell.get_normal_vectors();
  
  for (unsigned int component=0; component<data.n_components; ++component)
    for (unsigned int point=0; point<data.n_q_points; ++point)
      data.phi[this_thread][point][component] = data.psi[this_thread][point][component]*
					   data.normal_vectors[this_thread][point];
  
				   // if a coefficient was given: use that
				   // to scale the jump in the gradient
  if (data.coefficients != 0)
    {
				       // scalar coefficient
      if (data.coefficients->n_components == 1)
	{
	  
	  data.coefficients->value_list (fe_face_values_cell.get_quadrature_points(),
					 data.coefficient_values1[this_thread]);
	    for (unsigned int component=0; component<data.n_components; ++component)
	      for (unsigned int point=0; point<data.n_q_points; ++point)
		data.phi[this_thread][point][component] *=
		  data.coefficient_values1[this_thread][point];
	}
      else
					   // vector-valued coefficient
	{
	  data.coefficients->vector_value_list (fe_face_values_cell.get_quadrature_points(),
					   data.coefficient_values[this_thread]);
	  for (unsigned int component=0; component<data.n_components; ++component)
	    for (unsigned int point=0; point<data.n_q_points; ++point)
	      data.phi[this_thread][point][component] *=
		data.coefficient_values[this_thread][point](component);
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
      
      vector<Vector<double> > g(data.n_q_points, Vector<double>(data.n_components));
      data.neumann_bc.find(boundary_indicator)->second
	->vector_value_list (fe_face_values_cell.get_quadrature_points(),
			     g);
      
      for (unsigned int component=0; component<data.n_components; ++component)
	for (unsigned int point=0; point<data.n_q_points; ++point)
	  data.phi[this_thread][point][component] -= g[point](component);
    };


				   // now phi contains the following:
				   // - for an internal face, phi=[a du/dn]
				   // - for a neumann boundary face,
				   //   phi=a du/dn-g
				   // each component being the
				   // mentioned value at one of the
				   // quadrature points

  data.JxW_values[this_thread] = fe_face_values_cell.get_JxW_values();
  
				   // take the square of the phi[i]
				   // for integration, and sum up
  double face_integral = 0;
  for (unsigned int component=0; component<data.n_components; ++component)
    if (data.component_mask[component] == true)
      for (unsigned int p=0; p<data.n_q_points; ++p)
	face_integral += sqr(data.phi[this_thread][p][component]) *
			 data.JxW_values[this_thread][p];
  
  data.face_integrals[face] = face_integral;
};



template <int dim>
void KellyErrorEstimator<dim>::
integrate_over_irregular_face (Data                       &data,
			       const unsigned int          this_thread,
			       const active_cell_iterator &cell,
			       const unsigned int          face_no,
			       FEFaceValues<dim>          &fe_face_values,
			       FESubfaceValues<dim>       &fe_subface_values)
{
  const DoFHandler<dim>::cell_iterator neighbor = cell->neighbor(face_no);

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
      fe_subface_values.get_function_grads (data.solution, data.psi[this_thread]);

				       // restrict the finite element on the
				       // neighbor cell to the common @p{subface}.
				       // store the gradient in @p{neighbor_psi}
      
      fe_face_values.reinit (neighbor_child, neighbor_neighbor);
      fe_face_values.get_function_grads (data.solution, data.neighbor_psi[this_thread]);
      
				       // compute the jump in the gradients
      for (unsigned int component=0; component<data.n_components; ++component)
	for (unsigned int p=0; p<data.n_q_points; ++p)
	  data.psi[this_thread][p][component] -=
	    data.neighbor_psi[this_thread][p][component];

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
     
      data.normal_vectors[this_thread]=fe_face_values.get_normal_vectors();


      for (unsigned int component=0; component<data.n_components; ++component)
	for (unsigned int point=0; point<data.n_q_points; ++point)
	  data.phi[this_thread][point][component] =
	    data.psi[this_thread][point][component]*
	    data.normal_vectors[this_thread][point];
      
				       // if a coefficient was given: use that
				       // to scale the jump in the gradient
      if (data.coefficients != 0)
	{
					   // scalar coefficient
	  if (data.coefficients->n_components == 1)
	    {
	      data.coefficients->value_list (fe_face_values.get_quadrature_points(),
					     data.coefficient_values1[this_thread]);
	      for (unsigned int component=0; component<data.n_components; ++component)
		for (unsigned int point=0; point<data.n_q_points; ++point)
		  data.phi[this_thread][point][component] *=
		    data.coefficient_values1[this_thread][point];
	    }
	  else
					     // vector-valued coefficient
	    {
	      data.coefficients->vector_value_list (fe_face_values.get_quadrature_points(),
						    data.coefficient_values[this_thread]);
	      for (unsigned int component=0; component<data.n_components; ++component)
		for (unsigned int point=0; point<data.n_q_points; ++point)
		  data.phi[this_thread][point][component] *=
		    data.coefficient_values[this_thread][point](component);
	    };
	};
      
      data.JxW_values[this_thread] = fe_face_values.get_JxW_values();
      
				       // take the square of the phi[i]
				       // for integration, and sum up
      double face_integral = 0;
      for (unsigned int component=0; component<data.n_components; ++component)
	if (data.component_mask[component] == true)
	  for (unsigned int p=0; p<data.n_q_points; ++p)
	    face_integral += sqr(data.phi[this_thread][p][component]) *
			     data.JxW_values[this_thread][p];

      data.face_integrals[neighbor_child->face(neighbor_neighbor)] = face_integral;
    };


				   // finally loop over all subfaces to
				   // collect the contributions of the
				   // subfaces and store them with the
				   // mother face
  double sum=0;
  DoFHandler<dim>::face_iterator face = cell->face(face_no);
  for (unsigned int subface_no=0; subface_no<GeometryInfo<dim>::subfaces_per_face;
       ++subface_no) 
    {
      Assert (data.face_integrals.find(face->child(subface_no)) !=
	      data.face_integrals.end(),
	      ExcInternalError());
      Assert (data.face_integrals[face->child(subface_no)]>=0,
	      ExcInternalError());
      sum += data.face_integrals[face->child(subface_no)];
    };

  data.face_integrals[face] = sum;
};


// explicit instantiations

template class KellyErrorEstimator<deal_II_dimension>;
