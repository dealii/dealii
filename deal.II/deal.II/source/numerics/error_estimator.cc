//----------------------------  error_estimator.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  error_estimator.cc  ---------------------------


#include <base/thread_management.h>
#include <base/quadrature.h>
#include <base/quadrature_lib.h>
#include <lac/vector.h>
#include <lac/block_vector.h>
#include <lac/petsc_vector.h>
#include <grid/tria_iterator.h>
#include <grid/geometry_info.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <fe/fe.h>
#include <fe/fe_values.h>
#include <fe/fe_update_flags.h>
#include <fe/mapping_q1.h>
#include <numerics/error_estimator.h>

#include <numeric>
#include <algorithm>
#include <cmath>
#include <vector>



static
inline
double sqr (const double x)
{
  return x*x;
}


template <typename CellIterator>
inline
void advance_by_n (CellIterator &cell,
                   const unsigned int n)
{
  for (unsigned int t=0;
       ((t<n) && (cell!=cell->get_dof_handler().end()));
       ++t, ++cell);
}



#if deal_II_dimension == 1


template <typename InputVector>
void
KellyErrorEstimator<1>::
estimate (const Mapping<1>      &mapping,
          const DoFHandler<1>   &dof_handler,
          const Quadrature<0> &quadrature,
          const FunctionMap<1>::type &neumann_bc,
          const InputVector       &solution,
          Vector<float>           &error,
          const std::vector<bool> &component_mask,
          const Function<1>     *coefficients,
          const unsigned int       n_threads,
          const unsigned int       subdomain_id)
{
				   // just pass on to the other function
  const std::vector<const InputVector *> solutions (1, &solution);
  std::vector<Vector<float>*>              errors (1, &error);
  estimate (mapping, dof_handler, quadrature, neumann_bc, solutions, errors,
	    component_mask, coefficients, n_threads, subdomain_id);
}


template <typename InputVector>
void
KellyErrorEstimator<1>::
estimate (const DoFHandler<1>   &dof_handler,
          const Quadrature<0> &quadrature,
          const FunctionMap<1>::type &neumann_bc,
          const InputVector       &solution,
          Vector<float>           &error,
          const std::vector<bool> &component_mask,
          const Function<1>     *coefficients,
          const unsigned int       n_threads,
          const unsigned int       subdomain_id)
{
  Assert (DEAL_II_COMPAT_MAPPING, ExcCompatibility("mapping"));
  static const MappingQ1<1> mapping;
  estimate(mapping, dof_handler, quadrature, neumann_bc, solution,
	   error, component_mask, coefficients, n_threads, subdomain_id);
}



template <typename InputVector>
void
KellyErrorEstimator<1>::
estimate (const DoFHandler<1>   &dof_handler,
          const Quadrature<0> &quadrature,
          const FunctionMap<1>::type &neumann_bc,
          const std::vector<const InputVector*> &solutions,
          std::vector<Vector<float>*> &errors,
          const std::vector<bool> &component_mask,
          const Function<1>     *coefficients,
          const unsigned int       n_threads,
          const unsigned int       subdomain_id)
{
  Assert (DEAL_II_COMPAT_MAPPING, ExcCompatibility("mapping"));
  static const MappingQ1<1> mapping;
  estimate(mapping, dof_handler, quadrature, neumann_bc, solutions,
	   errors, component_mask, coefficients, n_threads, subdomain_id);
}



template <typename InputVector>
void KellyErrorEstimator<1>::
estimate (const Mapping<1>                    &mapping,
          const DoFHandler<1>                 &dof_handler,
          const Quadrature<0>                 &,
          const FunctionMap<1>::type          &neumann_bc,
          const std::vector<const InputVector *> &solutions,
          std::vector<Vector<float>*>              &errors,
          const std::vector<bool>                  &component_mask_,
          const Function<1>                   *coefficient,
          const unsigned int,
          const unsigned int                  subdomain_id)
{
  const unsigned int n_components       = dof_handler.get_fe().n_components();
  const unsigned int n_solution_vectors = solutions.size();

				   // sanity checks
  Assert (neumann_bc.find(255) == neumann_bc.end(),
	  ExcInvalidBoundaryIndicator());
  
  for (FunctionMap<1>::type::const_iterator i=neumann_bc.begin();
       i!=neumann_bc.end(); ++i)
    Assert (i->second->n_components == n_components, ExcInvalidBoundaryFunction());  
  
  Assert ((component_mask_.size() == 0) ||
          (component_mask_.size() == n_components), ExcInvalidComponentMask());
  Assert ((component_mask_.size() == 0) ||
          (std::count(component_mask_.begin(), component_mask_.end(),
                      true) > 0),
	  ExcInvalidComponentMask());

  Assert ((coefficient == 0) ||
	  (coefficient->n_components == n_components) ||
	  (coefficient->n_components == 1),
	  ExcInvalidCoefficient());  
  
  Assert (solutions.size() > 0,
	  ExcNoSolutions());
  Assert (solutions.size() == errors.size(),
	  ExcIncompatibleNumberOfElements(solutions.size(), errors.size()));
  for (unsigned int n=0; n<solutions.size(); ++n)
    Assert (solutions[n]->size() == dof_handler.n_dofs(),
	    ExcInvalidSolutionVector());
  
				   // if no mask given: treat all components
  std::vector<bool> component_mask ((component_mask_.size() == 0)    ?
				    std::vector<bool>(n_components, true) :
				    component_mask_);
  Assert (component_mask.size() == n_components, ExcInvalidComponentMask());
  Assert (count(component_mask.begin(), component_mask.end(), true) > 0,
	  ExcInvalidComponentMask());
  
  Assert ((coefficient == 0) ||
	  (coefficient->n_components == n_components) ||
	  (coefficient->n_components == 1),
	  ExcInvalidCoefficient());

  for (FunctionMap<1>::type::const_iterator i=neumann_bc.begin();
       i!=neumann_bc.end(); ++i)
    Assert (i->second->n_components == n_components,
            ExcInvalidBoundaryFunction());

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
  std::vector<std::vector<std::vector<Tensor<1,1> > > >
    gradients_here (n_solution_vectors,
		    std::vector<std::vector<Tensor<1,1> > >(2, std::vector<Tensor<1,1> >(n_components)));
  std::vector<std::vector<std::vector<Tensor<1,1> > > >
    gradients_neighbor (gradients_here);
  std::vector<Vector<double> >
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
  
  QTrapez<1> quadrature;
  FEValues<1> fe_values (mapping, dof_handler.get_fe(), quadrature,
                         update_gradients);
  
				   // loop over all cells and do something on
				   // the cells which we're told to work
				   // on. note that the error indicator is
				   // only a sum over the two contributions
				   // from the two vertices of each cell.
  DoFHandler<1>::active_cell_iterator cell = dof_handler.begin_active();
  for (unsigned int cell_index=0; cell != dof_handler.end();
       ++cell, ++cell_index)
    if ((subdomain_id == deal_II_numbers::invalid_unsigned_int)
        ||
        (cell->subdomain_id() == subdomain_id))
      {
        for (unsigned int n=0; n<n_solution_vectors; ++n)
          (*errors[n])(cell_index) = 0;
      
                                         // loop over the two points bounding
                                         // this line. n==0 is left point,
                                         // n==1 is right point
        for (unsigned int n=0; n<2; ++n)
          {
                                             // find left or right active
                                             // neighbor
            DoFHandler<1>::cell_iterator neighbor = cell->neighbor(n);
            if (neighbor.state() == IteratorState::valid)
              while (neighbor->has_children())
                neighbor = neighbor->child(n==0 ? 1 : 0);
      
                                             // now get the gradients on the
                                             // both sides of the point
            fe_values.reinit (cell);

            for (unsigned int s=0; s<n_solution_vectors; ++s)
              fe_values.get_function_grads (*solutions[s], gradients_here[s]);

            if (neighbor.state() == IteratorState::valid)
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
                    grad_neighbor[s](c)
                      = gradients_neighbor[s][n==0 ? 1 : 0][c][0];
              }
            else
              if (neumann_bc.find(n) != neumann_bc.end())
                                                 // if Neumann b.c., then fill
                                                 // the gradients field which
                                                 // will be used later on.
                {
                  if (n_components==1)
                    {
                      const double
                        v = neumann_bc.find(n)->second->value(cell->vertex(0));
		    
                      for (unsigned int s=0; s<n_solution_vectors; ++s)
                        grad_neighbor[s](0) = v;
                    }
                  else
                    {
                      Vector<double> v(n_components);
                      neumann_bc.find(n)->second->vector_value(cell->vertex(0), v);
		    
                      for (unsigned int s=0; s<n_solution_vectors; ++s)
                        grad_neighbor[s] = v;
                    };
                }
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
              }


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
                  }
          }
      
        for (unsigned int s=0; s<n_solution_vectors; ++s)
          (*errors[s])(cell_index) = std::sqrt((*errors[s])(cell_index));
      }
}


#else // #if deal_II_dimension !=1


template <int dim>
KellyErrorEstimator<dim>::PerThreadData::
PerThreadData (const unsigned int n_solution_vectors,
	       const unsigned int n_components,
	       const unsigned int n_q_points,
               const unsigned int subdomain_id)
                :
                subdomain_id (subdomain_id)
{
				   // Init the size of a lot of vectors
				   // needed in the calculations once
				   // per thread.
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
	}
    }

  for (unsigned int qp=0;qp<n_q_points;++qp)
    coefficient_values[qp].reinit(n_components);
}




// the following function is still independent of dimension, but it
// calls dimension dependent functions
template <int dim>
template <typename InputVector>
void
KellyErrorEstimator<dim>::
estimate (const Mapping<dim>      &mapping,
          const DoFHandler<dim>   &dof_handler,
          const Quadrature<dim-1> &quadrature,
          const typename FunctionMap<dim>::type &neumann_bc,
          const InputVector       &solution,
          Vector<float>           &error,
          const std::vector<bool> &component_mask,
          const Function<dim>     *coefficients,
          const unsigned int       n_threads,
          const unsigned int       subdomain_id)
{
				   // just pass on to the other function
  const std::vector<const InputVector *> solutions (1, &solution);
  std::vector<Vector<float>*>              errors (1, &error);
  estimate (mapping, dof_handler, quadrature, neumann_bc, solutions, errors,
	    component_mask, coefficients, n_threads, subdomain_id);
}


template <int dim>
template <typename InputVector>
void
KellyErrorEstimator<dim>::
estimate (const DoFHandler<dim>   &dof_handler,
          const Quadrature<dim-1> &quadrature,
          const typename FunctionMap<dim>::type &neumann_bc,
          const InputVector       &solution,
          Vector<float>           &error,
          const std::vector<bool> &component_mask,
          const Function<dim>     *coefficients,
          const unsigned int       n_threads,
          const unsigned int       subdomain_id)
{
  Assert (DEAL_II_COMPAT_MAPPING, ExcCompatibility("mapping"));
  static const MappingQ1<dim> mapping;
  estimate(mapping, dof_handler, quadrature, neumann_bc, solution,
	   error, component_mask, coefficients, n_threads, subdomain_id);
}




template <int dim>
template <typename InputVector>
void KellyErrorEstimator<dim>::
estimate_some (const Mapping<dim>                  &mapping,
               const DoFHandler<dim>               &dof_handler,
               const Quadrature<dim-1>             &quadrature,
               const typename FunctionMap<dim>::type &neumann_bc,
               const std::vector<const InputVector *> &solutions,
               const std::vector<bool>                  &component_mask,
               const Function<dim>                 *coefficients,
               const std::pair<unsigned int, unsigned int> this_thread,
               FaceIntegrals                       &face_integrals,
               PerThreadData                       &per_thread_data)
{
  const unsigned int n_solution_vectors = solutions.size();

  const unsigned int subdomain_id = per_thread_data.subdomain_id;
  
				   // make up a fe face values object
				   // for the restriction of the
				   // finite element function to a
				   // face, for the present cell and
				   // its neighbor. In principle we
				   // would only need one at a time,
				   // but this way we can have more
				   // fine grained access to what
				   // values really need to be
				   // computed (we need not compute
				   // all values on the neighbor
				   // cells, so using two objects
				   // gives us a performance gain).
                                   //
                                   // in debug mode, make sure that
                                   // some data matches, so compute
                                   // quadrature points always
  FEFaceValues<dim> fe_face_values_cell (mapping,
					 dof_handler.get_fe(),
					 quadrature,
                                         UpdateFlags(
                                           update_gradients      |
                                           update_JxW_values     |
                                           ((!neumann_bc.empty() ||
                                             (coefficients != 0))  ?
                                            update_q_points : 0) |
                                           update_normal_vectors));
  FEFaceValues<dim> fe_face_values_neighbor (mapping,
					     dof_handler.get_fe(),
					     quadrature,
					     update_gradients);
  FESubfaceValues<dim> fe_subface_values (mapping,
					  dof_handler.get_fe(),
					  quadrature,
                                          update_gradients);


  active_cell_iterator cell = dof_handler.begin_active();

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
  for (unsigned int t=0; (t<this_thread.first) && (cell!=dof_handler.end());
       ++t, ++cell);

  
				   // loop over all cells for this thread
  for (; cell!=dof_handler.end(); advance_by_n(cell,this_thread.second))
    {
      
				       // loop over all faces of this cell
      for (unsigned int face_no=0;
           face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
	{
					   // make sure we do work
					   // only once: this face
					   // may either be regular
					   // or irregular. if it is
					   // regular and has a
					   // neighbor, then we
					   // visit the face twice,
					   // once from every
					   // side. let the one with
					   // the lower index do the
					   // work. if it is at the
					   // boundary, or if the
					   // face is irregular,
					   // then do the work below
	  if ((cell->face(face_no)->has_children() == false) &&
	      !cell->at_boundary(face_no) &&
	      (cell->neighbor(face_no)->level() == cell->level()) &&
	      (cell->neighbor(face_no)->index() < cell->index()))
	    continue;
	  
					   // if we already visited
					   // this face: do
					   // nothing. only check
					   // component for first
					   // solution vector, as we
					   // treat them all at the
					   // same time
	  if (face_integrals[cell->face(face_no)][0] >=0)
	    continue;


					   // if the neighboring cell is less
					   // refined than the present one,
					   // then do nothing since we
					   // integrate over the subfaces when
					   // we visit the coarse cells.
	  if (cell->at_boundary(face_no) == false)
	    if (cell->neighbor(face_no)->level() < cell->level())
	      continue;
	  
					   // if this face is part of the
					   // boundary but not of the neumann
					   // boundary -> nothing to
					   // do. However, to make things
					   // easier when summing up the
					   // contributions of the faces of
					   // cells, we enter this face into
					   // the list of faces with
					   // contribution zero.
	  const unsigned char boundary_indicator
	    = cell->face(face_no)->boundary_indicator();
	  if (cell->face(face_no)->at_boundary()
	      &&
	      neumann_bc.find(boundary_indicator)==neumann_bc.end()) 
	    {
	      for (unsigned int n=0; n<n_solution_vectors; ++n)
		face_integrals[cell->face(face_no)][n] = 0;
	      continue;
	    }

                                           // finally: note that we only
                                           // have to do something if either
                                           // the present cell is on the
                                           // subdomain we care for, or if one
                                           // of the neighbors behind the face
                                           // is on the subdomain we care for
          if ( ! ((subdomain_id == deal_II_numbers::invalid_unsigned_int)
                  ||
                  (cell->subdomain_id() == subdomain_id)))
            {
                                               // ok, cell is unwanted, but
                                               // maybe its neighbor behind
                                               // the face we presently work
                                               // on? oh is there a face at
                                               // all?
              if (cell->at_boundary(face_no))
                continue;

              bool care_for_cell = false;
              if (cell->neighbor(face_no)->has_children() == false)
                care_for_cell |= (cell->neighbor(face_no)->subdomain_id()
                                  == subdomain_id);
              else
                {
                  for (unsigned int sf=0;
                       sf<GeometryInfo<dim>::subfaces_per_face; ++sf)
                    if (cell->neighbor_child_on_subface(face_no,sf)
                        ->subdomain_id() == subdomain_id)
                      {
                        care_for_cell = true;
                        break;
                      }
                }

                                               // so if none of the neighbors
                                               // cares for this subdomain
                                               // either, then try next face
              if (care_for_cell == false)
                continue;
            }
              

                                           // so now we know that we care for
                                           // this face, let's do something
                                           // about it:
	  if (cell->face(face_no)->has_children() == false)
					     // if the face is a regular one,
					     // i.e.  either on the other side
					     // there is nirvana (face is at
					     // boundary), or the other side's
					     // refinement level is the same
					     // as that of this side, then
					     // handle the integration of
					     // these both cases together
	    integrate_over_regular_face (dof_handler, quadrature,
                                         neumann_bc, solutions, component_mask,
                                         coefficients,
                                         face_integrals,
                                         per_thread_data,
					 cell, face_no,
					 fe_face_values_cell,
					 fe_face_values_neighbor);
	  
	  else
					     // otherwise we need to do some
					     // special computations which do
					     // not fit into the framework of
					     // the above function
	    integrate_over_irregular_face (dof_handler, quadrature,
                                           solutions, component_mask,
                                           coefficients,
                                           face_integrals,
                                           per_thread_data,
					   cell, face_no,
					   fe_face_values_cell,
					   fe_subface_values);
	}
    }
}



template <int dim>
template <typename InputVector>
void
KellyErrorEstimator<dim>::
estimate (const Mapping<dim>                  &mapping,
          const DoFHandler<dim>               &dof_handler,
          const Quadrature<dim-1>             &quadrature,
          const typename FunctionMap<dim>::type &neumann_bc,
          const std::vector<const InputVector *> &solutions,
          std::vector<Vector<float>*>              &errors,
          const std::vector<bool>                  &component_mask_,
          const Function<dim>                 *coefficients,
          const unsigned int                   n_threads_,
          const unsigned int                   subdomain_id)
{
  const unsigned int n_components = dof_handler.get_fe().n_components();

  				   // sanity checks
  Assert (solutions.size() > 0,
	  ExcNoSolutions());
  Assert (solutions.size() == errors.size(),
	  ExcIncompatibleNumberOfElements(solutions.size(), errors.size()));

  for (typename FunctionMap<dim>::type::const_iterator i=neumann_bc.begin();
       i!=neumann_bc.end(); ++i)
    Assert (i->second->n_components == n_components,
            ExcInvalidBoundaryFunction());  
  
  Assert ((component_mask_.size() == 0) ||
          (component_mask_.size() == n_components), ExcInvalidComponentMask());
  Assert ((component_mask_.size() == 0) ||
          (std::count(component_mask_.begin(), component_mask_.end(),
                      true) > 0),
	  ExcInvalidComponentMask());

  Assert ((coefficients == 0) ||
	  (coefficients->n_components == n_components) ||
	  (coefficients->n_components == 1),
	  ExcInvalidCoefficient());  

  for (unsigned int n=0; n<solutions.size(); ++n)
    Assert (solutions[n]->size() == dof_handler.n_dofs(),
	    ExcInvalidSolutionVector());

  Assert (n_threads_ > 0, ExcInternalError());
  
				   // if no mask given: treat all components
  std::vector<bool> component_mask ((component_mask_.size() == 0)    ?
				    std::vector<bool>(n_components, true) :
				    component_mask_);
  Assert (component_mask.size() == n_components, ExcInvalidComponentMask());
  Assert (std::count(component_mask.begin(), component_mask.end(), true) > 0,
	  ExcInvalidComponentMask());
	  
				   // if NOT multithreaded, set
				   // n_threads to one
  const unsigned int n_threads = (DEAL_II_USE_MT ? n_threads_ : 1);
  
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
  const double invalid_double = -10e20;
  std::vector<double> default_face_integrals (n_solution_vectors,
                                              invalid_double);
  FaceIntegrals face_integrals;
  for (active_cell_iterator cell=dof_handler.begin_active();
       cell!=dof_handler.end(); ++cell)
    for (unsigned int face_no=0;
         face_no<GeometryInfo<dim>::faces_per_cell;
         ++face_no)
      face_integrals[cell->face(face_no)] = default_face_integrals;


				   // all the data needed in the error
				   // estimator by each of the threads
				   // is gathered in the following
				   // stuctures
				   //
				   // note that if no component mask
				   // was given, then treat all
				   // components
  std::vector<PerThreadData*> data_structures (n_threads);
  for (unsigned int i=0; i<n_threads; ++i)
    data_structures[i]
      = new PerThreadData (solutions.size(),
                           dof_handler.get_fe().n_components(),
                           quadrature.n_quadrature_points,
                           subdomain_id);
  
				   // split all cells into threads if
				   // multithreading is used and run
				   // the whole thing. use the
				   // function pointer variable to
				   // work around another nasty bug in
				   // icc7
  Threads::ThreadGroup<> threads;
  void (*estimate_some_ptr) (const Mapping<dim>                  &,
			     const DoFHandler<dim>               &,
			     const Quadrature<dim-1>             &,
			     const typename FunctionMap<dim>::type &,
			     const std::vector<const InputVector *> &,
			     const std::vector<bool>               &,
			     const Function<dim>                 *,
			     const std::pair<unsigned int, unsigned int>,
			     FaceIntegrals                       &,
			     PerThreadData                       &)
    = &KellyErrorEstimator<dim>::template estimate_some<InputVector>;
  
  for (unsigned int i=0; i<n_threads; ++i)
    threads += Threads::spawn (estimate_some_ptr)
               (mapping, dof_handler,
                quadrature, neumann_bc, solutions,
                component_mask, coefficients,
                std::make_pair(i, n_threads),
                face_integrals,
                *data_structures[i]);
  threads.join_all();

				   // delete the structures for the
				   // different threads again. the
				   // results are in the
				   // face_integrals object
  for (unsigned int i=0; i<n_threads; ++i)
    {
      delete data_structures[i];
      data_structures[i] = 0;
    }
  
  
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

                                   // now walk over all cells and collect
                                   // information from the faces. only do
                                   // something if this is a cell we care for
                                   // based on the subdomain id
  for (active_cell_iterator cell=dof_handler.begin_active();
       cell!=dof_handler.end();
       ++cell, ++present_cell)
    if ((subdomain_id == deal_II_numbers::invalid_unsigned_int)
        ||
        (cell->subdomain_id() == subdomain_id))
      {
                                         // loop over all faces of this cell
        for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell;
             ++face_no)
          {
            Assert(face_integrals.find(cell->face(face_no))
                   != face_integrals.end(),
                   ExcInternalError());
            const double factor = cell->diameter() / 24;
	  
            for (unsigned int n=0; n<n_solution_vectors; ++n)
              {
                                                 // make sure that we have
                                                 // written a meaningful value
                                                 // into this slot
                Assert (face_integrals[cell->face(face_no)][n] >= 0,
                        ExcInternalError());
                
                (*errors[n])(present_cell)
                  += (face_integrals[cell->face(face_no)][n] * factor);
              }
          }

        for (unsigned int n=0; n<n_solution_vectors; ++n)
          (*errors[n])(present_cell) = std::sqrt((*errors[n])(present_cell));
      }
}


template <int dim>
template <typename InputVector>
void KellyErrorEstimator<dim>::estimate (const DoFHandler<dim>               &dof_handler,
					 const Quadrature<dim-1>             &quadrature,
					 const typename FunctionMap<dim>::type &neumann_bc,
					 const std::vector<const InputVector *> &solutions,
					 std::vector<Vector<float>*>              &errors,
					 const std::vector<bool>                  &component_mask,
					 const Function<dim>                 *coefficients,
					 const unsigned int                   n_threads,
                                         const unsigned int       subdomain_id)
{
  Assert (DEAL_II_COMPAT_MAPPING, ExcCompatibility("mapping"));  
  static const MappingQ1<dim> mapping;
  estimate(mapping, dof_handler, quadrature, neumann_bc, solutions,
	   errors, component_mask, coefficients, n_threads, subdomain_id);
}



template <int dim>
template <typename InputVector>
void KellyErrorEstimator<dim>::
integrate_over_regular_face (const DoFHandler<dim>               &dof_handler,
                             const Quadrature<dim-1>             &quadrature,
                             const typename FunctionMap<dim>::type &neumann_bc,
                             const std::vector<const InputVector *> &solutions,
                             const std::vector<bool>                  &component_mask,
                             const Function<dim>                 *coefficients,
                             FaceIntegrals                       &face_integrals,
                             PerThreadData              &per_thread_data,
			     const active_cell_iterator &cell,
			     const unsigned int          face_no,
			     FEFaceValues<dim>          &fe_face_values_cell,
			     FEFaceValues<dim>          &fe_face_values_neighbor)
{
  const typename DoFHandler<dim>::face_iterator face = cell->face(face_no);
  const unsigned int n_q_points         = quadrature.n_quadrature_points,
		     n_components       = dof_handler.get_fe().n_components(),
		     n_solution_vectors = solutions.size();
  
  
				   // initialize data of the restriction
				   // of this cell to the present face
  fe_face_values_cell.reinit (cell, face_no);
  
				   // get gradients of the finite element
				   // function on this cell
  for (unsigned int n=0; n<n_solution_vectors; ++n)
    fe_face_values_cell.get_function_grads (*solutions[n], per_thread_data.psi[n]);
  
				   // now compute over the other side of
				   // the face
  if (face->at_boundary() == false)
				     // internal face; integrate jump
				     // of gradient across this face
    {
      Assert (cell->neighbor(face_no).state() == IteratorState::valid,
	      ExcInternalError());      
      
      const active_cell_iterator neighbor = cell->neighbor(face_no);
      
				       // find which number the
				       // current face has relative to
				       // the neighboring cell
      const unsigned int neighbor_neighbor
        = cell->neighbor_of_neighbor (face_no);
      Assert (neighbor_neighbor<GeometryInfo<dim>::faces_per_cell,
              ExcInternalError());
      
				       // get restriction of finite element
				       // function of @p{neighbor} to the
				       // common face.
      fe_face_values_neighbor.reinit (neighbor, neighbor_neighbor);
      
				       // get gradients on neighbor cell
      for (unsigned int n=0; n<n_solution_vectors; ++n)
	{
	  fe_face_values_neighbor.get_function_grads (*solutions[n],
						      per_thread_data.neighbor_psi[n]);
      
					   // compute the jump in the gradients
	  for (unsigned int component=0; component<n_components; ++component)
	    for (unsigned int p=0; p<n_q_points; ++p)
	      per_thread_data.psi[n][p][component] -= per_thread_data.neighbor_psi[n][p][component];
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
  
  per_thread_data.normal_vectors=fe_face_values_cell.get_normal_vectors();
  
  for (unsigned int n=0; n<n_solution_vectors; ++n)
    for (unsigned int component=0; component<n_components; ++component)
      for (unsigned int point=0; point<n_q_points; ++point)
	per_thread_data.phi[n][point][component] = per_thread_data.psi[n][point][component]*
					per_thread_data.normal_vectors[point];
  
				   // if a coefficient was given: use that
				   // to scale the jump in the gradient
  if (coefficients != 0)
    {
				       // scalar coefficient
      if (coefficients->n_components == 1)
	{
	  
	  coefficients->value_list (fe_face_values_cell.get_quadrature_points(),
                                    per_thread_data.coefficient_values1);
	  for (unsigned int n=0; n<n_solution_vectors; ++n)
	    for (unsigned int component=0; component<n_components; ++component)
	      for (unsigned int point=0; point<n_q_points; ++point)
		per_thread_data.phi[n][point][component] *=
		  per_thread_data.coefficient_values1[point];
	}
      else
					   // vector-valued coefficient
	{
	  coefficients->vector_value_list (fe_face_values_cell.get_quadrature_points(),
                                           per_thread_data.coefficient_values);
	  for (unsigned int n=0; n<n_solution_vectors; ++n)
	    for (unsigned int component=0; component<n_components; ++component)
	      for (unsigned int point=0; point<n_q_points; ++point)
		per_thread_data.phi[n][point][component] *=
		  per_thread_data.coefficient_values[point](component);
	};
    };


  if (face->at_boundary() == true)
				     // neumann boundary face. compute
				     // difference between normal
				     // derivative and boundary function
    {
      const unsigned char boundary_indicator = face->boundary_indicator();
      
      Assert (neumann_bc.find(boundary_indicator) != neumann_bc.end(),
	      ExcInternalError ());
				       // get the values of the boundary
				       // function at the quadrature
				       // points
      if (n_components == 1)
	{
	  std::vector<double> g(n_q_points);
	  neumann_bc.find(boundary_indicator)->second
	    ->value_list (fe_face_values_cell.get_quadrature_points(), g);
	  
	  for (unsigned int n=0; n<n_solution_vectors; ++n)
	    for (unsigned int point=0; point<n_q_points; ++point)
	      per_thread_data.phi[n][point][0] -= g[point];
	}
      else
	{
	  std::vector<Vector<double> > g(n_q_points, Vector<double>(n_components));
	  neumann_bc.find(boundary_indicator)->second
	    ->vector_value_list (fe_face_values_cell.get_quadrature_points(),
				 g);
	  
	  for (unsigned int n=0; n<n_solution_vectors; ++n)
	    for (unsigned int component=0; component<n_components; ++component)
	      for (unsigned int point=0; point<n_q_points; ++point)
		per_thread_data.phi[n][point][component] -= g[point](component);
	};
    };


				   // now phi contains the following:
				   // - for an internal face, phi=[a du/dn]
				   // - for a neumann boundary face,
				   //   phi=a du/dn-g
				   // each component being the
				   // mentioned value at one of the
				   // quadrature points

  per_thread_data.JxW_values = fe_face_values_cell.get_JxW_values();
  
				   // take the square of the phi[i]
				   // for integration, and sum up
  std::vector<double> face_integral (n_solution_vectors, 0);
  for (unsigned int n=0; n<n_solution_vectors; ++n)
    for (unsigned int component=0; component<n_components; ++component)
      if (component_mask[component] == true)
	for (unsigned int p=0; p<n_q_points; ++p)
	  face_integral[n] += ::sqr(per_thread_data.phi[n][p][component]) *
			      per_thread_data.JxW_values[p];

				   // double check that the element
				   // already exists and that it was
				   // not already written to
  Assert (face_integrals.find (face) != face_integrals.end(),
	  ExcInternalError());
  Assert (face_integrals[face][0] < 0, ExcInternalError());
  
  face_integrals[face] = face_integral;
}



template <int dim>
template <typename InputVector>
void KellyErrorEstimator<dim>::
integrate_over_irregular_face (const DoFHandler<dim>               &dof_handler,
                               const Quadrature<dim-1>             &quadrature,
                               const std::vector<const InputVector *> &solutions,
                               const std::vector<bool>                  &component_mask,
                               const Function<dim>                 *coefficients,
                               FaceIntegrals                       &face_integrals,
                               PerThreadData              &per_thread_data,
			       const active_cell_iterator &cell,
			       const unsigned int          face_no,
			       FEFaceValues<dim>          &fe_face_values,
			       FESubfaceValues<dim>       &fe_subface_values)
{
  const typename DoFHandler<dim>::cell_iterator neighbor = cell->neighbor(face_no);
  const unsigned int n_q_points         = quadrature.n_quadrature_points,
		     n_components       = dof_handler.get_fe().n_components(),
		     n_solution_vectors = solutions.size();

  Assert (neighbor.state() == IteratorState::valid, ExcInternalError());
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
  
				   // store which number @p{cell} has
				   // in the list of neighbors of
				   // @p{neighbor}
  const unsigned int neighbor_neighbor
    = cell->neighbor_of_neighbor (face_no);
  Assert (neighbor_neighbor<GeometryInfo<dim>::faces_per_cell,
          ExcInternalError());
  
				   // loop over all subfaces
  for (unsigned int subface_no=0;
       subface_no<GeometryInfo<dim>::subfaces_per_face;
       ++subface_no)
    {
				       // get an iterator pointing to the
				       // cell behind the present subface
      const active_cell_iterator neighbor_child
        = cell->neighbor_child_on_subface (face_no, subface_no);
      Assert (!neighbor_child->has_children(),
	      ExcInternalError());
      
				       // restrict the finite element
				       // on the present cell to the
				       // subface
      fe_subface_values.reinit (cell, face_no, subface_no);

				       // restrict the finite element
				       // on the neighbor cell to the
				       // common @p{subface}.
      fe_face_values.reinit (neighbor_child, neighbor_neighbor);

                                       // store the gradient of the
				       // solution in psi
      for (unsigned int n=0; n<n_solution_vectors; ++n)
	fe_subface_values.get_function_grads (*solutions[n], per_thread_data.psi[n]);

                                       // store the gradient from the
				       // neighbor's side in
				       // @p{neighbor_psi}
      for (unsigned int n=0; n<n_solution_vectors; ++n)
	fe_face_values.get_function_grads (*solutions[n], per_thread_data.neighbor_psi[n]);
      
				       // compute the jump in the gradients
      for (unsigned int n=0; n<n_solution_vectors; ++n)
	for (unsigned int component=0; component<n_components; ++component)
	  for (unsigned int p=0; p<n_q_points; ++p)
	    per_thread_data.psi[n][p][component] -=
	      per_thread_data.neighbor_psi[n][p][component];

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
     
      per_thread_data.normal_vectors=fe_face_values.get_normal_vectors();


      for (unsigned int n=0; n<n_solution_vectors; ++n)
	for (unsigned int component=0; component<n_components; ++component)
	  for (unsigned int point=0; point<n_q_points; ++point)
	    per_thread_data.phi[n][point][component] = (per_thread_data.psi[n][point][component]*
					     per_thread_data.normal_vectors[point]);
      
				       // if a coefficient was given: use that
				       // to scale the jump in the gradient
      if (coefficients != 0)
	{
					   // scalar coefficient
	  if (coefficients->n_components == 1)
	    {
	      coefficients->value_list (fe_face_values.get_quadrature_points(),
                                        per_thread_data.coefficient_values1);
	      for (unsigned int n=0; n<n_solution_vectors; ++n)
		for (unsigned int component=0; component<n_components; ++component)
		  for (unsigned int point=0; point<n_q_points; ++point)
		    per_thread_data.phi[n][point][component] *=
		      per_thread_data.coefficient_values1[point];
	    }
	  else
					     // vector-valued coefficient
	    {
	      coefficients->vector_value_list (fe_face_values.get_quadrature_points(),
						    per_thread_data.coefficient_values);
	      for (unsigned int n=0; n<n_solution_vectors; ++n)
		for (unsigned int component=0; component<n_components; ++component)
		  for (unsigned int point=0; point<n_q_points; ++point)
		    per_thread_data.phi[n][point][component] *=
		      per_thread_data.coefficient_values[point](component);
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
      per_thread_data.JxW_values = fe_face_values.get_JxW_values();
      
				       // take the square of the phi[i]
				       // for integration, and sum up
      std::vector<double> face_integral (n_solution_vectors, 0);
      for (unsigned int n=0; n<n_solution_vectors; ++n)
	for (unsigned int component=0; component<n_components; ++component)
	  if (component_mask[component] == true)
	    for (unsigned int p=0; p<n_q_points; ++p)
	      face_integral[n] += ::sqr(per_thread_data.phi[n][p][component]) *
				  per_thread_data.JxW_values[p];

      face_integrals[neighbor_child->face(neighbor_neighbor)] = face_integral;
    };


				   // finally loop over all subfaces to
				   // collect the contributions of the
				   // subfaces and store them with the
				   // mother face
  std::vector<double> sum (n_solution_vectors, 0);
  typename DoFHandler<dim>::face_iterator face = cell->face(face_no);
  for (unsigned int subface_no=0; subface_no<GeometryInfo<dim>::subfaces_per_face;
       ++subface_no) 
    {
      Assert (face_integrals.find(face->child(subface_no)) !=
	      face_integrals.end(),
	      ExcInternalError());
      Assert (face_integrals[face->child(subface_no)][0] >= 0,
 	      ExcInternalError());
      
      for (unsigned int n=0; n<n_solution_vectors; ++n)
	sum[n] += face_integrals[face->child(subface_no)][n];
    };

  face_integrals[face] = sum;
}

#endif



// explicit instantiations
#if deal_II_dimension != 1
template class KellyErrorEstimator<deal_II_dimension>;
#endif

// instantiate the externally visible functions. define a list of functions
// for vectors, where the vector/matrix can be replaced using a preprocessor
// variable VectorType/MatrixType. note that we need a space between
// "VectorType" and ">" to disambiguate ">>" when VectorType trails in an
// angle bracket
#define INSTANTIATE(InputVector) \
template    \
void    \
KellyErrorEstimator<deal_II_dimension>::    \
estimate<InputVector > (const Mapping<deal_II_dimension>      &,    \
          const DoFHandler<deal_II_dimension>   &,    \
          const Quadrature<deal_II_dimension-1> &,    \
          const FunctionMap<deal_II_dimension>::type &,    \
          const InputVector       &,    \
          Vector<float>           &,    \
          const std::vector<bool> &,    \
          const Function<deal_II_dimension>     *,    \
          const unsigned int       , \
          const unsigned int);    \
    \
template    \
void    \
KellyErrorEstimator<deal_II_dimension>::    \
estimate<InputVector > (const DoFHandler<deal_II_dimension>   &,    \
          const Quadrature<deal_II_dimension-1> &,    \
          const FunctionMap<deal_II_dimension>::type &,    \
          const InputVector       &,    \
          Vector<float>           &,    \
          const std::vector<bool> &,    \
          const Function<deal_II_dimension>     *,    \
          const unsigned int       , \
          const unsigned int);    \
        \
template    \
void    \
KellyErrorEstimator<deal_II_dimension>::    \
estimate<InputVector > (const Mapping<deal_II_dimension>          &,    \
          const DoFHandler<deal_II_dimension>       &,    \
          const Quadrature<deal_II_dimension-1>     &,    \
          const FunctionMap<deal_II_dimension>::type &,    \
          const std::vector<const InputVector *> &,    \
          std::vector<Vector<float>*> &,    \
          const std::vector<bool>     &,    \
          const Function<deal_II_dimension>         *,    \
          const unsigned int           , \
          const unsigned int);    \
    \
template    \
void    \
KellyErrorEstimator<deal_II_dimension>::    \
estimate<InputVector > (const DoFHandler<deal_II_dimension>       &,    \
          const Quadrature<deal_II_dimension-1>     &,    \
          const FunctionMap<deal_II_dimension>::type &,    \
          const std::vector<const InputVector *> &,    \
          std::vector<Vector<float>*> &,    \
          const std::vector<bool>     &,    \
          const Function<deal_II_dimension>         *,    \
          const unsigned int           , \
          const unsigned int)

INSTANTIATE(Vector<double>);
INSTANTIATE(Vector<float>);
INSTANTIATE(BlockVector<double>);
INSTANTIATE(BlockVector<float>);

#ifdef DEAL_II_USE_PETSC
INSTANTIATE(PETScWrappers::Vector);
#endif
