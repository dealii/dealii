//----------------------------  gradient_estimator.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  gradient_estimator.cc  ---------------------------


#include <base/quadrature_lib.h>
#include <lac/vector.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_handler.h>
#include <fe/fe.h>
#include <fe/fe_values.h>
#include <numerics/gradient_estimator.h>

#ifdef DEAL_II_USE_MT
#  include <base/thread_management.h>
#  include <base/multithread_info.h>
#endif



template <int dim>
void 
GradientEstimator::estimate (const DoFHandler<dim> &dof_handler,
			     const Vector<double>  &solution,
			     Vector<float>         &error_per_cell)
{
  Assert (error_per_cell.size() == dof_handler.get_tria().n_active_cells(),
	  ExcInvalidVectorLength (error_per_cell.size(),
				  dof_handler.get_tria().n_active_cells()));
  Assert (dof_handler.get_fe().n_components() == 1,
	  ExcInternalError());

#ifdef DEAL_II_USE_MT
  const unsigned int n_threads = multithread_info.n_default_threads;
  vector<IndexInterval> index_intervals
    = Threads::split_interval (0, dof_handler.get_tria().n_active_cells(),
			       n_threads);
  Threads::ThreadManager thread_manager;
  for (unsigned int i=0; i<n_threads; ++i)
    Threads::spawn (thread_manager,
		    Threads::encapsulate (&GradientEstimator::
					  template estimate_threaded<dim>)
		    .collect_args (dof_handler, solution, index_intervals[i],
				   error_per_cell));
  thread_manager.wait ();
  
#else
  estimate_threaded (dof_handler, solution,
		     make_pair(0U, dof_handler.get_tria().n_active_cells()),
		     error_per_cell);
#endif
};



template <int dim>
void 
GradientEstimator::estimate_threaded (const DoFHandler<dim> &dof_handler,
				      const Vector<double>  &solution,
				      const IndexInterval   &index_interval,
				      Vector<float>         &error_per_cell)
{
  QMidpoint<dim> midpoint_rule;
  FEValues<dim>  fe_midpoint_value (dof_handler.get_fe(),
				    midpoint_rule,
				    UpdateFlags(update_values |
						update_q_points));
  
				   // matrix Y=sum_i y_i y_i^T
  Tensor<2,dim> Y;
  
				   // iterators over all cells and the
				   // respective entries in the output
				   // vector:
  Vector<float>::iterator
    error_on_this_cell = error_per_cell.begin() + index_interval.first;
  
  DoFHandler<dim>::active_cell_iterator cell, endc;
  cell = endc = dof_handler.begin_active();
				   // (static_cast to avoid warnings
				   // about unsigned always >=0)
  advance (cell, static_cast<int>(index_interval.first));
  advance (endc, static_cast<int>(index_interval.second));

				   // vector to hold iterators to all
				   // active neighbors of a cell
				   // reserve the maximal number of
				   // active neighbors
  vector<DoFHandler<dim>::active_cell_iterator> active_neighbors;
  active_neighbors.reserve (GeometryInfo<dim>::faces_per_cell *
			    GeometryInfo<dim>::subfaces_per_face);

  for (; cell!=endc; ++cell, ++error_on_this_cell)
    {
      Y.clear ();
				       // vector g=sum_i y_i (f(x+y_i)-f(x))/|y_i|
      Tensor<1,dim> projected_gradient;

				       // reinit fe values object...
      fe_midpoint_value.reinit (cell);

				       // ...and get the value of the
				       // solution...
      vector<double> this_midpoint_value(1);
      fe_midpoint_value.get_function_values (solution, this_midpoint_value);
				       // ...and the place where it lives
      Point<dim> this_center = fe_midpoint_value.quadrature_point(0);

      
				       // loop over all neighbors and
				       // accumulate the difference
				       // quotients from them. note
				       // that things get a bit more
				       // complicated if the neighbor
				       // is more refined than the
				       // present one
				       //
				       // to make processing simpler,
				       // first collect all neighbor
				       // cells in a vector, and then
				       // collect the data from them
      active_neighbors.clear ();
      for (unsigned int n=0; n<GeometryInfo<dim>::faces_per_cell; ++n)
	if (! cell->at_boundary(n))
	  {
	    DoFHandler<dim>::cell_iterator neighbor = cell->neighbor(n);
	    if (neighbor->active())
	      active_neighbors.push_back (neighbor);
	    else
	      {
						 // check children
						 // of
						 // neighbor. note
						 // that in 1d
						 // children of
						 // the neighbor
						 // may be further
						 // refined, while
						 // they can't in
						 // more than one
						 // dimension. however,
						 // in 1d the case
						 // is simpler
						 // since we know
						 // what children
						 // bound to the
						 // present cell
		if (dim == 1)
		  {
		    DoFHandler<dim>::cell_iterator neighbor_child = neighbor;
		    while (neighbor_child->has_children())
		      neighbor_child = neighbor_child->child (n==0 ? 1 : 0);
		    
		    Assert (neighbor_child->neighbor(n==0 ? 1 : 0)==cell,
			    ExcInternalError());
		    
		    active_neighbors.push_back (neighbor_child);
		  }
		else
						   // this neighbor has
						   // children. find out
						   // which border to the
						   // present cell
		  for (unsigned int c=0; c<GeometryInfo<dim>::children_per_cell; ++c)
		    for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
		      if (neighbor->child(c)->neighbor(f) == cell)
			active_neighbors.push_back (neighbor->child(c));
	      };
	  };

				       // now loop over all active
				       // neighbors and collect the
				       // data we need
      typename vector<DoFHandler<dim>::active_cell_iterator>::const_iterator
	neighbor_ptr = active_neighbors.begin();
      for (; neighbor_ptr!=active_neighbors.end(); ++neighbor_ptr)
	{
	  const DoFHandler<dim>::active_cell_iterator
	    neighbor = *neighbor_ptr;
	    
					   // reinit fe values object...
	  fe_midpoint_value.reinit (neighbor);
	  
					   // ...and get the value of the
					   // solution...
	  vector<double> neighbor_midpoint_value(1);
	  fe_midpoint_value.get_function_values (solution, this_midpoint_value);
					   // ...and the place where it lives
	  Point<dim> neighbor_center = fe_midpoint_value.quadrature_point(0);
	  
	  
					   // vector for the
					   // normalized
					   // direction between
					   // the centers of two
					   // cells
	  Point<dim>   y        = neighbor_center - this_center;
	  const double distance = sqrt(y.square());
					   // normalize y
	  y /= distance;
	  
					   // add up the
					   // contribution of
					   // this cell to Y
	  for (unsigned int i=0; i<dim; ++i)
	    for (unsigned int j=0; j<dim; ++j)
	      Y[i][j] += y[i] * y[j];
	  
					   // the update the sum
					   // of difference
					   // quotients
	  projected_gradient += (neighbor_midpoint_value[0] -
				 this_midpoint_value[0]) /
				distance *
				y;
	};

				       // can we determine an
				       // approximation of the
				       // gradient for the present
				       // cell? if so, then we need to
				       // have passed over vectors y_i
				       // which span the whole space,
				       // otherwise we would not have
				       // all components of the
				       // gradient
      AssertThrow (determinant(Y) != 0,
		   ExcInsufficientDirections());

                                       // compute Y^-1 g
      Point<dim> gradient;
      Tensor<2,dim> Y_inverse = invert(Y);
      
                                       // compute Y^-1 g
      contract (gradient, Y_inverse, projected_gradient);

      *error_on_this_cell = sqrt(gradient.square());
    };
};




// explicit instantiations
template
void 
GradientEstimator::estimate (const DoFHandler<deal_II_dimension> &dof_handler,
			     const Vector<double>  &solution,
			     Vector<float>         &error_per_cell);



