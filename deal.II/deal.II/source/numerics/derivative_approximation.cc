//----------------------------  gradient_estimator.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  derivative_approximation.cc  ---------------------------


#include <base/quadrature_lib.h>
#include <base/thread_management.h>
#include <base/multithread_info.h>
#include <lac/vector.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_handler.h>
#include <fe/fe.h>
#include <fe/fe_values.h>
#include <fe/mapping_q1.h>
#include <numerics/derivative_approximation.h>

#include <cmath>



template <typename T>
static T sqr (const T t)
{
  return t*t;
};





template <int dim>
inline
typename DerivativeApproximation::Gradient<dim>::ProjectedDerivative
DerivativeApproximation::Gradient<dim>::
get_projected_derivative (const FEValues<dim>  &fe_values,
			  const Vector<double> &solution,
			  const unsigned int    component)
{
  if (fe_values.get_fe().n_components() == 1)
    {
      std::vector<ProjectedDerivative> values (1);
      fe_values.get_function_values (solution, values);
      return values[0];
    }
  else
    {
      std::vector<Vector<double> > values
	(1, Vector<double>(fe_values.get_fe().n_components()));
      fe_values.get_function_values (solution, values);
      return values[0](component);
    };
};



template <int dim>
inline
double
DerivativeApproximation::Gradient<dim>::derivative_norm (const Derivative &d)
{
  double s = 0;
  for (unsigned int i=0; i<dim; ++i)
    s += d[i]*d[i];
  return sqrt(s);
};



template <int dim>
inline
void
DerivativeApproximation::Gradient<dim>::symmetrize (Derivative &)
{
				   // nothing to do here
};



template <int dim>
inline
typename DerivativeApproximation::SecondDerivative<dim>::ProjectedDerivative
DerivativeApproximation::SecondDerivative<dim>::
get_projected_derivative (const FEValues<dim>  &fe_values,
			  const Vector<double> &solution,
			  const unsigned int    component)
{
  if (fe_values.get_fe().n_components() == 1)
    {
      std::vector<ProjectedDerivative> values (1);
      fe_values.get_function_grads (solution, values);
      return values[0];
    }
  else
    {
      std::vector<std::vector<ProjectedDerivative> > values
	(1, std::vector<ProjectedDerivative>(fe_values.get_fe().n_components()));
      fe_values.get_function_grads (solution, values);
      return values[0][component];
    };
};


#if deal_II_dimension == 1

template <>
inline
double
DerivativeApproximation::SecondDerivative<1>::
derivative_norm (const Derivative &d)
{
  return fabs (d[0][0]);
};

#endif

#if deal_II_dimension == 2

template <>
inline
double
DerivativeApproximation::SecondDerivative<2>::
derivative_norm (const Derivative &d)
{
				   // note that d should be a
				   // symmetric 2x2 tensor, so the
				   // eigenvalues are:
				   //
				   // 1/2(a+b\pm\sqrt((a-b)^2+4c^2))
				   //
				   // if the d_11=a, d_22=b,
				   // d_12=d_21=c
  const double radicand = ::sqr(d[0][0] - d[1][1]) + 4*::sqr(d[0][1]);
  const double eigenvalues[2]
    = { 0.5*(d[0][0] + d[1][1] + sqrt(radicand)),
	0.5*(d[0][0] + d[1][1] - sqrt(radicand))  };
  
  return std::max (std::fabs (eigenvalues[0]),
		   std::fabs (eigenvalues[1]));
};

#endif


#if deal_II_dimension == 3

template <>
inline
double
DerivativeApproximation::SecondDerivative<3>::
derivative_norm (const Derivative &d)
{
				   // compute the three eigenvalues of
				   // the tensor @p{d} and take the
				   // largest:
  const double t1 = d[1][2]*d[1][2];
  const double t2 = d[0][0]*t1;
  const double t3 = d[0][1]*d[0][1];
  const double t4 = t3*d[2][2];
  const double t5 = d[0][2]*d[0][2];
  const double t6 = t5*d[1][1];
  const double t7 = d[0][0]*d[1][1];
  const double t8 = t7*d[2][2];
  const double t9 = d[0][1]*d[0][2];
  const double t10 = t9*d[1][2];
  const double t11 = t3*d[0][0];
  const double t12 = t3*d[1][1];
  const double t13 = d[0][0]*d[0][0];
  const double t14 = t13*d[1][1];
  const double t15 = d[1][1]*d[1][1];
  const double t16 = d[0][0]*t15;
  const double t17 = t13*d[2][2];
  const double t19 = d[2][2]*d[2][2];
  const double t24 = t5*d[2][2];
  const double t25 = t1*d[2][2];
  const double t27 = t13*d[0][0];
  const double t28 = t19*d[2][2];
  const double t29 = t15*d[1][1];
  const double t30 = t1*t1;
  const double t32 = d[1][1]*t28;
  const double t35 = t15*t19;
  const double t37 = t29*d[2][2];
  const double t39 = t13*t15;
  const double t41 = d[0][0]*d[2][2];
  const double t43 = t19*t19;
  const double t45 = t3*t1;
  const double t47 = t3*t13;
  const double t49 = t27*d[1][1];
  const double t51 = t27*d[2][2];
  const double t53 = t15*t15;
  const double t57 = -3.0*t30*t15-6.0*t32*t5+6.0*t11*t29-6.0*t35*t5+24.0*t37*t5-6.0*t39*
t5-24.0*t41*t30+6.0*t7*t43-6.0*t45*t15-24.0*t47*t15+24.0*t49*t1-6.0*t51*t15
-12.0*t5*t53-36.0*t5*t30-3.0*t30*t19;
  const double t59 = t13*t13;
  const double t61 = t5*t5;
  const double t64 = t3*t5;
  const double t71 = d[0][0]*t29;
  const double t78 = 24.0*t30*t13-12.0*t1*t59-3.0*t61*t19+6.0*t27*t28-6.0*t64*t13-60.0*
t5*t13*t1-60.0*t64*t19-60.0*t45*t19-6.0*t39*t1+24.0*t71*t5-6.0*t49*t19-3.0*t13*
t43+6.0*t29*t28-36.0*t61*t1-3.0*t61*t13;
  const double t81 = t3*t3;
  const double t82 = t81*d[0][0];
  const double t95 = d[0][0]*t28;
  const double t99 = 18.0*t39*t19-24.0*t82*d[2][2]-3.0*t59*t15+6.0*t27*t29+6.0*t59*d[1][1]*d[2][2]
-6.0*t14*t28-6.0*t71*t19-6.0*t16*t28-3.0*t13*t53+24.0*t81*t19-36.0*t81*t5-3.0*
t81*t13+6.0*t95*t5-6.0*t95*t1+30.0*t41*t61;
  const double t105 = t13*t19;
  const double t111 = d[1][1]*d[2][2];
  const double t114 = t5*t1;
  const double t119 = -3.0*t81*t15-12.0*t3*t43-6.0*t71*t1-24.0*t7*t61-24.0*t7*t30-6.0*
t105*t1+24.0*t51*t1-24.0*t105*t5-6.0*t17*t29+6.0*t51*t5-24.0*t111*t61+30.0*t111
*t30-6.0*t114*t19-60.0*t114*t15+6.0*t37*t1+6.0*t41*t114;
  const double t140 = 30.0*t7*t1*t19+114.0*t7*t114+30.0*t7*t5*t19+30.0*t16*t25+30.0*t14*
t24-60.0*t7*t3*t19+30.0*t16*t4-60.0*t14*t25+30.0*t14*t4-60.0*t16*t24+6.0*t12*
t25-12.0*t81*t3+6.0*t11*t24+6.0*t12*t2+114.0*t4*t2;
  const double t141 = t1*d[1][2];
  const double t163 = -216.0*d[0][0]*t141*t9+6.0*t11*t6+114.0*t4*t6+24.0*t61*t15+6.0*t32*t1
-3.0*t53*t19+6.0*d[0][0]*t53*d[2][2]-6.0*t49*t5-3.0*t15*t43-3.0*t59*t19-24.0*t81*d[1][1]*
d[2][2]-36.0*t3*t30-36.0*t3*t61-36.0*t81*t1+24.0*t9*d[1][2]*t27+24.0*t9*d[1][2]*t28;
  const double t169 = t5*d[0][2];
  const double t170 = d[0][1]*t169;
  const double t173 = d[1][2]*d[0][0];
  const double t182 = d[1][2]*d[1][1];
  const double t187 = d[1][2]*t13;
  const double t191 = t3*d[0][1];
  const double t192 = t191*d[0][2];
  const double t198 = 24.0*t9*d[1][2]*t29+108.0*t9*t141*d[1][1]+108.0*t170*d[1][2]*d[2][2]+108.0*t170*
t173+108.0*t9*t141*d[2][2]-36.0*t9*d[1][2]*t15*d[2][2]-36.0*t9*t173*t19-36.0*t9*t182*t19
-36.0*t9*t173*t15-36.0*t9*t187*d[2][2]-60.0*t47*t1+108.0*t192*t182-36.0*t9*t187*d[1][1]
+108.0*t192*t173+144.0*t8*t10;
  const double t209 = t3*t27;
  const double t222 = -216.0*t169*d[1][1]*d[0][1]*d[1][2]-216.0*t191*d[2][2]*d[0][2]*d[1][2]-12.0*t30*t1-12.0*
t61*t5+6.0*t111*t114-6.0*t47*t19+6.0*t209*d[1][1]-60.0*t64*t15+30.0*t82*d[1][1]-6.0*
t209*d[2][2]-6.0*t3*t15*t19-24.0*t35*t1+24.0*t11*t28+24.0*t12*t28-6.0*t3*t29*d[2][2]+
252.0*t64*t1;
  const double t226 = sqrt(t57+t78+t99+t119+t140+t163+t198+t222);
  const double t227 = -12.0*d[0][0]*t19-12.0*d[1][1]*t19-12.0*t15*d[2][2]+36.0*t5*d[0][0]+36.0*t24+36.0*
t25+36.0*t1*d[1][1]+8.0*t27+8.0*t28+8.0*t29+12.0*t226;
  const double t229 = pow(-72.0*t2-72.0*t4-72.0*t6+48.0*t8+216.0*t10+36.0*t11+36.0*t12
-12.0*t14-12.0*t16-12.0*t17+t227,1.0/3.0);
  const double t232 = (-t3/3+t7/9+t41/9+t111/9-t5/3-t1/3-t13/9-t19/9-t15/9)/t229;
  const double t234 = sqrt(3.0);
  const double t236 = t234*(t229/6+6.0*t232);

  const double eigenvalues[3]
  = { t229/6-6.0*t232+d[0][0]/3+d[2][2]/3+d[1][1]/3,
      -t229/12+3.0*t232+d[0][0]/3+d[2][2]/3+d[1][1]/3+sqrt(-1.0)*t236/2,
      -t229/12+3.0*t232+d[0][0]/3+d[2][2]/3+d[1][1]/3-sqrt(-1.0)*t236/2 };
  
  return std::max (std::fabs (eigenvalues[0]),
		   std::max (std::fabs (eigenvalues[1]),
			     std::fabs (eigenvalues[2])));
};

#endif


template <int dim>
inline
double
DerivativeApproximation::SecondDerivative<dim>::
derivative_norm (const Derivative &)
{
				   // computing the spectral norm is
				   // not so simple in general. it is
				   // feasible for dim==3 as shown
				   // above, since then there are
				   // still closed form expressions of
				   // the roots of the characteristic
				   // polynomial, and they can easily
				   // be computed using
				   // maple. however, for higher
				   // dimensions, some other method
				   // needs to be employed. maybe some
				   // steps of the power method would
				   // suffice?
  Assert (false, ExcNotImplemented());
  return 0;
};



template <int dim>
inline
void
DerivativeApproximation::SecondDerivative<dim>::symmetrize (Derivative &d)
{
				   // symmetrize non-diagonal entries
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=i+1; j<dim; ++j)
      {
	const double s = (d[i][j] + d[j][i]) / 2;
	d[i][j] = d[j][i] = s;
      };
};




template <int dim>
void 
DerivativeApproximation::
approximate_gradient (const Mapping<dim>    &mapping,
		      const DoFHandler<dim> &dof_handler,
		      const Vector<double>  &solution,
		      Vector<float>         &derivative_norm,
		      const unsigned int     component)
{
  approximate_derivative<Gradient<dim>,dim> (mapping,
					     dof_handler,
					     solution,
					     component,
					     derivative_norm);
};


template <int dim>
void 
DerivativeApproximation::
approximate_gradient (const DoFHandler<dim> &dof_handler,
		      const Vector<double>  &solution,
		      Vector<float>         &derivative_norm,
		      const unsigned int     component)
{
  static const MappingQ1<dim> mapping;
  approximate_derivative<Gradient<dim>,dim> (mapping,
					     dof_handler,
					     solution,
					     component,
					     derivative_norm);
};


template <int dim>
void 
DerivativeApproximation::
approximate_second_derivative (const Mapping<dim>    &mapping,
			       const DoFHandler<dim> &dof_handler,
			       const Vector<double>  &solution,
			       Vector<float>         &derivative_norm,
			       const unsigned int     component)
{
  approximate_derivative<SecondDerivative<dim>,dim> (mapping,
						     dof_handler,
						     solution,
						     component,
						     derivative_norm);
};


template <int dim>
void 
DerivativeApproximation::
approximate_second_derivative (const DoFHandler<dim> &dof_handler,
			       const Vector<double>  &solution,
			       Vector<float>         &derivative_norm,
			       const unsigned int     component)
{
  static const MappingQ1<dim> mapping;
  approximate_derivative<SecondDerivative<dim>,dim> (mapping,
						     dof_handler,
						     solution,
						     component,
						     derivative_norm);
};


template <class DerivativeDescription, int dim>
void 
DerivativeApproximation::
approximate_derivative (const Mapping<dim>    &mapping,
			const DoFHandler<dim> &dof_handler,
			const Vector<double>  &solution,
			const unsigned int     component,
			Vector<float>         &derivative_norm)
{
  Assert (derivative_norm.size() == dof_handler.get_tria().n_active_cells(),
	  ExcInvalidVectorLength (derivative_norm.size(),
				  dof_handler.get_tria().n_active_cells()));

  const unsigned int n_threads = multithread_info.n_default_threads;
  std::vector<IndexInterval> index_intervals
    = Threads::split_interval (0, dof_handler.get_tria().n_active_cells(),
			       n_threads);
  Threads::ThreadManager thread_manager;
  for (unsigned int i=0; i<n_threads; ++i)
    Threads::spawn (thread_manager,
		    Threads::encapsulate
		    (&DerivativeApproximation::
		     template approximate<DerivativeDescription,dim>)
		    .collect_args (mapping, dof_handler, solution, component,
				   index_intervals[i],
				   derivative_norm));
  thread_manager.wait ();
};



template <class DerivativeDescription, int dim>
void 
DerivativeApproximation::approximate (const Mapping<dim>    &mapping,
				      const DoFHandler<dim> &dof_handler,
				      const Vector<double>  &solution,
				      const unsigned int     component,
				      const IndexInterval   &index_interval,
				      Vector<float>         &derivative_norm)
{
  QMidpoint<dim> midpoint_rule;
  FEValues<dim>  fe_midpoint_value (mapping, dof_handler.get_fe(),
				    midpoint_rule,
				    UpdateFlags(DerivativeDescription::update_flags |
						update_q_points));
  
				   // matrix Y=sum_i y_i y_i^T
  Tensor<2,dim> Y;
  
				   // iterators over all cells and the
				   // respective entries in the output
				   // vector:
  Vector<float>::iterator
    derivative_norm_on_this_cell
    = derivative_norm.begin() + index_interval.first;
  
  typename DoFHandler<dim>::active_cell_iterator cell, endc;
  cell = endc = dof_handler.begin_active();
				   // (static_cast to avoid warnings
				   // about unsigned always >=0)
  advance (cell, static_cast<int>(index_interval.first));
  advance (endc, static_cast<int>(index_interval.second));

				   // vector to hold iterators to all
				   // active neighbors of a cell
				   // reserve the maximal number of
				   // active neighbors
  std::vector<typename DoFHandler<dim>::active_cell_iterator> active_neighbors;
  active_neighbors.reserve (GeometryInfo<dim>::faces_per_cell *
			    GeometryInfo<dim>::subfaces_per_face);

  for (; cell!=endc; ++cell, ++derivative_norm_on_this_cell)
    {
      Y.clear ();
				       // vector
				       // g=sum_i y_i (f(x+y_i)-f(x))/|y_i|
				       // or related type for higher
				       // derivatives
      typename DerivativeDescription::Derivative projected_derivative;

				       // reinit fe values object...
      fe_midpoint_value.reinit (cell);

				       // ...and get the value of the
				       // projected derivative...
      const typename DerivativeDescription::ProjectedDerivative
	this_midpoint_value
	= DerivativeDescription::get_projected_derivative (fe_midpoint_value,
							   solution,
							   component);
      				       // ...and the place where it lives
      const Point<dim> this_center = fe_midpoint_value.quadrature_point(0);

      
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
	    typename DoFHandler<dim>::cell_iterator
	      neighbor = cell->neighbor(n);
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
		    typename DoFHandler<dim>::cell_iterator
		      neighbor_child = neighbor;
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
      typename std::vector<typename DoFHandler<dim>::active_cell_iterator>::const_iterator
	neighbor_ptr = active_neighbors.begin();
      for (; neighbor_ptr!=active_neighbors.end(); ++neighbor_ptr)
	{
	  const typename DoFHandler<dim>::active_cell_iterator
	    neighbor = *neighbor_ptr;
	    
					   // reinit fe values object...
	  fe_midpoint_value.reinit (neighbor);
	  
					   // ...and get the value of the
					   // solution...
	  const typename DerivativeDescription::ProjectedDerivative
	    neighbor_midpoint_value
	    = DerivativeDescription::get_projected_derivative (fe_midpoint_value,
							       solution, component);
	  
					   // ...and the place where it lives
	  const Point<dim>
	    neighbor_center = fe_midpoint_value.quadrature_point(0);
	  
	  
					   // vector for the
					   // normalized
					   // direction between
					   // the centers of two
					   // cells
	  Point<dim>   y        = neighbor_center - this_center;
	  const double distance = sqrt(y.square());
					   // normalize y
	  y /= distance;
					   // *** note that unlike in
					   // the docs, y denotes the
					   // normalized vector
					   // connecting the centers
					   // of the two cells, rather
					   // than the normal
					   // difference! ***
	  
					   // add up the
					   // contribution of
					   // this cell to Y
	  for (unsigned int i=0; i<dim; ++i)
	    for (unsigned int j=0; j<dim; ++j)
	      Y[i][j] += y[i] * y[j];
	  
					   // then update the sum
					   // of difference
					   // quotients
	  typename DerivativeDescription::ProjectedDerivative
	    projected_finite_difference
	    = (neighbor_midpoint_value -
	       this_midpoint_value);
	  projected_finite_difference /= distance;
	  
	  typename DerivativeDescription::Derivative projected_derivative_update;
	  outer_product (projected_derivative_update,
			 y,
			 projected_finite_difference);
	  projected_derivative += projected_derivative_update;
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

				       // first symmetrize g
      DerivativeDescription::symmetrize (projected_derivative);
      
                                       // compute Y^-1 g
      typename DerivativeDescription::Derivative derivative;
      Tensor<2,dim> Y_inverse = invert(Y);
      
      contract (derivative, Y_inverse, projected_derivative);

      *derivative_norm_on_this_cell
	= DerivativeDescription::derivative_norm (derivative);
    };
};




// explicit instantiations
template
void 
DerivativeApproximation::
approximate_gradient (const Mapping<deal_II_dimension> &mapping,
		      const DoFHandler<deal_II_dimension> &dof_handler,
		      const Vector<double>  &solution,
		      Vector<float>         &derivative_norm,
		      const unsigned int     component);
template
void 
DerivativeApproximation::
approximate_gradient (const DoFHandler<deal_II_dimension> &dof_handler,
		      const Vector<double>  &solution,
		      Vector<float>         &derivative_norm,
		      const unsigned int     component);

template
void 
DerivativeApproximation::
approximate_second_derivative (const Mapping<deal_II_dimension>    &mapping,
			       const DoFHandler<deal_II_dimension> &dof_handler,
			       const Vector<double>  &solution,
			       Vector<float>         &derivative_norm,
			       const unsigned int     component);
template
void 
DerivativeApproximation::
approximate_second_derivative (const DoFHandler<deal_II_dimension> &dof_handler,
			       const Vector<double>  &solution,
			       Vector<float>         &derivative_norm,
			       const unsigned int     component);



