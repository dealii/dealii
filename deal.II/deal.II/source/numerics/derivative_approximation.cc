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



// static variables
template <int dim>
const UpdateFlags DerivativeApproximation::Gradient<dim>::update_flags;

template <int dim>
const UpdateFlags DerivativeApproximation::SecondDerivative<dim>::update_flags;



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
/*
  compute the three eigenvalues of the tensor @p{d} and take the
  largest. one could use the following maple script to generate C
  code:
  
  with(linalg);
  readlib(C);
  A:=matrix(3,3,[[a00,a01,a02],[a01,a11,a12],[a02,a12,a22]]);
  E:=eigenvals(A);
  EE:=vector(3,[E[1],E[2],E[3]]);
  C(EE);

  Unfortunately, with both optimized and non-optimized output, at some
  places the code `sqrt(-1.0)' is emitted, and I don't know what
  Maple intends to do with it. This happens both with Maple4 and
  Maple5.

  Fortunately, Roger Young provided the following Fortran code, which
  is transcribed below to C. The code uses an algorithm that uses the
  invariants of a symmetric matrix. (The translated algorithm is
  augmented by a test for R>0, since R==0 indicates that all three
  eigenvalues are equal.)

  
      PROGRAM MAIN

C FIND EIGENVALUES OF REAL SYMMETRIC MATRIX
C (ROGER YOUNG, 2001)

      IMPLICIT NONE

      REAL*8 A11, A12, A13, A22, A23, A33
      REAL*8 I1, J2, J3, AM
      REAL*8 S11, S12, S13, S22, S23, S33
      REAL*8 SS12, SS23, SS13
      REAL*8 R,R3, XX,YY, THETA
      REAL*8 A1,A2,A3
      REAL*8 PI
      PARAMETER (PI=3.141592653587932384D0)
      REAL*8 A,B,C, TOL
      PARAMETER (TOL=1.D-14)

C DEFINE A TEST MATRIX

      A11 = -1.D0
      A12 = 5.D0
      A13 = 3.D0
      A22 = -2.D0
      A23 = 0.5D0
      A33 = 4.D0
      

      I1 = A11 + A22 + A33
      AM = I1/3.D0

      S11 = A11 - AM
      S22 = A22 - AM
      S33 = A33 - AM
      S12 = A12
      S13 = A13
      S23 = A23

      SS12 = S12*S12
      SS23 = S23*S23
      SS13 = S13*S13

      J2 = S11*S11 + S22*S22 + S33*S33
      J2 = J2 + 2.D0*(SS12 + SS23 + SS13)
      J2 = J2/2.D0

      J3 = S11**3 + S22**3 + S33**3
      J3 = J3 + 3.D0*S11*(SS12 + SS13)
      J3 = J3 + 3.D0*S22*(SS12 + SS23)
      J3 = J3 + 3.D0*S33*(SS13 + SS23)
      J3 = J3 + 6.D0*S12*S23*S13
      J3 = J3/3.D0

      R = SQRT(4.D0*J2/3.D0)
      R3 = R*R*R
      XX = 4.D0*J3/R3

      YY = 1.D0 - DABS(XX)
      IF(YY.LE.0.D0)THEN
         IF(YY.GT.(-TOL))THEN
            WRITE(6,*)'Equal roots: XX= ',XX
            A = -(XX/DABS(XX))*SQRT(J2/3.D0)
            B = AM + A
            C = AM - 2.D0*A
            WRITE(6,*)B,' (twice) ',C
            STOP
         ELSE
            WRITE(6,*)'Error: XX= ',XX
            STOP
         ENDIF
      ENDIF

      THETA = (ACOS(XX))/3.D0
      
      A1 = AM + R*COS(THETA)
      A2 = AM + R*COS(THETA + 2.D0*PI/3.D0)
      A3 = AM + R*COS(THETA + 4.D0*PI/3.D0)

      WRITE(6,*)A1,A2,A3

      STOP
      END
  
*/

  const double pi = 3.141592653587932384;
  const double am = trace(d) / 3.;

				   // s := d - trace(d) I
  Tensor<2,3> s = d;
  for (unsigned int i=0; i<3; ++i)
    s[i][i] -= am;
  
  const double ss01 = s[0][1] * s[0][1],
	       ss12 = s[1][2] * s[1][2],
	       ss02 = s[0][2] * s[0][2];

  const double J2 = (s[0][0]*s[0][0] + s[1][1]*s[1][1] + s[2][2]*s[2][2]
		     + 2 * (ss01 + ss02 + ss12))  / 2.;
  const double J3 = (std::pow(s[0][0],3) + std::pow(s[1][1],3) + std::pow(s[2][2],3)
		     + 3. * s[0][0] * (ss01 + ss02)
		     + 3. * s[1][1] * (ss01 + ss12)
		     + 3. * s[2][2] * (ss02 + ss12)
		     + 6. * s[0][1] * s[0][2] * s[1][2]) / 3.;
  
  const double R  = std::sqrt (4. * J2 / 3.);

  double EE[3] = { 0, 0, 0 };
				   // the eigenvalues are away from
				   // @p{am} in the order of R. thus,
				   // if R<<AM, then we have
				   // degenerate case with three
				   // identical eigenvalues. check
				   // this first
  if (R < 1e-14*am)
    EE[0] = EE[1] = EE[2] = am;
  else
    {
				       // at least two eigenvalues are
				       // distinct
      const double R3 = R*R*R;
      const double XX = 4. * J3 / R3;
      const double YY = 1. - std::fabs(XX);
      
      Assert (YY > -1e-14, ExcInternalError());
      
      if (YY < 0)
	{
					   // two roots are equal
	  const double a = (XX>0 ? -1. : 1.) * R / 2;
	  EE[0] = EE[1] = am + a;
	  EE[2] = am - 2.*a;
	}
      else
	{
	  const double theta = std::acos(XX) / 3.;
	  EE[0] = am + R*std::cos(theta);
	  EE[1] = am + R*std::cos(theta + 2./3.*pi);
	  EE[2] = am + R*std::cos(theta + 4./3.*pi);
	};
    };

  cout << "EE=" << EE[0] << ' ' << EE[1] << ' ' << EE[2] << endl;
  
  return std::max (std::fabs (EE[0]),
		   std::max (std::fabs (EE[1]),
			     std::fabs (EE[2])));
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
  Assert (DEAL_II_COMPAT_MAPPING, ExcCompatibility("mapping"));
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
  Assert (DEAL_II_COMPAT_MAPPING, ExcCompatibility("mapping"));
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



