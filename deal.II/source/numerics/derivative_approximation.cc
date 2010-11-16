//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <base/quadrature_lib.h>
#include <base/thread_management.h>
#include <base/multithread_info.h>
#include <lac/vector.h>
#include <lac/block_vector.h>
#include <lac/petsc_vector.h>
#include <lac/petsc_block_vector.h>
#include <grid/tria_iterator.h>
#include <grid/grid_tools.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_handler.h>
#include <fe/fe.h>
#include <fe/fe_values.h>
#include <hp/fe_values.h>
#include <fe/mapping_q1.h>
#include <hp/q_collection.h>
#include <hp/fe_collection.h>
#include <hp/mapping_collection.h>
#include <numerics/derivative_approximation.h>

#include <cmath>

DEAL_II_NAMESPACE_OPEN



template <typename T>
static inline T sqr (const T t)
{
  return t*t;
}



// static variables
template <int dim>
const UpdateFlags DerivativeApproximation::Gradient<dim>::update_flags = update_values;

template <int dim>
const UpdateFlags DerivativeApproximation::SecondDerivative<dim>::update_flags = update_gradients;

template <int dim>
const UpdateFlags DerivativeApproximation::ThirdDerivative<dim>::update_flags = update_hessians;



template <int dim>
template <class InputVector, int spacedim>
inline
typename DerivativeApproximation::Gradient<dim>::ProjectedDerivative
DerivativeApproximation::Gradient<dim>::
get_projected_derivative (const FEValues<dim,spacedim>  &fe_values,
			  const InputVector    &solution,
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
    }
}



template <int dim>
inline
double
DerivativeApproximation::Gradient<dim>::derivative_norm (const Derivative &d)
{
  double s = 0;
  for (unsigned int i=0; i<dim; ++i)
    s += d[i]*d[i];
  return std::sqrt(s);
}



template <int dim>
inline
void
DerivativeApproximation::Gradient<dim>::symmetrize (Derivative &)
{
				   // nothing to do here
}



template <int dim>
template <class InputVector, int spacedim>
inline
typename DerivativeApproximation::SecondDerivative<dim>::ProjectedDerivative
DerivativeApproximation::SecondDerivative<dim>::
get_projected_derivative (const FEValues<dim,spacedim>  &fe_values,
			  const InputVector    &solution,
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
}



template <>
inline
double
DerivativeApproximation::SecondDerivative<1>::
derivative_norm (const Derivative &d)
{
  return std::fabs (d[0][0]);
}



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
  const double radicand = dealii::sqr(d[0][0] - d[1][1]) +
			  4*dealii::sqr(d[0][1]);
  const double eigenvalues[2]
    = { 0.5*(d[0][0] + d[1][1] + std::sqrt(radicand)),
	0.5*(d[0][0] + d[1][1] - std::sqrt(radicand))  };

  return std::max (std::fabs (eigenvalues[0]),
		   std::fabs (eigenvalues[1]));
}



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
				   // if R<<AM, then we have the
				   // degenerate case with three
				   // identical eigenvalues. check
				   // this first
  if (R <= 1e-14*std::fabs(am))
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
	  EE[1] = am + R*std::cos(theta + 2./3.*numbers::PI);
	  EE[2] = am + R*std::cos(theta + 4./3.*numbers::PI);
	};
    };

  return std::max (std::fabs (EE[0]),
		   std::max (std::fabs (EE[1]),
			     std::fabs (EE[2])));
}



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
}



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
}


template <int dim>
template <class InputVector, int spacedim>
inline
typename DerivativeApproximation::ThirdDerivative<dim>::ProjectedDerivative
DerivativeApproximation::ThirdDerivative<dim>::
get_projected_derivative (const FEValues<dim,spacedim>  &fe_values,
			  const InputVector    &solution,
			  const unsigned int    component)
{
  if (fe_values.get_fe().n_components() == 1)
    {
      std::vector<ProjectedDerivative> values (1);
      fe_values.get_function_hessians (solution, values);
      return values[0];
    }
  else
    {
      std::vector<std::vector<ProjectedDerivative> > values
	(1, std::vector<ProjectedDerivative>(fe_values.get_fe().n_components()));
      fe_values.get_function_hessians (solution, values);
      return values[0][component];
    };
}



template <>
inline
double
DerivativeApproximation::ThirdDerivative<1>::
derivative_norm (const Derivative &d)
{
  return std::fabs (d[0][0][0]);
}



template <int dim>
inline
double
DerivativeApproximation::ThirdDerivative<dim>::
derivative_norm (const Derivative &d)
{
				   // return the Frobenius-norm. this is a
				   // member function of Tensor<rank_,dim>
  return d.norm();
}


template <int dim>
inline
void
DerivativeApproximation::ThirdDerivative<dim>::symmetrize (Derivative &d)
{
				   // symmetrize non-diagonal entries

				   // first do it in the case, that i,j,k are
				   // pairwise different (which can onlky happen
				   // in dim >= 3)
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=i+1; j<dim; ++j)
      for (unsigned int k=j+1; k<dim; ++k)
	{
	  const double s = (d[i][j][k] +
			    d[i][k][j] +
			    d[j][i][k] +
			    d[j][k][i] +
			    d[k][i][j] +
			    d[k][j][i]) / 6;
	  d[i][j][k]
	    = d[i][k][j]
	    = d[j][i][k]
	    = d[j][k][i]
	    = d[k][i][j]
	    = d[k][j][i]
	    = s;
	};
				   // now do the case, where two indices are
				   // equal
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=i+1; j<dim; ++j)
      {
					 // case 1: index i (lower one) is
					 // double
	const double s = (d[i][i][j] +
			  d[i][j][i] +
			  d[j][i][i] ) / 3;
	d[i][i][j]
	  = d[i][j][i]
	  = d[j][i][i]
	  = s;

					 // case 2: index j (higher one) is
					 // double
	const double t = (d[i][j][j] +
			  d[j][i][j] +
			  d[j][j][i] ) / 3;
	d[i][j][j]
	  = d[j][i][j]
	  = d[j][j][i]
	  = t;
      };

}



template <int dim, template <int, int> class DH, class InputVector, int spacedim>
void
DerivativeApproximation::
approximate_gradient (const Mapping<dim,spacedim>    &mapping,
		      const DH<dim,spacedim>         &dof_handler,
		      const InputVector     &solution,
		      Vector<float>         &derivative_norm,
		      const unsigned int     component)
{
  approximate_derivative<Gradient<dim>,dim> (mapping,
                                             dof_handler,
                                             solution,
                                             component,
                                             derivative_norm);
}


template <int dim, template <int, int> class DH, class InputVector, int spacedim>
void
DerivativeApproximation::
approximate_gradient (const DH<dim,spacedim>         &dof_handler,
		      const InputVector     &solution,
		      Vector<float>         &derivative_norm,
		      const unsigned int     component)
{
  Assert (DEAL_II_COMPAT_MAPPING, ExcCompatibility("mapping"));
  approximate_derivative<Gradient<dim>,dim> (StaticMappingQ1<dim>::mapping,
                                             dof_handler,
                                             solution,
                                             component,
                                             derivative_norm);
}


template <int dim, template <int, int> class DH, class InputVector, int spacedim>
void
DerivativeApproximation::
approximate_second_derivative (const Mapping<dim,spacedim>    &mapping,
			       const DH<dim,spacedim>         &dof_handler,
			       const InputVector     &solution,
			       Vector<float>         &derivative_norm,
			       const unsigned int     component)
{
  approximate_derivative<SecondDerivative<dim>,dim> (mapping,
						     dof_handler,
						     solution,
						     component,
						     derivative_norm);
}


template <int dim, template <int, int> class DH, class InputVector, int spacedim>
void
DerivativeApproximation::
approximate_second_derivative (const DH<dim,spacedim>         &dof_handler,
			       const InputVector     &solution,
			       Vector<float>         &derivative_norm,
			       const unsigned int     component)
{
  Assert (DEAL_II_COMPAT_MAPPING, ExcCompatibility("mapping"));
  approximate_derivative<SecondDerivative<dim>,dim> (StaticMappingQ1<dim>::mapping,
						     dof_handler,
						     solution,
						     component,
						     derivative_norm);
}


template <int dim, template <int, int> class DH, class InputVector, int order, int spacedim>
void
DerivativeApproximation::
approximate_derivative_tensor (const Mapping<dim,spacedim>                           &mapping,
			       const DH<dim,spacedim>                                &dof,
			       const InputVector                            &solution,
			       const typename DH<dim,spacedim>::active_cell_iterator &cell,
			       Tensor<order,dim>                            &derivative,
			       const unsigned int                            component)
{
  approximate_cell<typename DerivativeSelector<order,dim>::DerivDescr,dim,DH,InputVector>
    (mapping,
     dof,
     solution,
     component,
     cell,
     derivative);
}



template <int dim, template <int, int> class DH, class InputVector, int order, int spacedim>
void
DerivativeApproximation::
approximate_derivative_tensor (const DH<dim,spacedim>                                &dof,
			       const InputVector                            &solution,
			       const typename DH<dim,spacedim>::active_cell_iterator &cell,
			       Tensor<order,dim>                            &derivative,
			       const unsigned int                            component)
{
				   // just call the respective function with Q1
				   // mapping
  Assert (DEAL_II_COMPAT_MAPPING, ExcCompatibility("mapping"));
  approximate_derivative_tensor (StaticMappingQ1<dim>::mapping,
				 dof,
				 solution,
				 cell,
				 derivative,
				 component);
}



template <class DerivativeDescription, int dim,
          template <int, int> class DH, class InputVector, int spacedim>
void
DerivativeApproximation::
approximate_derivative (const Mapping<dim,spacedim>    &mapping,
			const DH<dim,spacedim>         &dof_handler,
			const InputVector     &solution,
			const unsigned int     component,
			Vector<float>         &derivative_norm)
{
  Assert (derivative_norm.size() == dof_handler.get_tria().n_active_cells(),
	  ExcInvalidVectorLength (derivative_norm.size(),
				  dof_handler.get_tria().n_active_cells()));
  Assert (component < dof_handler.get_fe().n_components(),
	  ExcIndexRange (component, 0, dof_handler.get_fe().n_components()));

  const unsigned int n_threads = multithread_info.n_default_threads;
  std::vector<IndexInterval> index_intervals
    = Threads::split_interval (0, dof_handler.get_tria().n_active_cells(),
			       n_threads);

  typedef void (*FunPtr) (const Mapping<dim,spacedim>    &,
                          const DH<dim,spacedim>         &,
                          const InputVector     &,
                          const unsigned int     ,
                          const IndexInterval   &,
                          Vector<float>         &);
  FunPtr fun_ptr = &DerivativeApproximation::
		   template approximate<DerivativeDescription,dim,
		   DH,InputVector>;

//TODO: Use WorkStream here
  Threads::TaskGroup<> tasks;
  for (unsigned int i=0; i<n_threads; ++i)
    tasks += Threads::new_task (fun_ptr,
				mapping, dof_handler,
				solution, component,
				index_intervals[i],
				derivative_norm);
  tasks.join_all ();
}



template <class DerivativeDescription, int dim,
          template <int, int> class DH, class InputVector, int spacedim>
void
DerivativeApproximation::approximate (const Mapping<dim,spacedim>    &mapping,
				      const DH<dim,spacedim>         &dof_handler,
				      const InputVector     &solution,
				      const unsigned int     component,
				      const IndexInterval   &index_interval,
				      Vector<float>         &derivative_norm)
{
				   // iterators over all cells and the
				   // respective entries in the output
				   // vector:
  Vector<float>::iterator
    derivative_norm_on_this_cell
    = derivative_norm.begin() + index_interval.first;

  typename DH<dim,spacedim>::active_cell_iterator cell, endc;
  cell = endc = dof_handler.begin_active();
				   // (static_cast to avoid warnings
				   // about unsigned always >=0)
  std::advance (cell, static_cast<int>(index_interval.first));
  std::advance (endc, static_cast<int>(index_interval.second));

  for (; cell!=endc; ++cell, ++derivative_norm_on_this_cell)
    {
      typename DerivativeDescription::Derivative derivative;
				       // call the function doing the actual
				       // work on this cell
      DerivativeApproximation::
	template approximate_cell<DerivativeDescription,dim,DH,InputVector>
	(mapping,
	 dof_handler,
	 solution,
	 component,
	 cell,
	 derivative);
				       // evaluate the norm and fill the vector
      *derivative_norm_on_this_cell
	= DerivativeDescription::derivative_norm (derivative);
    }
}


template <class DerivativeDescription, int dim,
          template <int, int> class DH, class InputVector, int spacedim>
void
DerivativeApproximation::
approximate_cell (const Mapping<dim,spacedim>                   &mapping,
		  const DH<dim,spacedim>                        &dof_handler,
		  const InputVector                             &solution,
		  const unsigned int                             component,
		  const typename DH<dim,spacedim>::active_cell_iterator  &cell,
		  typename DerivativeDescription::Derivative    &derivative)
{
  QMidpoint<dim> midpoint_rule;

				   // create collection objects from
				   // single quadratures, mappings,
				   // and finite elements. if we have
				   // an hp DoFHandler,
				   // dof_handler.get_fe() returns a
				   // collection of which we do a
				   // shallow copy instead
  const hp::QCollection<dim>       q_collection (midpoint_rule);
  const hp::FECollection<dim>      fe_collection(dof_handler.get_fe());
  const hp::MappingCollection<dim> mapping_collection (mapping);

  hp::FEValues<dim> x_fe_midpoint_value (mapping_collection, fe_collection,
					 q_collection,
					 DerivativeDescription::update_flags |
					 update_quadrature_points);

				   // matrix Y=sum_i y_i y_i^T
  Tensor<2,dim> Y;


				   // vector to hold iterators to all
				   // active neighbors of a cell
				   // reserve the maximal number of
				   // active neighbors
  std::vector<typename DH<dim,spacedim>::active_cell_iterator> active_neighbors;
  active_neighbors.reserve (GeometryInfo<dim>::faces_per_cell *
			    GeometryInfo<dim>::max_children_per_face);

				       // vector
				       // g=sum_i y_i (f(x+y_i)-f(x))/|y_i|
				       // or related type for higher
				       // derivatives
  typename DerivativeDescription::Derivative projected_derivative;

				       // reinit fe values object...
  x_fe_midpoint_value.reinit (cell);
  const FEValues<dim> &fe_midpoint_value
    = x_fe_midpoint_value.get_present_fe_values();

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
  GridTools::template get_active_neighbors<DH<dim,spacedim> >(cell, active_neighbors);

				       // now loop over all active
				       // neighbors and collect the
				       // data we need
  typename std::vector<typename DH<dim,spacedim>::active_cell_iterator>::const_iterator
    neighbor_ptr = active_neighbors.begin();
  for (; neighbor_ptr!=active_neighbors.end(); ++neighbor_ptr)
    {
      const typename DH<dim,spacedim>::active_cell_iterator
	neighbor = *neighbor_ptr;

					   // reinit fe values object...
      x_fe_midpoint_value.reinit (neighbor);
      const FEValues<dim> &fe_midpoint_value
	= x_fe_midpoint_value.get_present_fe_values();

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
      const double distance = std::sqrt(y.square());
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

                                       // compute Y^-1 g
  const Tensor<2,dim> Y_inverse = invert(Y);

  contract (derivative, Y_inverse, projected_derivative);

				       // finally symmetrize the derivative
  DerivativeDescription::symmetrize (derivative);
}


template <int dim, int order>
double
DerivativeApproximation::
derivative_norm(const Tensor<order,dim> &derivative)
{
  return DerivativeSelector<order,dim>::DerivDescr::derivative_norm(derivative);
}



// --------------------------- explicit instantiations ---------------------
#include "derivative_approximation.inst"


DEAL_II_NAMESPACE_CLOSE
