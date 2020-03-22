/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2010 Sandia National Laboratories.  Developed at the
    University of Wisconsin--Madison under SNL contract number
    624796.  The U.S. Government and the University of Wisconsin
    retain certain rights to this software.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License 
    (lgpl.txt) along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    (2010) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file TMetric.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "TMetric.hpp"
#include "TMetricBarrier.hpp"
#include "MsqMatrix.hpp"
#include "MsqError.hpp"
#include <limits>

namespace MESQUITE_NS {

template <unsigned Dim>
static inline double
do_finite_difference( int r, int c, TMetric* metric, 
                      MsqMatrix<Dim, Dim> A,
                      double value, MsqError& err )
{
  const double INITIAL_STEP = std::max( 1e-6, fabs(1e-14*value) );
  const double init = A(r,c);
  bool valid;
  double diff_value;
  for (double step = INITIAL_STEP; step > std::numeric_limits<double>::epsilon(); step *= 0.1) {
    A(r,c) = init + step;
    valid = metric->evaluate( A, diff_value, err );
    MSQ_ERRZERO(err);
    if (valid)
      return (diff_value - value) / step;
  }
  
    // If we couldn't find a valid step, try stepping in the other
    // direciton
  for (double step = INITIAL_STEP; step > std::numeric_limits<double>::epsilon(); step *= 0.1) {
    A(r,c) = init - step;
    valid = metric->evaluate( A, diff_value, err );
    MSQ_ERRZERO(err);
    if (valid)
      return (value - diff_value) / step;
  }  
    // If that didn't work either, then give up.
  MSQ_SETERR(err)("No valid step size for finite difference of 2D target metric.",
                  MsqError::INTERNAL_ERROR);
  return 0.0;
}


template <unsigned Dim>
static inline bool
do_numerical_gradient( TMetric* mu,
                       MsqMatrix<Dim, Dim> A,
                       double& result,
                       MsqMatrix<Dim,Dim>& wrt_A,
                       MsqError& err )
{
  bool valid;
  valid = mu->evaluate( A, result, err );
  MSQ_ERRZERO(err);
  if (MSQ_CHKERR(err) || !valid)
    return valid;
  
  switch (Dim) {
    case 3:
  wrt_A(0,2) = do_finite_difference( 0, 2, mu, A, result, err ); MSQ_ERRZERO(err);
  wrt_A(1,2) = do_finite_difference( 1, 2, mu, A, result, err ); MSQ_ERRZERO(err);
  wrt_A(2,0) = do_finite_difference( 2, 0, mu, A, result, err ); MSQ_ERRZERO(err);
  wrt_A(2,1) = do_finite_difference( 2, 1, mu, A, result, err ); MSQ_ERRZERO(err);
  wrt_A(2,2) = do_finite_difference( 2, 2, mu, A, result, err ); MSQ_ERRZERO(err);
    case 2:
  wrt_A(0,1) = do_finite_difference( 0, 1, mu, A, result, err ); MSQ_ERRZERO(err);
  wrt_A(1,0) = do_finite_difference( 1, 0, mu, A, result, err ); MSQ_ERRZERO(err);
  wrt_A(1,1) = do_finite_difference( 1, 1, mu, A, result, err ); MSQ_ERRZERO(err);
    case 1:
  wrt_A(0,0) = do_finite_difference( 0, 0, mu, A, result, err ); MSQ_ERRZERO(err);
    break;
    default:
     assert(false);
  }
  return true;
}

template <unsigned Dim>
static inline bool
do_numerical_hessian( TMetric* metric, 
                      MsqMatrix<Dim, Dim> A,
                      double& value,
                      MsqMatrix<Dim, Dim>& grad, 
                      MsqMatrix<Dim, Dim>* Hess, 
                      MsqError& err )
{
    // zero hessian data
  const int num_block = Dim * (Dim + 1) / 2;
  for (int i = 0; i < num_block; ++i)
    Hess[i].zero();

    // evaluate gradient for input values
  bool valid;
  valid = metric->evaluate_with_grad( A, value, grad, err );
  if (MSQ_CHKERR(err) || !valid)
    return false;
  
    // do finite difference for each term of A
  const double INITAL_STEP = std::max( 1e-6, fabs(1e-14*value) );
  double value2;
  MsqMatrix<Dim,Dim> grad2;
  for (unsigned r = 0; r < Dim; ++r) {  // for each row of A
    for (unsigned c = 0; c < Dim; ++c) {  // for each column of A
      const double in_val = A(r,c);
      double step;
      for (step = INITAL_STEP; step > std::numeric_limits<double>::epsilon(); step *= 0.1) {
        A(r,c) = in_val + step;
        valid = metric->evaluate_with_grad( A, value2, grad2, err );
        MSQ_ERRZERO(err);
        if (valid)
          break;
      }
      
        // if no valid step size, try step in other direction
      if (!valid) {
        for (step = -INITAL_STEP; step < -std::numeric_limits<double>::epsilon(); step *= 0.1) {
          A(r,c) = in_val + step;
          valid = metric->evaluate_with_grad( A, value2, grad2, err );
          MSQ_ERRZERO(err);
          if (valid)
            break;
        }
        
          // if still no valid step size, give up.
        if (!valid) {
          MSQ_SETERR(err)("No valid step size for finite difference of 2D target metric.",
                          MsqError::INTERNAL_ERROR);
          return false;
        }
      }
      
        // restore A.
      A(r,c) = in_val;
      
        // add values into result matrix
        // values of grad2, in row-major order, are a single 9-value row of the Hessian
      grad2 -= grad;
      grad2 /= step;
      for (unsigned b = 0; b < r; ++b) {
        const int idx = Dim*b - b*(b+1)/2 + r;
        Hess[idx].add_column( c, transpose( grad2.row(b) ) );
      }
      for (unsigned b = r; b < Dim; ++b) {
        const int idx = Dim*r - r*(r+1)/2 + b;
        Hess[idx].add_row( c, grad2.row(b) );
      }
    } // for (c)
  } // for (r)
  
    // Values in non-diagonal blocks were added twice.
  for (unsigned r = 0, h = 1; r < Dim-1; ++r, ++h)
    for (unsigned c = r + 1; c < Dim; ++c, ++h)
      Hess[h] *= 0.5;
  
  return true;
}


TMetric::~TMetric() {}

bool TMetric::evaluate( const MsqMatrix<2,2>& T, 
               double& result, 
               MsqError& err )
{
  return false;
}

bool TMetric::evaluate( const MsqMatrix<3,3>& T, 
               double& result, 
               MsqError& err )
{
  return false;
}

bool TMetric::evaluate_with_grad( const MsqMatrix<2,2>& T,
                                  double& result,
                                  MsqMatrix<2,2>& wrt_T,
                                  MsqError& err )
{
  return do_numerical_gradient( this, T, result, wrt_T, err );
}

bool TMetric::evaluate_with_grad( const MsqMatrix<3,3>& T,
                                  double& result,
                                  MsqMatrix<3,3>& wrt_T,
                                  MsqError& err )
{
  return do_numerical_gradient( this, T, result, wrt_T, err );
}

bool TMetric::evaluate_with_hess( const MsqMatrix<2,2>& T,
                                  double& result,
                                  MsqMatrix<2,2>& deriv_wrt_T,
                                  MsqMatrix<2,2> hess_wrt_T[3],
                                  MsqError& err )
{
  return do_numerical_hessian( this, T, result, deriv_wrt_T, hess_wrt_T, err );
}

bool TMetric::evaluate_with_hess( const MsqMatrix<3,3>& T,
                                  double& result,
                                  MsqMatrix<3,3>& deriv_wrt_T,
                                  MsqMatrix<3,3> hess_wrt_T[6],
                                  MsqError& err )
{
  return do_numerical_hessian( this, T, result, deriv_wrt_T, hess_wrt_T, err );
}

TMetric2D::~TMetric2D() {}
TMetric3D::~TMetric3D() {}

bool TMetric2D::evaluate( const MsqMatrix<3,3>&, double&, MsqError& err )
{
  MSQ_SETERR(err)("2D target metric cannot be evaluated for volume elements",
                  MsqError::UNSUPPORTED_ELEMENT);
  return false;
}

bool TMetric3D::evaluate( const MsqMatrix<2,2>&, double&, MsqError& err )
{
  MSQ_SETERR(err)("2D target metric cannot be evaluated for volume elements",
                  MsqError::UNSUPPORTED_ELEMENT);
  return false;
}


} // namespace Mesquite

