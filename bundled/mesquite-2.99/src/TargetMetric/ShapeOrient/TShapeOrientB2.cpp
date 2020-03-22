/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Sandia National Laboratories.  Developed at the
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
 
    (2009) kraftche@cae.wisc.edu
   
  ***************************************************************** */


/** \file TShapeOrientB2.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "TShapeOrientB2.hpp"
#include "MsqMatrix.hpp"
#include "MsqError.hpp"
#include "TMPDerivs.hpp"
#include "TMPCommon.hpp"

namespace MESQUITE_NS {

std::string TShapeOrientB2::get_name() const
  { return "TShapeOrientB2"; }

TShapeOrientB2::~TShapeOrientB2() {}

bool TShapeOrientB2::evaluate( const MsqMatrix<2,2>& T, 
                               double& result, 
                               MsqError& err )
{
  double tau = det(T);
  if (TMetric::invalid_determinant(tau)) {
    MSQ_SETERR(err)( barrier_violated_msg, MsqError::BARRIER_VIOLATED );
    return false;
  }
  const double tr = trace(T);
  result = (0.5/tau) * (sqr_Frobenius( T ) - 0.5 * tr * fabs(tr));
  return true;
}

bool TShapeOrientB2::evaluate_with_grad( const MsqMatrix<2,2>& T, 
                                         double& result, 
                                         MsqMatrix<2,2>& deriv_wrt_T,
                                         MsqError& err )
{
  double tau = det(T);
  if (TMetric::invalid_determinant(tau)) {
    MSQ_SETERR(err)( barrier_violated_msg, MsqError::BARRIER_VIOLATED );
    return false;
  }
  const double b = 0.5/tau;

  const double tr = trace(T);
  const double f = 0.5 * fabs(tr);
  result = sqr_Frobenius( T ) - f * tr;
  
  deriv_wrt_T = T;
  pluseq_scaled_I( deriv_wrt_T, -f );
  deriv_wrt_T *= 2*tau;
  deriv_wrt_T -= result * transpose_adj(T);
  
  result *= b;
  deriv_wrt_T *= b/tau;
  
  return true;
}

bool TShapeOrientB2::evaluate_with_hess( const MsqMatrix<2,2>& T, 
                                         double& result, 
                                         MsqMatrix<2,2>& deriv_wrt_T,
                                         MsqMatrix<2,2> second_wrt_T[3],
                                         MsqError& err )
{
  double tau = det(T);
  if (TMetric::invalid_determinant(tau)) {
    MSQ_SETERR(err)( barrier_violated_msg, MsqError::BARRIER_VIOLATED );
    return false;
  }
  const double b = 0.5/tau;

    // calculate non-barrier value (ShapeOrientAlt1)
  const double tr = trace(T);
  const double f = 0.5 * fabs(tr);
  result = sqr_Frobenius( T ) - f * tr;
  
    // calculate non-barrier first derivatives
  deriv_wrt_T = T;
  pluseq_scaled_I( deriv_wrt_T, -f );
  deriv_wrt_T *= 2;
  
    // calculate barrier second derivs
  const MsqMatrix<2,2> adjt = transpose_adj(T);
  set_scaled_sum_outer_product( second_wrt_T, -b/tau, deriv_wrt_T, adjt );
  pluseq_scaled_outer_product( second_wrt_T, result/(tau*tau*tau), adjt );
  pluseq_scaled_2nd_deriv_of_det( second_wrt_T, -result * b / tau, T );
    // calculate non-barrier barrier portion of second derivs
  pluseq_scaled_I( second_wrt_T, 1/tau );
  pluseq_scaled_outer_product_I_I( second_wrt_T, 0.5/tau * (tr < 0 ? 1 : -1) );
  
    // calculate barrier derivs from non-barrier
  deriv_wrt_T *= tau;
  deriv_wrt_T -= result * adjt;
  deriv_wrt_T *= b/tau;
  
    // barrier value from non-barrier
  result *= b;
  return true;
}

bool TShapeOrientB2::evaluate( const MsqMatrix<3,3>& T, 
                               double& result, 
                               MsqError& err )
{
  double tau = det(T);
  if (TMetric::invalid_determinant(tau)) {
    MSQ_SETERR(err)( barrier_violated_msg, MsqError::BARRIER_VIOLATED );
    return false;
  }
  const double tr = trace(T);
  result = (0.5/tau) * (sqr_Frobenius( T ) - MSQ_ONE_THIRD * tr * fabs(tr));
  return true;
}

bool TShapeOrientB2::evaluate_with_grad( const MsqMatrix<3,3>& T, 
                                         double& result, 
                                         MsqMatrix<3,3>& deriv_wrt_T,
                                         MsqError& err )
{
  double tau = det(T);
  if (TMetric::invalid_determinant(tau)) {
    MSQ_SETERR(err)( barrier_violated_msg, MsqError::BARRIER_VIOLATED );
    return false;
  }
  const double b = 0.5/tau;

  const double tr = trace(T);
  const double f = MSQ_ONE_THIRD * fabs(tr);
  result = sqr_Frobenius( T ) - f * tr;
  
  deriv_wrt_T = T;
  pluseq_scaled_I( deriv_wrt_T, -f );
  deriv_wrt_T *= 2*tau;
  deriv_wrt_T -= result * transpose_adj(T);
  
  result *= b;
  deriv_wrt_T *= b/tau;
  
  return true;
}

bool TShapeOrientB2::evaluate_with_hess( const MsqMatrix<3,3>& T, 
                                         double& result, 
                                         MsqMatrix<3,3>& deriv_wrt_T,
                                         MsqMatrix<3,3> second_wrt_T[6],
                                         MsqError& err )
{
  double tau = det(T);
  if (TMetric::invalid_determinant(tau)) {
    MSQ_SETERR(err)( barrier_violated_msg, MsqError::BARRIER_VIOLATED );
    result = 0.0;
    return false;
  }
  const double b = 0.5/tau;

    // calculate non-barrier value (ShapeOrientAlt1)
  const double tr = trace(T);
  const double f = MSQ_ONE_THIRD * fabs(tr);
  result = sqr_Frobenius( T ) - f * tr;
  
    // calculate non-barrier first derivatives
  deriv_wrt_T = T;
  pluseq_scaled_I( deriv_wrt_T, -f );
  deriv_wrt_T *= 2;
  
    // calculate barrier second derivs
  const MsqMatrix<3,3> adjt = transpose_adj(T);
  set_scaled_sum_outer_product( second_wrt_T, -b/tau, deriv_wrt_T, adjt );
  pluseq_scaled_outer_product( second_wrt_T, result/(tau*tau*tau), adjt );
  pluseq_scaled_2nd_deriv_of_det( second_wrt_T, -result * b / tau, T );
    // calculate non-barrier barrier portion of second derivs
  pluseq_scaled_I( second_wrt_T, 1/tau );
  pluseq_scaled_outer_product_I_I( second_wrt_T, MSQ_ONE_THIRD/tau * (tr < 0 ? 1 : -1) );
  
    // calculate barrier derivs from non-barrier
  deriv_wrt_T *= tau;
  deriv_wrt_T -= result * adjt;
  deriv_wrt_T *= b/tau;
  
    // barrier value from non-barrier
  result *= b;
  return true;
}

} // namespace Mesquite
