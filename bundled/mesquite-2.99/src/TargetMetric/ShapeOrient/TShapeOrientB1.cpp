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
 
    (2006) kraftche@cae.wisc.edu
    (2010) kraftche@cae.wisc.edu
   
  ***************************************************************** */


/** \file TShapeOrientB1.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "TShapeOrientB1.hpp"
#include "MsqMatrix.hpp"
#include "MsqError.hpp"
#include "TMPDerivs.hpp"
#include "TMPCommon.hpp"

namespace MESQUITE_NS {

std::string TShapeOrientB1::get_name() const
  { return "TShapeOrientB1"; }

TShapeOrientB1::~TShapeOrientB1() {}

bool TShapeOrientB1::evaluate( const MsqMatrix<2,2>& T, 
                               double& result, 
                               MsqError& err )
{
  const double tau = det(T);
  if (TMetric::invalid_determinant(tau)) { // barrier
    MSQ_SETERR(err)( barrier_violated_msg, MsqError::BARRIER_VIOLATED );
    return false;
  }
  result = 0.5/tau * (Frobenius( T ) - trace(T)/MSQ_SQRT_TWO);
  return true;
}

bool TShapeOrientB1::evaluate_with_grad( const MsqMatrix<2,2>& T,
                                         double& result,
                                         MsqMatrix<2,2>& deriv_wrt_T,
                                         MsqError& err )
{
  const double norm = Frobenius(T);
  const double invroot = 1.0/MSQ_SQRT_TWO;
  const double tau = det(T);
  if (TMetric::invalid_determinant(tau)) { // barrier
    MSQ_SETERR(err)( barrier_violated_msg, MsqError::BARRIER_VIOLATED );
    return false;
  }
  const double inv_tau = 1.0/tau;
  const double invnorm = 1.0/norm;
  
  result = 0.5*inv_tau*(norm - invroot * trace(T));

  deriv_wrt_T = invnorm * T;
  pluseq_scaled_I( deriv_wrt_T, -invroot );
  deriv_wrt_T *= 0.5;
  deriv_wrt_T -= result * transpose_adj(T);
  deriv_wrt_T *= inv_tau;
  return true;
}

bool TShapeOrientB1::evaluate_with_hess( const MsqMatrix<2,2>& T,
                                         double& result,
                                         MsqMatrix<2,2>& deriv_wrt_T,
                                         MsqMatrix<2,2> second_wrt_T[3],
                                         MsqError& err )
{
  const double norm = Frobenius(T);
  const double invroot = 1.0/MSQ_SQRT_TWO;
  const double tau = det(T);
  if (TMetric::invalid_determinant(tau)) { // barrier
    MSQ_SETERR(err)( barrier_violated_msg, MsqError::BARRIER_VIOLATED );
    return false;
  }
  const double inv_tau = 1.0/tau;
  const double invnorm = 1.0/norm;
  
  const double f = norm - invroot * trace(T);
  result = 0.5 * inv_tau * f;

  const MsqMatrix<2,2> adjt = transpose_adj(T);
  deriv_wrt_T = invnorm * T;
  pluseq_scaled_I( deriv_wrt_T, -invroot );
  deriv_wrt_T *= 0.5;
  deriv_wrt_T -= result * adjt;
  deriv_wrt_T *= inv_tau;
  
  const double a = 0.5 * inv_tau * invnorm;
  set_scaled_outer_product( second_wrt_T, -a*invnorm*invnorm, T );
  pluseq_scaled_I( second_wrt_T, a );
  pluseq_scaled_outer_product( second_wrt_T, f*inv_tau*inv_tau*inv_tau, adjt );
  pluseq_scaled_2nd_deriv_of_det( second_wrt_T, -0.5*f*inv_tau*inv_tau, T );
  pluseq_scaled_sum_outer_product( second_wrt_T, -0.5*inv_tau*inv_tau*invnorm, T, adjt );
  pluseq_scaled_sum_outer_product_I( second_wrt_T, 0.5*inv_tau*inv_tau*invroot, adjt );
  return true;
}

bool TShapeOrientB1::evaluate( const MsqMatrix<3,3>& T, 
                               double& result, 
                               MsqError& err )
{
  const double tau = det(T);
  if (TMetric::invalid_determinant(tau)) { // barrier
    MSQ_SETERR(err)( barrier_violated_msg, MsqError::BARRIER_VIOLATED );
    return false;
  }
  result = 0.5/tau * (Frobenius( T ) - trace(T)/MSQ_SQRT_THREE);
  return true;
}


  bool TShapeOrientB1::evaluate_with_grad( const MsqMatrix<3,3>& T,
                                           double& result,
                                           MsqMatrix<3,3>& deriv_wrt_T,
                                           MsqError& err )
{
  const double norm = Frobenius(T);
  const double invroot = 1.0/MSQ_SQRT_THREE;
  const double tau = det(T);
  if (TMetric::invalid_determinant(tau)) { // barrier
    MSQ_SETERR(err)( barrier_violated_msg, MsqError::BARRIER_VIOLATED );
    return false;
  }
  const double inv_tau = 1.0/tau;
  const double invnorm = 1.0/norm;
  
  result = 0.5*inv_tau*(norm - invroot * trace(T));

  deriv_wrt_T = invnorm * T;
  pluseq_scaled_I( deriv_wrt_T, -invroot );
  deriv_wrt_T *= 0.5;
  deriv_wrt_T -= result * transpose_adj(T);
  deriv_wrt_T *= inv_tau;
  return true;
}

  bool TShapeOrientB1::evaluate_with_hess( const MsqMatrix<3,3>& T,
                                           double& result,
                                           MsqMatrix<3,3>& deriv_wrt_T,
                                           MsqMatrix<3,3> second_wrt_T[6],
                                           MsqError& err )
{
  const double norm = Frobenius(T);
  const double invroot = 1.0/MSQ_SQRT_THREE;
  const double tau = det(T);
  if (TMetric::invalid_determinant(tau)) { // barrier
    MSQ_SETERR(err)( barrier_violated_msg, MsqError::BARRIER_VIOLATED );
    return false;
  }
  const double inv_tau = 1.0/tau;
  const double invnorm = 1.0/norm;
  
  const double f = norm - invroot * trace(T);
  result = 0.5 * inv_tau * f;

  const MsqMatrix<3,3> adjt = transpose_adj(T);
  deriv_wrt_T = invnorm * T;
  pluseq_scaled_I( deriv_wrt_T, -invroot );
  deriv_wrt_T *= 0.5;
  deriv_wrt_T -= result * adjt;
  deriv_wrt_T *= inv_tau;
  
  const double a = 0.5 * inv_tau * invnorm;
  set_scaled_outer_product( second_wrt_T, -a*invnorm*invnorm, T );
  pluseq_scaled_I( second_wrt_T, a );
  pluseq_scaled_outer_product( second_wrt_T, f*inv_tau*inv_tau*inv_tau, adjt );
  pluseq_scaled_2nd_deriv_of_det( second_wrt_T, -0.5*f*inv_tau*inv_tau, T );
  pluseq_scaled_sum_outer_product( second_wrt_T, -0.5*inv_tau*inv_tau*invnorm, T, adjt );
  pluseq_scaled_sum_outer_product_I( second_wrt_T, 0.5*inv_tau*inv_tau*invroot, adjt );
  return true;
}

} // namespace Mesquite
