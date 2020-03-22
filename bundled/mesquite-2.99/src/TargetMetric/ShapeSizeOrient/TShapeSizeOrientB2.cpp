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
   
  ***************************************************************** */


/** \file TRel2DShapeSizeOrientBarrierAlt1.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "TShapeSizeOrientB2.hpp"
#include "MsqMatrix.hpp"
#include "MsqError.hpp"
#include "TMPDerivs.hpp"

namespace MESQUITE_NS {

std::string TShapeSizeOrientB2::get_name() const
  { return "TShapeSizeOrientB2"; }

TShapeSizeOrientB2::~TShapeSizeOrientB2() {}

bool TShapeSizeOrientB2::evaluate( const MsqMatrix<2,2>& T, 
                                   double& result,
                                   MsqError& err )
{
  double d = det(T);
  if (TMetric::invalid_determinant(d)) {
    MSQ_SETERR(err)( barrier_violated_msg, MsqError::BARRIER_VIOLATED );
    return false;
  }
  MsqMatrix<2,2> T_inv = 1/d * adj(T);
  pluseq_scaled_I( T_inv, -1.0 );
  result = sqr_Frobenius(T_inv);
  return true;
}

bool TShapeSizeOrientB2::evaluate( const MsqMatrix<3,3>& T, 
                                   double& result,
                                   MsqError& err)
{
  double d = det(T);
  if (TMetric::invalid_determinant(d)) {
    MSQ_SETERR(err)( barrier_violated_msg, MsqError::BARRIER_VIOLATED );
    return false;
  }
  MsqMatrix<3,3> T_inv = 1/d * adj(T);
  pluseq_scaled_I( T_inv, -1.0 );
  result = sqr_Frobenius(T_inv);
  return true;
}

/** \f$ \frac{1}{\tau^2}|T|^2 - \frac{2}{\tau}tr(adj T) + 2 \f$ */
bool TShapeSizeOrientB2::evaluate_with_grad( const MsqMatrix<2,2>& T,
                                             double& result,
                                             MsqMatrix<2,2>& deriv_wrt_T,
                                             MsqError& err )
{
  const double tau = det(T);
  if (invalid_determinant(tau)) {
    MSQ_SETERR(err)( barrier_violated_msg, MsqError::BARRIER_VIOLATED );
    return false;
  }
  
  const MsqMatrix<2,2> adjt = transpose_adj(T);
  const double it = 1.0/tau;
  result = it*(it*sqr_Frobenius(T) - 2.0*trace(T)) + 2.0;
  deriv_wrt_T = T;
  deriv_wrt_T *= it*it;
  deriv_wrt_T(0,0) -= it;
  deriv_wrt_T(1,1) -= it;
  deriv_wrt_T += it*it*(trace(T)-it*sqr_Frobenius(T))*adjt;
  deriv_wrt_T *= 2.0;
  return true;
}

/** \f$ \frac{1}{\tau^2}|adj T|^2 - \frac{2}{\tau}tr(adj T) + 3 \f$ */
bool TShapeSizeOrientB2::evaluate_with_grad( const MsqMatrix<3,3>& T,
                                             double& result,
                                             MsqMatrix<3,3>& deriv_wrt_T,
                                             MsqError& err )
{
  const double tau = det(T);
  if (invalid_determinant(tau)) {
    MSQ_SETERR(err)( barrier_violated_msg, MsqError::BARRIER_VIOLATED );
    return false;
  }
  
  const MsqMatrix<3,3> adjt = adj(T);
  const double it = 1.0/tau;
  result = it*(it*sqr_Frobenius(adjt) - 2.0*trace(adjt)) + 3.0;
  
  deriv_wrt_T = T;
  deriv_wrt_T *= sqr_Frobenius(T);
  deriv_wrt_T -= T * transpose(T) * T;
  deriv_wrt_T *= it*it;
  
  deriv_wrt_T += it*it*(trace(adjt)-it*sqr_Frobenius(adjt))*transpose(adjt);

  double f = trace(T) * it;
  deriv_wrt_T(0,0) -= f;
  deriv_wrt_T(1,1) -= f;
  deriv_wrt_T(2,2) -= f;
  
  deriv_wrt_T += it*transpose(T);

  deriv_wrt_T *= 2.0;
  return true;
}

bool TShapeSizeOrientB2::evaluate_with_hess( const MsqMatrix<2,2>& T,
                                             double& result,
                                             MsqMatrix<2,2>& deriv_wrt_T,
                                             MsqMatrix<2,2> second[3],
                                             MsqError& err )
{
  const double tau = det(T);
  if (invalid_determinant(tau)) {
    MSQ_SETERR(err)( barrier_violated_msg, MsqError::BARRIER_VIOLATED );
    return false;
  }
  
  const MsqMatrix<2,2> adjt = transpose_adj(T);
  const double it = 1.0/tau;
  result = it*(it*sqr_Frobenius(T) - 2.0*trace(T)) + 2.0;
  deriv_wrt_T = T;
  deriv_wrt_T *= it*it;
  deriv_wrt_T(0,0) -= it;
  deriv_wrt_T(1,1) -= it;
  deriv_wrt_T += it*it*(trace(T)-it*sqr_Frobenius(T))*adjt;
  deriv_wrt_T *= 2.0;
  
  set_scaled_outer_product( second, it*it*it*(6*it*sqr_Frobenius(T) - 4*trace(T)), adjt );
  pluseq_scaled_I( second, 2*it*it );
  pluseq_scaled_2nd_deriv_of_det( second, 2*it*it*(trace(T) - it*sqr_Frobenius(T)) );
  pluseq_scaled_sum_outer_product( second, -4*it*it*it, T, adjt );
  pluseq_scaled_sum_outer_product_I( second, 2*it*it, adjt );
  
  return true;
}

bool TShapeSizeOrientB2::evaluate_with_hess( const MsqMatrix<3,3>& T,
                                             double& result,
                                             MsqMatrix<3,3>& deriv_wrt_T,
                                             MsqMatrix<3,3> second[6],
                                             MsqError& err )
{
  const double tau = det(T);
  if (invalid_determinant(tau)) {
    MSQ_SETERR(err)( barrier_violated_msg, MsqError::BARRIER_VIOLATED );
    return false;
  }
  
  const MsqMatrix<3,3> adjt = adj(T);
  const double it = 1.0/tau;
  const double nadjt = sqr_Frobenius(adjt);
  const double nT = sqr_Frobenius(T);
  const double tadjT = trace(adjt);
  result = it*(it*nadjt - 2.0*tadjT) + 3.0;
  
  const MsqMatrix<3,3> TTtT = T * transpose(T) * T;
  deriv_wrt_T = T;
  deriv_wrt_T *= nT;
  deriv_wrt_T -= TTtT;
  deriv_wrt_T *= it*it;
 
  deriv_wrt_T += it*it*(tadjT-it*nadjt)*transpose(adjt);

  const double tT = trace(T);
  double f = tT * it;
  deriv_wrt_T(0,0) -= f;
  deriv_wrt_T(1,1) -= f;
  deriv_wrt_T(2,2) -= f;
  
  deriv_wrt_T += it*transpose(T);

  deriv_wrt_T *= 2.0;
  
  set_scaled_2nd_deriv_norm_sqr_adj( second, it*it, T );

  const double yf = -it*it*it*it;
  const double sf = -2;
  const double zf = -it*it*sf;

  pluseq_scaled_2nd_deriv_of_det( second, yf*2*nadjt*tau + zf*tadjT, T );
  pluseq_scaled_outer_product( second, yf*-6*nadjt - 2*zf*tadjT*it, transpose(adjt) );
  MsqMatrix<3,3> dnadj_dT = 2 * (nT * T - TTtT);
  pluseq_scaled_sum_outer_product( second, yf * 2 * tau, dnadj_dT, transpose(adjt) );
  pluseq_scaled_2nd_deriv_tr_adj( second, sf * it );
  MsqMatrix<3,3> dtradj_dT = -transpose(T);
  dtradj_dT(0,0) += tT;
  dtradj_dT(1,1) += tT;
  dtradj_dT(2,2) += tT;
  pluseq_scaled_sum_outer_product( second, zf, dtradj_dT, transpose(adjt) );

  return true;
}


} // namespace Mesquite
