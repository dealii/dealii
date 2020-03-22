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


/** \file TShapeSizeB1.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "TShapeSizeB1.hpp"
#include "TMPDerivs.hpp"
#include "MsqError.hpp"

namespace MESQUITE_NS {

std::string TShapeSizeB1::get_name() const
  { return "TShapeSizeB1"; }

TShapeSizeB1::~TShapeSizeB1() {}

bool TShapeSizeB1::evaluate( const MsqMatrix<2,2>& T, 
                             double& result, 
                             MsqError& err )
{
  const double tau = det(T);
  if (invalid_determinant(tau)) { // barrier
    return false;
  }
  
  const double nT = sqr_Frobenius(T);
  const double f = 1/(tau*tau);
  result = (1 + f) * nT - 4;
  return true;
}

bool TShapeSizeB1::evaluate( const MsqMatrix<3,3>& T, 
                             double& result, 
                             MsqError& err )
{
  const double tau = det(T);
  if (invalid_determinant(tau)) { // barrier
    return false;
  }
  
  const double nT = sqr_Frobenius(T);
  const double nadj = sqr_Frobenius(transpose_adj(T));
  const double f = 1/(tau*tau);
  result = nT + f*nadj - 6;
  return true;
}

bool TShapeSizeB1::evaluate_with_grad( const MsqMatrix<2,2>& T,
                                       double& result,
                                       MsqMatrix<2,2>& deriv_wrt_T,
                                       MsqError& err )
{
  const double tau = det(T);
  if (invalid_determinant(tau)) { // barrier
    MSQ_SETERR(err)( barrier_violated_msg, MsqError::BARRIER_VIOLATED );
    return false;
  }
  
  const MsqMatrix<2,2> adjt = transpose_adj(T);
  const double nT = sqr_Frobenius(T);
  const double f = 1/(tau*tau);
  result = (1 + f) * nT - 4;
  
  deriv_wrt_T = T;
  deriv_wrt_T *= 2 + 2*f;
  deriv_wrt_T -= 2 * f/tau * nT * adjt;
  
  return true;
}

bool TShapeSizeB1::evaluate_with_grad( const MsqMatrix<3,3>& T,
                                       double& result,
                                       MsqMatrix<3,3>& deriv_wrt_T,
                                       MsqError& err )
{
  const double tau = det(T);
  if (invalid_determinant(tau)) { // barrier
    MSQ_SETERR(err)( barrier_violated_msg, MsqError::BARRIER_VIOLATED );
    return false;
  }
  
  const MsqMatrix<3,3> adjt = transpose_adj(T);
  const double nT = sqr_Frobenius(T);
  const double nadj = sqr_Frobenius(adjt);
  const double f = 1/(tau*tau);
  result = nT + f*nadj - 6;
  
  deriv_wrt_T = T;
  deriv_wrt_T *= (1+f*nT);
  deriv_wrt_T -= f * T * transpose(T) * T;
  deriv_wrt_T -= f/tau * nadj * adjt;
  deriv_wrt_T *= 2;

  return true;
}

bool TShapeSizeB1::evaluate_with_hess( const MsqMatrix<2,2>& T,
                                       double& result,
                                       MsqMatrix<2,2>& deriv_wrt_T,
                                       MsqMatrix<2,2> second_wrt_T[3],
                                       MsqError& err )
{
  const double tau = det(T);
  if (invalid_determinant(tau)) { // barrier
    MSQ_SETERR(err)( barrier_violated_msg, MsqError::BARRIER_VIOLATED );
    return false;
  }
  
  const MsqMatrix<2,2> adjt = transpose_adj(T);
  const double nT = sqr_Frobenius(T);
  const double f = 1/(tau*tau);
  result = (1 + f) * nT - 4;
  
  deriv_wrt_T = T;
  deriv_wrt_T *= 2 + 2*f;
  deriv_wrt_T -= 2 * f/tau * nT * adjt;

  set_scaled_sum_outer_product( second_wrt_T, -4*f/tau, T, adjt );
  pluseq_scaled_I( second_wrt_T, 2 + 2*f );
  pluseq_scaled_outer_product( second_wrt_T, 6*nT*f*f, adjt );
  pluseq_scaled_2nd_deriv_of_det( second_wrt_T, -2*nT*f/tau );

  return true;
}
bool TShapeSizeB1::evaluate_with_hess( const MsqMatrix<3,3>& T,
                                       double& result,
                                       MsqMatrix<3,3>& deriv_wrt_T,
                                       MsqMatrix<3,3> second_wrt_T[6],
                                       MsqError& err )
{
  const double tau = det(T);
  if (invalid_determinant(tau)) { // barrier
    MSQ_SETERR(err)( barrier_violated_msg, MsqError::BARRIER_VIOLATED );
    return false;
  }
  
  const MsqMatrix<3,3> adjt = transpose_adj(T);
  const double nT = sqr_Frobenius(T);
  const double nadj = sqr_Frobenius(adjt);
  const double f = 1/(tau*tau);
  result = nT + f*nadj - 6;
  
  //! \f$ \frac{\partial}{\partial T} |adj T|^2 \f$
  const MsqMatrix<3,3> dNadj_dT = 2 * (nT * T - T * transpose(T) * T);
  deriv_wrt_T = T;
  deriv_wrt_T -= f/tau * nadj * adjt;
  deriv_wrt_T *= 2;
  deriv_wrt_T += f * dNadj_dT;
 
    // calculate negative of 2nd wrt T of (|adj T|^2 / tau^2) (sec 3.2.2)
  set_scaled_2nd_deriv_norm_sqr_adj( second_wrt_T,    f,            T );
  pluseq_scaled_2nd_deriv_of_det(    second_wrt_T, -2*f*f*nadj*tau, T );
  pluseq_scaled_outer_product(       second_wrt_T,  6*f*f*nadj,     adjt );
  pluseq_scaled_sum_outer_product(   second_wrt_T, -2*f*f     *tau, adjt, dNadj_dT );
    // calculate 2nd wrt T of this metric
  pluseq_scaled_I( second_wrt_T, 2.0 );

  return true;
}

} // namespace Mesquite
