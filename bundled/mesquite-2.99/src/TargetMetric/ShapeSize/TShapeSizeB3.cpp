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


/** \file TShapeSizeB3.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "TShapeSizeB3.hpp"
#include "TMPDerivs.hpp"
#include "MsqError.hpp"

#include <iostream>

namespace MESQUITE_NS {

std::string TShapeSizeB3::get_name() const
  { return "TShapeSizeB3"; }

TShapeSizeB3::~TShapeSizeB3() {}

bool TShapeSizeB3::evaluate( const MsqMatrix<2,2>& T, 
                                           double& result, 
                                           MsqError& err )
{
  const double tau = det(T);
  if (invalid_determinant(tau)) { // barrier
    MSQ_SETERR(err)( barrier_violated_msg, MsqError::BARRIER_VIOLATED );
    return false;
  }
  
  result = sqr_Frobenius(T) - 2.0 * std::log(tau) - 2;
  return true;
}

bool TShapeSizeB3::evaluate_with_grad( const MsqMatrix<2,2>& T,
                                                     double& result,
                                                     MsqMatrix<2,2>& deriv_wrt_T,
                                                     MsqError& err )
{
  const double tau = det(T);
  if (invalid_determinant(tau)) { // barrier
    MSQ_SETERR(err)( barrier_violated_msg, MsqError::BARRIER_VIOLATED );
    return false;
  }
  
  result = sqr_Frobenius(T) - 2.0 * std::log(tau) - 2;
  deriv_wrt_T = T;
  deriv_wrt_T -= 1/tau * transpose_adj(T);
  deriv_wrt_T *= 2;
  
  return true;
}

bool TShapeSizeB3::evaluate_with_hess( const MsqMatrix<2,2>& T,
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
  
  result = sqr_Frobenius(T) - 2.0 * std::log(tau) - 2;

  const MsqMatrix<2,2> adjt = transpose_adj(T);
  const double it = 1/tau;
  deriv_wrt_T = T;
  deriv_wrt_T -= it * adjt;
  deriv_wrt_T *= 2;
  
  set_scaled_outer_product( second_wrt_T, 2*it*it, adjt );
  pluseq_scaled_2nd_deriv_of_det( second_wrt_T, -2*it );
  pluseq_scaled_I( second_wrt_T, 2.0 );

  return true;
}


bool TShapeSizeB3::evaluate( const MsqMatrix<3,3>& T, 
                             double& result, 
                             MsqError& err )
{
  const double tau = det(T);
  if (invalid_determinant(tau)) { // barrier
    MSQ_SETERR(err)( barrier_violated_msg, MsqError::BARRIER_VIOLATED );
    return false;
  }
  
  double n = Frobenius(T);
  result = n*n*n - 3*MSQ_SQRT_THREE*( log(tau) + 1 );
  return true;
}

bool TShapeSizeB3::evaluate_with_grad( const MsqMatrix<3,3>& T,
                                       double& result,
                                       MsqMatrix<3,3>& deriv_wrt_T,
                                       MsqError& err )
{
  const double tau = det(T);
  if (invalid_determinant(tau)) { // barrier
    MSQ_SETERR(err)( barrier_violated_msg, MsqError::BARRIER_VIOLATED );
    return false;
  }
  
  double n = Frobenius(T);
  result = n*n*n - 3*MSQ_SQRT_THREE*( log(tau) + 1 );
  
  const MsqMatrix<3,3> adjt = transpose_adj(T);
  deriv_wrt_T = T;
  deriv_wrt_T *= 3*n;
  deriv_wrt_T -= 3*MSQ_SQRT_THREE/tau * adjt;
  
  return true;
}

bool TShapeSizeB3::evaluate_with_hess( const MsqMatrix<3,3>& T,
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
  
  double n = Frobenius(T);
  result = n*n*n - 3*MSQ_SQRT_THREE*( log(tau) + 1 );
  
  const MsqMatrix<3,3> adjt = transpose_adj(T);
  const double it = 1/tau;
  deriv_wrt_T = T;
  deriv_wrt_T *= 3*n;
  deriv_wrt_T -= 3*MSQ_SQRT_THREE*it * adjt;

  if (n > 1e-50) 
  {
    set_scaled_outer_product( second_wrt_T, 3/n, T );
    pluseq_scaled_I( second_wrt_T, 3*n );
    pluseq_scaled_2nd_deriv_of_det( second_wrt_T, -3*MSQ_SQRT_THREE*it, T );
    pluseq_scaled_outer_product( second_wrt_T, 3*MSQ_SQRT_THREE*it*it, adjt );
  }
  else
  {
    std::cout << "Warning: Division by zero avoided in TShapeSizeB3::evaluate_with_hess()" << std::endl;
  }
  

  return true;
}

} // namespace Mesquite
