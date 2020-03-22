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


/** \file TShapeSize2DB2.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "TShapeSize2DB2.hpp"
#include "MsqMatrix.hpp"
#include "MsqError.hpp"
#include "TMPDerivs.hpp"

#include <iostream>

namespace MESQUITE_NS {

std::string TShapeSize2DB2::get_name() const
  { return "TShapeSize2DB2"; }

TShapeSize2DB2::~TShapeSize2DB2() {}

bool TShapeSize2DB2::evaluate( const MsqMatrix<2,2>& T, 
                               double& result, 
                               MsqError& err )
{
  const double two_det = 2.0 * det(T);
  if (invalid_determinant(two_det)) { // barrier
    MSQ_SETERR(err)( barrier_violated_msg, MsqError::BARRIER_VIOLATED );
    return false;
  }
    
  const double frob_sqr = sqr_Frobenius(T);
  result = (frob_sqr - 2.0 * sqrt( frob_sqr + two_det ) + 2.0)/two_det;
  return true;
}

bool TShapeSize2DB2::evaluate_with_grad( const MsqMatrix<2,2>& T, 
                                         double& result, 
                                         MsqMatrix<2,2>& deriv_wrt_T,
                                         MsqError& err )
{
  const double d = det(T);
  if (invalid_determinant(d)) { // barrier
    MSQ_SETERR(err)( barrier_violated_msg, MsqError::BARRIER_VIOLATED );
    return false;
  }
  const double frob_sqr = sqr_Frobenius(T);
  const double psi = sqrt( frob_sqr + 2.0*det(T) );
  const double v = frob_sqr - 2.0 * psi + 2.0;
  result = v / (2*d);

    // deriv of V wrt T
  MsqMatrix<2,2> adjt = transpose_adj(T);
  MsqMatrix<2,2> v_wrt_T(T);
  if (psi > 1e-50)
  {
    v_wrt_T *= (1.0 - 1.0/psi);
    v_wrt_T -= 1.0/psi * adjt;
    v_wrt_T *= 2;
  }
  else
  {
    std::cout << "Warning: Division by zero avoided in TShapeSize2DB2::evaluate_with_grad()" << std::endl;
  }
  
    // deriv of mu wrt T
  deriv_wrt_T = v_wrt_T;
  deriv_wrt_T *= 0.5/d;
  deriv_wrt_T -= v / (2*d*d) * adjt;
  
  return true;
}

bool TShapeSize2DB2::evaluate_with_hess( const MsqMatrix<2,2>& T, 
                                         double& result, 
                                         MsqMatrix<2,2>& deriv_wrt_T,
                                         MsqMatrix<2,2> second[3],
                                         MsqError& err )
{
  const double d = det(T);
  if (invalid_determinant(d)) { // barrier
    MSQ_SETERR(err)( barrier_violated_msg, MsqError::BARRIER_VIOLATED );
    return false;
  }
  const double frob_sqr = sqr_Frobenius(T);
  const double psi = sqrt( frob_sqr + 2.0*det(T) );
  const double v = frob_sqr - 2.0 * psi + 2.0;
  result = v / (2*d);

    // deriv of V wrt T
  MsqMatrix<2,2> adjt = transpose_adj(T);
  MsqMatrix<2,2> v_wrt_T(T);
  if (psi > 1e-50) 
  {
    v_wrt_T *= (1.0 - 1.0/psi);
    v_wrt_T -= 1.0/psi * adjt;
    v_wrt_T *= 2;
  }
  else
  {
    std::cout << "Warning: Division by zero avoided in TShapeSize2DB2::evaluate_with_hess()" << std::endl;
  }  
  
    // deriv of mu wrt T
  deriv_wrt_T = v_wrt_T;
  deriv_wrt_T *= 0.5/d;
  deriv_wrt_T -= v / (2*d*d) * adjt;
  
    // second of V wrt T 
  const double s = T(0,1) - T(1,0);
  const double tr = trace(T);
  const double f = -2.0/(psi*psi*psi);
  second[0](0,0) = second[1](0,1) = second[2](1,1) =  f*s*s;
  second[0](0,1) = second[0](1,0) = second[1](1,1) = -f*s*tr;
  second[1](0,0) = second[2](0,1) = second[2](1,0) =  f*s*tr;
  second[0](1,1) = second[2](0,0) = -(second[1](1,0) = -f*tr*tr);
  pluseq_scaled_I( second, 2 );
  
    // second of mu wrt T 
  const double x = 1.0/(2*d);
  second[0] *= x;
  second[1] *= x;
  second[2] *= x;
  pluseq_scaled_2nd_deriv_of_det( second, v/(-2*d*d) );
  pluseq_scaled_outer_product( second, v/(d*d*d), adjt );
  pluseq_scaled_sum_outer_product( second, -1/(2*d*d), v_wrt_T, adjt );
  
  return true;
}

} // namespace Mesquite
