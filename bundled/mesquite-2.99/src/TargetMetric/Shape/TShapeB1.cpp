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


/** \file TShapeB1.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "TShapeB1.hpp"
#include "MsqMatrix.hpp"
#include "MsqError.hpp"
#include "TMPDerivs.hpp"

namespace MESQUITE_NS {

std::string TShapeB1::get_name() const
  { return "TShapeB1"; }

TShapeB1::~TShapeB1() {}

bool TShapeB1::evaluate( const MsqMatrix<2,2>& T, 
                         double& result, 
                         MsqError& err)
{
  const double d = det(T);
  if (invalid_determinant(d)) { // barrier
    MSQ_SETERR(err)( barrier_violated_msg, MsqError::BARRIER_VIOLATED );
    return false;
  }
    
  result = 0.5 * sqr_Frobenius(T) / d - 1;
  return true;
}

bool TShapeB1::evaluate_with_grad( const MsqMatrix<2,2>& T,
                                   double& result,
                                   MsqMatrix<2,2>& deriv_wrt_T,
                                   MsqError& err )
{
  const double d = det(T);
  if (invalid_determinant(d)) { // barrier
    MSQ_SETERR(err)( barrier_violated_msg, MsqError::BARRIER_VIOLATED );
    return false;
  }
  
  double inv_d = 1.0/d;
  result = 0.5 * sqr_Frobenius(T) * inv_d;
  deriv_wrt_T = T;
  deriv_wrt_T -= result * transpose_adj(T);
  deriv_wrt_T *= inv_d;

  result -= 1.0;
  return true;
}

/** \f$ \frac{\partial^2 \mu}{\partial T^2} 
      = \frac{1}{\tau} I_4 
      - \frac{1}{\tau^2} \left( T \otimes \frac{\partial \tau}{\partial T} 
                          + \frac{\partial \tau}{\partial T} \otimes T \right) 
      + \frac{|T|^2}{\tau^3} \left( \frac{\partial \tau}{\partial T} \otimes
                               \frac{\partial \tau}{\partial T} \right) 
      - \frac{|T|^2}{2 \tau^3} \frac{\partial^2 \tau}{\partial T^2} \f$
  */
bool TShapeB1::evaluate_with_hess( const MsqMatrix<2,2>& T,
                                   double& result,
                                   MsqMatrix<2,2>& deriv_wrt_T,
                                   MsqMatrix<2,2> second_wrt_T[3],
                                   MsqError& err )
{
  const double d = det(T);
  if (invalid_determinant(d)) { // barrier
    MSQ_SETERR(err)( barrier_violated_msg, MsqError::BARRIER_VIOLATED );
    return false;
  }
    
  double inv_d = 1.0/d;
  double f1 = sqr_Frobenius(T) * inv_d;
  result = 0.5 * f1;
  const MsqMatrix<2,2> adjt = transpose_adj(T);
  deriv_wrt_T = T;
  deriv_wrt_T -= result * adjt;
  deriv_wrt_T *= inv_d;
  
  set_scaled_outer_product( second_wrt_T, f1 * inv_d * inv_d, adjt );
  pluseq_scaled_sum_outer_product( second_wrt_T, -inv_d*inv_d, T, adjt );
  pluseq_scaled_I( second_wrt_T, inv_d );
  pluseq_scaled_2nd_deriv_of_det( second_wrt_T, -result * inv_d );

  result -= 1.0;
  return true;
}


bool TShapeB1::evaluate( const MsqMatrix<3,3>& T, 
                         double& result, 
                         MsqError& err)
{
  double f = Frobenius(T);
  double d = det(T);
  double den = 3 * MSQ_SQRT_THREE * d;
  
  if (invalid_determinant(d)) {
    MSQ_SETERR(err)( barrier_violated_msg, MsqError::BARRIER_VIOLATED );
    return false;
  }
  result = (f*f*f)/den - 1.0;
  return true;
}


bool TShapeB1::evaluate_with_grad( const MsqMatrix<3,3>& T, 
                                   double& result, 
                                   MsqMatrix<3,3>& wrt_T,
                                   MsqError& err )
{
  double d = det(T);
  if (invalid_determinant(d)) {
    MSQ_SETERR(err)( barrier_violated_msg, MsqError::BARRIER_VIOLATED );
    return false;
  }
    
  double norm = Frobenius(T);
  double den = 1.0/(3 * MSQ_SQRT_THREE * d);
  double norm_cube = norm*norm*norm;
  result = norm_cube * den - 1.0;
  wrt_T = T;
  wrt_T *= 3 * norm * den;
  wrt_T -= norm_cube * den/d * transpose_adj(T);
   return true;
}

bool TShapeB1::evaluate_with_hess( const MsqMatrix<3,3>& T,
                                   double& result,
                                   MsqMatrix<3,3>& deriv_wrt_T,
                                   MsqMatrix<3,3> second_wrt_T[6],
                                   MsqError& err )
{
  double d = det(T);
  if (invalid_determinant(d)) {
    MSQ_SETERR(err)( barrier_violated_msg, MsqError::BARRIER_VIOLATED );
    return false;
  }
  
  double id = 1.0/d;
  double norm = Frobenius(T);
  double den = 1.0/(3 * MSQ_SQRT_THREE * d);
  double norm_cube = norm*norm*norm;
  result = norm_cube * den - 1.0;
  MsqMatrix<3,3> adjt = transpose_adj(T);
  deriv_wrt_T = T;
  deriv_wrt_T *= 3 * norm * den;
  deriv_wrt_T -= norm_cube * den * id * transpose_adj(T);
 
  set_scaled_outer_product( second_wrt_T, 3 * den / norm, T );
  pluseq_scaled_I( second_wrt_T, 3 * norm * den );
  pluseq_scaled_2nd_deriv_of_det( second_wrt_T, -den * norm_cube * id, T );
  pluseq_scaled_outer_product( second_wrt_T, 2 * den * norm_cube * id * id , adjt );
  pluseq_scaled_sum_outer_product( second_wrt_T, -3 * norm * den * id, T, adjt );

  return true;
}

} // namespace Mesquite
