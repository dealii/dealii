/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2009 Sandia National Laboratories.  Developed at the
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


/** \file TShapeSize3DB2.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "TShapeSize3DB2.hpp"
#include "TMPDerivs.hpp"
#include "MsqError.hpp"

namespace MESQUITE_NS {

std::string TShapeSize3DB2::get_name() const
  { return "TShapeSize3DB2"; }

TShapeSize3DB2::~TShapeSize3DB2() {}

bool TShapeSize3DB2::evaluate( const MsqMatrix<3,3>& T, 
                               double& result, 
                               MsqError& err )
{
  const double tau = det(T);
  if (invalid_determinant(tau)) { // barrier
    MSQ_SETERR(err)( barrier_violated_msg, MsqError::BARRIER_VIOLATED );
    return false;
  }
  
  const double f = sqr_Frobenius(T);
  const double g = sqr_Frobenius(adj(T));
  result = (f + g)/(6 * tau) - 1;
  return true;
}

bool TShapeSize3DB2::evaluate_with_grad( const MsqMatrix<3,3>& T,
                                         double& result,
                                         MsqMatrix<3,3>& deriv_wrt_T,
                                         MsqError& err )
{
  const double tau = det(T);
  if (invalid_determinant(tau)) { // barrier
    MSQ_SETERR(err)( barrier_violated_msg, MsqError::BARRIER_VIOLATED );
    return false;
  }
  
  const double f = sqr_Frobenius(T);
  const double g = sqr_Frobenius(adj(T));
  result = (f + g)/(6 * tau);
  
  deriv_wrt_T = -transpose(T) * T;
  deriv_wrt_T(0,0) += 1+f;
  deriv_wrt_T(1,1) += 1+f;
  deriv_wrt_T(2,2) += 1+f;
  deriv_wrt_T = T * deriv_wrt_T;
  deriv_wrt_T -= 3*result * transpose_adj(T);
  deriv_wrt_T *= 1.0/(3*tau);
  
  result -= 1.0;
  return true;
}


bool TShapeSize3DB2::evaluate_with_hess( const MsqMatrix<3,3>& T,
                                         double& result,
                                         MsqMatrix<3,3>& wrt_T,
                                         MsqMatrix<3,3> second[6],
                                         MsqError& err )
{
  const double tau = det(T);
  if (invalid_determinant(tau)) { // barrier
    MSQ_SETERR(err)( barrier_violated_msg, MsqError::BARRIER_VIOLATED );
    return false;
  }
  
  const double f = sqr_Frobenius(T);
  const double g = sqr_Frobenius(adj(T));
  result = (f + g)/(6 * tau);
  
  MsqMatrix<3,3> dtau = transpose_adj(T);
  MsqMatrix<3,3> dg = -transpose(T) * T;
  dg(0,0) += f;
  dg(1,1) += f;
  dg(2,2) += f;
  dg = T * dg;
  dg *= 2;
  
  wrt_T = T;
  wrt_T += 0.5*dg;
  wrt_T *= 1.0/3.0;
  wrt_T -= result * dtau;
  wrt_T *= 1.0/tau;
  
  set_scaled_2nd_deriv_norm_sqr_adj( second, 1.0/6.0, T );
  pluseq_scaled_I( second, 1.0/3.0 );
  pluseq_scaled_sum_outer_product( second, -1./3./tau, T, dtau );
  pluseq_scaled_sum_outer_product( second, -1./6./tau, dg, dtau );
  pluseq_scaled_outer_product( second, 2*result/tau, dtau );
  pluseq_scaled_2nd_deriv_of_det( second, -result, T );
  hess_scale( second, 1.0/tau );
  
  result -= 1.0;
  return true;
}

} // namespace Mesquite
