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


/** \file TShape3DB2.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "TShape3DB2.hpp"
#include "MsqMatrix.hpp"
#include "MsqError.hpp"
#include "TMPDerivs.hpp"

namespace MESQUITE_NS {

TShape3DB2::~TShape3DB2() {}

std::string TShape3DB2::get_name() const
  { return "TShape3DB2"; }

// \mu_3(T) = \frac{ |T|^2 |adj(T)|^2 } {9 \tau^2} - 1
bool TShape3DB2::evaluate( const MsqMatrix<3,3>& T, 
                           double& result, 
                           MsqError& err )
{
  double f = sqr_Frobenius(T);
  double g = sqr_Frobenius(adj(T));
  double d = det(T);
  if (invalid_determinant(d)) {
    MSQ_SETERR(err)( barrier_violated_msg, MsqError::BARRIER_VIOLATED );
    return false;
  }
  result = (f*g) / (9*d*d) - 1;
  return true;
}


bool TShape3DB2::evaluate_with_grad( const MsqMatrix<3,3>& T, 
                                     double& result, 
                                     MsqMatrix<3,3>& wrt_T,
                                     MsqError& err )
{
  double f = sqr_Frobenius(T);
  double g = sqr_Frobenius(adj(T));
  double d = det(T);
  if (invalid_determinant(d)) {
    MSQ_SETERR(err)( barrier_violated_msg, MsqError::BARRIER_VIOLATED );
    return false;
  }
  result = (f*g) / (9*d*d) - 1;
  
  wrt_T = T;
  wrt_T *= (g + f*f);
  wrt_T -= f * (T * transpose(T) * T);
  wrt_T -= f * g / d * transpose_adj(T);
  wrt_T *= 2 / (9*d*d);
  
  return true;
}


bool TShape3DB2::evaluate_with_hess( const MsqMatrix<3,3>& T,
                                     double& result,
                                     MsqMatrix<3,3>& wrt_T,
                                     MsqMatrix<3,3> second[6],
                                     MsqError& err )
{
  double f = sqr_Frobenius(T);
  double g = sqr_Frobenius(adj(T));
  double d = det(T);
  if (invalid_determinant(d)) {
    MSQ_SETERR(err)( barrier_violated_msg, MsqError::BARRIER_VIOLATED );
    return false;
  }
  const double den = 1.0/(9*d*d);
  result = f*g*den- 1;
  
  MsqMatrix<3,3> dg = 2 * (f * T - T * transpose(T) * T);
  MsqMatrix<3,3> df = 2 * T;
  MsqMatrix<3,3> dtau = transpose_adj(T);
  
  wrt_T = g*df + f*dg - 2*f*g/d * transpose_adj(T);
  wrt_T *= den;
  
  set_scaled_2nd_deriv_norm_sqr_adj( second, den*f, T );
  pluseq_scaled_I( second, 2*den*g );
  pluseq_scaled_sum_outer_product( second, den, dg, df );
  pluseq_scaled_sum_outer_product( second, -2*den*g/d, df, dtau );
  pluseq_scaled_sum_outer_product( second, -2*den*f/d, dg, dtau );
  pluseq_scaled_outer_product( second, 6*den*f*g/(d*d), dtau );
  pluseq_scaled_2nd_deriv_of_det( second, -2*den*f*g/d, T );
  
  return true;
}


} // namespace Mesquite
