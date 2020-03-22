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


/** \file TShapeSize2DNB2.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "TShapeSize2DNB2.hpp"
#include "MsqMatrix.hpp"
#include "TMPDerivs.hpp"

namespace MESQUITE_NS {

std::string TShapeSize2DNB2::get_name() const
  { return "TShapeSize2DNB2"; }

TShapeSize2DNB2::~TShapeSize2DNB2() {}

/** \f$ \mu(T) = \frac{|T|^2+2}{2\psi(T)} - 1 \f$
 *  \f$ \psi(T) = \sqrt{|T|^2 + 2 \tau}\f$
 *  \f$ \tau = det(T) \f$
 */
bool TShapeSize2DNB2::evaluate( const MsqMatrix<2,2>& T, 
                                double& result, 
                                MsqError&  )
{
  double frob_sqr = sqr_Frobenius(T);
  double psi = sqrt( frob_sqr + 2.0*det(T) );

  double a = 1e-12;
  while (!Mesquite::divide(frob_sqr+2,2*psi,result)) {
    MsqMatrix<2,2> Tdelta(T);
    Tdelta(0,0) += a;
    Tdelta(1,1) += a;
    a *= 2.0;
    frob_sqr = sqr_Frobenius(Tdelta);
    psi = sqrt( frob_sqr + 2.0*det(Tdelta) );
    if (psi > 1e-50) 
     result = (frob_sqr+2) / 2*psi;
  }

  result -= 1.0;
  return true;
}


bool TShapeSize2DNB2::evaluate_with_grad( const MsqMatrix<2,2>& T, 
                                          double& result, 
                                          MsqMatrix<2,2>& deriv_wrt_T,
                                          MsqError& err )
{
  double frob_sqr = sqr_Frobenius(T);
  double psi = sqrt( frob_sqr + 2.0*det(T) );

  double a = 1e-12;
  while (!Mesquite::divide(frob_sqr+2,2*psi,result)) {
    MsqMatrix<2,2> Tdelta(T);
    Tdelta(0,0) += a;
    Tdelta(1,1) += a;
    a *= 2.0;
    frob_sqr = sqr_Frobenius(Tdelta);
    psi = sqrt( frob_sqr + 2.0*det(Tdelta) );
    if (psi > 1e-50) 
      result = (frob_sqr+2) / 2*psi;
  }
  
  //MsqMatrix<2,2> d_psi = 1.0/psi * (T + transpose_adj(T));
  //deriv_wrt_T = (1.0/psi) * (T - result * d_psi);
  if (psi > 1e-50) 
  {
    const double f = result/(psi*psi);
    deriv_wrt_T = transpose_adj(T);
    deriv_wrt_T *= -f;
    deriv_wrt_T += (1.0/psi - f) * T;
  }
  
  result -= 1.0;
  
  return true;
}


bool TShapeSize2DNB2::evaluate_with_hess( const MsqMatrix<2,2>& T, 
                                          double& result, 
                                          MsqMatrix<2,2>& deriv_wrt_T,
                                          MsqMatrix<2,2> second[3],
                                          MsqError& err )
{
  double frob_sqr = sqr_Frobenius(T);
  double psi = sqrt( frob_sqr + 2.0*det(T) );

  double a = 1e-12;
  while (!Mesquite::divide(frob_sqr+2,2*psi,result)) {
    MsqMatrix<2,2> Tdelta(T);
    Tdelta(0,0) += a;
    Tdelta(1,1) += a;
    a *= 2.0;
    frob_sqr = sqr_Frobenius(Tdelta);
    psi = sqrt( frob_sqr + 2.0*det(Tdelta) );
    if (psi > 1e-50) 
      result = (frob_sqr+2) / 2*psi;
  }
  
  const double inv_psi = 1.0/psi;
  MsqMatrix<2,2> d_psi = T + transpose_adj(T);
  d_psi *= inv_psi;

  deriv_wrt_T = d_psi;
  deriv_wrt_T *= -result;
  deriv_wrt_T += T;
  deriv_wrt_T *= inv_psi; 
  
  set_scaled_2nd_deriv_wrt_psi( second, -result*inv_psi, psi, T );
  pluseq_scaled_outer_product( second,  2*result*inv_psi*inv_psi, d_psi );
  pluseq_scaled_sum_outer_product( second, -inv_psi*inv_psi, d_psi, T );
  pluseq_scaled_I( second, inv_psi );
  
  result -= 1.0;
  
  return true;
}


} // namespace Mesquite
