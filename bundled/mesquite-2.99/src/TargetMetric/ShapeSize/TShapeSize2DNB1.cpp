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


/** \file TShapeSize2DNB1.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "TShapeSize2DNB1.hpp"
#include "MsqMatrix.hpp"
#include "TMPDerivs.hpp"

#include <iostream>

namespace MESQUITE_NS {

std::string TShapeSize2DNB1::get_name() const
  { return "TShapeSize2DNB1"; }

TShapeSize2DNB1::~TShapeSize2DNB1() {}

/** \f$ \mu(T) = |T|^2 - 2 \psi(T) + 2 \f$
 *  \f$ \psi(T) = \sqrt{|T|^2 + 2 \tau} \f$
 *  \f$ \tau = det(T) \f$
 */
bool TShapeSize2DNB1::evaluate( const MsqMatrix<2,2>& T, 
                                double& result, 
                                MsqError&  )
{
  double frob_sqr = sqr_Frobenius(T);
  double psi = sqrt( frob_sqr + 2.0*det(T) );

  MsqMatrix<2,2> Tdelta(T);
  while (fabs(psi) < DBL_EPSILON) {
    Tdelta(0,0) += 1e-12;
    Tdelta(1,1) += 1e-12;
    frob_sqr = sqr_Frobenius(Tdelta);
    psi = sqrt( frob_sqr + 2.0*det(Tdelta) );
  }

  result = frob_sqr - 2.0 * psi + 2.0;
  return true;
}


bool TShapeSize2DNB1::evaluate_with_grad( const MsqMatrix<2,2>& T, 
                                          double& result, 
                                          MsqMatrix<2,2>& deriv_wrt_T,
                                          MsqError& err )
{
  double frob_sqr = sqr_Frobenius(T);
  double psi = sqrt( frob_sqr + 2.0*det(T) );

  MsqMatrix<2,2> Tdelta(T);
  while (fabs(psi) < DBL_EPSILON) {
    Tdelta(0,0) += 1e-12;
    Tdelta(1,1) += 1e-12;
    frob_sqr = sqr_Frobenius(Tdelta);
    psi = sqrt( frob_sqr + 2.0*det(Tdelta) );
  }

  result = frob_sqr - 2.0 * psi + 2.0;

  deriv_wrt_T = T;
  if (psi > 1e-50)
  {
    deriv_wrt_T *= (1.0 - 1.0/psi);
    deriv_wrt_T -= 1.0/psi * transpose_adj(T);
    deriv_wrt_T *= 2;
  }
  else
  {
    std::cout << "Warning: Division by zero avoided in TShapeSize2DNB2::evaluate_with_grad()" << std::endl;
  }

  return true;
}

bool TShapeSize2DNB1::evaluate_with_hess( const MsqMatrix<2,2>& T, 
                                          double& result, 
                                          MsqMatrix<2,2>& deriv_wrt_T,
                                          MsqMatrix<2,2> second[3],
                                          MsqError& err )
{
  double frob_sqr = sqr_Frobenius(T);
  double psi = sqrt( frob_sqr + 2.0*det(T) );

  MsqMatrix<2,2> Tdelta(T);
  while (fabs(psi) < DBL_EPSILON) {
    Tdelta(0,0) += 1e-12;
    Tdelta(1,1) += 1e-12;
    frob_sqr = sqr_Frobenius(Tdelta);
    psi = sqrt( frob_sqr + 2.0*det(Tdelta) );
  }

  result = frob_sqr - 2.0 * psi + 2.0;

  deriv_wrt_T = T;
  if (psi > 1e-50)
  {
    deriv_wrt_T *= (1.0 - 1.0/psi);
    deriv_wrt_T -= 1.0/psi * transpose_adj(T);
    deriv_wrt_T *= 2;
  }
  else
  {
    std::cout << "Warning: Division by zero avoided in TShapeSize2DNB2::evaluate_with_hess()" << std::endl;
  }

  set_scaled_2nd_deriv_wrt_psi( second, -2.0, psi, T );
  pluseq_scaled_I( second, 2 );
  
  return true;
}

} // namespace Mesquite
