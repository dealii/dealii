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


/** \file TShapeNB1.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "TShapeNB1.hpp"
#include "MsqMatrix.hpp"
#include "TMPDerivs.hpp"

#include <iostream>

namespace MESQUITE_NS {

std::string TShapeNB1::get_name() const
  { return "TShapeNB1"; }

TShapeNB1::~TShapeNB1() {}

bool TShapeNB1::evaluate( const MsqMatrix<2,2>& T, 
                          double& result, 
                          MsqError&  )
{
  result = sqr_Frobenius(T) - 2.0*det(T);
  return true;  
}

bool TShapeNB1::evaluate_with_grad( const MsqMatrix<2,2>& T,
                                    double& result,
                                    MsqMatrix<2,2>& deriv_wrt_T,
                                    MsqError& err )
{
  result = sqr_Frobenius(T) - 2.0*det(T);
  deriv_wrt_T = T;
  deriv_wrt_T -= transpose_adj(T);
  deriv_wrt_T *= 2;
  return true;  
  
}

bool TShapeNB1::evaluate_with_hess( const MsqMatrix<2,2>& T,
                                    double& result,
                                    MsqMatrix<2,2>& deriv_wrt_T,
                                    MsqMatrix<2,2> second_wrt_T[3],
                                    MsqError& err )
{
  result = sqr_Frobenius(T) - 2.0*det(T);
  deriv_wrt_T = T;
  deriv_wrt_T -= transpose_adj(T);
  deriv_wrt_T *= 2;
  set_scaled_I( second_wrt_T, 2.0 );
  pluseq_scaled_2nd_deriv_of_det( second_wrt_T, -2.0 );
  return true;
}


bool TShapeNB1::evaluate( const MsqMatrix<3,3>& T, 
                          double& result, 
                          MsqError& )
{
  double f = Frobenius(T);
  double d = det(T);
  result = f*f*f - 3*MSQ_SQRT_THREE*d;
  return true;
}


bool TShapeNB1::evaluate_with_grad( const MsqMatrix<3,3>& T, 
                                    double& result, 
                                    MsqMatrix<3,3>& deriv_wrt_T,
                                    MsqError& err )
{
  double f = Frobenius(T);
  double d = det(T);
  result = f*f*f - 3*MSQ_SQRT_THREE*d;

  deriv_wrt_T = T;
  deriv_wrt_T *= f;
  deriv_wrt_T -= MSQ_SQRT_THREE*transpose_adj(T);
  deriv_wrt_T *= 3;
  return true;
}

bool TShapeNB1::evaluate_with_hess( const MsqMatrix<3,3>& T, 
                                    double& result, 
                                    MsqMatrix<3,3>& deriv_wrt_T,
                                    MsqMatrix<3,3> second_wrt_T[6],
                                    MsqError& err )
{
  double f = Frobenius(T);
  double d = det(T);
  result = f*f*f - 3*MSQ_SQRT_THREE*d;

  deriv_wrt_T = T;
  deriv_wrt_T *= f;
  deriv_wrt_T -= MSQ_SQRT_THREE*transpose_adj(T);
  deriv_wrt_T *= 3;
  
  set_scaled_2nd_deriv_of_det( second_wrt_T, -3 * MSQ_SQRT_THREE, T );
  if (f > 1e-50)
    pluseq_scaled_outer_product( second_wrt_T, 3.0/f, T );
  else
    std::cout << "Warning: Division by zero avoided in TShapeNB1::evaluate_with_hess()" << std::endl;
  pluseq_scaled_I( second_wrt_T, 3.0*f );
  return true;
}

} // namespace Mesquite
