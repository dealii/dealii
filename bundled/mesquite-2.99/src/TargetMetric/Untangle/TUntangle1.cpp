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


/** \file TUntangle1.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "TUntangle1.hpp"
#include "TMPDerivs.hpp"
#include "TMPCommon.hpp"

#include <iostream>

namespace MESQUITE_NS {


TUntangle1::~TUntangle1()
{}

std::string TUntangle1::get_name() const
  { return "Untangle1"; }


template <unsigned DIM> inline
bool TUntangle1::eval( const MsqMatrix<DIM,DIM>& T, 
                       double& result )
{
  double tau = det(T);
  result = 0.5 * (sqrt(tau*tau + mFactor) - tau);
  return true;
}

template <unsigned DIM> inline
bool TUntangle1::grad( const MsqMatrix<DIM,DIM>& T, 
                       double& result, 
                       MsqMatrix<DIM,DIM>& deriv_wrt_T )
{
  double tau = det(T);
  double g = sqrt(tau*tau + mFactor);
  if (g > 1e-50)
  {
    double f = tau/g - 1;
    result = 0.5 * (g - tau);
    deriv_wrt_T = transpose_adj(T);
    deriv_wrt_T *= 0.5 * f;
    return true;
  }
  else
  {
    std::cout << "Warning: Division by zero avoided in TUntangle1::grad()" << std::endl;
    return false;
  }
}

template <unsigned DIM> inline
bool TUntangle1::hess( const MsqMatrix<DIM,DIM>& T, 
                       double& result, 
                       MsqMatrix<DIM,DIM>& deriv_wrt_T, 
                       MsqMatrix<DIM,DIM>* second_wrt_T )
{
  const MsqMatrix<DIM,DIM> adjt = transpose_adj(T);
  double tau = det(T);
  double g = sqrt(tau*tau + mFactor);
  if (g == 0)
    return false;
  double f = 0.5 * (tau/g - 1);
  result = 0.5 * (g - tau);
  
  deriv_wrt_T = adjt;
  deriv_wrt_T *= f;
  
  set_scaled_outer_product( second_wrt_T, 0.5*mFactor/(g*g*g), adjt );
  pluseq_scaled_2nd_deriv_of_det( second_wrt_T, f, T );
  
  return true;
}


TMP_T_TEMPL_IMPL_COMMON(TUntangle1)


} // namespace MESQUITE_NS
