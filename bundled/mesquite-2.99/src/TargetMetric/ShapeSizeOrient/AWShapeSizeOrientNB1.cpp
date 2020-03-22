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


/** \file AWShapeSizeOrientNB1.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "AWShapeSizeOrientNB1.hpp"
#include "MsqMatrix.hpp"
#include "TMPDerivs.hpp"
#include "TMPCommon.hpp"

namespace MESQUITE_NS {

std::string AWShapeSizeOrientNB1::get_name() const
  { return "AWShapeSizeOrientNB1"; }

AWShapeSizeOrientNB1::~AWShapeSizeOrientNB1() {}

template <unsigned DIM> static inline
bool eval( const MsqMatrix<DIM,DIM>& A, 
           const MsqMatrix<DIM,DIM>& W, 
           double& result)
{
  result = sqr_Frobenius( A - W );
  return true;
}

template <unsigned DIM> static inline
bool grad( const MsqMatrix<DIM,DIM>& A, 
           const MsqMatrix<DIM,DIM>& W, 
           double& result, 
           MsqMatrix<DIM,DIM>& deriv )
{
  deriv = A - W;
  result = sqr_Frobenius( deriv );
  deriv *= 2.0;
  return true;
}

template <unsigned DIM> static inline
bool hess( const MsqMatrix<DIM,DIM>& A, 
           const MsqMatrix<DIM,DIM>& W, 
           double& result, 
           MsqMatrix<DIM,DIM>& deriv, 
           MsqMatrix<DIM,DIM>* second )
{
  deriv = A - W;
  result = sqr_Frobenius( deriv );
  deriv *= 2.0;
  set_scaled_I( second, 2.0 );
  return true;
}

TMP_AW_TEMPL_IMPL_COMMON(AWShapeSizeOrientNB1)

} // namespace Mesquite
