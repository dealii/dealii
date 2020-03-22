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


/** \file TPower2.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "TMetricBarrier.hpp"
#include "TPower2.hpp"
#include "MsqMatrix.hpp"
#include "TMPDerivs.hpp"
#include "MsqError.hpp"
#include "TMPCommon.hpp"

namespace MESQUITE_NS {

std::string TPower2::get_name() const
  { return "sqr(" + mMetric->get_name() + ')'; }

TPower2::~TPower2() {}

template <unsigned DIM> inline
bool TPower2::eval( const MsqMatrix<DIM,DIM>& T, 
                    double& result,
                    MsqError& err )
{
  bool rval = mMetric->evaluate( T, result, err );
  MSQ_ERRZERO(err);
  result *= result;
  return rval;
}

template <unsigned DIM> inline
bool TPower2::grad( const MsqMatrix<DIM,DIM>& T, 
                    double& result, 
                    MsqMatrix<DIM,DIM>& deriv_wrt_T,
                    MsqError& err )
{
  bool rval = mMetric->evaluate_with_grad( T, result, deriv_wrt_T, err );
  MSQ_ERRZERO(err);
  deriv_wrt_T *= 2 * result;
  result *= result;
  return rval;
}

template <unsigned DIM> inline
bool TPower2::hess( const MsqMatrix<DIM,DIM>& T, 
                    double& result, 
                    MsqMatrix<DIM,DIM>& deriv_wrt_T, 
                    MsqMatrix<DIM,DIM>* second_wrt_T,
                    MsqError& err )
{
  bool rval = mMetric->evaluate_with_hess( T, result, deriv_wrt_T, second_wrt_T, err );//
  MSQ_ERRZERO(err);
  hess_scale( second_wrt_T, 2*result );
  pluseq_scaled_outer_product( second_wrt_T, 2.0, deriv_wrt_T );
  deriv_wrt_T *= 2 * result;
  result *= result;
  return rval;
}


TMP_T_TEMPL_IMPL_COMMON_ERR(TPower2)


} // namespace MESQUITE_NS
