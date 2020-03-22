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


/** \file InvTransBarrier.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "InvTransBarrier.hpp"
#include "MsqMatrix.hpp"
#include "MsqError.hpp"

namespace MESQUITE_NS {

std::string InvTransBarrier::get_name() const
  { return "InvTransBarrier"; }

InvTransBarrier::~InvTransBarrier() {}

bool InvTransBarrier::evaluate( const MsqMatrix<2,2>& T, 
                                double& result, MsqError& err )
{
  double tau = det(T);
  if (invalid_determinant(tau)) {
    MSQ_SETERR(err)( barrier_violated_msg, MsqError::BARRIER_VIOLATED );
    return false;
  }
  MsqMatrix<2,2> Tp = transpose_adj(T);
  Tp *= 1.0/tau;
  bool rval = metricPtr->evaluate( Tp, result, err );
  return !MSQ_CHKERR(err) && rval;
}

bool InvTransBarrier::evaluate( const MsqMatrix<3,3>& T, 
                                double& result, 
                                MsqError& err )
{
  double tau = det(T);
  if (invalid_determinant(tau)) {
    MSQ_SETERR(err)( barrier_violated_msg, MsqError::BARRIER_VIOLATED );
    return false;
  }
  MsqMatrix<3,3> Tp = transpose_adj(T);
  Tp *= 1.0/tau;
  bool rval = metricPtr->evaluate( Tp, result, err );
  return !MSQ_CHKERR(err) && rval;
}


} // namespace Mesquite
