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


/** \file ElementPMeanP.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "ElementPMeanP.hpp"
#include "ElemSampleQM.hpp"
#include "MsqError.hpp"
#include "PatchData.hpp"

namespace MESQUITE_NS {

ElementPMeanP::ElementPMeanP( double p, ElemSampleQM* metric )
    : PMeanPMetric( p ), 
      mMetric(metric)
    {}
    
ElementPMeanP::~ElementPMeanP() {}

std::string ElementPMeanP::get_name() const
{
  std::string result("ElementPMeanP(");
  result += mMetric->get_name();
  result += ")";
  return result;
}


int ElementPMeanP::get_negate_flag() const
{
  return get_quality_metric()->get_negate_flag();
}

bool ElementPMeanP::evaluate( PatchData& pd, 
                              size_t handle, 
                              double& value, 
                              MsqError& err )
{
  ElemSampleQM* qm = get_quality_metric();
  mHandles.clear();
  qm->get_element_evaluations( pd, handle, mHandles, err ); MSQ_ERRFALSE(err);
  bool result = average( pd, qm, mHandles, value, err );
  return !MSQ_CHKERR(err) && result;
}

bool ElementPMeanP::evaluate_with_gradient( PatchData& pd, 
                                            size_t handle, 
                                            double& value, 
                                            std::vector<size_t>& indices,
                                            std::vector<Vector3D>& gradient,
                                            MsqError& err )
{
  ElemSampleQM* qm = get_quality_metric();
  mHandles.clear();
  qm->get_element_evaluations( pd, handle, mHandles, err ); MSQ_ERRFALSE(err);
  bool result = average_with_gradient( pd, qm, mHandles, value, indices, gradient, err );
  return !MSQ_CHKERR(err) && result;
}

bool ElementPMeanP::evaluate_with_Hessian( PatchData& pd, 
                          size_t handle, 
                          double& value, 
                          std::vector<size_t>& indices,
                          std::vector<Vector3D>& gradient,
                          std::vector<Matrix3D>& Hessian,
                          MsqError& err )
{
  ElemSampleQM* qm = get_quality_metric();
  mHandles.clear();
  qm->get_element_evaluations( pd, handle, mHandles, err ); MSQ_ERRFALSE(err);
  bool result = average_with_Hessian( pd, qm, mHandles, value, indices, gradient, Hessian, err );
  return !MSQ_CHKERR(err) && result;
}

bool ElementPMeanP::evaluate_with_Hessian_diagonal( PatchData& pd, 
                          size_t handle, 
                          double& value, 
                          std::vector<size_t>& indices,
                          std::vector<Vector3D>& gradient,
                          std::vector<SymMatrix3D>& diagonal,
                          MsqError& err )
{
  ElemSampleQM* qm = get_quality_metric();
  mHandles.clear();
  qm->get_element_evaluations( pd, handle, mHandles, err ); MSQ_ERRFALSE(err);
  bool result = average_with_Hessian_diagonal( pd, qm, mHandles, value, indices, gradient, diagonal, err );
  return !MSQ_CHKERR(err) && result;
}


} // namespace Mesquite
