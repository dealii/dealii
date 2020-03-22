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


/** \file StdDevTemplate.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "StdDevTemplate.hpp"
#include "QualityMetric.hpp"
#include "MsqError.hpp"
#include "MsqHessian.hpp"
#include "PatchData.hpp"

namespace MESQUITE_NS {

ObjectiveFunction* StdDevTemplate::clone() const
  { return new StdDevTemplate(*this); }

bool StdDevTemplate::evaluate( EvalType type, 
                               PatchData& pd,
                               double& value_out,
                               bool free,
                               MsqError& err )
{
  bool result = VarianceTemplate::evaluate( type, pd, value_out, free, err );
  if (MSQ_CHKERR(err) || !result)
    return false;
  
  const double neg = get_quality_metric()->get_negate_flag();
  value_out = neg * sqrt( neg * value_out );
  return true;
}

bool StdDevTemplate::evaluate_with_gradient( EvalType type, 
                                             PatchData& pd,
                                             double& value_out,
                                             std::vector<Vector3D>& grad_out,
                                             MsqError& err )
{
  bool result = VarianceTemplate::evaluate_with_gradient( type, pd, value_out, grad_out, err );
  if (MSQ_CHKERR(err) || !result)
    return false;
  
  const double neg = get_quality_metric()->get_negate_flag();
  value_out *= neg; // undo any negation done by VarianceTemplate
  value_out = sqrt( value_out ); // standard deviation
  const double factor = 1.0/(2.0 * value_out);  
  for (std::vector<Vector3D>::iterator i = grad_out.begin(); i != grad_out.end(); ++i)
    *i *= factor;
  value_out *= neg; // redo any negation done by VariandeTemplate
  
  return true;
}


bool StdDevTemplate::evaluate_with_Hessian_diagonal( EvalType type, 
                                        PatchData& pd,
                                        double& value_out,
                                        std::vector<Vector3D>& grad_out,
                                        std::vector<SymMatrix3D>& hess_diag_out,
                                        MsqError& err )
{
  bool result = VarianceTemplate::evaluate_with_Hessian_diagonal( type, pd, value_out, grad_out, hess_diag_out, err );
  if (MSQ_CHKERR(err) || !result)
    return false;
  
  const double neg = get_quality_metric()->get_negate_flag();
  value_out *= neg; // undo any negation done by VarianceTemplate
  value_out = sqrt( value_out ); // standard deviation
  const double f1 = 1.0/(2.0 * value_out);
  const double f2 = neg * -0.25 / (value_out * value_out * value_out);
  for (size_t i = 0; i < grad_out.size(); ++i) {
    hess_diag_out[i] *= f1;
    hess_diag_out[i] += f2 * outer( grad_out[i] );
    grad_out[i] *= f1;
  }
  
  value_out *= neg; // redo any negation done by VariandeTemplate
  return true;
}


} // namespace Mesquite
