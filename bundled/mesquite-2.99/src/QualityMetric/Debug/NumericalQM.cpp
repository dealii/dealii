/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2009 Sandia National Laboratories.  Developed at the
University NumericalQM::of Wisconsin--Madison under SNL contract number
    624796.  The U.S. Government and the University of Wisconsin
retain NumericalQM::certain rights to this software.

This NumericalQM::library is free software; you can redistribute it and/or
modify NumericalQM::it under the terms of the GNU Lesser General Public
License NumericalQM::as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

This NumericalQM::library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY NumericalQM::or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

You NumericalQM::should have received a copy of the GNU Lesser General Public License 
    (lgpl.txt) along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    (2009) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file NumericalQM.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "NumericalQM.hpp"

namespace MESQUITE_NS {


NumericalQM::NumericalQM( QualityMetric* real_metric,
                          bool numerical_gradient,
                          bool numerical_hessian )
  : realMetric(real_metric),
    numericGrad(numerical_gradient),
    numericHess(numerical_hessian)
{}


QualityMetric::MetricType NumericalQM::get_metric_type() const
  { return realMetric->get_metric_type(); }

std::string NumericalQM::get_name() const
  { return realMetric->get_name(); }

int NumericalQM::get_negate_flag() const
  { return realMetric->get_negate_flag(); }

void NumericalQM::get_evaluations( PatchData& pd, 
                                   std::vector<size_t>& handles, 
                                   bool free,
                                   MsqError& err )
  { return realMetric->get_evaluations( pd, handles, free, err ); }

bool NumericalQM::evaluate( PatchData& pd, 
                            size_t handle, 
                            double& value, 
                            MsqError& err )
  { return realMetric->evaluate( pd, handle, value, err ); }

bool NumericalQM::evaluate_with_indices( PatchData& pd,
                                         size_t handle,
                                         double& value,
                                         std::vector<size_t>& indices,
                                         MsqError& err )
  { return realMetric->evaluate_with_indices( pd, handle, value, indices, err ); }

bool NumericalQM::evaluate_with_gradient( PatchData& pd,
                                          size_t handle,
                                          double& value,
                                          std::vector<size_t>& indices,
                                          std::vector<Vector3D>& gradient,
                                          MsqError& err )
{
  if (numericGrad) 
    return realMetric->QualityMetric::evaluate_with_gradient( pd, 
                                                              handle,
                                                              value,
                                                              indices,
                                                              gradient, 
                                                              err );
  else
    return realMetric->evaluate_with_gradient( pd, 
                                               handle,
                                               value,
                                               indices,
                                               gradient, 
                                               err );
}

bool NumericalQM::evaluate_with_Hessian_diagonal( PatchData& pd,
                                      size_t handle,
                                      double& value,
                                      std::vector<size_t>& indices,
                                      std::vector<Vector3D>& gradient,
                                      std::vector<SymMatrix3D>& hess,
                                      MsqError& err )
{
  if (numericHess) 
    return realMetric->QualityMetric::evaluate_with_Hessian_diagonal( pd, 
                                                               handle,
                                                               value,
                                                               indices,
                                                               gradient, 
                                                               hess,
                                                               err );
  else
    return realMetric->evaluate_with_Hessian_diagonal( pd, 
                                                       handle,
                                                       value,
                                                       indices,
                                                       gradient, 
                                                       hess,
                                                       err );
}

bool NumericalQM::evaluate_with_Hessian( PatchData& pd,
                                         size_t handle,
                                         double& value,
                                         std::vector<size_t>& indices,
                                         std::vector<Vector3D>& gradient,
                                         std::vector<Matrix3D>& Hessian,
                                         MsqError& err )

{
  if (numericHess) 
    return realMetric->QualityMetric::evaluate_with_Hessian( pd, 
                                                             handle,
                                                             value,
                                                             indices,
                                                             gradient, 
                                                             Hessian,
                                                             err );
  else
    return realMetric->evaluate_with_Hessian( pd, 
                                              handle,
                                              value,
                                              indices,
                                              gradient, 
                                              Hessian,
                                              err );
}

} // namespace MESQUITE_NS
