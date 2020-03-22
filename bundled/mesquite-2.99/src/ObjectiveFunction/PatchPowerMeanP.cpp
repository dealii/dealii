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


/** \file PatchPowerMeanP.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "PatchPowerMeanP.hpp"
#include "QualityMetric.hpp"
#include "MsqError.hpp"
#include "MsqHessian.hpp"
#include "PatchData.hpp"
#include "PatchIterator.hpp"

namespace MESQUITE_NS {

ObjectiveFunction* PatchPowerMeanP::clone() const
  { return new PatchPowerMeanP(*this); }


bool PatchPowerMeanP::initialize_block_coordinate_descent( Mesh* mesh, 
                                                      MeshDomain* domain, 
                                                      const Settings* settings,
                                                      PatchSet* patch_set,
                                                      MsqError& err )
{
  clear();
  PatchIterator patches( patch_set );
  
  PatchData pd;
  pd.set_mesh( mesh );
  pd.set_domain( domain );
  if (settings)
    pd.attach_settings(settings);
  
  bool result = true;
  while (patches.get_next_patch( pd, err ) && !MSQ_CHKERR(err))
  {
    double value;
    bool b = evaluate( ObjectiveFunction::ACCUMULATE, pd, value, false, err ); 
    MSQ_ERRZERO(err);
    result = result && b;
  }
  return result;
}

bool PatchPowerMeanP::evaluate( EvalType type, 
                                PatchData& pd,
                                double& value_out,
                                bool free,
                                MsqError& err )
{
  QualityMetric* qm = get_quality_metric();
  if (type == ObjectiveFunction::ACCUMULATE)
    qm->get_single_pass( pd, qmHandles, free, err );
  else
    qm->get_evaluations( pd, qmHandles, free, err );  
  MSQ_ERRFALSE(err);
  
    // calculate OF value for just the patch
  std::vector<size_t>::const_iterator i;
  double value, working_sum = 0.0;
  for (i = qmHandles.begin(); i != qmHandles.end(); ++i)
  {
    bool result = qm->evaluate( pd, *i, value, err );
    if (MSQ_CHKERR(err) || !result)
      return false;
    
    working_sum += mPower.raise( value );
  }
  working_sum /= qmHandles.size();
  
    // get overall OF value, update member data, etc.
  size_t global_count;
  value_out = qm->get_negate_flag() 
            * get_value( working_sum, 1, type, global_count );
  return true;
}

bool PatchPowerMeanP::evaluate_with_gradient( EvalType type, 
                                              PatchData& pd,
                                              double& value_out,
                                              std::vector<Vector3D>& grad_out,
                                              MsqError& err )
{
  QualityMetric* qm = get_quality_metric();
  qm->get_evaluations( pd, qmHandles, OF_FREE_EVALS_ONLY, err );  MSQ_ERRFALSE(err);
  
    // zero gradient
  grad_out.clear();
  grad_out.resize( pd.num_free_vertices(), Vector3D(0.0,0.0,0.0) );
  
    // calculate OF value and gradient for just the patch
  std::vector<size_t>::const_iterator i;
  double value, working_sum = 0.0;
  const double f = qm->get_negate_flag() * mPower.value() / qmHandles.size();
  for (i = qmHandles.begin(); i != qmHandles.end(); ++i)
  {
    bool result = qm->evaluate_with_gradient( pd, *i, value, mIndices, mGradient, err );
    if (MSQ_CHKERR(err) || !result)
      return false;
    if (fabs(value) < DBL_EPSILON)
      continue;
    
    const double qmp = mPower.raise( value );
    working_sum += mPower.raise( value );
    value = f * qmp / value;

    for (size_t j = 0; j < mIndices.size(); ++j) {
      mGradient[j] *= value;
      grad_out[mIndices[j]] += mGradient[j];
    }
  }
  
    // get overall OF value, update member data, etc.
  size_t global_count;
  value_out = qm->get_negate_flag() 
            * get_value( working_sum, 1, type, global_count );
  const double inv_n = 1.0 / global_count;
  std::vector<Vector3D>::iterator g;
  for (g = grad_out.begin(); g != grad_out.end(); ++g)
    *g *= inv_n;
  return true;
}

bool PatchPowerMeanP::evaluate_with_Hessian( EvalType type, 
                                             PatchData& pd,
                                             double& value_out,
                                             std::vector<Vector3D>& grad_out,
                                             MsqHessian& Hessian_out,
                                             MsqError& err )
{
  QualityMetric* qm = get_quality_metric();
  qm->get_evaluations( pd, qmHandles, OF_FREE_EVALS_ONLY, err );  MSQ_ERRFALSE(err);
  
    // zero gradient and hessian
  grad_out.clear();
  grad_out.resize( pd.num_free_vertices(), 0.0 );
  Hessian_out.zero_out();
  
    // calculate OF value and gradient for just the patch
  std::vector<size_t>::const_iterator i;
  size_t j, k, n;
  double value, working_sum = 0.0;
  const double f1 = qm->get_negate_flag() * mPower.value() / qmHandles.size();
  const double f2 = f1 * (mPower.value() - 1);
  Matrix3D m;
  for (i = qmHandles.begin(); i != qmHandles.end(); ++i)
  {
    bool result = qm->evaluate_with_Hessian( pd, *i, value, mIndices, mGradient, mHessian, err );
    if (MSQ_CHKERR(err) || !result)
      return false;
    if (fabs(value) < DBL_EPSILON)
      continue;
    
    const double qmp = mPower.raise( value );
    const double hf = f2 * qmp / (value*value);
    const double gf = f1 * qmp / value;
    working_sum += qmp;

    const size_t nfree = mIndices.size();
    n = 0;
    for (j = 0; j < nfree; ++j) {
      for (k = j; k < nfree; ++k) {
        m.outer_product( mGradient[j], mGradient[k] );
        m *= hf;
        mHessian[n] *= gf;
        m += mHessian[n];
        ++n;
        Hessian_out.add( mIndices[j], mIndices[k], m, err );  MSQ_ERRFALSE(err);
      }
    }
    
    for (j = 0; j < nfree; ++j) {
      mGradient[j] *= gf;
      grad_out[mIndices[j]] += mGradient[j];
    }
  }
  
    // get overall OF value, update member data, etc.
  size_t global_count;
  value_out = qm->get_negate_flag() 
            * get_value( working_sum, 1, type, global_count );
  const double inv_n = 1.0 / global_count;
  std::vector<Vector3D>::iterator g;
  for (g = grad_out.begin(); g != grad_out.end(); ++g)
    *g *= inv_n;
  Hessian_out.scale( inv_n );
  return true;
}



} // namespace Mesquite
