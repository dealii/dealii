/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Lawrence Livermore National Laboratory.  Under 
    the terms of Contract B545069 with the University of Wisconsin -- 
    Madison, Lawrence Livermore National Laboratory retains certain
    rights in this software.

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


/** \file PMeanPTemplate.cpp
 *  \brief previous name: PowerMeanP.cpp
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "PMeanPTemplate.hpp"
#include "QualityMetric.hpp"
#include "MsqError.hpp"
#include "MsqHessian.hpp"
#include "PatchData.hpp"

namespace MESQUITE_NS {

ObjectiveFunction* PMeanPTemplate::clone() const
  { return new PMeanPTemplate(*this); }
  
void PMeanPTemplate::clear()
{
  mCount = 0;
  mPowSum = 0;
  saveCount = 0;
  savePowSum = 0;
}

double PMeanPTemplate::get_value( double power_sum, 
                              size_t count, 
                              EvalType type,
                              size_t& global_count )
{
  double result = 0;
  switch (type) 
  {
    case CALCULATE:
      result = power_sum;
      global_count = count;
      break;
    
    case ACCUMULATE:
      mPowSum += power_sum;
      mCount += count;
      result = mPowSum;
      global_count = mCount;
      break;
    
    case SAVE:
      savePowSum = power_sum;
      saveCount = count;
      result = mPowSum;
      global_count = mCount;
      break;
    
    case UPDATE:
      mPowSum -= savePowSum;
      mCount -= saveCount;
      savePowSum = power_sum;
      saveCount = count;
      mPowSum += savePowSum;
      mCount += saveCount;
      result = mPowSum;
      global_count = mCount;
      break;
    
    case TEMPORARY:
      result = mPowSum - savePowSum + power_sum;
      global_count = mCount + count - saveCount;
      break;
  }
  
  return global_count ? result/global_count : 0.0;
}

bool PMeanPTemplate::evaluate( EvalType type, 
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
  
    // get overall OF value, update member data, etc.
  size_t global_count = 0;
  value_out = qm->get_negate_flag() 
            * get_value( working_sum, qmHandles.size(), type, global_count );
  return true;
}

bool PMeanPTemplate::evaluate_with_gradient( EvalType type, 
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
  const double f = qm->get_negate_flag() * mPower.value();
  for (i = qmHandles.begin(); i != qmHandles.end(); ++i)
  {
    bool result = qm->evaluate_with_gradient( pd, *i, value, mIndices, mGradient, err );
    if (MSQ_CHKERR(err) || !result)
      return false;
    if (fabs(value) < DBL_EPSILON)
      continue;
    
    const double r1 = mPowerMinus1.raise( value );
    const double qmp = r1 * value;
    working_sum += qmp;
    value = f * r1;

    for (size_t j = 0; j < mIndices.size(); ++j) {
      mGradient[j] *= value;
      grad_out[mIndices[j]] += mGradient[j];
    }
  }
  
    // get overall OF value, update member data, etc.
  size_t global_count = 0;
  value_out = qm->get_negate_flag() 
            * get_value( working_sum, qmHandles.size(), type, global_count );
  if (global_count)
  {
    const double inv_n = 1.0 / global_count;
    std::vector<Vector3D>::iterator g;
    for (g = grad_out.begin(); g != grad_out.end(); ++g)
      *g *= inv_n;
  }
  return true;
}

bool PMeanPTemplate::evaluate_with_Hessian_diagonal( EvalType type, 
                                        PatchData& pd,
                                        double& value_out,
                                        std::vector<Vector3D>& grad_out,
                                        std::vector<SymMatrix3D>& hess_diag_out,
                                        MsqError& err )
{
  QualityMetric* qm = get_quality_metric();
  qm->get_evaluations( pd, qmHandles, OF_FREE_EVALS_ONLY, err );  MSQ_ERRFALSE(err);
  
    // zero gradient and hessian
  const size_t s = pd.num_free_vertices();
  grad_out.clear();
  grad_out.resize( s, 0.0 );
  hess_diag_out.clear();
  hess_diag_out.resize( s, 0.0 );
  
    // calculate OF value and gradient for just the patch
  std::vector<size_t>::const_iterator i;
  size_t j;
  double value, working_sum = 0.0;
  const double f1 = qm->get_negate_flag() * mPower.value();
  const double f2 = f1 * (mPower.value() - 1);
  for (i = qmHandles.begin(); i != qmHandles.end(); ++i)
  {
    bool result = qm->evaluate_with_Hessian_diagonal( pd, *i, value, mIndices, mGradient, mDiag, err );
    if (MSQ_CHKERR(err) || !result)
      return false;
    if (fabs(value) < DBL_EPSILON)
      continue;
    
    const size_t nfree = mIndices.size();
    if (mPower.value() == 1.0) {
      working_sum += mPower.raise( value );
      for (j = 0; j < nfree; ++j) {
        const size_t idx = mIndices[j];
        hess_diag_out[idx] += f1 * mDiag[j];
        mGradient[j] *= f1;
        grad_out[idx] += mGradient[j];
      }
    }
    else {
      const double r2 = mPowerMinus2.raise( value );
      const double r1 = r2 * value;
      working_sum += r1 * value;
      const double hf = f2 * r2;
      const double gf = f1 * r1;
      for (j = 0; j < nfree; ++j) {
        const size_t idx = mIndices[j];

        hess_diag_out[idx] += hf * outer( mGradient[j] );
        hess_diag_out[idx] += gf * mDiag[j];

        mGradient[j] *= gf;
        grad_out[idx] += mGradient[j];
      }
    }
  }
  
    // get overall OF value, update member data, etc.
  size_t global_count = 0;
  value_out = qm->get_negate_flag() 
            * get_value( working_sum, qmHandles.size(), type, global_count );
  if (global_count)
  {
    const double inv_n = 1.0 / global_count;
    for (j = 0; j < s; ++j) {
      grad_out[j] *= inv_n;
      hess_diag_out[j] *= inv_n;
    }
  }  
  return true;
}

bool PMeanPTemplate::evaluate_with_Hessian( EvalType type, 
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
  const double f1 = qm->get_negate_flag() * mPower.value();
  const double f2 = f1 * (mPower.value() - 1);
  Matrix3D m;
  for (i = qmHandles.begin(); i != qmHandles.end(); ++i)
  {
    bool result = qm->evaluate_with_Hessian( pd, *i, value, mIndices, mGradient, mHessian, err );
    if (MSQ_CHKERR(err) || !result)
      return false;
    if (fabs(value) < DBL_EPSILON)
      continue;
    
    const size_t nfree = mIndices.size();
    n = 0;
    if (mPower.value() == 1.0) {
      working_sum += mPower.raise( value );
      for (j = 0; j < nfree; ++j) {
        mGradient[j] *= f1;
        grad_out[mIndices[j]] += mGradient[j];
        for (k = j; k < nfree; ++k) {
          mHessian[n] *= f1;
          Hessian_out.add( mIndices[j], mIndices[k], mHessian[n], err );  MSQ_ERRFALSE(err);
          ++n;
        }
      }
    }
    else {
      const double r2 = mPowerMinus2.raise( value );
      const double r1 = r2 * value;
      working_sum += r1 * value;
      const double hf = f2 * r2;
      const double gf = f1 * r1;
      for (j = 0; j < nfree; ++j) {
        for (k = j; k < nfree; ++k) {
          m.outer_product( mGradient[j], mGradient[k] );
          m *= hf;
          mHessian[n] *= gf;
          m += mHessian[n];
          Hessian_out.add( mIndices[j], mIndices[k], m, err );  MSQ_ERRFALSE(err);
          ++n;
        }
      }
      for (j = 0; j < nfree; ++j) {
        mGradient[j] *= gf;
        grad_out[mIndices[j]] += mGradient[j];
      }
    }
  }
  
    // get overall OF value, update member data, etc.
  size_t global_count = 0;
  value_out = qm->get_negate_flag() 
            * get_value( working_sum, qmHandles.size(), type, global_count );
  if (global_count)
  {
    const double inv_n = 1.0 / global_count;
    std::vector<Vector3D>::iterator g;
    for (g = grad_out.begin(); g != grad_out.end(); ++g)
      *g *= inv_n;
    Hessian_out.scale( inv_n );
  }
  return true;
}



} // namespace Mesquite
