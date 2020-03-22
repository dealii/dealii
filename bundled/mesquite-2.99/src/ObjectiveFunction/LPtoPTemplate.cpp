/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Sandia Corporation and Argonne National
    Laboratory.  Under the terms of Contract DE-AC04-94AL85000 
    with Sandia Corporation, the U.S. Government retains certain 
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
 
    diachin2@llnl.gov, djmelan@sandia.gov, mbrewer@sandia.gov, 
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov      
   
  ***************************************************************** */
/*!
  \file   LPtoPTemplate.cpp
  \brief  

  This Objective Function is evaluated using an L P norm to the pth power.
  total=(sum (x_i)^pVal)
  \author Michael Brewer
  \author Thomas Leurent
  \date   2002-01-23
*/
#include <math.h>
#include "LPtoPTemplate.hpp"
#include "MsqFreeVertexIndexIterator.hpp"
#include "MsqTimer.hpp"
#include "MsqHessian.hpp"
#include "MsqDebug.hpp"
#include "QualityMetric.hpp"

using  namespace Mesquite;  

LPtoPTemplate::LPtoPTemplate(QualityMetric *qualitymetric, short Pinput, MsqError &err)
  : ObjectiveFunctionTemplate(qualitymetric)
{
  pVal=Pinput;
  if(pVal<1){
    MSQ_SETERR(err)("P_VALUE must be greater than 0.", MsqError::INVALID_ARG);
    return;
  }

  dividingByN=false;

  clear();
}

LPtoPTemplate::LPtoPTemplate( short P, QualityMetric* qm )
  : ObjectiveFunctionTemplate(qm), pVal(P), dividingByN(false)
  { clear(); }
  
void LPtoPTemplate::clear()
{
  mCount = 0;
  mPowSum = 0;
  saveCount = 0;
  savePowSum = 0;
}

//Michael:  need to clean up here
LPtoPTemplate::~LPtoPTemplate(){

}

ObjectiveFunction* LPtoPTemplate::clone() const
  { return new LPtoPTemplate(*this); }

double LPtoPTemplate::get_value( double power_sum, size_t count, EvalType type,
                                 size_t& global_count, MsqError& err )
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
  
//  if (!global_count)
//    {
//      MSQ_SETERR(err)(" global_count is zero, possibly due to an invalid mesh.", MsqError::INVALID_MESH);
//      return -1;  // result is invalid
//    }   
  if (dividingByN && global_count)
    result /= global_count;
  return result;
}

bool LPtoPTemplate::evaluate( EvalType type, 
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
    
    double tmp_val = value;
    for (short j = 1; j < pVal; ++j)
      tmp_val *= value;
    working_sum += fabs(tmp_val);
  }
  
    // get overall OF value, update member data, etc.
  size_t global_count;
  value_out = qm->get_negate_flag() * get_value( working_sum, qmHandles.size(), type, global_count, err );
//  if (!global_count)
//    return false;  // invalid mesh
//  else
  return true;
}

bool LPtoPTemplate::evaluate_with_gradient( EvalType type, 
                                            PatchData& pd,
                                            double& OF_val,
                                            std::vector<Vector3D>& grad_out,
                                            MsqError& err )
{
  QualityMetric* qm = get_quality_metric();
  qm->get_evaluations( pd, qmHandles, OF_FREE_EVALS_ONLY, err );  MSQ_ERRFALSE(err);
  
    // zero gradient
  grad_out.clear();
  grad_out.resize( pd.num_free_vertices(), Vector3D(0.0,0.0,0.0) );
  bool qm_bool=true;
  double QM_val;
  OF_val = 0.;
  int p1;
  
    // calculate OF value and gradient for just the patch
  std::vector<size_t>::const_iterator i;
  for (i = qmHandles.begin(); i != qmHandles.end(); ++i)
  {
    qm_bool = qm->evaluate_with_gradient( pd, *i, QM_val, mIndices, mGradient, err );
    if (MSQ_CHKERR(err) || !qm_bool)
      return false;
    
    QM_val = fabs(QM_val);
    double QM_pow = 1.0;
    double factor = qm->get_negate_flag();
    if (pVal == 1) 
      QM_pow = 1.0;
    else {
      QM_pow = QM_val;
      for (p1 = 2; p1 < pVal; ++p1)
        QM_pow *= QM_val;
      factor *= QM_pow * pVal;
    }
    
    OF_val += QM_pow * QM_val;
    for (size_t j = 0; j < mIndices.size(); ++j) {
      mGradient[j] *= factor;
      grad_out[mIndices[j]] += mGradient[j];
    }
  }
  
    // get overall OF value, update member data, etc.
  size_t global_count;
  OF_val = qm->get_negate_flag() * get_value( OF_val, qmHandles.size(), type, global_count, err );
//  if (!global_count)
//    return false;  // invalid mesh

  if (dividingByN && global_count) {
    const double inv_n = 1.0/global_count;
    std::vector<Vector3D>::iterator g;
    for (g = grad_out.begin(); g != grad_out.end(); ++g)
      *g *= inv_n;
  }  
  
  return true;
}
  
bool LPtoPTemplate::evaluate_with_Hessian_diagonal( EvalType type, 
                                        PatchData& pd,
                                        double& OF_val,
                                        std::vector<Vector3D>& grad,
                                        std::vector<SymMatrix3D>& hess_diag,
                                        MsqError& err )
{
  QualityMetric* qm = get_quality_metric();
  qm->get_evaluations( pd, qmHandles, OF_FREE_EVALS_ONLY, err );  MSQ_ERRFALSE(err);
  
    // zero gradient and hessian
  grad.clear();
  grad.resize( pd.num_free_vertices(), 0.0 );
  hess_diag.clear();
  hess_diag.resize( pd.num_free_vertices(), 0.0 );
  
  double QM_val, QM_pow = 1.0;
  double fac1, fac2;
  const double negate_flag = qm->get_negate_flag();
  bool qm_bool;
  size_t i;
  short p;
   
  // Loops over all elements in the patch.
  OF_val = 0.0;
  std::vector<size_t>::const_iterator k;
  for (k = qmHandles.begin(); k != qmHandles.end(); ++k)
  {
    // Computes \nabla^2 Q(e). Only the free vertices will have non-zero entries. 
    qm_bool = qm->evaluate_with_Hessian_diagonal( pd, *k, QM_val, mIndices, mGradient, mDiag, err );
    if (MSQ_CHKERR(err) || !qm_bool) return false;
    QM_val = fabs(QM_val);

    // **** Computes Hessian ****
    const size_t nve = mIndices.size();
    if (pVal == 1) {
      QM_pow = 1.0;
      for (i=0; i<nve; ++i) {
        mDiag[i] *= negate_flag;
        hess_diag[mIndices[i]] += mDiag[i];
      }
      fac1 = 1;
    }
    else if (pVal >= 2) {
       // Computes the coefficients:
      QM_pow = 1.0;
      for (p=0; p<pVal-2; ++p)
        QM_pow *= QM_val;
      // 1 - computes p(p-1)Q(e)^{p-2}
      fac2 = pVal* (pVal-1) * QM_pow;
      // 2 - computes  pQ(e)^{p-1}
      QM_pow *= QM_val;
      fac1 = pVal * QM_pow;

        //fac1 *= qm->get_negate_flag();
        //fac2 *= qm->get_negate_flag();

      for (i=0; i<nve; ++i) {
        SymMatrix3D op(mGradient[i]);
        op *= fac2;
        mDiag[i] *= fac1;
        op += mDiag[i];
        op *= negate_flag;
        hess_diag[mIndices[i]] += op;
      }
    } else {
      MSQ_SETERR(err)(" invalid P value.", MsqError::INVALID_STATE);
      return false;
    }


    // **** Computes Gradient ****

    // For each vertex in the element ... 
    for (i=0; i<nve; ++i) {
      // ... computes p*q^{p-1}*grad(q) ...
      mGradient[i] *= fac1*negate_flag;
      // ... and accumulates it in the objective function gradient.
        //also scale the gradient by the scaling factor
      assert (mIndices[i] < pd.num_free_vertices());
      grad[mIndices[i]] += mGradient[i];
    }

    // **** computes Objective Function value \sum_{i=1}^{N_e} |q_i|^P ****
    OF_val += QM_pow * QM_val;
  }

  size_t global_count;
  OF_val = negate_flag
         * get_value( OF_val, qmHandles.size(), type, global_count, err );
//  if (!global_count)
//    return false;  // invalid mesh

  if (dividingByN && global_count) {
    const double inv_n = 1.0 / global_count;
    for (i = 0; i < pd.num_free_vertices(); ++i) {
      grad[i] *= inv_n;
      hess_diag[i] *= inv_n;
    }
  }
  
  return true;
}
	
/*\ For each element, each entry to be accumulated in the Hessian for
    this objective function (\f$ \sum_{e \in E} Q(e)^p \f$ where \f$ E \f$
    is the set of all elements in the patch) has the form:
    \f$ pQ(e)^{p-1} \nabla^2 Q(e) + p(p-1)Q(e)^{p-2} \nabla Q(e) [\nabla Q(e)]^T \f$.

    For \f$ p=2 \f$, this simplifies to
    \f$ 2Q(e) \nabla^2 Q(e) + 2 \nabla Q(e) [\nabla Q(e)]^T \f$.

    For \f$ p=1 \f$, this simplifies to \f$ \nabla^2 Q(e) \f$.

    The \f$ p=1 \f$ simplified version is implemented directly
    to speed up computation. 

    This function does not support vertex-based metrics.
    
    \param pd The PatchData object for which the objective function
           hessian is computed.
    \param hessian this object must have been previously initialized.
*/
bool LPtoPTemplate::evaluate_with_Hessian( EvalType type, 
                                           PatchData& pd,
                                           double& OF_val,
                                           std::vector<Vector3D>& grad,
                                           MsqHessian& hessian,
                                           MsqError& err )
{
  QualityMetric* qm = get_quality_metric();
  qm->get_evaluations( pd, qmHandles, OF_FREE_EVALS_ONLY, err );  MSQ_ERRFALSE(err);
  double negate_flag = qm->get_negate_flag();
  
    // zero gradient and hessian
  grad.clear();
  grad.resize( pd.num_free_vertices(), 0.0 );
  hessian.zero_out();
  
  double QM_val, QM_pow = 1.0;
  double fac1, fac2;
  Matrix3D elem_outer_product;
  bool qm_bool;
  size_t i, j, n;
  short p;
   
  // Loops over all elements in the patch.
  OF_val = 0.0;
  std::vector<size_t>::const_iterator k;
  for (k = qmHandles.begin(); k != qmHandles.end(); ++k)
  {
    // Computes \nabla^2 Q(e). Only the free vertices will have non-zero entries. 
    qm_bool = qm->evaluate_with_Hessian( pd, *k, QM_val, mIndices, mGradient, mHessian, err );
    if (MSQ_CHKERR(err) || !qm_bool) 
      return false;
    QM_val = fabs(QM_val);

    // **** Computes Hessian ****
    const size_t nve = mIndices.size();
    if (pVal == 1) {
      QM_pow = 1.0;
      n=0;
      for (i=0; i<nve; ++i) {
        for (j=i; j<nve; ++j) {
            //negate if necessary
          mHessian[n] *= negate_flag;
          hessian.add( mIndices[i], mIndices[j], mHessian[n], err ); MSQ_ERRFALSE(err);
          ++n;
        }
      }
      fac1 = 1;
    }
    else if (pVal >= 2) {
       // Computes the coefficients:
      QM_pow = 1.0;
      for (p=0; p<pVal-2; ++p)
        QM_pow *= QM_val;
      // 1 - computes p(p-1)Q(e)^{p-2}
      fac2 = pVal* (pVal-1) * QM_pow;
      // 2 - computes  pQ(e)^{p-1}
      QM_pow *= QM_val;
      fac1 = pVal * QM_pow;

        //fac1 *= qm->get_negate_flag();
        //fac2 *= qm->get_negate_flag();

      n=0;
      for (i=0; i<nve; ++i) {
        for (j=i; j<nve; ++j) {
          elem_outer_product.outer_product(mGradient[i], mGradient[j]);
          elem_outer_product *= fac2;
          mHessian[n] *= fac1;
          mHessian[n] += elem_outer_product;
          mHessian[n] *= negate_flag;
          hessian.add( mIndices[i], mIndices[j], mHessian[n], err ); MSQ_ERRFALSE(err);
          ++n;
        }
      }
    } else {
      MSQ_SETERR(err)(" invalid P value.", MsqError::INVALID_STATE);
      return false;
    }


    // **** Computes Gradient ****

    // For each vertex in the element ... 
    for (i=0; i<nve; ++i) {
      // ... computes p*q^{p-1}*grad(q) ...
      mGradient[i] *= fac1*negate_flag;
      // ... and accumulates it in the objective function gradient.
        //also scale the gradient by the scaling factor
      assert (mIndices[i] < pd.num_free_vertices());
      grad[mIndices[i]] += mGradient[i];
    }

    // **** computes Objective Function value \sum_{i=1}^{N_e} |q_i|^P ****
    OF_val += QM_pow * QM_val;
  }

  size_t global_count;
  OF_val = negate_flag
         * get_value( OF_val, qmHandles.size(), type, global_count, err );
//  if (!global_count)
//    return false;  // invalid mesh

  if (dividingByN && global_count) {
    const double inv_n = 1.0 / global_count;
    std::vector<Vector3D>::iterator g;
    for (g = grad.begin(); g != grad.end(); ++g)
      *g *= inv_n;
    hessian.scale( inv_n );
  }
  
  return true;
}
