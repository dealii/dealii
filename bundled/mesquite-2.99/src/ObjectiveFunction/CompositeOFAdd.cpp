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
  \file    CompositeOFAdd.cpp
  \brief  

  This Objective Function combines two Objective Functions by addition
  \author Michael Brewer
  \date   2002-06-24
*/
#include <math.h>
#include "ObjectiveFunction.hpp"
#include "CompositeOFAdd.hpp"
#include "MsqTimer.hpp"
#include "PatchData.hpp"

namespace MESQUITE_NS {


/*!
Sets the QualityMetric pointer to the metric associated with Obj1 and Obj2
if Obj1 and Obj2 are associated with the same metric.  Otherwise, it sets
the QualityMetric pointer to NULL.  The new
ObjectiveFunction's negateFlag is always set to one, because the values
produced by obj1 and obj2 have already been multiplied by negative one
if it was needed.  Defaults to the analytical gradient.
  \param Obj1 (ObjectiveFunction*)
  \param Obj2 (ObjectiveFunction*)
 */
CompositeOFAdd::CompositeOFAdd(ObjectiveFunction* Obj1,
                               ObjectiveFunction* Obj2,
                               bool delete_OFs)
  : deleteObjFuncs(delete_OFs)
{
  objFunc1=Obj1;
  objFunc2=Obj2;
}

ObjectiveFunction* CompositeOFAdd::clone() const
  { return new CompositeOFAdd( objFunc1->clone(), objFunc2->clone(), true ); }
  
void CompositeOFAdd::clear()
{
  objFunc1->clear();
  objFunc2->clear();
}

//Michael:  need to clean up here
CompositeOFAdd::~CompositeOFAdd(){
  if (deleteObjFuncs) {
    delete objFunc1;
    delete objFunc2;
  }
}

void CompositeOFAdd::initialize_queue( MeshDomainAssoc* mesh_and_domain,
                                       const Settings* settings,
                                       MsqError& err )
{
  objFunc1->initialize_queue( mesh_and_domain, settings, err ); MSQ_ERRRTN(err);
  objFunc2->initialize_queue( mesh_and_domain, settings, err ); MSQ_ERRRTN(err);
}


bool CompositeOFAdd::initialize_block_coordinate_descent( 
                                                       MeshDomainAssoc* mesh_and_domain,
                                                       const Settings* settings,
                                                       PatchSet* user_set,
                                                       MsqError& err )
{
  bool rval1, rval2;
  rval1 = objFunc1->initialize_block_coordinate_descent( mesh_and_domain, settings, user_set, err );
  MSQ_ERRZERO(err);
  rval2 = objFunc2->initialize_block_coordinate_descent( mesh_and_domain, settings, user_set, err );
  return !MSQ_CHKERR(err) && rval1 && rval2;
}

bool CompositeOFAdd::evaluate( EvalType type, 
                               PatchData& pd,
                               double& value_out,
                               bool free,
                               MsqError& err )
{
  double value_2;
  bool ok;
  
  ok = objFunc1->evaluate( type, pd, value_out, free, err );
  if (MSQ_CHKERR(err) || !ok) return false;
  ok = objFunc2->evaluate( type, pd, value_2, free, err );
  if (MSQ_CHKERR(err) || !ok) return false;
  
  value_out += value_2;
  return true;
}

bool CompositeOFAdd::evaluate_with_gradient( EvalType type, 
                                             PatchData& pd,
                                             double& value_out,
                                             std::vector<Vector3D>& grad_out,
                                             MsqError& err )
{
  double value_2;
  bool ok;
  
  ok = objFunc1->evaluate_with_gradient( type, pd, value_out, grad_out, err );
  if (MSQ_CHKERR(err) || !ok) return false;
  ok = objFunc2->evaluate_with_gradient( type, pd, value_2, mGradient, err );
  if (MSQ_CHKERR(err) || !ok) return false;
  
  assert( grad_out.size() == pd.num_free_vertices() );
  assert( mGradient.size() == pd.num_free_vertices() );
  
  std::vector<Vector3D>::iterator i = grad_out.begin(), j = mGradient.begin();
  while (i != grad_out.end()) {
    *i += *j;
    ++i;
    ++j;
  }
  value_out += value_2;
  return true;
}

bool CompositeOFAdd::evaluate_with_Hessian_diagonal( EvalType type, 
                                            PatchData& pd,
                                            double& value_out,
                                            std::vector<Vector3D>& grad_out,
                                            std::vector<SymMatrix3D>& diag_out,
                                            MsqError& err )
{
  double value_2;
  bool valid;

  valid = objFunc1->evaluate_with_Hessian_diagonal( type, pd, value_out, grad_out, diag_out, err );
  if (MSQ_CHKERR(err) || !valid) return false;
  valid = objFunc2->evaluate_with_Hessian_diagonal( type, pd, value_2, mGradient, mDiagonal, err );
  if (MSQ_CHKERR(err) || !valid) return false;

  for (size_t i = 0; i < pd.num_free_vertices(); ++i) {
    grad_out[i] += mGradient[i];
    diag_out[i] += mDiagonal[i];
  }
  
  value_out += value_2;
  return true;
}

bool CompositeOFAdd::evaluate_with_Hessian( EvalType type, 
                                            PatchData& pd,
                                            double& value_out,
                                            std::vector<Vector3D>& grad_out,
                                            MsqHessian& Hessian_out,
                                            MsqError& err )
{
  double value_2;
  bool ok;
  
  mHessian.initialize( Hessian_out );
  
  ok = objFunc1->evaluate_with_Hessian( type, pd, value_out, grad_out, Hessian_out, err );
  if (MSQ_CHKERR(err) || !ok) return false;
  ok = objFunc2->evaluate_with_Hessian( type, pd, value_2, mGradient, mHessian, err );
  if (MSQ_CHKERR(err) || !ok) return false;
  
  value_out += value_2;
  
  assert( grad_out.size() == pd.num_free_vertices() );
  assert( mGradient.size() == pd.num_free_vertices() );
  
  for (size_t i = 0; i < pd.num_free_vertices(); ++i) 
    grad_out[i] += mGradient[i];
  Hessian_out.add( mHessian );
  return true;
}

int CompositeOFAdd::min_patch_layers() const
{
  const int v1 = objFunc1->min_patch_layers();
  const int v2 = objFunc2->min_patch_layers();
  return v1 > v2 ? v1 : v2;
}

} //namespace Mesquite
