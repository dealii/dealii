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

/** \file OFEvaluator.cpp
 *  \brief 
 *  \author Jason Kraftcheck
 */

#include "OFEvaluator.hpp"
#include "MsqError.hpp"
#include "ObjectiveFunction.hpp"

namespace MESQUITE_NS {


OFEvaluator::OFEvaluator( ObjectiveFunction* of ) : OF(of), doBCD(false)
  { }


bool OFEvaluator::initialize( MeshDomainAssoc* mesh_and_domain,
                              const Settings* settings,
                              PatchSet* user_set,
                              MsqError& err )
{
  if (doBCD) {
    tempType = ObjectiveFunction::TEMPORARY;
    firstType = ObjectiveFunction::SAVE;
    updateType = ObjectiveFunction::UPDATE;
  }
  else {
    tempType = firstType = updateType = ObjectiveFunction::CALCULATE;
  }
  reset();

  if (!doBCD) // Nash
    return true;
  
  if (!have_objective_function()) {
    MSQ_SETERR(err)("Cannot perform block coordinate descent algorithm"
                    " without an ObjectiveFunction", MsqError::INVALID_STATE);
    return false;
  }
  
  bool result = get_objective_function()->
    initialize_block_coordinate_descent( mesh_and_domain, settings, user_set, err );
  return !MSQ_CHKERR(err) && result;
}

void OFEvaluator::initialize_queue( MeshDomainAssoc* mesh_and_domain,
                                    const Settings* settings,
                                    MsqError& err )
{
  if (get_objective_function())
    get_objective_function()->initialize_queue( mesh_and_domain, settings, err );
  MSQ_ERRRTN(err);
}

bool OFEvaluator::reset() 
{
  currUpdateType = firstType;
  return true;
}

bool OFEvaluator::update( PatchData& pd, double& value, MsqError& err )
{
  if (!have_objective_function()) 
    { MSQ_SETERR(err)("No ObjectiveFunction",MsqError::INVALID_STATE); return false; }
  bool b = get_objective_function()
    ->evaluate( currUpdateType, pd, value, OF_FREE_EVALS_ONLY, err );
  currUpdateType = updateType;
  return !MSQ_CHKERR(err) && b;
}


bool OFEvaluator::update( PatchData& pd, double& value, 
                          std::vector<Vector3D>& grad,
                          MsqError& err )
{
  if (!have_objective_function()) 
    { MSQ_SETERR(err)("No ObjectiveFunction",MsqError::INVALID_STATE); return false; }
  bool b = get_objective_function()
    ->evaluate_with_gradient( currUpdateType, pd, value, grad, err );
  currUpdateType = updateType;
  return !MSQ_CHKERR(err) && b;
}


bool OFEvaluator::update( PatchData& pd, double& value, 
                          std::vector<Vector3D>& grad,
                          std::vector<SymMatrix3D>& hess_diag,
                          MsqError& err )
{
  if (!have_objective_function()) 
    { MSQ_SETERR(err)("No ObjectiveFunction",MsqError::INVALID_STATE); return false; }
  bool b = get_objective_function()
    ->evaluate_with_Hessian_diagonal( currUpdateType, pd, value, grad, hess_diag, err );
  currUpdateType = updateType;
  return !MSQ_CHKERR(err) && b;
}


bool OFEvaluator::update( PatchData& pd, double& value, 
                          std::vector<Vector3D>& grad,
                          MsqHessian& Hessian,
                          MsqError& err )
{
  if (!have_objective_function()) 
    { MSQ_SETERR(err)("No ObjectiveFunction",MsqError::INVALID_STATE); return false; }
  bool b = get_objective_function()
    ->evaluate_with_Hessian( currUpdateType, pd, value, grad, Hessian, err );
  currUpdateType = updateType;
  return !MSQ_CHKERR(err) && b;
}


bool OFEvaluator::evaluate( PatchData& pd, double& value, MsqError& err ) const
{ 
  if (!have_objective_function()) 
    { MSQ_SETERR(err)("No ObjectiveFunction",MsqError::INVALID_STATE); return false; }
  bool b = get_objective_function()->evaluate( tempType, pd, value, OF_FREE_EVALS_ONLY, err ); 
  return !MSQ_CHKERR(err) && b;
}


bool OFEvaluator::evaluate( PatchData& pd, double& value, 
                            std::vector<Vector3D>& grad,
                            MsqError& err ) const
{ 
  if (!have_objective_function()) 
    { MSQ_SETERR(err)("No ObjectiveFunction",MsqError::INVALID_STATE); return false; }
  bool b = get_objective_function()
    ->evaluate_with_gradient( tempType, pd, value, grad, err ); 
  return !MSQ_CHKERR(err) && b;
}


bool OFEvaluator::evaluate( PatchData& pd, double& value, 
                            std::vector<Vector3D>& grad,
                            std::vector<SymMatrix3D>& hess_diag,
                            MsqError& err ) const
{ 
  if (!have_objective_function()) 
    { MSQ_SETERR(err)("No ObjectiveFunction",MsqError::INVALID_STATE); return false; }
  bool b = get_objective_function()
    ->evaluate_with_Hessian_diagonal( tempType, pd, value, grad, hess_diag, err ); 
  return !MSQ_CHKERR(err) && b;
}


bool OFEvaluator::evaluate( PatchData& pd, double& value, 
                            std::vector<Vector3D>& grad,
                            MsqHessian& Hessian,
                            MsqError& err ) const
{ 
  if (!have_objective_function()) 
    { MSQ_SETERR(err)("No ObjectiveFunction",MsqError::INVALID_STATE); return false; }
  bool b = get_objective_function()
    ->evaluate_with_Hessian( tempType, pd, value, grad, Hessian, err ); 
  return !MSQ_CHKERR(err) && b;
}

} // namespace Mesquite

