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
  \file    CompositeOFScalarMultiply.cpp
  \brief  

  This Objective Function combines two Objective Functions by mulitplication
  \author Michael Brewer
  \date   2002-01-23
*/
#include <math.h>
#include "CompositeOFScalarMultiply.hpp"
#include "MsqTimer.hpp"
#include "MsqError.hpp"
#include "MsqDebug.hpp"
#include "MsqHessian.hpp"

namespace MESQUITE_NS {


/*!
Sets the QualityMetric pointer to the metric associated with Obj. 
However, if alp is less than zero, the new
ObjectiveFunction's negateFlag is the opposite of Obj's.  This objective
function defaults to the analytical gradient.
  \param alp (double)
  \param Obj (ObjectiveFunction*)
 */
CompositeOFScalarMultiply::CompositeOFScalarMultiply(double alp, 
                                           ObjectiveFunction* Obj,
                                           bool delete_of)
  : deleteObjFunc(delete_of)
{
  objFunc=Obj;
  mAlpha=alp;
}

//Michael:  need to clean up here
CompositeOFScalarMultiply::~CompositeOFScalarMultiply(){
  if (deleteObjFunc)
    delete objFunc;
}

ObjectiveFunction* CompositeOFScalarMultiply::clone() const
  { return new CompositeOFScalarMultiply( mAlpha, objFunc->clone(), true ); }
  
void CompositeOFScalarMultiply::clear()
{
  objFunc->clear();
}

void CompositeOFScalarMultiply::initialize_queue( MeshDomainAssoc* mesh_and_domain,
                                                  const Settings* settings,
                                                  MsqError& err )
{
  objFunc->initialize_queue( mesh_and_domain, settings, err ); MSQ_ERRRTN(err);
}

bool CompositeOFScalarMultiply::initialize_block_coordinate_descent( 
                                                       MeshDomainAssoc* mesh_and_domain,
                                                       const Settings* settings,
                                                       PatchSet* user_set,
                                                       MsqError& err )
{
  bool rval = objFunc->initialize_block_coordinate_descent( mesh_and_domain, settings, user_set, err );
  return !MSQ_CHKERR(err) && rval;
}

bool CompositeOFScalarMultiply::evaluate( EvalType type, 
                                          PatchData& pd,
                                          double& value_out,
                                          bool free,
                                          MsqError& err )
{
  bool ok = objFunc->evaluate( type, pd, value_out, free, err );
  value_out *= mAlpha;
  return !MSQ_CHKERR(err) && ok;
}

bool CompositeOFScalarMultiply::evaluate_with_gradient( EvalType type, 
                                             PatchData& pd,
                                             double& value_out,
                                             std::vector<Vector3D>& grad_out,
                                             MsqError& err )
{
  bool ok = objFunc->evaluate_with_gradient( type, pd, value_out, grad_out, err );
  value_out *= mAlpha;
  for (std::vector<Vector3D>::iterator i = grad_out.begin(); i != grad_out.end(); ++i)
    *i *= mAlpha;
  return !MSQ_CHKERR(err) && ok;
}

bool CompositeOFScalarMultiply::evaluate_with_Hessian_diagonal( EvalType type, 
                                            PatchData& pd,
                                            double& value_out,
                                            std::vector<Vector3D>& grad_out,
                                            std::vector<SymMatrix3D>& diag_out,
                                            MsqError& err )
{
  bool ok = objFunc->evaluate_with_Hessian_diagonal( type, pd, value_out, grad_out, diag_out, err );
  value_out *= mAlpha;
  for (size_t i = 0; i < pd.num_free_vertices(); ++i) {
    grad_out[i] *= mAlpha;
    diag_out[i] *= mAlpha;
  }
  return !MSQ_CHKERR(err) && ok;
}

bool CompositeOFScalarMultiply::evaluate_with_Hessian( EvalType type, 
                                            PatchData& pd,
                                            double& value_out,
                                            std::vector<Vector3D>& grad_out,
                                            MsqHessian& Hessian_out,
                                            MsqError& err )
{
  bool ok = objFunc->evaluate_with_Hessian( type, pd, value_out, grad_out, Hessian_out, err );
  value_out *= mAlpha;
  for (std::vector<Vector3D>::iterator i = grad_out.begin(); i != grad_out.end(); ++i)
    *i *= mAlpha;
  Hessian_out.scale( mAlpha );
  return !MSQ_CHKERR(err) && ok;
}

int CompositeOFScalarMultiply::min_patch_layers() const
  { return objFunc->min_patch_layers(); }

} // namespace Mesquite
