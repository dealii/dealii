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
  \file    CompositeOFScalarAdd.cpp
  \brief  

  This Objective Function combines two Objective Functions by addition
  \author Michael Brewer
  \date   2002-06-24
*/
#include <math.h>
#include "CompositeOFScalarAdd.hpp"
#include "PatchData.hpp"
#include "MsqTimer.hpp"
using namespace Mesquite;


/*!
Sets the QualityMetric pointer to the metric associated with Obj1.  
The new ObjectiveFunction's negateFlag is also the
same as that of Obj1.  This objective function defaults to the analytical
gradient which essentially just calls Obj1's gradient function.
  \param alp (double)
  \param Obj1 (ObjectiveFunction*)
 */
CompositeOFScalarAdd::CompositeOFScalarAdd( double alp, 
                                            ObjectiveFunction* Obj1,
                                            bool delete_OF)
  : deleteObjFunc(delete_OF)
{
  objFunc=Obj1;
  mAlpha=alp;
}


//Michael:  need to clean up here
CompositeOFScalarAdd::~CompositeOFScalarAdd(){
  if (deleteObjFunc)
    delete objFunc;
}

ObjectiveFunction* CompositeOFScalarAdd::clone() const
  { return new CompositeOFScalarAdd( mAlpha, objFunc->clone(), true ); }
  
void CompositeOFScalarAdd::clear()
{
  objFunc->clear();
}


void CompositeOFScalarAdd::initialize_queue( MeshDomainAssoc* mesh_and_domain,
                                             const Settings* settings,
                                             MsqError& err )
{
  objFunc->initialize_queue( mesh_and_domain, settings, err ); MSQ_ERRRTN(err);
}

bool CompositeOFScalarAdd::initialize_block_coordinate_descent( 
                                                       MeshDomainAssoc* mesh_and_domain,
                                                       const Settings* settings,
                                                       PatchSet* user_set,
                                                       MsqError& err )
{
  bool rval = objFunc->initialize_block_coordinate_descent( mesh_and_domain, settings, user_set, err );
  return !MSQ_CHKERR(err) && rval;
}

bool CompositeOFScalarAdd::evaluate( EvalType type, 
                                     PatchData& pd,
                                     double& value_out,
                                     bool free,
                                     MsqError& err )
{
  bool ok = objFunc->evaluate( type, pd, value_out, free, err );
  value_out += mAlpha;
  return !MSQ_CHKERR(err) && ok;
}

bool CompositeOFScalarAdd::evaluate_with_gradient( EvalType type, 
                                             PatchData& pd,
                                             double& value_out,
                                             std::vector<Vector3D>& grad_out,
                                             MsqError& err )
{
  bool ok = objFunc->evaluate_with_gradient( type, pd, value_out, grad_out, err );
  value_out += mAlpha;
  return !MSQ_CHKERR(err) && ok;
}

bool CompositeOFScalarAdd::evaluate_with_Hessian_diagonal( EvalType type, 
                                            PatchData& pd,
                                            double& value_out,
                                            std::vector<Vector3D>& grad_out,
                                            std::vector<SymMatrix3D>& diag_out,
                                            MsqError& err )
{
  bool ok = objFunc->evaluate_with_Hessian_diagonal( type, pd, value_out, grad_out, diag_out, err );
  value_out += mAlpha;
  return !MSQ_CHKERR(err) && ok;
}

bool CompositeOFScalarAdd::evaluate_with_Hessian( EvalType type, 
                                            PatchData& pd,
                                            double& value_out,
                                            std::vector<Vector3D>& grad_out,
                                            MsqHessian& Hessian_out,
                                            MsqError& err )
{
  bool ok = objFunc->evaluate_with_Hessian( type, pd, value_out, grad_out, Hessian_out, err );
  value_out += mAlpha;
  return !MSQ_CHKERR(err) && ok;
}

int CompositeOFScalarAdd::min_patch_layers() const
  { return objFunc->min_patch_layers(); }
