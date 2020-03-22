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
  \file    CompositeOFMultiply.cpp
  \brief  

  This Objective Function combines two Objective Functions by mulitplication
  \author Michael Brewer
  \date   2002-01-23
*/
#include <math.h>
#include "CompositeOFMultiply.hpp"
#include "MsqTimer.hpp"
#include "PatchData.hpp"

namespace MESQUITE_NS {


/*!
Sets the QualityMetric pointer to the metric associated with Obj1 and Obj2
if Obj1 and Obj2 are associated with the same metric.  Otherwise, it sets
the QualityMetric pointer to NULL.  The new
ObjectiveFunction's negateFlag is set to negative one only if both Obj1 and
Obj2's negateFlag are negative one (because obj1 and obj2's evaluate function
multiply their return values by negative one if their respective
function needs to be maximized.  If both of these functions needed to
be maximized, then the negative ones will have cancelled out).  Otherwise,
the negateFlag is set to one.  Defaults to the analytical gradient.
  \param Obj1 (ObjectiveFunction*)
  \param Obj2 (ObjectiveFunction*)
 */
CompositeOFMultiply::CompositeOFMultiply( ObjectiveFunction* Obj1, 
                                          ObjectiveFunction* Obj2,
                                          bool delete_OFs)
  : deleteObjFuncs(delete_OFs)
{
  objFunc1=Obj1;
  objFunc2=Obj2;
}

//Michael:  need to clean up here
CompositeOFMultiply::~CompositeOFMultiply(){
  if (deleteObjFuncs) {
    delete objFunc1;
    delete objFunc2;
  }
}


ObjectiveFunction* CompositeOFMultiply::clone() const
  { return new CompositeOFMultiply( objFunc1->clone(), objFunc2->clone(), true ); }
  
void CompositeOFMultiply::clear()
{
  objFunc1->clear();
  objFunc2->clear();
}

void CompositeOFMultiply::initialize_queue( MeshDomainAssoc* mesh_and_domain,
                                            const Settings* settings,
                                            MsqError& err )
{
  objFunc1->initialize_queue( mesh_and_domain, settings, err ); MSQ_ERRRTN(err);
  objFunc2->initialize_queue( mesh_and_domain, settings, err ); MSQ_ERRRTN(err);
}

bool CompositeOFMultiply::initialize_block_coordinate_descent( 
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

bool CompositeOFMultiply::evaluate( EvalType type, 
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
  
  value_out *= value_2;
  return true;
}

bool CompositeOFMultiply::evaluate_with_gradient( EvalType type, 
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
    *i *= value_2;
    *j *= value_out;
    *i += *j;
    ++i;
    ++j;
  }
  value_out *= value_2;
  return true;
}

bool CompositeOFMultiply::evaluate_with_Hessian_diagonal( EvalType type, 
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
    diag_out[i] *= value_2;
    mDiagonal[i] *= value_out;
    diag_out[i] += mDiagonal[i];
    diag_out[i] += outer_plus_transpose( grad_out[i], mGradient[i] );

    grad_out[i] *= value_2;
    mGradient[i] *= value_out;
    grad_out[i] += mGradient[i];
  }
  
  value_out *= value_2;
  return true;
}

bool CompositeOFMultiply::evaluate_with_Hessian( EvalType , 
                                            PatchData& ,
                                            double& ,
                                            std::vector<Vector3D>& ,
                                            MsqHessian& ,
                                            MsqError& err )
{
  MSQ_SETERR(err)("Mesquite is not capable of representing the dense "
                  "Hessian of the product of two objective fuctions. "
                  "Either choose a solver that does not require the "
                  "Hessian of the objective function or do not use the "
                  "CompositeOFMultiple objective function .",
                  MsqError::INVALID_STATE );
  return false;
}
	

int CompositeOFMultiply::min_patch_layers() const
{
  const int v1 = objFunc1->min_patch_layers();
  const int v2 = objFunc2->min_patch_layers();
  return v1 > v2 ? v1 : v2;
}
	
} // namespace Mesquite
