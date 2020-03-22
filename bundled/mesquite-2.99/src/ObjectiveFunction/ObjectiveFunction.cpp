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
  \file   ObjectiveFunction.cpp
  \brief  

  \author Michael Brewer
  \author Thomas Leurent
  
  \date   2002-08-02
*/

#include "ObjectiveFunction.hpp"
#include "MsqVertex.hpp"
#include "MsqDebug.hpp"
#include "PatchData.hpp"
#include "MsqError.hpp"
#include "MsqHessian.hpp"
#include "SymMatrix3D.hpp"
#include <memory>  // auto_ptr

namespace MESQUITE_NS {

ObjectiveFunction::~ObjectiveFunction() {}

/*!Returns an appropiate value (eps) to use as a delta step for
  MsqVertex vertex in dimension k (i.e. k=0 -> x, k=1 -> y, k=2 -> z).
  The objective function value at the perturbed vertex position is given
  in local_val.
*/
double ObjectiveFunction::get_eps( PatchData &pd, 
                                   EvalType type,
                                   double &local_val,
                                   int dim,
                                   size_t vertex_index, 
                                   MsqError& err)
{
  double eps = 1.e-07;
  const double rho = 0.5;
  const int imax = 20;
  bool feasible = false;
  double tmp_var = 0.0;
  Vector3D delta(0,0,0);
  for (int i = 0; i < imax; ++i)
  {
    //perturb kth coord val and check feas if needed
    tmp_var = pd.vertex_by_index(vertex_index)[dim];
    delta[dim] = eps;
    pd.move_vertex( delta, vertex_index, err );
    feasible = evaluate( type, pd, local_val, OF_FREE_EVALS_ONLY, err ); MSQ_ERRZERO(err);
    //revert kth coord val
    delta = pd.vertex_by_index(vertex_index);
    delta[dim] = tmp_var;
    pd.set_vertex_coordinates( delta, vertex_index, err );
    //if step was too big, shorten it and go again
    if (feasible)
      return eps;
    else
      eps *= rho;
  }//end while looking for feasible eps
  return 0.0;
}//end function get_eps

bool ObjectiveFunction::compute_subpatch_numerical_gradient(
                                 EvalType type,
                                 EvalType subtype,
                                 PatchData& pd,
                                 double& flocal,
                                 Vector3D& grad,
                                 MsqError& err )
{
  assert( pd.num_free_vertices() == 1 );
  
  double flocald=0;
  double eps=0;
  
  bool b = evaluate( type, pd, flocal, OF_FREE_EVALS_ONLY, err );
  if(MSQ_CHKERR(err) || !b) {
    return false;
  }

    //loop over the three coords x,y,z
  for(int j=0; j<3; ++j) {
    eps = get_eps( pd, subtype, flocald, j, 0, err ); MSQ_ERRZERO(err);
    if (eps==0) {
      MSQ_SETERR(err)("Dividing by zero in Objective Functions numerical grad",
                      MsqError::INVALID_STATE);
      return false;
    }
    grad[j]=(flocald-flocal)/eps;
  }
  return true;
}

bool ObjectiveFunction::compute_patch_numerical_gradient(  EvalType type,
                                                           EvalType subtype,
                                                           PatchData& pd,
                                                           double& flocal,
                                                           std::vector<Vector3D>& grad,
                                                           MsqError& err )
{
  double flocald=0;
  double eps=0;
  
  bool b = evaluate( type, pd, flocal, OF_FREE_EVALS_ONLY, err );
  if(MSQ_CHKERR(err) || !b) {
    return false;
  }

  for (size_t i = 0; i < pd.num_free_vertices(); ++i) {
      //loop over the three coords x,y,z
    for(int j=0; j<3; ++j) {
      eps = get_eps( pd, subtype, flocald, j, i, err ); MSQ_ERRZERO(err);
      if (eps==0) {
        MSQ_SETERR(err)("Dividing by zero in Objective Functions numerical grad",
                        MsqError::INVALID_STATE);
        return false;
      }
      grad[i][j]=(flocald-flocal)/eps;
    }
  }
  
  return true;
}
  


/*! 
  Numerically Calculates the gradient of the ObjectiveFunction for the
  free vertices in the patch.  Returns 'false' if the patch is outside
  of a required feasible region, returns 'ture' otherwise.
  The behavior of the function depends on the value of the boolean
  useLocalGradient.  If useLocalGradient is set to
  'true', compute_numerical_gradient creates a sub-patch around a free
  vertex, and then perturbs that vertex in one of the coordinate directions.
  Only the ObjectiveFunction value on the local sub-patch is used in the
  computation of the gradient.  Therefore, useLocalGradient should only
  be set to 'true' for ObjectiveFunctions which can use this method.  Unless
  the concrete ObjectiveFunction sets useLocalGradient to 'true' in its
  constructor, the value will be 'false'.  In this case, the objective
  function value for the entire patch is used in the calculation of the
  gradient.  This is computationally expensive, but it is numerically
  correct for all (C_1) functions.
  \param pd  PatchData on which the gradient is taken.
  \param grad  Array of Vector3D of length the number of vertices used to store gradient.
  \param OF_val will be set to the objective function value.
 */
bool ObjectiveFunction::evaluate_with_gradient( EvalType eval_type,
                                                PatchData &pd,
                                                double& OF_val,
                                                std::vector<Vector3D>& grad,
                                                MsqError &err )
{
  bool b;
  grad.resize( pd.num_free_vertices() );
  
    // Fast path for single-free-vertex patch
  if (pd.num_free_vertices() == 1) {
    const EvalType sub_type = (eval_type == CALCULATE) ? CALCULATE : TEMPORARY;
    b = compute_subpatch_numerical_gradient( eval_type, sub_type, pd, OF_val, grad[0], err );
    return !MSQ_CHKERR(err) && b;
  }
  
  ObjectiveFunction* of = this;
  std::auto_ptr<ObjectiveFunction> deleter;
  if (eval_type == CALCULATE) {
    of->clear();
    b = of->evaluate( ACCUMULATE, pd, OF_val, OF_FREE_EVALS_ONLY, err );
    if (err) { // OF doesn't support BCD type evals, try slow method
      err.clear();
      of->clear();
      b = compute_patch_numerical_gradient( CALCULATE, CALCULATE, pd, OF_val, grad, err );
      return !MSQ_CHKERR(err) && b;
    }
    else if (!b)
      return b;
  } 
  else {
    b = this->evaluate( eval_type, pd, OF_val, OF_FREE_EVALS_ONLY, err );
    if (MSQ_CHKERR(err) || !b)
      return false;
    of = this->clone();
    deleter = std::auto_ptr<ObjectiveFunction>(of);
  }

    // Determine number of layers of adjacent elements based on metric type.
  unsigned layers = min_patch_layers();
  
    // Create a subpatch for each free vertex and use it to evaluate the
    // gradient for that vertex.
  double flocal;
  PatchData subpatch;
  for (size_t i = 0; i < pd.num_free_vertices(); ++i)
  {
    pd.get_subpatch( i, layers, subpatch, err ); MSQ_ERRZERO(err);
    b = of->compute_subpatch_numerical_gradient( SAVE, TEMPORARY, subpatch, flocal, grad[i], err );
    if (MSQ_CHKERR(err) || !b) {
      of->clear();
      return false;
    }
  }
  
  of->clear();
  return true;
}

bool ObjectiveFunction::evaluate_with_Hessian_diagonal( EvalType type, 
                                        PatchData& pd,
                                        double& value_out,
                                        std::vector<Vector3D>& grad_out,
                                        std::vector<SymMatrix3D>& hess_diag_out,
                                        MsqError& err )
{
  MsqHessian hess;
  hess.initialize( pd, err ); MSQ_ERRZERO(err);
  bool val = evaluate_with_Hessian( type, pd, value_out, grad_out, hess, err );
  MSQ_ERRZERO(err);
  hess_diag_out.resize( hess.size() );
  for (size_t i = 0; i < hess.size(); ++i)
    hess_diag_out[i] = hess.get_block(i,i)->upper();
  return val;
}

bool ObjectiveFunction::evaluate_with_Hessian( EvalType type, 
                                               PatchData& pd,
                                               double& value_out,
                                               std::vector<Vector3D>& grad_out,
                                               MsqHessian& Hessian_out,
                                               MsqError& err ) 
{
      MSQ_SETERR(err)("No Hessian available for this objective function.\n"
                      "Choose either a different objective function or a "
                      "different solver.\n",
                    MsqError::INVALID_STATE);
      return false;
}



} // namespace Mesquite

