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
  \file   SteepestDescent.cpp
  \brief  

  Implements the SteepestDescent class member functions.
  
  \author Thomas Leurent
  \date   2002-06-13
*/

#include "SteepestDescent.hpp"
#include "MsqFreeVertexIndexIterator.hpp"
#include "MsqTimer.hpp"
#include "MsqDebug.hpp"

#include <memory>

namespace MESQUITE_NS {

std::string SteepestDescent::get_name() const
  { return "SteepestDescent"; }
  
PatchSet* SteepestDescent::get_patch_set()
  { return PatchSetUser::get_patch_set(); }

SteepestDescent::SteepestDescent(ObjectiveFunction* of) 
  : VertexMover(of),
    PatchSetUser(true),
    projectGradient(false) //,
    //cosineStep(false)
{
}  
  

void SteepestDescent::initialize(PatchData &/*pd*/, MsqError &/*err*/)
{
}

void SteepestDescent::initialize_mesh_iteration(PatchData &/*pd*/, MsqError &/*err*/)
{
}

void SteepestDescent::optimize_vertex_positions(PatchData &pd, 
                                                MsqError &err)
{
  MSQ_FUNCTION_TIMER( "SteepestDescent::optimize_vertex_positions" );

  const int SEARCH_MAX = 100;
  const double c1 = 1e-4;
  //std::vector<Vector3D> unprojected(pd.num_free_vertices()); 
  std::vector<Vector3D> gradient(pd.num_free_vertices()); 
  bool feasible=true;//bool for OF values
  double min_edge_len, max_edge_len;
  double step_size=0, original_value=0, new_value=0;
  double norm_squared=0;
  PatchDataVerticesMemento* pd_previous_coords;
  TerminationCriterion* term_crit=get_inner_termination_criterion();
  OFEvaluator& obj_func = get_objective_function_evaluator();
  
    // get vertex memento so we can restore vertex coordinates for bad steps.
  pd_previous_coords = pd.create_vertices_memento( err ); MSQ_ERRRTN(err);
    // use auto_ptr to automatically delete memento when we exit this function
  std::auto_ptr<PatchDataVerticesMemento> memento_deleter( pd_previous_coords );

    // Evaluate objective function.
    //
    // Always use 'update' version when beginning optimization so that
    // if doing block coordinate descent the OF code knows the set of
    // vertices we are modifying during the optimziation (the subset
    // of the mesh contained in the current patch.)  This has to be
    // done up-front because typically an OF will just store the portion
    // of the OF value (e.g. the numeric contribution to the sum for an
    // averaging OF) for the initial patch.
  feasible = obj_func.update( pd, original_value, gradient, err ); MSQ_ERRRTN(err);
    // calculate gradient dotted with itself
  norm_squared = length_squared( gradient );
  
    //set an error if initial patch is invalid.
  if(!feasible){
    MSQ_SETERR(err)("SteepestDescent passed invalid initial patch.",
                    MsqError::INVALID_ARG);
    return;
  }

    // use edge length as an initial guess for for step size
  pd.get_minmax_edge_length( min_edge_len, max_edge_len );
  //step_size = max_edge_len / std::sqrt(norm_squared);
  //if (!finite(step_size))  // zero-length gradient
  //  return;
//  if (norm_squared < DBL_EPSILON)
//    return;
  if (norm_squared >= DBL_EPSILON)
    step_size = max_edge_len / std::sqrt(norm_squared) * pd.num_free_vertices();

    // The steepest descent loop...
    // We loop until the user-specified termination criteria are met.
  while (!term_crit->terminate()) {
    MSQ_DBGOUT(3) << "Iteration " << term_crit->get_iteration_count() << std::endl;
    MSQ_DBGOUT(3) << "  o  original_value: " << original_value << std::endl;
    MSQ_DBGOUT(3) << "  o  grad norm suqared: " << norm_squared << std::endl;

      // Save current vertex coords so that they can be restored if
      // the step was bad.
    pd.recreate_vertices_memento( pd_previous_coords, err ); MSQ_ERRRTN(err);

      // Reduce step size until it satisfies Armijo condition
    int counter = 0;
    for (;;) {
      if (++counter > SEARCH_MAX || step_size < DBL_EPSILON) {
        MSQ_DBGOUT(3) << "    o  No valid step found.  Giving Up." << std::endl;
        return;
      }
      
      // Move vertices to new positions.
      // Note: step direction is -gradient so we pass +gradient and 
      //       -step_size to achieve the same thing.
      pd.move_free_vertices_constrained( arrptr(gradient), gradient.size(), -step_size, err ); MSQ_ERRRTN(err);
      // Evaluate objective function for new vertices.  We call the
      // 'evaluate' form here because we aren't sure yet if we want to
      // keep these vertices.  Until we call 'update', we have the option
      // of reverting a block coordinate decent objective function's state
      // to that of the initial vertex coordinates.  However, for block
      // coordinate decent to work correctly, we will need to call an
      // 'update' form if we decide to keep the new vertex coordinates.
      feasible = obj_func.evaluate( pd, new_value, err ); 
      if (err.error_code() == err.BARRIER_VIOLATED) 
        err.clear();  // barrier violated does not represent an actual error here
      MSQ_ERRRTN(err);
      MSQ_DBGOUT(3) << "    o  step_size: " << step_size << std::endl;
      MSQ_DBGOUT(3) << "    o  new_value: " << new_value << std::endl;

      if (!feasible) {
        // OF value is invalid, decrease step_size a lot
        step_size *= 0.2;
      }
      else if (new_value > original_value - c1 * step_size * norm_squared) {
        // Armijo condition not met.
        step_size *= 0.5;
      }
      else {
        // Armijo condition met, stop
        break;
      }
      
      // undo previous step : restore vertex coordinates
      pd.set_to_vertices_memento( pd_previous_coords, err );  MSQ_ERRRTN(err);
    }
   
      // Re-evaluate objective function to get gradient.
      // Calling the 'update' form here incorporates the new vertex 
      // positions into the 'accumulated' value if we are doing a 
      // block coordinate descent optimization.
    obj_func.update(pd, original_value, gradient, err ); MSQ_ERRRTN(err);
    if (projectGradient) {
      //if (cosineStep) {
      //  unprojected = gradient;
      //  pd.project_gradient( gradient, err ); MSQ_ERRRTN(err);
      //  double dot = inner_product( arrptr(gradient), arrptr(unprojected), gradient.size() );
      //  double lensqr1 = length_squared( gradient );
      //  double lensqr2 = length_squared( unprojected );
      //  double cossqr = dot * dot / lensqr1 / lensqr2;
      //  step_size *= sqrt(cossqr);
      //}
      //else {
        pd.project_gradient( gradient, err ); MSQ_ERRRTN(err);
      //}      
    }
      
      // Update terination criterion for next iteration.
      // This is necessary for efficiency.  Some values can be adjusted
      // for each iteration so we don't need to re-caculate the value
      // over the entire mesh.
    term_crit->accumulate_patch( pd, err );  MSQ_ERRRTN(err);
    term_crit->accumulate_inner( pd, original_value, arrptr(gradient), err ); MSQ_ERRRTN(err); 
      
      // Calculate initial step size for next iteration using step size 
      // from this iteration
    step_size *= norm_squared;
    norm_squared = length_squared( gradient );
//    if (norm_squared < DBL_EPSILON)
//      break;
  if (norm_squared >= DBL_EPSILON)
    step_size /= norm_squared;
  }
}

void SteepestDescent::terminate_mesh_iteration(PatchData &/*pd*/, MsqError &/*err*/)
{
  //  cout << "- Executing SteepestDescent::iteration_complete()\n";
}
  
void SteepestDescent::cleanup()
{
  //  cout << "- Executing SteepestDescent::iteration_end()\n";
}
  
} // namespace Mesquite
