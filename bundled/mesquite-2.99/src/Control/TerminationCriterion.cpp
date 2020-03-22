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
// -*- Mode : c++; tab-width: 2; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 2 -*-
//
// DESCRIPTION:
// ============
/*! \file TerminationCriterion.cpp
  
    \brief  Member functions of the Mesquite::TerminationCriterion class

    \author Michael Brewer
    \author Thomas Leurent
    \date   Feb. 14, 2003
 */

#include "TerminationCriterion.hpp"
#include "MsqVertex.hpp"
#include "MsqInterrupt.hpp"
#include "OFEvaluator.hpp"
#include "MsqError.hpp"
#include "MsqDebug.hpp"
#include "PatchData.hpp"
#include "MeshWriter.hpp"
#include "MeshUtil.hpp"
#include "SimpleStats.hpp"

#include <sstream>
#include <set>

namespace MESQUITE_NS {

extern int get_parallel_rank();
extern int get_parallel_size();
extern double reduce_parallel_max(double value);

#define MSQ_DBGOUT_P0_ONLY(flag) if (!get_parallel_rank()) MSQ_DBGOUT(flag)

#define RPM(val) val

// this causes race conditions - don't use it
#define RPM1(val) reduce_parallel_max(val)

 /*! \enum TCType  defines the termination criterion */
enum TCType {
   NONE    = 0,
   //! checks the gradient \f$\nabla f \f$ of objective function 
   //! \f$f : I\!\!R^{3N} \rightarrow I\!\!R \f$ against a double \f$d\f$  
   //! and stops when \f$\sqrt{\sum_{i=1}^{3N}\nabla f_i^2}<d\f$  
   GRADIENT_L2_NORM_ABSOLUTE = 1<<0,  
   //! checks the gradient \f$\nabla f \f$ of objective function 
   //! \f$f : I\!\!R^{3N} \rightarrow I\!\!R \f$ against a double \f$d\f$  
   //! and stops when \f$ \max_{i=1}^{3N} \nabla f_i < d \f$  
   GRADIENT_INF_NORM_ABSOLUTE = 1<<1,
     //!terminates on the j_th iteration when
     //! \f$\sqrt{\sum_{i=1}^{3N}\nabla f_{i,j}^2}<d\sqrt{\sum_{i=1}^{3N}\nabla f_{i,0}^2}\f$
     //!  That is, terminates when the norm of the gradient is small
     //! than some scaling factor times the norm of the original gradient. 
   GRADIENT_L2_NORM_RELATIVE = 1<<2,
   //!terminates on the j_th iteration when
     //! \f$\max_{i=1 \cdots 3N}\nabla f_{i,j}<d \max_{i=1 \cdots 3N}\nabla f_{i,0}\f$
     //!  That is, terminates when the norm of the gradient is small
     //! than some scaling factor times the norm of the original gradient.
     //! (Using the infinity norm.)
   GRADIENT_INF_NORM_RELATIVE = 1<<3,
     //! Not yet implemented.
   KKT  = 1<<4,
     //!Terminates when the objective function value is smaller than
     //! the given scalar value.
   QUALITY_IMPROVEMENT_ABSOLUTE = 1<<5,
     //!Terminates when the objective function value is smaller than
     //! the given scalar value times the original objective function
     //! value.
   QUALITY_IMPROVEMENT_RELATIVE = 1<<6,
     //!Terminates when the number of iterations exceeds a given integer.
   NUMBER_OF_ITERATES = 1<<7,
     //!Terminates when the algorithm exceeds an allotted time limit
     //! (given in seconds).
   CPU_TIME  = 1<<8,
     //!Terminates when a the maximum distance moved by any vertex
     //! during the previous iteration is below the given value.
   VERTEX_MOVEMENT_ABSOLUTE  = 1<<9,
     //!Terminates when a the maximum distance moved by any vertex
     //! during the previous iteration is below a value calculated
     //! from the edge lengths of the initial mesh.
   VERTEX_MOVEMENT_ABS_EDGE_LENGTH  = 1<<10,
     //!Terminates when a the maximum distance moved by any vertex
     //! during the previous iteration is below the given value
     //! times the maximum distance moved by any vertex over the
     //! entire course of the optimization.
   VERTEX_MOVEMENT_RELATIVE  = 1<<11,
     //!Terminates when the decrease in the objective function value since
     //! the previous iteration is below the given value.
   SUCCESSIVE_IMPROVEMENTS_ABSOLUTE = 1<<12,
     //!Terminates when the decrease in the objective function value since
     //! the previous iteration is below the given value times the
     //! decrease in the objective function value since the beginning
     //! of this optimization process.
   SUCCESSIVE_IMPROVEMENTS_RELATIVE = 1<<13,
     //!Terminates when any vertex leaves the bounding box, defined
     //! by the given value, d.  That is, when the absolute value of
     //! a single coordinate of vertex's position exceeds d.
   BOUNDED_VERTEX_MOVEMENT = 1<<14,
    //! Terminate when no elements are inverted
   UNTANGLED_MESH = 1<<15
};


const unsigned long GRAD_FLAGS = GRADIENT_L2_NORM_ABSOLUTE |
                                 GRADIENT_INF_NORM_ABSOLUTE |
                                 GRADIENT_L2_NORM_RELATIVE |
                                 GRADIENT_INF_NORM_RELATIVE;
const unsigned long OF_FLAGS   = QUALITY_IMPROVEMENT_ABSOLUTE |
                                 QUALITY_IMPROVEMENT_RELATIVE |
                                 SUCCESSIVE_IMPROVEMENTS_ABSOLUTE |
                                 SUCCESSIVE_IMPROVEMENTS_RELATIVE;

const unsigned long MOVEMENT_FLAGS = VERTEX_MOVEMENT_ABSOLUTE |
                                     VERTEX_MOVEMENT_ABS_EDGE_LENGTH |
                                     VERTEX_MOVEMENT_RELATIVE;

/*!Constructor initializes all of the data members which are not
  necessarily automatically initialized in their constructors.*/
  TerminationCriterion::TerminationCriterion(std::string name, InnerOuterType innerOuterType)
  : mGrad(8),
    initialVerticesMemento(0),
    previousVerticesMemento(0),
    debugLevel(2),
    timeStepFileType(NOTYPE), moniker(name), innerOuterType(innerOuterType)
{
  terminationCriterionFlag=NONE;
  cullingMethodFlag=NONE;
  cullingEps=0.0;
  cullingGlobalPatch=false;
  initialOFValue=0.0;
  previousOFValue=0.0;
  currentOFValue = 0.0;
  lowerOFBound=0.0;
  initialGradL2NormSquared=0.0;
  initialGradInfNorm=0.0;
    //initial size of the gradient array
  gradL2NormAbsoluteEpsSquared=0.0;
  gradL2NormRelativeEpsSquared=0.0;
  gradInfNormAbsoluteEps=0.0;
  gradInfNormRelativeEps=0.0;
  qualityImprovementAbsoluteEps=0.0;
  qualityImprovementRelativeEps=0.0;
  iterationBound=0;
  iterationCounter=0;
  timeBound=0.0;
  vertexMovementAbsoluteEps=0.0;
  vertexMovementRelativeEps=0.0;
  vertexMovementAvgBeta = -1;
  vertexMovementAbsoluteAvgEdge = -1; 
  successiveImprovementsAbsoluteEps=0.0;
  successiveImprovementsRelativeEps=0.0;
  boundedVertexMovementEps=0.0;
  globalInvertedCount = 0;
  patchInvertedCount = 0;
  
}

std::string TerminationCriterion::par_string() 
{
  if (get_parallel_size())
    {
      std::ostringstream str;
      str << "P[" << get_parallel_rank() << "] " + moniker + " ";
      return str.str();
    }
  return moniker + " ";
}

void TerminationCriterion::add_absolute_gradient_L2_norm( double eps )
{
       terminationCriterionFlag|=GRADIENT_L2_NORM_ABSOLUTE;
       gradL2NormAbsoluteEpsSquared=eps*eps;
}  

void TerminationCriterion::add_absolute_gradient_inf_norm( double eps )
{
       terminationCriterionFlag|=GRADIENT_INF_NORM_ABSOLUTE;
       gradInfNormAbsoluteEps=eps;
}
    
void TerminationCriterion::add_relative_gradient_L2_norm( double eps )
{
       terminationCriterionFlag|=GRADIENT_L2_NORM_RELATIVE;
       gradL2NormRelativeEpsSquared=eps*eps;
}

void TerminationCriterion::add_relative_gradient_inf_norm( double eps )
{
       terminationCriterionFlag|=GRADIENT_INF_NORM_RELATIVE;
       gradInfNormRelativeEps=eps;
}

void TerminationCriterion::add_absolute_quality_improvement( double eps )
{
       terminationCriterionFlag|=QUALITY_IMPROVEMENT_ABSOLUTE;
       qualityImprovementAbsoluteEps=eps;
}

void TerminationCriterion::add_relative_quality_improvement( double eps )
{
       terminationCriterionFlag|=QUALITY_IMPROVEMENT_RELATIVE;
       qualityImprovementRelativeEps=eps;
}

void TerminationCriterion::add_absolute_vertex_movement( double eps )
{
       terminationCriterionFlag|=VERTEX_MOVEMENT_ABSOLUTE;
         //we actually compare squared movement to squared epsilon
       vertexMovementAbsoluteEps=(eps*eps);
}

void TerminationCriterion::add_absolute_vertex_movement_edge_length( double beta )
{
       terminationCriterionFlag|=VERTEX_MOVEMENT_ABS_EDGE_LENGTH;
       vertexMovementAvgBeta=beta;
}

void TerminationCriterion::add_relative_vertex_movement( double eps )
{
       terminationCriterionFlag|=VERTEX_MOVEMENT_RELATIVE;
         //we actually compare squared movement to squared epsilon
       vertexMovementRelativeEps=(eps*eps);
}

void TerminationCriterion::add_absolute_successive_improvement( double eps )
{
       terminationCriterionFlag|=SUCCESSIVE_IMPROVEMENTS_ABSOLUTE;
       successiveImprovementsAbsoluteEps=eps;
}

void TerminationCriterion::add_relative_successive_improvement( double eps )
{
       terminationCriterionFlag|=SUCCESSIVE_IMPROVEMENTS_RELATIVE;
       successiveImprovementsRelativeEps=eps;
}
    
void TerminationCriterion::add_cpu_time( double seconds )
{
       terminationCriterionFlag|=CPU_TIME;
       timeBound=seconds;
}
    
void TerminationCriterion::add_iteration_limit( unsigned int max_iterations )
{
       terminationCriterionFlag|=NUMBER_OF_ITERATES;
       iterationBound=max_iterations;
}
    
void TerminationCriterion::add_bounded_vertex_movement( double eps )
{
       terminationCriterionFlag|=BOUNDED_VERTEX_MOVEMENT;
       boundedVertexMovementEps=eps;
}

void TerminationCriterion::add_untangled_mesh()
{
  terminationCriterionFlag |= UNTANGLED_MESH;
}
    
void TerminationCriterion::remove_all_criteria()
{
  terminationCriterionFlag=0;
}
    
void TerminationCriterion::cull_on_absolute_quality_improvement( double limit )
{
       cullingMethodFlag=QUALITY_IMPROVEMENT_ABSOLUTE;
       cullingEps=limit;
}
void TerminationCriterion::cull_on_relative_quality_improvement( double limit )
{
       cullingMethodFlag=QUALITY_IMPROVEMENT_RELATIVE;
       cullingEps=limit;
}
void TerminationCriterion::cull_on_absolute_vertex_movement( double limit )
{
       cullingMethodFlag=VERTEX_MOVEMENT_ABSOLUTE;
       cullingEps=limit;
}
void TerminationCriterion::cull_on_relative_vertex_movement( double limit )
{
       cullingMethodFlag=VERTEX_MOVEMENT_RELATIVE;
       cullingEps=limit;
}
void TerminationCriterion::cull_on_absolute_vertex_movement_edge_length( double limit )
{
       cullingMethodFlag=VERTEX_MOVEMENT_ABS_EDGE_LENGTH;
       vertexMovementAvgBeta=limit;
}
void TerminationCriterion::cull_on_absolute_successive_improvement( double limit )
{
       cullingMethodFlag=SUCCESSIVE_IMPROVEMENTS_ABSOLUTE;
       cullingEps=limit;
}
void TerminationCriterion::cull_on_relative_successive_improvement( double limit )
{
       cullingMethodFlag=SUCCESSIVE_IMPROVEMENTS_RELATIVE;
       cullingEps=limit;
}
void TerminationCriterion::cull_untangled_mesh( )
{
  cullingMethodFlag = UNTANGLED_MESH;
}
    
    
void TerminationCriterion::remove_culling()
{
  cullingMethodFlag=NONE;
}

void TerminationCriterion::cull_for_global_patch(bool val)
{
  cullingGlobalPatch = val;
}



/*!This version of reset is called using a MeshSet, which implies
  it is only called when this criterion is used as the 'outer' termination
  criterion.  
 */
void TerminationCriterion::reset_outer(Mesh* mesh, 
                                       MeshDomain* domain,
                                       OFEvaluator& obj_eval,
                                       const Settings* settings,
                                       MsqError &err)
{
  const unsigned long totalFlag = terminationCriterionFlag | cullingMethodFlag;
  PatchData global_patch;
  if (settings)
    global_patch.attach_settings( settings );
  
    //if we need to fill out the global patch data object.
  if ((totalFlag & (GRAD_FLAGS | OF_FLAGS | VERTEX_MOVEMENT_RELATIVE | UNTANGLED_MESH))
     || timeStepFileType)
  {
    global_patch.set_mesh( mesh );
    global_patch.set_domain( domain );
    global_patch.fill_global_patch( err );
    MSQ_ERRRTN(err);
  }

    //now call the other reset
  reset_inner( global_patch, obj_eval, err ); MSQ_ERRRTN(err);
}
    
/*!Reset function using using a PatchData object.  This function is
  called for the inner-stopping criterion directly from the
  loop over mesh function in VertexMover.  For outer criterion,
  it is called from the reset function which takes a MeshSet object.
  This function prepares the object to be used by setting the initial
  values of some of the data members.  As examples, if needed, it resets
  the cpu timer to zero, the iteration counter to zero, and the
  initial and previous objective function values to the current
  objective function value for this patch.
  The return value for this function is similar to that of terminate().
  The function returns false if the checked criteria have not been
  satisfied, and true if they have been.  reset() only checks the
  GRADIENT_INF_NORM_ABSOLUTE, GRADIENT_L2_NORM_ABSOLUTE, and the
  QUALITY_IMPROVEMENT_ABSOLUTE criteria.  Checking these criteria
  allows the QualityImprover to skip the entire optimization if
  the initial mesh satisfies the appropriate conditions.
 */
void TerminationCriterion::reset_inner(PatchData &pd, OFEvaluator& obj_eval,
                                    MsqError &err)
{
  const unsigned long totalFlag = terminationCriterionFlag | cullingMethodFlag;
  
    // clear flag for BOUNDED_VERTEX_MOVEMENT
  vertexMovementExceedsBound = 0;
  
    // Use -1 to denote that this isn't initialized yet.
    // As all valid values must be >= 0.0, a negative
    // value indicates that it is uninitialized and is
    // always less than any valid value.
  maxSquaredMovement = -1;
  
    // Clear the iteration count.
  iterationCounter = 0;
  
    //reset the inner timer if needed
  if(totalFlag & CPU_TIME){
    mTimer.reset();
  }
   
    //GRADIENT
  currentGradInfNorm = initialGradInfNorm = 0.0;
  currentGradL2NormSquared = initialGradL2NormSquared = 0.0;
  if(totalFlag & GRAD_FLAGS)
  {
    if (!obj_eval.have_objective_function()) {
      MSQ_SETERR(err)("Error termination criteria set which uses objective "
                      "functions, but no objective function is available.",
                      MsqError::INVALID_STATE);   
      return;
    } 
    int num_vertices=pd.num_free_vertices();
    mGrad.resize( num_vertices );

      //get gradient and make sure it is valid
    bool b = obj_eval.evaluate(pd, currentOFValue, mGrad, err); MSQ_ERRRTN(err);
    if (!b) {
      MSQ_SETERR(err)("Initial patch is invalid for gradient computation.", 
                      MsqError::INVALID_STATE);
      return;
    } 

      //get the gradient norms
    if (totalFlag & (GRADIENT_INF_NORM_ABSOLUTE|GRADIENT_INF_NORM_RELATIVE))
    {
      currentGradInfNorm = initialGradInfNorm = Linf(mGrad);
      MSQ_DBGOUT_P0_ONLY(debugLevel) << par_string() << "  o Initial gradient Inf norm: " << " "
                                     << RPM(initialGradInfNorm) << std::endl;
    }  
      
    if (totalFlag & (GRADIENT_L2_NORM_ABSOLUTE|GRADIENT_L2_NORM_RELATIVE))
    {
      currentGradL2NormSquared = initialGradL2NormSquared = length_squared(mGrad);
      MSQ_DBGOUT_P0_ONLY(debugLevel) << par_string() << "  o Initial gradient L2 norm: " << " "
                                     << RPM(std::sqrt(initialGradL2NormSquared)) << std::endl;
    }  

      //the OFvalue comes for free, so save it
    previousOFValue=currentOFValue;
    initialOFValue=currentOFValue;
  }
  //find the initial objective function value if needed and not already
  //computed.  If we needed the gradient, we have the OF value for free.
  // Also, if possible, get initial OF value if writing plot file.  Solvers
  // often supply the OF value for subsequent iterations so by calculating
  // the initial value we can generate OF value plots.
  else if ((totalFlag & OF_FLAGS) || 
           (plotFile.is_open() && pd.num_free_vertices() && obj_eval.have_objective_function()))
  {
      //ensure the obj_ptr is not null
    if(!obj_eval.have_objective_function()){
      MSQ_SETERR(err)("Error termination criteria set which uses objective "
                      "functions, but no objective function is available.",
                      MsqError::INVALID_STATE);
      return;
    }
    
    bool b = obj_eval.evaluate(pd, currentOFValue, err); MSQ_ERRRTN(err);
    if (!b){
      MSQ_SETERR(err)("Initial patch is invalid for evaluation.",MsqError::INVALID_STATE);
      return;
    }
      //std::cout<<"\nReseting initial of value = "<<initialOFValue;
    previousOFValue=currentOFValue;
    initialOFValue=currentOFValue;
  }
  
  if (totalFlag & (GRAD_FLAGS|OF_FLAGS))
    MSQ_DBGOUT_P0_ONLY(debugLevel) << par_string() << "  o Initial OF value: " << " " << RPM(initialOFValue) << std::endl;
  
    // Store current vertex locations now, because we'll
    // need them later to compare the current movement with.
  if (totalFlag & VERTEX_MOVEMENT_RELATIVE)
  {
    if (initialVerticesMemento)
    {
      pd.recreate_vertices_memento( initialVerticesMemento, err );
    }
    else
    {
      initialVerticesMemento = pd.create_vertices_memento( err );
    }
    MSQ_ERRRTN(err);
    maxSquaredInitialMovement = DBL_MAX;
  }
  else {
    maxSquaredInitialMovement = 0;
  }
  
  if (terminationCriterionFlag & UNTANGLED_MESH) {
    globalInvertedCount = count_inverted( pd, err );
    //if (innerOuterType==TYPE_OUTER) MSQ_DBGOUT_P0_ONLY(debugLevel) << par_string() << "  o Num Inverted: " << " " << globalInvertedCount << std::endl;
    patchInvertedCount = 0;
    MSQ_ERRRTN(err);
  }

  if (timeStepFileType) {
      // If didn't already calculate gradient abive, calculate it now.
    if (!(totalFlag & GRAD_FLAGS)) {
      mGrad.resize( pd.num_free_vertices() );
      obj_eval.evaluate(pd, currentOFValue, mGrad, err);
      err.clear();
    }
    write_timestep( pd, mGrad.empty() ? 0 : arrptr(mGrad), err);
  }
    
  if (plotFile.is_open()) {
      // two newlines so GNU plot knows that we are starting a new data set
    plotFile << std::endl << std::endl;
      // write column headings as comment in data file
    plotFile << "#Iter\tCPU\tObjFunc\tGradL2\tGradInf\tMovement\tInverted" << std::endl;
      // write initial values
    plotFile << 0 
     << '\t' << mTimer.since_birth() 
     << '\t' << initialOFValue 
     << '\t' << std::sqrt( currentGradL2NormSquared ) 
     << '\t' << currentGradInfNorm 
     << '\t' << 0.0
     << '\t' << globalInvertedCount
     << std::endl;
  }
}

void TerminationCriterion::reset_patch(PatchData &pd, MsqError &err)
{
  const unsigned long totalFlag = terminationCriterionFlag | cullingMethodFlag;
  if (totalFlag & MOVEMENT_FLAGS)
  {
    if (previousVerticesMemento)
      pd.recreate_vertices_memento(previousVerticesMemento,err); 
    else
      previousVerticesMemento = pd.create_vertices_memento(err);
    MSQ_ERRRTN(err);
  }

  if (totalFlag & UNTANGLED_MESH) {
    patchInvertedCount = count_inverted( pd, err );
    //MSQ_DBGOUT_P0_ONLY(debugLevel) << par_string() << "  o Num Patch Inverted: " << " " << patchInvertedCount << std::endl;
    MSQ_ERRRTN(err);
  }
}

void TerminationCriterion::accumulate_inner( PatchData& pd, 
                                             OFEvaluator& of_eval,
                                             MsqError& err )
{
  double of_value = 0;
  
  if (terminationCriterionFlag & GRAD_FLAGS)
  {
    mGrad.resize( pd.num_free_vertices() );
    bool b = of_eval.evaluate(pd, of_value, mGrad, err);
    MSQ_ERRRTN(err);
    if (!b) {
      MSQ_SETERR(err)("Initial patch is invalid for gradient compuation.",
                      MsqError::INVALID_MESH);
      return;
    }
  }
  else if (terminationCriterionFlag & OF_FLAGS)
  {
    bool b = of_eval.evaluate(pd, of_value, err); MSQ_ERRRTN(err);
    if (!b) {
      MSQ_SETERR(err)("Invalid patch passed to TerminationCriterion.",
                      MsqError::INVALID_MESH);
      return;
    }
  }

  accumulate_inner( pd, of_value, mGrad.empty() ? 0 : arrptr(mGrad), err );  MSQ_CHKERR(err);
}


void TerminationCriterion::accumulate_inner( PatchData& pd, 
                                             double of_value,
                                             Vector3D* grad_array,
                                             MsqError& err )
{
  //if terminating on the norm of the gradient
  //currentGradL2NormSquared = HUGE_VAL;
  if (terminationCriterionFlag & (GRADIENT_L2_NORM_ABSOLUTE | GRADIENT_L2_NORM_RELATIVE)) 
  {
    currentGradL2NormSquared = length_squared(grad_array, pd.num_free_vertices()); // get the L2 norm
    MSQ_DBGOUT_P0_ONLY(debugLevel) << par_string() << "  o Info -- gradient L2 norm: " << " "
                                   << RPM(std::sqrt(currentGradL2NormSquared)) << std::endl;
  }
  //currentGradInfNorm = 10e6;
  if (terminationCriterionFlag & (GRADIENT_INF_NORM_ABSOLUTE | GRADIENT_INF_NORM_RELATIVE)) 
  {
    currentGradInfNorm = Linf(grad_array, pd.num_free_vertices()); // get the Linf norm
    MSQ_DBGOUT_P0_ONLY(debugLevel) << par_string() << "  o Info -- gradient Inf norm: " << " "
                                   << RPM(currentGradInfNorm) << std::endl;
  } 
  
  if (terminationCriterionFlag & VERTEX_MOVEMENT_RELATIVE)
  {
    maxSquaredInitialMovement = pd.get_max_vertex_movement_squared(
                               initialVerticesMemento, err );  MSQ_ERRRTN(err);
    MSQ_DBGOUT_P0_ONLY(debugLevel) << par_string() << "  o Info -- max initial vertex movement: " << " "
                                   << RPM(maxSquaredInitialMovement) << std::endl;
  }
  
  previousOFValue = currentOFValue;
  currentOFValue = of_value;
  if (terminationCriterionFlag & OF_FLAGS) {
    MSQ_DBGOUT_P0_ONLY(debugLevel) << par_string() << "  o Info -- OF Value: " << " " << RPM(of_value) << " iterationCounter= " << iterationCounter << std::endl;
  }
  else if (grad_array) {
    MSQ_DBGOUT_P0_ONLY(debugLevel) << par_string() << "  o OF Value: " << " " << RPM(of_value) << " iterationCounter= " << iterationCounter 
      //<< " terminationCriterionFlag= " << terminationCriterionFlag << " OF_FLAGS = " << OF_FLAGS
                                   << std::endl;
  }
  
  ++iterationCounter;
  if (timeStepFileType)
    write_timestep( pd, grad_array, err);
    
  if (plotFile.is_open()) 
    plotFile << iterationCounter 
     << '\t' << mTimer.since_birth() 
     << '\t' << of_value 
     << '\t' << std::sqrt( currentGradL2NormSquared ) 
     << '\t' << currentGradInfNorm 
     << '\t' << (maxSquaredMovement > 0.0 ? std::sqrt( maxSquaredMovement ) : 0.0)
     << '\t' << globalInvertedCount
     << std::endl;
}


void TerminationCriterion::accumulate_outer(Mesh* mesh, 
                                            MeshDomain* domain, 
                                            OFEvaluator& of_eval,
                                            const Settings* settings,
                                            MsqError &err)
{
  PatchData global_patch;
  if (settings)
    global_patch.attach_settings( settings );
  
    //if we need to fill out the global patch data object.
  if ((terminationCriterionFlag & (GRAD_FLAGS|OF_FLAGS|VERTEX_MOVEMENT_RELATIVE))
      || timeStepFileType)
  {
    global_patch.set_mesh( mesh );
    global_patch.set_domain( domain );
    global_patch.fill_global_patch( err ); MSQ_ERRRTN(err);
  }
  
  accumulate_inner( global_patch, of_eval, err ); MSQ_ERRRTN(err);
}


void TerminationCriterion::accumulate_patch( PatchData& pd, MsqError& err )
{
  if (terminationCriterionFlag & MOVEMENT_FLAGS)
  {
    double patch_max_dist = pd.get_max_vertex_movement_squared( previousVerticesMemento, err );
    if (patch_max_dist > maxSquaredMovement)
      maxSquaredMovement = patch_max_dist;
    pd.recreate_vertices_memento( previousVerticesMemento, err );  MSQ_ERRRTN(err);
  }
    
    //if terminating on bounded vertex movement (a bounding box for the mesh)
  if(terminationCriterionFlag & BOUNDED_VERTEX_MOVEMENT)
  {
    const MsqVertex* vert = pd.get_vertex_array(err);
    int num_vert = pd.num_free_vertices();
    int i=0;
      //for each vertex
    for(i=0;i<num_vert;++i)
    {
        //if any of the coordinates are greater than eps
      if( (vert[i][0]>boundedVertexMovementEps) ||
          (vert[i][1]>boundedVertexMovementEps) ||
          (vert[i][2]>boundedVertexMovementEps) )
      {
        ++vertexMovementExceedsBound;
      }
    }
  }


  if ((terminationCriterionFlag|cullingMethodFlag) & UNTANGLED_MESH) {
    size_t new_count = count_inverted( pd, err );
      // be careful here because size_t is unsigned
    globalInvertedCount += new_count;
    globalInvertedCount -= patchInvertedCount;
    patchInvertedCount = new_count;
    //if (innerOuterType==TYPE_OUTER) 
    //  MSQ_DBGOUT_P0_ONLY(debugLevel) << par_string() << "  o Num Patch Inverted: " << " " << patchInvertedCount << " globalInvertedCount= " << globalInvertedCount << std::endl;
      
    MSQ_ERRRTN(err);
  }
}


/*!  This function evaluates the needed information and then evaluates
  the termination criteria.  If any of the selected criteria are satisfied,
  the function returns true.  Otherwise, the function returns false.
 */
bool TerminationCriterion::terminate( )
{
  bool return_flag = false;
  //std::cout<<"\nInside terminate(pd,of,err):  flag = "<<terminationCriterionFlag << std::endl;
  int type=0;

    //First check for an interrupt signal
  if (MsqInterrupt::interrupt())
  {
     MSQ_DBGOUT_P0_ONLY(debugLevel) << par_string() << "  o TermCrit -- INTERRUPTED" << " " << std::endl;
    return true;
  }
  
    //if terminating on numbering of inner iterations
  if (NUMBER_OF_ITERATES & terminationCriterionFlag
    && iterationCounter >= iterationBound)
  {
    return_flag = true;
    type=1;
    MSQ_DBGOUT_P0_ONLY(debugLevel) << par_string() << "  o TermCrit -- Reached " << " " << iterationBound << " iterations." << std::endl;
  }
  
  if (CPU_TIME & terminationCriterionFlag && mTimer.since_birth()>=timeBound)
  {
    return_flag=true;
    type=2;
    MSQ_DBGOUT_P0_ONLY(debugLevel) << par_string() << "  o TermCrit -- Exceeded CPU time. " << " " << std::endl;
  }
  
  
  if (MOVEMENT_FLAGS & terminationCriterionFlag
      && maxSquaredMovement >= 0.0)
  {
    MSQ_DBGOUT_P0_ONLY(debugLevel) << par_string() << "  o Info -- Maximum vertex movement: " << " "
                                   << RPM(sqrt(maxSquaredMovement)) << std::endl;

    if (VERTEX_MOVEMENT_ABSOLUTE & terminationCriterionFlag 
        && maxSquaredMovement <= vertexMovementAbsoluteEps)
    {
      MSQ_DBGOUT_P0_ONLY(debugLevel) << par_string() << "  o TermCrit -- VERTEX_MOVEMENT_ABSOLUTE: " << " "
                                     << RPM(sqrt(maxSquaredMovement)) << std::endl;
      return_flag = true;
    type=3;
    }
  
    if (VERTEX_MOVEMENT_RELATIVE & terminationCriterionFlag
        && maxSquaredMovement <= vertexMovementRelativeEps*maxSquaredInitialMovement)
    {
      MSQ_DBGOUT_P0_ONLY(debugLevel) << par_string() << "  o TermCrit -- VERTEX_MOVEMENT_RELATIVE: " << " "
                                     << RPM(sqrt(maxSquaredMovement)) << std::endl;
      return_flag = true;
    type=4;
    }

    if (VERTEX_MOVEMENT_ABS_EDGE_LENGTH & terminationCriterionFlag)
    {
      assert( vertexMovementAbsoluteAvgEdge > -1e-12 ); // make sure value actually got calculated
      if (maxSquaredMovement <= vertexMovementAbsoluteAvgEdge)
      {
        MSQ_DBGOUT_P0_ONLY(debugLevel) << par_string() << "  o TermCrit -- VERTEX_MOVEMENT_ABS_EDGE_LENGTH: " << " "
                                       << RPM(sqrt(maxSquaredMovement)) << std::endl;
        return_flag = true;
        type=5;
      }
    }
  }

  if (GRADIENT_L2_NORM_ABSOLUTE & terminationCriterionFlag &&
      currentGradL2NormSquared <= gradL2NormAbsoluteEpsSquared)
  {
    MSQ_DBGOUT_P0_ONLY(debugLevel) << par_string() << "  o TermCrit -- GRADIENT_L2_NORM_ABSOLUTE: " << " " << RPM(currentGradL2NormSquared) << std::endl;
    return_flag = true;
    type=6;
  }
  
  if (GRADIENT_INF_NORM_ABSOLUTE & terminationCriterionFlag &&
      currentGradInfNorm <= gradInfNormAbsoluteEps)
  {
    MSQ_DBGOUT_P0_ONLY(debugLevel) << par_string() << "  o TermCrit -- GRADIENT_INF_NORM_ABSOLUTE: " << " " << RPM(currentGradInfNorm) << std::endl;
    return_flag = true;
    type=7;
  }
  
  if (GRADIENT_L2_NORM_RELATIVE & terminationCriterionFlag &&
      currentGradL2NormSquared <= (gradL2NormRelativeEpsSquared * initialGradL2NormSquared))
  {
    MSQ_DBGOUT_P0_ONLY(debugLevel) << par_string() << "  o TermCrit -- GRADIENT_L2_NORM_RELATIVE: " << " " << RPM(currentGradL2NormSquared) << std::endl;
    return_flag = true;
    type=8;
  }
  
  if (GRADIENT_INF_NORM_RELATIVE & terminationCriterionFlag &&
      currentGradInfNorm <= (gradInfNormRelativeEps * initialGradInfNorm))
  {
    MSQ_DBGOUT_P0_ONLY(debugLevel) << par_string() << "  o TermCrit -- GRADIENT_INF_NORM_RELATIVE: " << " " << RPM(currentGradInfNorm) << std::endl;
    return_flag = true;
    type=9;
  }
    //Quality Improvement and Successive Improvements are below.
    // The relative forms are only valid after the first iteration.
  if ((QUALITY_IMPROVEMENT_ABSOLUTE & terminationCriterionFlag) &&
      currentOFValue <= qualityImprovementAbsoluteEps)
  {
    MSQ_DBGOUT_P0_ONLY(debugLevel) << par_string() << "  o TermCrit -- QUALITY_IMPROVEMENT_ABSOLUTE: " << " " << RPM(currentOFValue) << std::endl;
    return_flag = true;
    type=10;
  }
  
    //only valid if an iteration has occurred, see above.
  if(iterationCounter > 0){
    if (SUCCESSIVE_IMPROVEMENTS_ABSOLUTE & terminationCriterionFlag &&
        (previousOFValue - currentOFValue) <= successiveImprovementsAbsoluteEps)
    {  
      MSQ_DBGOUT_P0_ONLY(debugLevel) << par_string() << "  o TermCrit -- SUCCESSIVE_IMPROVEMENTS_ABSOLUTE: previousOFValue= " << " " << previousOFValue << " currentOFValue= " << RPM(currentOFValue) 
                             << " successiveImprovementsAbsoluteEps= " << successiveImprovementsAbsoluteEps
                             << std::endl;
      return_flag = true;
    type=11;
    }
    if (QUALITY_IMPROVEMENT_RELATIVE & terminationCriterionFlag &&
        (currentOFValue - lowerOFBound) <= 
        qualityImprovementRelativeEps * (initialOFValue - lowerOFBound))
    {
      MSQ_DBGOUT_P0_ONLY(debugLevel) << par_string() << "  o TermCrit -- QUALITY_IMPROVEMENT_RELATIVE: " << " " << std::endl;
      return_flag = true;
    type=12;
    }
    if (SUCCESSIVE_IMPROVEMENTS_RELATIVE & terminationCriterionFlag &&
        (previousOFValue - currentOFValue) <= 
        successiveImprovementsRelativeEps * (initialOFValue - currentOFValue))
    {
      MSQ_DBGOUT_P0_ONLY(debugLevel) << par_string() << "  o TermCrit -- SUCCESSIVE_IMPROVEMENTS_RELATIVE: previousOFValue= " << " " << previousOFValue << " currentOFValue= " << RPM(currentOFValue) 
                             << " successiveImprovementsRelativeEps= " << successiveImprovementsRelativeEps
                             << std::endl;
      return_flag = true;
    type=13;
    }
  }
  
  if (BOUNDED_VERTEX_MOVEMENT & terminationCriterionFlag && vertexMovementExceedsBound)
  {
    return_flag = true;
    type=14;
    MSQ_DBGOUT_P0_ONLY(debugLevel) << par_string() << "  o TermCrit -- " << " " << vertexMovementExceedsBound
                           << " vertices out of bounds." << std::endl;
  }
  
  if (UNTANGLED_MESH & terminationCriterionFlag)
  {
    if (innerOuterType==TYPE_OUTER) {
      //MSQ_DBGOUT_P0_ONLY(debugLevel) 
      MSQ_DBGOUT(debugLevel) 
        << par_string() << "  o Num Inverted: " << " " << globalInvertedCount << std::endl;
    }
    if (!globalInvertedCount)
      {
        MSQ_DBGOUT_P0_ONLY(debugLevel) << par_string() << "  o TermCrit -- UNTANGLED_MESH: "<< " " << std::endl;
        return_flag = true;
        type=15;
      }
  }
  
    // clear this value at the end of each iteration
  vertexMovementExceedsBound = 0;
  maxSquaredMovement = -1.0;

  if (timeStepFileType == GNUPLOT && return_flag) {
    MsqError err;
    MeshWriter::write_gnuplot_overlay( iterationCounter, timeStepFileName.c_str(), err );
  }

  if (0 && return_flag && MSQ_DBG(2))
    std::cout << "P[" << get_parallel_rank() << "] tmp TerminationCriterion::terminate: " << moniker << " return_flag= " << return_flag << " type= " << type
              << " terminationCriterionFlag= " << terminationCriterionFlag << " debugLevel= " << debugLevel << std::endl;
  
    //if none of the criteria were satisfied
  return return_flag;
}

bool TerminationCriterion::criterion_is_set()
{
  if (!terminationCriterionFlag)
    return false;
  else
    return true;
}


/*!This function checks the culling method criterion supplied to the object
  by the user.  If the user does not supply a culling method criterion,
  the default criterion is NONE, and in that case, no culling is performed.
  If the culling method criterion is satisfied, the interior vertices
  of the given patch are flagged as soft_fixed.  Otherwise, the soft_fixed
  flag is removed from each of the vertices in the patch (interior and
  boundary vertices).  Also, if the criterion was satisfied, then the
  function returns true.  Otherwise, the function returns false.
 */
bool TerminationCriterion::cull_vertices(PatchData &pd,
                                      OFEvaluator& of_eval,
                                      MsqError &err)
{
    //PRINT_INFO("CULLING_METHOD FLAG = %i",cullingMethodFlag);
  
    //cull_bool will be changed to true if the criterion is satisfied
  bool b, cull_bool=false;
  double prev_m, init_m;
  switch(cullingMethodFlag){
      //if no culling is requested, always return false
    case NONE:
       return cull_bool;
         //if culling on quality improvement absolute
    case QUALITY_IMPROVEMENT_ABSOLUTE:
         //get objective function value
       b = of_eval.evaluate(pd, currentOFValue, err);
       if (MSQ_CHKERR(err)) return false;
       if (!b) {
         MSQ_SETERR(err)(MsqError::INVALID_MESH);
         return false;
       }
         //if the improvement was enough, cull
       if(currentOFValue <= cullingEps)
       {
         cull_bool=true;  
       }
         //PRINT_INFO("\ncurrentOFValue = %f, bool = %i\n",currentOFValue,cull_bool);
       
       break;
         //if culing on quality improvement relative
    case QUALITY_IMPROVEMENT_RELATIVE:
         //get objective function value
       b = of_eval.evaluate(pd, currentOFValue, err);
       if (MSQ_CHKERR(err)) return false;
       if(!b){
         MSQ_SETERR(err)(MsqError::INVALID_MESH);
         return false;
       }
         //if the improvement was enough, cull
       if((currentOFValue-lowerOFBound)<=
          (cullingEps*(initialOFValue-lowerOFBound)))
       {
         cull_bool=true;  
       }
       break;
         //if culling on vertex movement absolute
    case VERTEX_MOVEMENT_ABSOLUTE:
    case VERTEX_MOVEMENT_ABS_EDGE_LENGTH:
         //if movement was enough, cull
       prev_m = pd.get_max_vertex_movement_squared(previousVerticesMemento,err);
       MSQ_ERRZERO(err);
       if(prev_m <= cullingEps*cullingEps){
         cull_bool=true;  
       }
       
       break;
         //if culling on vertex movement relative
    case VERTEX_MOVEMENT_RELATIVE:
         //if movement was small enough, cull
       prev_m = pd.get_max_vertex_movement_squared(previousVerticesMemento,err);
       MSQ_ERRZERO(err);
       init_m = pd.get_max_vertex_movement_squared(initialVerticesMemento,err);
       MSQ_ERRZERO(err);
       if(prev_m <= (cullingEps*cullingEps * init_m)){
         cull_bool=true;  
       }
       break;
    case UNTANGLED_MESH:
      if (!patchInvertedCount)
        cull_bool = true;
      break;
    default:
       MSQ_SETERR(err)("Requested culling method not yet implemented.",
                       MsqError::NOT_IMPLEMENTED);
       return false;
  };
    //Now actually have patch data cull vertices
  if(cull_bool)
  {
    pd.set_free_vertices_soft_fixed(err); MSQ_ERRZERO(err);
  }
  else
  {
    pd.set_all_vertices_soft_free(err); MSQ_ERRZERO(err);
  }
  return cull_bool;
}

/*!This function is activated when cullingGlobalPatch is true.  It supplies
  cull_vertices with a single vertex-based patch at a time.  If the patch
  satisfies the culling criterion, it's free vertices are then soft-fixed.
 */
bool TerminationCriterion::cull_vertices_global(PatchData &global_patch,
                                                Mesh *mesh, MeshDomain *domain, const Settings *settings,
                                                OFEvaluator& of_eval,
                                                MsqError &err)
{
  if (!cullingGlobalPatch) return false;

    //PRINT_INFO("CULLING_METHOD FLAG = %i",cullingMethodFlag);
  
    //cull_bool will be changed to true if the criterion is satisfied
  bool cull_bool=false;

  std::vector<Mesh::VertexHandle> mesh_vertices;
  //std::vector<Mesh::VertexHandle> patch_vertices;
  //std::vector<Mesh::ElementHandle> patch_elements;
  //std::vector<Mesh::VertexHandle> fixed_vertices;
  //std::vector<Mesh::VertexHandle> free_vertices;

  // FIXME, verify global_patch is a global patch... how, is this right?
  mesh->get_all_vertices(mesh_vertices, err);
  size_t mesh_num_nodes = mesh_vertices.size();
  size_t global_patch_num_nodes = global_patch.num_nodes() ;
  if (0)  std::cout << "tmp srk mesh_num_nodes= " << mesh_num_nodes << " global_patch_num_nodes= " 
            << global_patch_num_nodes << std::endl;
  if (mesh_num_nodes != global_patch_num_nodes)
    {
      std::cout << "tmp srk cull_vertices_global found non global patch" << std::endl;
      exit(123);
      return false;
    }
  PatchData patch;
  patch.set_mesh( (Mesh*) mesh );
  patch.set_domain( domain );
  patch.attach_settings( settings );

  const MsqVertex* global_patch_vertex_array = global_patch.get_vertex_array( err );
  Mesh::VertexHandle* global_patch_vertex_handles = global_patch.get_vertex_handles_array();

  int num_culled = 0;
  for (unsigned iv=0; iv < global_patch_num_nodes; iv++)
    {
      // form a patch for this vertex; if it is culled, set it to be soft fixed
      Mesh::VertexHandle vert = global_patch_vertex_handles[iv];
      std::vector<Mesh::ElementHandle> elements;
      std::vector<size_t> offsets;
      mesh->vertices_get_attached_elements(&vert, 1, elements, offsets, err);
      
      std::set<Mesh::VertexHandle> patch_free_vertices_set;

      for (unsigned ie=0; ie < elements.size(); ie++)
        {
          std::vector<Mesh::VertexHandle> vert_handles;
          std::vector<size_t> v_offsets;
          mesh->elements_get_attached_vertices(&elements[ie], 1, vert_handles, v_offsets, err);
          for (unsigned jv=0; jv < vert_handles.size(); jv++)
            {
              unsigned char bt;
              mesh->vertex_get_byte(vert_handles[jv], &bt, err);
              MsqVertex v;
              v.set_flags(bt);
              if (v.is_free_vertex())
                patch_free_vertices_set.insert(vert_handles[jv]);
            }
        }

      std::vector<Mesh::VertexHandle> patch_free_vertices_vector(patch_free_vertices_set.begin(), patch_free_vertices_set.end());
      //std::vector<unsigned char> byte_vector(patch_vertices_vector.size());
      //mesh->vertices_get_byte(&vert_handles[0], &byte_vector[0], vert_handles.size(), err);

      patch.set_mesh_entities( elements, patch_free_vertices_vector, err );
      if (cull_vertices(patch, of_eval, err))
        {
          //std::cout << "tmp srk cull_vertices_global found culled patch" << std::endl;
          Mesh::VertexHandle* patch_vertex_handles = patch.get_vertex_handles_array();
          const MsqVertex* patch_vertex_array = patch.get_vertex_array( err );
          for (unsigned jv=0; jv < patch.num_nodes(); jv++)
            {
              if (patch_vertex_handles[jv] == global_patch_vertex_handles[iv])
                {
                  if (patch_vertex_array[jv].is_flag_set(MsqVertex::MSQ_CULLED))
                    {
                      global_patch.set_vertex_culled(iv);
                      ++num_culled;
                      cull_bool = true;
                      //std::cout << "tmp srk cull_vertices_global found culled vertex" << std::endl;
                    }
                }
            }
        }
    }
    if (0)  std::cout << "tmp srk cull_vertices_global found " << num_culled << " culled vertices out of " << global_patch_num_nodes << std::endl;

  return cull_bool;
}

size_t TerminationCriterion::count_inverted( PatchData& pd, MsqError& err )
{
  size_t num_elem = pd.num_elements();
  size_t count=0;
  int inverted, samples;
  for (size_t i = 0; i < num_elem; i++) {
    pd.element_by_index(i).check_element_orientation(pd, inverted, samples, err);
    if (inverted)
      ++count;
  }
  return count;
}

/*!
  Currently this only deletes the memento of the vertex positions and the
  mGrad vector if neccessary.
  When culling, we remove the soft fixed flags from all of the vertices.
 */
void TerminationCriterion::cleanup(Mesh* mesh, MeshDomain*, MsqError &err)
{
  delete previousVerticesMemento;
  delete initialVerticesMemento;
  previousVerticesMemento = 0;
  initialVerticesMemento = 0;
}

void TerminationCriterion::write_timestep( PatchData& pd, 
                                           const Vector3D* gradient,
                                           MsqError& err )
{
  std::ostringstream str;
  if (timeStepFileType == VTK) {
    str << timeStepFileName << '_' << iterationCounter << ".vtk";
    MeshWriter::write_vtk( pd, str.str().c_str(), err, gradient );
  }
  else if (timeStepFileType == GNUPLOT) {
    str << timeStepFileName << '.' << iterationCounter;
    MeshWriter::write_gnuplot( pd.get_mesh(), str.str().c_str(), err );
  }
}

void TerminationCriterion::write_iterations( const char* filename, MsqError& err )
{
  if (filename) {
    plotFile.open( filename, std::ios::trunc );
    if (!plotFile)
      MSQ_SETERR(err)( MsqError::FILE_ACCESS, "Failed to open plot data file: '%s'", filename );
  }
  else {
    plotFile.close();
  }
}
    
    
void TerminationCriterion::initialize_queue( MeshDomainAssoc* mesh_and_domain,
                                             const Settings* ,
                                             MsqError& err )
{
  if (VERTEX_MOVEMENT_ABS_EDGE_LENGTH & (terminationCriterionFlag|cullingMethodFlag)) 
  {
    Mesh* mesh = mesh_and_domain->get_mesh();
    MeshUtil tool(mesh);
    SimpleStats stats;
    tool.edge_length_distribution( stats, err );
    MSQ_ERRRTN(err);
    double limit = vertexMovementAvgBeta * (stats.average() - stats.standard_deviation());
      // we actually calculate the square of the length
    vertexMovementAbsoluteAvgEdge = limit * limit;
    if (VERTEX_MOVEMENT_ABS_EDGE_LENGTH & cullingMethodFlag)
      cullingEps = limit;
  }
}

} //namespace Mesquite
