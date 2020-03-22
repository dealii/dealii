/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2010 Sandia National Laboratories.  Developed at the
    University of Wisconsin--Madison under SNL contract number
    624796.  The U.S. Government and the University of Wisconsin
    retain certain rights to this software.

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

    (2010) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file LaplaceWrapper.cpp
 *  \brief Implement LaplaceWrapper class
 *  \author Jason Kraftcheck 
 */

#include "LaplaceWrapper.hpp"
#include "IdealWeightInverseMeanRatio.hpp" 
#include "LaplacianSmoother.hpp"
#include "QualityAssessor.hpp"
#include "InstructionQueue.hpp"
#include "TerminationCriterion.hpp"
#include "MsqError.hpp"

namespace MESQUITE_NS {

const double DEFAULT_MOVEMENT_FACTOR = 0.001;
const bool CULLING_DEFAULT = true;
const int DEFAULT_ITERATION_LIMIT = 100;

LaplaceWrapper::LaplaceWrapper() 
  : maxTime(-1.0),
    movementFactor(DEFAULT_MOVEMENT_FACTOR),
    iterationLimit(DEFAULT_ITERATION_LIMIT),
    doCulling(CULLING_DEFAULT)
{}

LaplaceWrapper::~LaplaceWrapper()
{}

void LaplaceWrapper::run_wrapper( MeshDomainAssoc* mesh_and_domain,
                                  ParallelMesh* pmesh,
                                  Settings* settings,
                                  QualityAssessor* qa,
                                  MsqError& err )
{
  if (maxTime <= 0.0 && movementFactor <= 0.0 && iterationLimit <= 0) {
    MSQ_SETERR(err)("No termination criterion set.  "
                    "LaplaceWrapper will run forever.", 
                    MsqError::INVALID_STATE);
    return;
  }
  
  IdealWeightInverseMeanRatio qa_metric;
  qa->add_quality_assessment( &qa_metric );
  
  LaplacianSmoother smoother;
  TerminationCriterion outer("<type:laplace_outer>"), inner("<type:laplace_inner>");
  if (maxTime > 0.0)
    outer.add_cpu_time( maxTime );
  if (iterationLimit > 0)
    outer.add_iteration_limit( iterationLimit );
  if (doCulling && movementFactor > 0.0) {
    inner.cull_on_absolute_vertex_movement_edge_length( movementFactor );
    smoother.set_inner_termination_criterion( &inner );
  }
  else if (movementFactor > 0.0) {
    outer.add_absolute_vertex_movement_edge_length( movementFactor );
  }
  smoother.set_outer_termination_criterion( &outer );
  
  InstructionQueue q;
  q.add_quality_assessor( qa, err ); MSQ_ERRRTN(err);
  q.set_master_quality_improver( &smoother, err ); MSQ_ERRRTN(err);
  q.add_quality_assessor( qa, err ); MSQ_ERRRTN(err);
  q.run_common( mesh_and_domain, pmesh, settings, err ); MSQ_ERRRTN(err);
}

} // namespace MESQUITE_NS
