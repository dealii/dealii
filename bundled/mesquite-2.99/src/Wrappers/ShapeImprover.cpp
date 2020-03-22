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


/** \file ShapeImprover.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 *  \author Boyd Tidwell
 *  \date   10-23-12
 */


#include "ShapeImprover.hpp"
#include "MsqTimer.hpp"
#include "MsqDebug.hpp"
#include "PMeanPTemplate.hpp"
#include "ElementPMeanP.hpp"
#include "ConjugateGradient.hpp"
#include "TerminationCriterion.hpp"
#include "InstructionQueue.hpp"
#include "QualityAssessor.hpp"
#include "ElementPatches.hpp"

#include "TShapeB1.hpp"
#include "TShapeNB1.hpp"
#include "TQualityMetric.hpp"
#include "IdealShapeTarget.hpp"

namespace MESQUITE_NS {

const double DEFAULT_BETA = 0.0;
const int DEFUALT_PARALLEL_ITERATIONS = 10;

ShapeImprover::ShapeImprover() 
 : maxTime(300.0), 
   mBeta(DEFAULT_BETA),
   parallelIterations(DEFUALT_PARALLEL_ITERATIONS)
{}

void ShapeImprover::set_vertex_movement_limit_factor( double beta )
{
  assert(beta > 0.0);
  assert(beta < 1.0);
  mBeta = beta;
}

void ShapeImprover::set_cpu_time_limit( double limit )
{
  assert(limit >= 0.0);
  maxTime = limit;
}

void ShapeImprover::set_parallel_iterations( int count )
{
  if (count < 1) {
    assert(false);
    count = 1;
  }
  parallelIterations = count;
}

void ShapeImprover::run_wrapper( MeshDomainAssoc* mesh_and_domain,
                                 ParallelMesh* pmesh,
                                 Settings* settings,
                                 QualityAssessor* qa,
                                 MsqError& err )
{
    // Quality Metrics
  IdealShapeTarget target;
  Mesh* mesh = mesh_and_domain->get_mesh();
  MeshDomain* domain = mesh_and_domain->get_domain();

  // only calc min edge length if user hasn't set the vertex_movement_limit_factor
  if (!mBeta)
  {  
    // Calculate minimum edge length in mesh

      // create a temp global patch to get min edge len from
    PatchData patch;
    patch.set_mesh(mesh);
    patch.set_domain(domain);
    if (settings)
      patch.attach_settings(settings);
    std::vector<Mesh::ElementHandle> patch_elems;
    std::vector<Mesh::VertexHandle> patch_verts;
    mesh->get_all_elements(patch_elems, err);
    mesh->get_all_vertices(patch_verts, err);
    patch.set_mesh_entities(patch_elems, patch_verts, err);

      // get min edge length from temp global patch
    double min_edge_len = 0.0, max_edge_len = 0.0;
    patch.get_minmax_edge_length(min_edge_len, max_edge_len);
    mBeta = min_edge_len / 10.0;
  }

  TerminationCriterion inner_b;

  // check for inverted elements 

    // create QualityAssessor instance without a Quality Metric
  QualityAssessor check_inverted(false, false);
  InstructionQueue q_invert_check;
    // a QuallityAssessor without a metric will just check for inverted elements and samples
  q_invert_check.add_quality_assessor(&check_inverted, err);
  q_invert_check.run_common( mesh_and_domain, pmesh, settings, err ); MSQ_ERRRTN(err);
  int inverted_elems = 0, inverted_samples = 0;
  check_inverted.get_inverted_element_count(inverted_elems, inverted_samples, err);
  if (inverted_elems || inverted_samples)
  {
    // use non barrier shape improvement on tangled mesh
    TShapeNB1 mu_no;
    TQualityMetric metric_no_0( &target, &mu_no );
    ElementPMeanP metric_no( 1.0, &metric_no_0 );
    
    // QualityAssessor
    qa->add_quality_assessment( &metric_no );
   
    PMeanPTemplate obj_func_no( 1.0, &metric_no );
    ConjugateGradient improver_no( &obj_func_no );
    improver_no.use_global_patch();
    TerminationCriterion inner_no;
    if (maxTime > 0.0)
      inner_no.add_cpu_time( maxTime ); 
    inner_no.add_absolute_vertex_movement_edge_length( mBeta );
    improver_no.set_inner_termination_criterion( &inner_no );
    InstructionQueue q_no;
    q_no.add_quality_assessor( qa, err ); MSQ_ERRRTN(err);
    q_no.set_master_quality_improver( &improver_no, err ); MSQ_ERRRTN(err);
    q_no.add_quality_assessor( qa, err ); MSQ_ERRRTN(err);
    q_no.run_common( mesh_and_domain, pmesh, settings, err ); MSQ_ERRRTN(err);
  }
  else
  {
    // use barrer metric on non-tangled mesh

    Timer totalTimer;
    TShapeB1 mu_b;
    TQualityMetric metric_b_0( &target, &mu_b );
    ElementPMeanP metric_b( 1.0, &metric_b_0 );
    qa->add_quality_assessment( &metric_b );
    PMeanPTemplate obj_func_b( 1.0, &metric_b );
    ConjugateGradient improver_b( &obj_func_b );
    improver_b.use_global_patch();

    if (maxTime > 0.0)
      inner_b.add_cpu_time( maxTime ); 
    inner_b.add_absolute_vertex_movement_edge_length( mBeta );

    improver_b.set_inner_termination_criterion( &inner_b );
    InstructionQueue q_b;
    q_b.add_quality_assessor( qa, err ); MSQ_ERRRTN(err);
    q_b.set_master_quality_improver( &improver_b, err ); MSQ_ERRRTN(err);
    q_b.add_quality_assessor( qa, err ); MSQ_ERRRTN(err);
    q_b.run_common( mesh_and_domain, pmesh, settings, err ); MSQ_ERRRTN(err);
  }
}

} // namespace MESQUITE_NS
