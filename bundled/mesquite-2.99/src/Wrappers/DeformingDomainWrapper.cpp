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


/** \file DeformingDomainWrapper.cpp
 *  \brief Implement DeformingDomainWrapper class
 *  \author Jason Kraftcheck 
 */

#include "DeformingDomainWrapper.hpp"
#include "TagVertexMesh.hpp"
#include "ReferenceMesh.hpp"
#include "RefMeshTargetCalculator.hpp"
#include "InstructionQueue.hpp"
#include "QualityImprover.hpp"
#include "SteepestDescent.hpp"
#include "TerminationCriterion.hpp"
#include "PMeanPTemplate.hpp"
#include "ElementPMeanP.hpp"
#include "TQualityMetric.hpp"
#include "TMixed.hpp"
#include "TShapeNB1.hpp"
#include "TShapeSize2DNB1.hpp"
#include "TShapeSize3DNB1.hpp"
#include "TShapeSizeOrientNB1.hpp"
#include "MeshInterface.hpp"
#include "MsqError.hpp"
#include "MsqVertex.hpp"
#include "QualityAssessor.hpp"
#include "CurveDomain.hpp"

#include <numeric> // for std::accumulate


namespace MESQUITE_NS {

const DeformingDomainWrapper::MeshCharacteristic DEFAULT_METRIC_TYPE = 
  DeformingDomainWrapper::SHAPE;

const bool DEFAULT_CULLING = true;
const double DEFAULT_CPU_TIME = 0.0;
const double DEFAULT_MOVEMENT_FACTOR = 0.01;
const int DEFAULT_INNER_ITERATIONS = 2;
const char DEFAULT_CURVE_TAG[] = "MesquiteCurveFraction";
const DeformingCurveSmoother::Scheme DEFAULT_CURVE_TYPE = 
  DeformingCurveSmoother::PROPORTIONAL;

DeformingDomainWrapper::DeformingDomainWrapper()
  : metricType( DEFAULT_METRIC_TYPE ),
    initVertexCoords(""),
    doCulling(DEFAULT_CULLING),
    cpuTime(DEFAULT_CPU_TIME),
    movementFactor(DEFAULT_MOVEMENT_FACTOR)
  {}

DeformingDomainWrapper::~DeformingDomainWrapper() {}

void DeformingDomainWrapper::store_initial_mesh( Mesh* mesh, MsqError& err )
{
  TagVertexMesh tool( err, mesh, false, initVertexCoords ); MSQ_ERRRTN(err);
  InstructionQueue q;
  q.add_tag_vertex_mesh( &tool, err ); MSQ_ERRRTN(err);
  q.run_instructions( mesh, err ); MSQ_ERRRTN(err);
  if (initVertexCoords.empty())
    initVertexCoords = tool.get_tag_name();
}

void DeformingDomainWrapper::set_vertex_movement_limit_factor( double f )
{
  assert(f > 0.0);
  movementFactor = f;
}

void DeformingDomainWrapper::run_wrapper( MeshDomainAssoc* mesh_and_domain,
                                          ParallelMesh* pmesh,
                                          Settings* settings,
                                          QualityAssessor* qa,
                                          MsqError& err )
{
  if (movementFactor <= 0) {
    MSQ_SETERR(err)(MsqError::INVALID_STATE,
      "Optimization will not terminate with non-positive movement factor %f",
      movementFactor);
    return;
  }

  Mesh* mesh = mesh_and_domain->get_mesh();
  MeshDomain* geom = mesh_and_domain->get_domain();

    // Move initial mesh to domain in case caller did not
  move_to_domain( mesh, geom, err ); MSQ_ERRRTN(err);

  const double P = 1.0;
  TagVertexMesh init_mesh( err, mesh ); MSQ_ERRRTN(err);
  ReferenceMesh ref_mesh( &init_mesh );
  RefMeshTargetCalculator W( &ref_mesh );
  
  TShapeNB1 mu_s;
  TShapeSize2DNB1 mu_2d_ss;
  TShapeSize3DNB1 mu_3d_ss;
  TMixed mu_ss( &mu_2d_ss, &mu_3d_ss );
  TShapeSizeOrientNB1 mu_sso;
  TMetric* mu = 0;
  switch (metricType) {
    case SHAPE:             mu = &mu_s  ; break;
    case SHAPE_SIZE:        mu = &mu_ss ; break;
    case SHAPE_SIZE_ORIENT: mu = &mu_sso; break;
  }

  TQualityMetric sample_metric( &W, mu );
  ElementPMeanP elem_metric( P, &sample_metric );
  PMeanPTemplate obj_func( P, &elem_metric );
  
  TerminationCriterion inner, outer;
  SteepestDescent improver( &obj_func );
  improver.use_element_on_vertex_patch(); // Nash optimization
  inner.add_iteration_limit( DEFAULT_INNER_ITERATIONS );
  if (doCulling)
    inner.cull_on_absolute_vertex_movement_edge_length( movementFactor );
  else
    outer.add_absolute_vertex_movement_edge_length( movementFactor );
  if (cpuTime > 0.0)
    outer.add_cpu_time( cpuTime );
  improver.set_inner_termination_criterion( &inner );
  improver.set_outer_termination_criterion( &outer );
  
  qa->add_quality_assessment(&elem_metric);
  InstructionQueue q;
  q.add_quality_assessor( qa, err ); MSQ_ERRRTN(err);
  q.set_master_quality_improver( &improver, err ); MSQ_ERRRTN(err);
  q.add_quality_assessor( qa, err ); MSQ_ERRRTN(err);
  q.run_common( mesh_and_domain, pmesh, settings, err ); MSQ_ERRRTN(err);
}

void DeformingDomainWrapper::move_to_domain( Mesh* mesh, 
                                             MeshDomain* geom, 
                                             MsqError& err )
{
  std::vector<Mesh::VertexHandle> verts;
  mesh->get_all_vertices( verts, err ); MSQ_ERRRTN(err);
  
  MsqVertex coords;
  std::vector<Mesh::VertexHandle>::const_iterator i;
  for (i = verts.begin(); i != verts.end(); ++i) {
    mesh->vertices_get_coordinates( &*i, &coords, 1, err ); MSQ_ERRRTN(err);
    geom->snap_to( *i, coords );
    mesh->vertex_set_coordinates( *i, coords, err ); MSQ_ERRRTN(err);
  }
}

DeformingCurveSmoother::DeformingCurveSmoother()
  : metricType(DEFAULT_CURVE_TYPE), initFractTag(DEFAULT_CURVE_TAG)
  {}

DeformingCurveSmoother::~DeformingCurveSmoother()
  {}

void DeformingCurveSmoother::store_initial_mesh( Mesh* mesh,
                                                 const Mesh::VertexHandle* verts,
                                                 int nverts,
                                                 CurveDomain* geom,
                                                 MsqError& err )
{
  if (nverts < 2) {
    MSQ_SETERR(err)("Invalid curve mesh.  Need at least two end vertices.",
                    MsqError::INVALID_MESH);
    return;
  }

  // get edge lengths
  std::vector<double> vals(nverts - 1);
  MsqVertex coords[2];
  int prev = 0;
  mesh->vertices_get_coordinates( verts, coords+prev, 1, err ); MSQ_ERRRTN(err);
  for (int i = 1; i < nverts; ++i) {
    int next = 1 - prev;
    mesh->vertices_get_coordinates( verts+i, coords+next, 1, err ); MSQ_ERRRTN(err);
    vals[i-1] = (coords[0] - coords[1]).length();
    prev = next;
  }
  
  // convert to length fraction before each iterior vertex
  // (sum of lengths of adjacent edges over total curve length)
  const double total = std::accumulate( vals.begin(), vals.end(), 0.0 );
  for (int i = 1; i < nverts-1; ++i)
    vals[i-1] = vals[i-1]/total;
  vals.resize(nverts-2);
  
  // create tag
  TagHandle tag = mesh->tag_create( initFractTag, Mesh::DOUBLE, 1, 0, err );
  if (err.error_code() == MsqError::TAG_ALREADY_EXISTS) {
    err.clear();
    tag = get_tag( mesh, err ); 
  }
  MSQ_ERRRTN(err);
  
  // store tag data on interior vertices
  mesh->tag_set_vertex_data( tag, nverts-2, verts+1, &vals[0], err ); MSQ_ERRRTN(err);
}

void DeformingCurveSmoother::smooth_curve( Mesh* mesh,
                                           const Mesh::VertexHandle* verts,
                                           int nverts,
                                           CurveDomain* geom,
                                           Scheme type,
                                           MsqError& err )
{
  if (nverts < 2) {
    MSQ_SETERR(err)("Invalid curve mesh.  Need at least two end vertices.",
                    MsqError::INVALID_MESH);
    return;
  }

  // Verify that end vertices are on curve.  We cannot snap end vertices
  // to ends of curve for application because we don't know where the
  // ends of the curve are.
  MsqVertex coords[2], coords2;
  Mesh::VertexHandle ends[2] = { verts[0], verts[nverts-1] };
  for (int i = 0; i < 2; ++i) {
    mesh->vertices_get_coordinates( ends+i, coords+i, 1, err ); MSQ_ERRRTN(err);
    geom->position_from_length( coords[i].to_array(), 0, coords2.to_array(), err );
    MSQ_ERRRTN(err);
    if ((coords[i] - coords2).length_squared() > DBL_EPSILON) {
      MSQ_SETERR(err)("Curve end vertices do not line on curve.  Move ends to curve end points first",
                      MsqError::INVALID_MESH);
      return;
    }
  }
  const double total = geom->arc_length( coords[0].to_array(), coords[1].to_array(), err ); MSQ_ERRRTN(err);
  
  std::vector<double> vals(nverts-1);
  if (metricType == EQUAL) {
    std::fill( vals.begin(), vals.end(), 1.0/(nverts-1) );
    //fracsum = 1.0;
  }
  else { // metricType == PROPORTIONAL
    TagHandle tag = get_tag( mesh, err ); MSQ_ERRRTN(err);
    mesh->tag_get_vertex_data( tag, nverts-2, verts+1, &vals[0], err ); MSQ_ERRRTN(err);
    double sum = std::accumulate( vals.begin(), vals.end()-1, 0.0 ); 
    if (1.0 - sum > 1e-8)
      vals.back() = 1.0 - sum;
    else {
      vals.back() = *std::min_element( vals.begin(), vals.end()-1 );
      sum += vals.back();
      for (size_t i = 0; i < vals.size(); ++i)
        vals[i] /= sum;
    }
  }
  
  double frac_sum = 0.0;
  for (int i = 1; i < nverts-1; ++i) {
    frac_sum += vals[i-1];
    geom->position_from_length( coords[0].to_array(), total * frac_sum, coords[1].to_array(), err ); MSQ_ERRRTN(err);
    mesh->vertex_set_coordinates( verts[i], coords[1], err ); MSQ_ERRRTN(err);
  }
}

TagHandle DeformingCurveSmoother::get_tag( Mesh* mesh, MsqError& err )
{
  TagHandle h = mesh->tag_get( initFractTag, err ); MSQ_ERRZERO(err);
  std::string name;
  Mesh::TagType type;
  unsigned len;
  mesh->tag_properties( h, name, type, len, err ); MSQ_ERRZERO(err);
  if (type != Mesh::DOUBLE || len != 1) {
    MSQ_SETERR(err)(MsqError::INVALID_MESH,
                    "Tag \"%s\" exists but is not 1 double value per vertex.",
                    initFractTag.c_str());
  }
  return h;
}

} // namespace MESQUITE_NS
