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
/*! \file EdgeLengthQualityMetric.cpp
  \author Michael Brewer
  \date 2002-05-14
  Evaluates the lengths of the edges attached to the given vertex.
  By default, the averaging method is set to SUM.
*/


#include "EdgeLengthQualityMetric.hpp"
#include "Vector3D.hpp"
#include "QualityMetric.hpp"
#include "MsqVertex.hpp"
#include "PatchData.hpp"
#include "MsqError.hpp"

#include <math.h>

using namespace Mesquite;

std::string EdgeLengthQualityMetric::get_name() const
  { return "Edge Length"; }

int EdgeLengthQualityMetric::get_negate_flag() const
  { return 1; }

bool EdgeLengthQualityMetric::evaluate_common(PatchData &pd, 
                                              size_t this_vert,
                                              double &fval, 
                                              std::vector<size_t>& adj_verts,
                                              MsqError &err)
{
  fval=0.0;
  pd.get_adjacent_vertex_indices(this_vert,adj_verts,err);  MSQ_ERRZERO(err);
  const unsigned num_sample_points=adj_verts.size();
  double *metric_values=new double[num_sample_points];
  const MsqVertex* verts = pd.get_vertex_array(err);  MSQ_ERRZERO(err);
  for (unsigned i = 0; i < num_sample_points; ++i) 
    metric_values[i] = (verts[this_vert] - verts[adj_verts[i]]).length();
  fval=average_metrics(metric_values,num_sample_points,err);  
  delete[] metric_values;
  return !MSQ_CHKERR(err);
  
}

bool EdgeLengthQualityMetric::evaluate( PatchData& pd, 
                                        size_t vertex, 
                                        double& value, 
                                        MsqError& err )
{
  std::vector<size_t> verts;
  bool rval = evaluate_common( pd, vertex, value, verts, err );
  return !MSQ_CHKERR(err) && rval;
}

bool EdgeLengthQualityMetric::evaluate_with_indices( PatchData& pd,
                                                     size_t vertex,
                                                     double& value,
                                                     std::vector<size_t>& indices,
                                                     MsqError& err )
{
  indices.clear();
  bool rval = evaluate_common( pd, vertex, value, indices, err );
  
  std::vector<size_t>::iterator r, w;
  for (r = w = indices.begin(); r != indices.end(); ++r) {
    if (*r < pd.num_free_vertices()) {
      *w = *r;
      ++w;
    }
  }
  indices.erase( w, indices.end() );
  if (vertex < pd.num_free_vertices())
    indices.push_back( vertex );
  
  return !MSQ_CHKERR(err) && rval;
}
