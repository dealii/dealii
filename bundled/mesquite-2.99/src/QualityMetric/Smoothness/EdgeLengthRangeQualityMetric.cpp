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
/*! \file EdgeLengthRangeQualityMetric.cpp
  \author Michael Brewer
  \date 2002-05-14
  Evaluates the lengths of the edges attached to the given vertex.
  By default, the averaging method is set to SUM.
*/


#include "EdgeLengthRangeQualityMetric.hpp"
#include "Vector3D.hpp"
#include "QualityMetric.hpp"
#include "MsqVertex.hpp"
#include "PatchData.hpp"
#include "MsqDebug.hpp"
#include "MsqError.hpp"

#include <vector>
using std::vector;

using namespace Mesquite;

EdgeLengthRangeQualityMetric::EdgeLengthRangeQualityMetric( double low_a, double high_a )
  : AveragingQM(SUM),
    highVal(high_a),
    lowVal(low_a)
{
  if (lowVal > highVal)
    std::swap( lowVal, highVal );
}

EdgeLengthRangeQualityMetric::~EdgeLengthRangeQualityMetric()
{}

std::string EdgeLengthRangeQualityMetric::get_name() const
  { return "Edge Length Range Metric"; }
  
int EdgeLengthRangeQualityMetric::get_negate_flag() const
  { return 1; }

/*!For the given vertex, vert, with connected edges of lengths l_j for
  j=1...k, the metric value is the average (where the default average
  type is SUM) of
        u_j = ( | l_j - lowVal | - (l_j - lowVal) )^2 +
              ( | highVal - l_j | - (highVal - l_j) )^2.
*/
bool EdgeLengthRangeQualityMetric::evaluate_common(PatchData &pd, 
                                             size_t this_vert,
                                             double &fval, 
                                             std::vector<size_t>& adj_verts,
                                             MsqError &err)
{
  fval=0.0;
  Vector3D edg;
  pd.get_adjacent_vertex_indices(this_vert,adj_verts,err);  MSQ_ERRZERO(err);
  int num_sample_points=adj_verts.size();
  double *metric_values=new double[num_sample_points];
  const MsqVertex* verts = pd.get_vertex_array(err);  MSQ_ERRZERO(err);
    //store the length of the edge, and the first and second component of
    //metric values, respectively.
  double temp_length=0.0;
  double temp_first=0.0;
  double temp_second=0.0;
    //PRINT_INFO("INSIDE ELR, vertex = %f,%f,%f\n",verts[this_vert][0],verts[this_vert][1],verts[this_vert][2]);
    //loop while there are still more adjacent vertices.
  for (unsigned i = 0; i < adj_verts.size(); ++i) 
  {
    edg = verts[this_vert] - verts[adj_verts[i]];
      //compute the edge length
    temp_length=edg.length();
      //get the first component
    temp_first = temp_length - lowVal;
    temp_first = fabs(temp_first) - (temp_first);
    temp_first*=temp_first;
      //get the second component
    temp_second = highVal - temp_length;
    temp_second = fabs(temp_second) - (temp_second);
    temp_second*=temp_second;
      //combine the two components
    metric_values[i]=temp_first+temp_second;
  }
    //average the metric values of the edges
  fval=average_metrics(metric_values,num_sample_points,err);
    //clean up
  delete[] metric_values;
    //always return true because mesh is always valid wrt this metric.
  return !MSQ_CHKERR(err);
  
}


bool EdgeLengthRangeQualityMetric::evaluate( PatchData& pd, 
                                        size_t vertex, 
                                        double& value, 
                                        MsqError& err )
{
  std::vector<size_t> verts;
  bool rval = evaluate_common( pd, vertex, value, verts, err );
  return !MSQ_CHKERR(err) && rval;
}

bool EdgeLengthRangeQualityMetric::evaluate_with_indices( PatchData& pd,
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
