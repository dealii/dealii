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
  \file   IdealWeightMeanRatio.cpp
  \brief  

  \author Michael Brewer
  \date   2002-11-11
*/
#include "IdealWeightMeanRatio.hpp"
#include "MeanRatioFunctions.hpp"
#include "Vector3D.hpp"
#include "MsqTimer.hpp"
#include "MsqDebug.hpp"
#include "MsqError.hpp"
#include "PatchData.hpp"

#include <math.h>
#include <vector>
using std::vector;

using namespace Mesquite;

std::string IdealWeightMeanRatio::get_name() const
  { return "Mean Ratio"; }

int IdealWeightMeanRatio::get_negate_flag() const
  { return -1; }

bool IdealWeightMeanRatio::evaluate( PatchData& pd, 
                                     size_t handle, 
                                     double& m, 
                                     MsqError& err )
{
  const MsqMeshEntity* e = &pd.element_by_index(handle);
  EntityTopology topo = e->get_element_type();

  const MsqVertex *vertices = pd.get_vertex_array(err);
  const size_t *v_i = e->get_vertex_index_array();

  Vector3D n;			// Surface normal for 2D objects

  // Prism and Hex element descriptions
  static const int locs_pri[6][4] = {{0, 1, 2, 3}, {1, 2, 0, 4},
				     {2, 0, 1, 5}, {3, 5, 4, 0},
				     {4, 3, 5, 1}, {5, 4, 3, 2}};
  static const int locs_hex[8][4] = {{0, 1, 3, 4}, {1, 2, 0, 5},
				     {2, 3, 1, 6}, {3, 0, 2, 7},
				     {4, 7, 5, 0}, {5, 4, 6, 1},
				     {6, 5, 7, 2}, {7, 6, 4, 3}};

  const Vector3D d_con(1.0, 1.0, 1.0);

  int i;

  m = 0.0;
  bool metric_valid = false;
  switch(topo) {
  case TRIANGLE:
    pd.get_domain_normal_at_element(e, n, err); MSQ_ERRZERO(err);
    n = n / n.length();		// Need unit normal
    mCoords[0] = vertices[v_i[0]];
    mCoords[1] = vertices[v_i[1]];
    mCoords[2] = vertices[v_i[2]];
    metric_valid = m_fcn_2e(m, mCoords, n, a2Con, b2Con, c2Con);
    if (!metric_valid) return false;
    break;
    
  case QUADRILATERAL:
    pd.get_domain_normal_at_element(e, n, err); MSQ_ERRZERO(err);
    for (i = 0; i < 4; ++i) {
      n = n / n.length();	// Need unit normal
      mCoords[0] = vertices[v_i[locs_hex[i][0]]];
      mCoords[1] = vertices[v_i[locs_hex[i][1]]];
      mCoords[2] = vertices[v_i[locs_hex[i][2]]];
      metric_valid = m_fcn_2i(mMetrics[i], mCoords, n, 
			      a2Con, b2Con, c2Con, d_con);
      if (!metric_valid) return false;
    }
    m = average_metrics(mMetrics, 4, err); MSQ_ERRZERO(err);
    break;

  case TETRAHEDRON:
    mCoords[0] = vertices[v_i[0]];
    mCoords[1] = vertices[v_i[1]];
    mCoords[2] = vertices[v_i[2]];
    mCoords[3] = vertices[v_i[3]];
    metric_valid = m_fcn_3e(m, mCoords, a3Con, b3Con, c3Con);
    if (!metric_valid) return false;
    break;

  case PYRAMID:
    for (i = 0; i < 4; ++i) {
      mCoords[0] = vertices[v_i[ i     ]];
      mCoords[1] = vertices[v_i[(i+1)%4]];
      mCoords[2] = vertices[v_i[(i+3)%4]];
      mCoords[3] = vertices[v_i[ 4     ]];
      metric_valid = m_fcn_3p(mMetrics[i], mCoords, a3Con, b3Con, c3Con);
      if (!metric_valid) return false;
    }
    m = average_metrics(mMetrics, 4, err); MSQ_ERRZERO(err);
    break;

  case PRISM:
    for (i = 0; i < 6; ++i) {
      mCoords[0] = vertices[v_i[locs_pri[i][0]]];
      mCoords[1] = vertices[v_i[locs_pri[i][1]]];
      mCoords[2] = vertices[v_i[locs_pri[i][2]]];
      mCoords[3] = vertices[v_i[locs_pri[i][3]]];
      metric_valid = m_fcn_3w(mMetrics[i], mCoords, a3Con, b3Con, c3Con);
      if (!metric_valid) return false;
    }
    m = average_metrics(mMetrics, 6, err); MSQ_ERRZERO(err);
    break;

  case HEXAHEDRON:
    for (i = 0; i < 8; ++i) {
      mCoords[0] = vertices[v_i[locs_hex[i][0]]];
      mCoords[1] = vertices[v_i[locs_hex[i][1]]];
      mCoords[2] = vertices[v_i[locs_hex[i][2]]];
      mCoords[3] = vertices[v_i[locs_hex[i][3]]];
      metric_valid = m_fcn_3i(mMetrics[i], mCoords, 
			      a3Con, b3Con, c3Con, d_con);
      if (!metric_valid) return false;
    }
    m = average_metrics(mMetrics, 8, err); MSQ_ERRZERO(err);
    break;

  default:
    MSQ_SETERR(err)(MsqError::UNSUPPORTED_ELEMENT,
                    "Element type (%d) not supported in IdealWeightMeanRatio",
                    (int)topo);
    return false;
  } // end switch over element type
  return true;
}

bool IdealWeightMeanRatio::evaluate_with_gradient( PatchData& pd,
                    size_t handle,
                    double& m,
                    std::vector<size_t>& indices,
                    std::vector<Vector3D>& g,
                    MsqError& err )
{
//  FUNCTION_TIMER_START(__FUNC__);
  const MsqMeshEntity* e = &pd.element_by_index(handle);
  EntityTopology topo = e->get_element_type();

  if (!analytical_average_gradient() &&
      topo != TRIANGLE &&
      topo != TETRAHEDRON) {
    static bool print = true;
    if (print) {
      MSQ_DBGOUT(1) << "Analyical gradient not available for selected averaging scheme. "
                    << "Using (possibly much slower) numerical approximation of gradient"
                    << " of quality metric. " << std::endl;
      print = false;
    }
    return QualityMetric::evaluate_with_gradient( pd, handle, m, indices, g, err );
  }

  const MsqVertex *vertices = pd.get_vertex_array(err);
  const size_t *v_i = e->get_vertex_index_array();

  Vector3D n;			// Surface normal for 2D objects

  //double   nm, t=0;

  // Prism and Hex element descriptions
  static const int locs_pri[6][4] = {{0, 1, 2, 3}, {1, 2, 0, 4},
				     {2, 0, 1, 5}, {3, 5, 4, 0},
				     {4, 3, 5, 1}, {5, 4, 3, 2}};
  static const int locs_hex[8][4] = {{0, 1, 3, 4}, {1, 2, 0, 5},
				     {2, 3, 1, 6}, {3, 0, 2, 7},
				     {4, 7, 5, 0}, {5, 4, 6, 1},
				     {6, 5, 7, 2}, {7, 6, 4, 3}};

  const Vector3D d_con(1.0, 1.0, 1.0);

  int i;

  bool metric_valid = false;
  const uint32_t fm = fixed_vertex_bitmap( pd, e, indices );

  m = 0.0;

  switch(topo) {
  case TRIANGLE:
    pd.get_domain_normal_at_element(e, n, err); MSQ_ERRZERO(err);
    n = n / n.length();		// Need unit normal
    mCoords[0] = vertices[v_i[0]];
    mCoords[1] = vertices[v_i[1]];
    mCoords[2] = vertices[v_i[2]];
    g.resize(3);
    if (!g_fcn_2e(m, arrptr(g), mCoords, n, a2Con, b2Con, c2Con)) return false;
    break;

  case QUADRILATERAL:
    pd.get_domain_normal_at_element(e, n, err); MSQ_ERRZERO(err);
    n /= n.length();	// Need unit normal
    for (i = 0; i < 4; ++i) {
      mCoords[0] = vertices[v_i[locs_hex[i][0]]];
      mCoords[1] = vertices[v_i[locs_hex[i][1]]];
      mCoords[2] = vertices[v_i[locs_hex[i][2]]];
      if (!g_fcn_2i(mMetrics[i], mGradients+3*i, mCoords, n,
		    a2Con, b2Con, c2Con, d_con)) return false;
    }
    
    g.resize(4);
    m = average_corner_gradients( QUADRILATERAL, fm, 4,
                                  mMetrics, mGradients, 
                                  arrptr(g), err ); MSQ_ERRZERO(err);
    break;

  case TETRAHEDRON:
    mCoords[0] = vertices[v_i[0]];
    mCoords[1] = vertices[v_i[1]];
    mCoords[2] = vertices[v_i[2]];
    mCoords[3] = vertices[v_i[3]];
    g.resize(4);
    metric_valid = g_fcn_3e(m, arrptr(g), mCoords, a3Con, b3Con, c3Con);
    if (!metric_valid) return false;
    break;

  case PYRAMID:
    for (i = 0; i < 4; ++i) {
      mCoords[0] = vertices[v_i[ i     ]];
      mCoords[1] = vertices[v_i[(i+1)%4]];
      mCoords[2] = vertices[v_i[(i+3)%4]];
      mCoords[3] = vertices[v_i[ 4     ]];
      metric_valid = g_fcn_3p(mMetrics[i], mGradients+4*i, mCoords, a3Con, b3Con, c3Con);
      if (!metric_valid) return false;
    }
    g.resize(5);
    m = average_corner_gradients( PYRAMID, fm, 4,
                                  mMetrics, mGradients,
                                  arrptr(g), err ); MSQ_ERRZERO(err);
    break;
    
  case PRISM:
    for (i = 0; i < 6; ++i) {
      mCoords[0] = vertices[v_i[locs_pri[i][0]]];
      mCoords[1] = vertices[v_i[locs_pri[i][1]]];
      mCoords[2] = vertices[v_i[locs_pri[i][2]]];
      mCoords[3] = vertices[v_i[locs_pri[i][3]]];
      if (!g_fcn_3w(mMetrics[i], mGradients+4*i, mCoords, 
		    a3Con, b3Con, c3Con)) return false;
    }
    g.resize(6);
    m = average_corner_gradients( PRISM, fm, 6,
                                  mMetrics, mGradients,
                                  arrptr(g), err ); MSQ_ERRZERO(err);
    break;

  case HEXAHEDRON:
    for (i = 0; i < 8; ++i) {
      mCoords[0] = vertices[v_i[locs_hex[i][0]]];
      mCoords[1] = vertices[v_i[locs_hex[i][1]]];
      mCoords[2] = vertices[v_i[locs_hex[i][2]]];
      mCoords[3] = vertices[v_i[locs_hex[i][3]]];
      if (!g_fcn_3i(mMetrics[i], mGradients+4*i, mCoords, 
		    a3Con, b3Con, c3Con, d_con)) return false;
    }
    g.resize(8);
    m = average_corner_gradients( HEXAHEDRON, fm, 8,
                                  mMetrics, mGradients,
                                  arrptr(g), err ); MSQ_ERRZERO(err);
     break;

  default:
    MSQ_SETERR(err)(MsqError::UNSUPPORTED_ELEMENT,
                    "Element type (%d) not supported in IdealWeightMeanRatio",
                    (int)topo);
    return false;
  }

  remove_fixed_gradients( topo, fm, g );
  return true;
}

bool IdealWeightMeanRatio::evaluate_with_Hessian_diagonal( PatchData& pd,
                    size_t handle,
                    double& m,
                    std::vector<size_t>& indices,
                    std::vector<Vector3D>& g,
                    std::vector<SymMatrix3D>& h,
                    MsqError& err )
{
  const MsqMeshEntity* e = &pd.element_by_index(handle);
  EntityTopology topo = e->get_element_type();

  if (!analytical_average_hessian() &&
      topo != TRIANGLE &&
      topo != TETRAHEDRON) {
    static bool print = true;
    if (print) {
      MSQ_DBGOUT(1) << "Analyical gradient not available for selected averaging scheme. "
                    << "Using (possibly much slower) numerical approximation of gradient"
                    << " of quality metric. " << std::endl;
      print = false;
    }
    return QualityMetric::evaluate_with_Hessian_diagonal( pd, handle, m, indices, g, h, err );
  }

  const MsqVertex *vertices = pd.get_vertex_array(err);
  const size_t *v_i = e->get_vertex_index_array();


  Vector3D n;			// Surface normal for 2D objects

  // Prism and Hex element descriptions
  static const int locs_pri[6][4] = {{0, 1, 2, 3}, {1, 2, 0, 4},
				     {2, 0, 1, 5}, {3, 5, 4, 0},
				     {4, 3, 5, 1}, {5, 4, 3, 2}};
  static const int locs_hex[8][4] = {{0, 1, 3, 4}, {1, 2, 0, 5},
				     {2, 3, 1, 6}, {3, 0, 2, 7},
				     {4, 7, 5, 0}, {5, 4, 6, 1},
				     {6, 5, 7, 2}, {7, 6, 4, 3}};

  const Vector3D d_con(1.0, 1.0, 1.0);

  int i;

  bool metric_valid = false;
  const uint32_t fm = fixed_vertex_bitmap( pd, e, indices );

  m = 0.0;

  switch(topo) {
  case TRIANGLE:
    pd.get_domain_normal_at_element(e, n, err); MSQ_ERRZERO(err);
    n = n / n.length();		// Need unit normal
    mCoords[0] = vertices[v_i[0]];
    mCoords[1] = vertices[v_i[1]];
    mCoords[2] = vertices[v_i[2]];
    g.resize(3), h.resize(3);
    if (!h_fcn_2e(m, arrptr(g), mHessians, mCoords, n, a2Con, b2Con, c2Con)) return false;
    h[0] = mHessians[0].upper();
    h[1] = mHessians[3].upper();
    h[2] = mHessians[5].upper();
    break;

  case QUADRILATERAL:
    pd.get_domain_normal_at_element(e, n, err); MSQ_ERRZERO(err);
    n = n / n.length();	// Need unit normal
    for (i = 0; i < 4; ++i) {
      mCoords[0] = vertices[v_i[locs_hex[i][0]]];
      mCoords[1] = vertices[v_i[locs_hex[i][1]]];
      mCoords[2] = vertices[v_i[locs_hex[i][2]]];
      if (!h_fcn_2i(mMetrics[i], mGradients+3*i, mHessians+6*i, mCoords, n,
		    a2Con, b2Con, c2Con, d_con)) return false;
    }

    g.resize(4), h.resize(4);
    m = average_corner_hessian_diagonals( QUADRILATERAL, fm, 4,
                                 mMetrics, mGradients, mHessians,
                                 arrptr(g), arrptr(h), err );
    MSQ_ERRZERO( err );
    break;

  case TETRAHEDRON:
    mCoords[0] = vertices[v_i[0]];
    mCoords[1] = vertices[v_i[1]];
    mCoords[2] = vertices[v_i[2]];
    mCoords[3] = vertices[v_i[3]];
    g.resize(4), h.resize(4);
    metric_valid = h_fcn_3e(m, arrptr(g), mHessians, mCoords, a3Con, b3Con, c3Con);
    if (!metric_valid) return false;
    h[0] = mHessians[0].upper();
    h[1] = mHessians[4].upper();
    h[2] = mHessians[7].upper();
    h[3] = mHessians[9].upper();
    break;

  case PYRAMID:
    for (i = 0; i < 4; ++i) {
      mCoords[0] = vertices[v_i[ i     ]];
      mCoords[1] = vertices[v_i[(i+1)%4]];
      mCoords[2] = vertices[v_i[(i+3)%4]];
      mCoords[3] = vertices[v_i[ 4     ]];
      metric_valid = h_fcn_3p(mMetrics[i], mGradients+4*i, 
                              mHessians+10*i, mCoords, a3Con, b3Con, c3Con);
      if (!metric_valid) return false;
    }

    g.resize(5), h.resize(5);
    m = average_corner_hessian_diagonals( PYRAMID, fm, 4,
                                 mMetrics, mGradients, mHessians,
                                 arrptr(g), arrptr(h), err );
    MSQ_ERRZERO( err );
    break;

  case PRISM:
    for (i = 0; i < 6; ++i) {
      mCoords[0] = vertices[v_i[locs_pri[i][0]]];
      mCoords[1] = vertices[v_i[locs_pri[i][1]]];
      mCoords[2] = vertices[v_i[locs_pri[i][2]]];
      mCoords[3] = vertices[v_i[locs_pri[i][3]]];
      if (!h_fcn_3w(mMetrics[i], mGradients+4*i, mHessians+10*i, mCoords,
		    a3Con, b3Con, c3Con)) return false;
    }

    g.resize(6), h.resize(6);
    m = average_corner_hessian_diagonals( PRISM, fm, 6,
                                 mMetrics, mGradients, mHessians,
                                 arrptr(g), arrptr(h), err );
    MSQ_ERRZERO( err );
    break;

  case HEXAHEDRON:
    for (i = 0; i < 8; ++i) {
      mCoords[0] = vertices[v_i[locs_hex[i][0]]];
      mCoords[1] = vertices[v_i[locs_hex[i][1]]];
      mCoords[2] = vertices[v_i[locs_hex[i][2]]];
      mCoords[3] = vertices[v_i[locs_hex[i][3]]];
      if (!h_fcn_3i(mMetrics[i], mGradients+4*i, mHessians+10*i, mCoords,
		    a3Con, b3Con, c3Con, d_con)) return false;
    }

    g.resize(8), h.resize(8);
    m = average_corner_hessian_diagonals( HEXAHEDRON, fm, 8,
                                 mMetrics, mGradients, mHessians,
                                 arrptr(g), arrptr(h), err );
    MSQ_ERRZERO( err );
    break;

  default:
    MSQ_SETERR(err)(MsqError::UNSUPPORTED_ELEMENT,
                    "Element type (%d) not supported in IdealWeightMeanRatio",
                    (int)topo);
    return false;
  } // end switch over element type

  remove_fixed_diagonals( topo, fm, g, h );
  return true;
}


bool IdealWeightMeanRatio::evaluate_with_Hessian( PatchData& pd,
                    size_t handle,
                    double& m,
                    std::vector<size_t>& indices,
                    std::vector<Vector3D>& g,
                    std::vector<Matrix3D>& h,
                    MsqError& err )
{
//  FUNCTION_TIMER_START(__FUNC__);
  const MsqMeshEntity* e = &pd.element_by_index(handle);
  EntityTopology topo = e->get_element_type();

  if (!analytical_average_hessian() &&
      topo != TRIANGLE &&
      topo != TETRAHEDRON) {
    static bool print = true;
    if (print) {
      MSQ_DBGOUT(1) << "Analyical gradient not available for selected averaging scheme. "
                    << "Using (possibly much slower) numerical approximation of gradient"
                    << " of quality metric. " << std::endl;
      print = false;
    }
    return QualityMetric::evaluate_with_Hessian( pd, handle, m, indices, g, h, err );
  }

  const MsqVertex *vertices = pd.get_vertex_array(err);
  const size_t *v_i = e->get_vertex_index_array();


  Vector3D n;			// Surface normal for 2D objects

  // Prism and Hex element descriptions
  static const int locs_pri[6][4] = {{0, 1, 2, 3}, {1, 2, 0, 4},
				     {2, 0, 1, 5}, {3, 5, 4, 0},
				     {4, 3, 5, 1}, {5, 4, 3, 2}};
  static const int locs_hex[8][4] = {{0, 1, 3, 4}, {1, 2, 0, 5},
				     {2, 3, 1, 6}, {3, 0, 2, 7},
				     {4, 7, 5, 0}, {5, 4, 6, 1},
				     {6, 5, 7, 2}, {7, 6, 4, 3}};

  const Vector3D d_con(1.0, 1.0, 1.0);

  int i;

  bool metric_valid = false;
  const uint32_t fm = fixed_vertex_bitmap( pd, e, indices );

  m = 0.0;

  switch(topo) {
  case TRIANGLE:
    pd.get_domain_normal_at_element(e, n, err); MSQ_ERRZERO(err);
    n = n / n.length();		// Need unit normal
    mCoords[0] = vertices[v_i[0]];
    mCoords[1] = vertices[v_i[1]];
    mCoords[2] = vertices[v_i[2]];
    g.resize(3), h.resize(6);
    if (!h_fcn_2e(m, arrptr(g), arrptr(h), mCoords, n, a2Con, b2Con, c2Con)) return false;
    break;

  case QUADRILATERAL:
    pd.get_domain_normal_at_element(e, n, err); MSQ_ERRZERO(err);
    n = n / n.length();	// Need unit normal
    for (i = 0; i < 4; ++i) {
      mCoords[0] = vertices[v_i[locs_hex[i][0]]];
      mCoords[1] = vertices[v_i[locs_hex[i][1]]];
      mCoords[2] = vertices[v_i[locs_hex[i][2]]];
      if (!h_fcn_2i(mMetrics[i], mGradients+3*i, mHessians+6*i, mCoords, n,
		    a2Con, b2Con, c2Con, d_con)) return false;
    }

    g.resize(4), h.resize(10);
    m = average_corner_hessians( QUADRILATERAL, fm, 4,
                                 mMetrics, mGradients, mHessians,
                                 arrptr(g), arrptr(h), err );
    MSQ_ERRZERO( err );
    break;

  case TETRAHEDRON:
    mCoords[0] = vertices[v_i[0]];
    mCoords[1] = vertices[v_i[1]];
    mCoords[2] = vertices[v_i[2]];
    mCoords[3] = vertices[v_i[3]];
    g.resize(4), h.resize(10);
    metric_valid = h_fcn_3e(m, arrptr(g), arrptr(h), mCoords, a3Con, b3Con, c3Con);
    if (!metric_valid) return false;
    break;

  case PYRAMID:
    for (i = 0; i < 4; ++i) {
      mCoords[0] = vertices[v_i[ i     ]];
      mCoords[1] = vertices[v_i[(i+1)%4]];
      mCoords[2] = vertices[v_i[(i+3)%4]];
      mCoords[3] = vertices[v_i[ 4     ]];
      metric_valid = h_fcn_3p(mMetrics[i], mGradients+4*i, 
                              mHessians+10*i, mCoords, a3Con, b3Con, c3Con);
      if (!metric_valid) return false;
    }

    g.resize(5), h.resize(15);
    m = average_corner_hessians( PYRAMID, fm, 4,
                                 mMetrics, mGradients, mHessians,
                                 arrptr(g), arrptr(h), err );
    MSQ_ERRZERO( err );
    break;

  case PRISM:
    for (i = 0; i < 6; ++i) {
      mCoords[0] = vertices[v_i[locs_pri[i][0]]];
      mCoords[1] = vertices[v_i[locs_pri[i][1]]];
      mCoords[2] = vertices[v_i[locs_pri[i][2]]];
      mCoords[3] = vertices[v_i[locs_pri[i][3]]];
      if (!h_fcn_3w(mMetrics[i], mGradients+4*i, mHessians+10*i, mCoords,
		    a3Con, b3Con, c3Con)) return false;
    }

    g.resize(6), h.resize(21);
    m = average_corner_hessians( PRISM, fm, 6,
                                 mMetrics, mGradients, mHessians,
                                 arrptr(g), arrptr(h), err );
    MSQ_ERRZERO( err );
    break;

  case HEXAHEDRON:
    for (i = 0; i < 8; ++i) {
      mCoords[0] = vertices[v_i[locs_hex[i][0]]];
      mCoords[1] = vertices[v_i[locs_hex[i][1]]];
      mCoords[2] = vertices[v_i[locs_hex[i][2]]];
      mCoords[3] = vertices[v_i[locs_hex[i][3]]];
      if (!h_fcn_3i(mMetrics[i], mGradients+4*i, mHessians+10*i, mCoords,
		    a3Con, b3Con, c3Con, d_con)) return false;
    }

    g.resize(8), h.resize(36);
    m = average_corner_hessians( HEXAHEDRON, fm, 8,
                                 mMetrics, mGradients, mHessians,
                                 arrptr(g), arrptr(h), err );
    MSQ_ERRZERO( err );
    break;

  default:
    MSQ_SETERR(err)(MsqError::UNSUPPORTED_ELEMENT,
                    "Element type (%d) not supported in IdealWeightMeanRatio",
                    (int)topo);
    return false;
  } // end switch over element type

  remove_fixed_gradients( topo, fm, g );
  remove_fixed_hessians( topo, fm, h );
//  FUNCTION_TIMER_END();
  return true;
}
