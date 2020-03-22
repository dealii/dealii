/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Sandia National Laboratories.  Developed at the
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
 
    (2006) kraftche@cae.wisc.edu
   
  ***************************************************************** */


/** \file TriLagrangeShape.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "TriLagrangeShape.hpp"
#include "MsqError.hpp"
#include <assert.h>

namespace MESQUITE_NS {

EntityTopology TriLagrangeShape::element_topology() const
  { return TRIANGLE; }
  
int TriLagrangeShape::num_nodes() const
  { return 6; }

NodeSet TriLagrangeShape::sample_points( NodeSet ns ) const
{
  if (ns.have_any_mid_node()) {
    ns.set_all_corner_nodes(TRIANGLE);
  }
  else {
    ns.clear();
    ns.set_mid_face_node(0);
  }
  return ns;
}

void TriLagrangeShape::coefficients( Sample loc,
                                     NodeSet nodeset,
                                     double* coeff_out,
                                     size_t* indices_out,
                                     size_t& num_coeff,
                                     MsqError& err ) const
{
  if (nodeset.have_any_mid_face_node()) {
    MSQ_SETERR(err)("TriLagrangeShape does not support mid-element nodes",
                    MsqError::UNSUPPORTED_ELEMENT);
    return;
  }
  
  switch (loc.dimension) {
    case 0:
      num_coeff = 1;
      indices_out[0] = loc.number;
      coeff_out[0] = 1.0;
      break;
    case 1:
      if (nodeset.mid_edge_node(loc.number)) { // if mid-edge node is present
        num_coeff = 1;
        indices_out[0] = 3+loc.number;
        coeff_out[0] = 1.0;
      }
      else { // no mid node on edge
        num_coeff = 2;
        indices_out[0] = loc.number;
        indices_out[1] = (loc.number+1)%3;
        coeff_out[0] = 0.5;
        coeff_out[1] = 0.5;
      }
      break;
    case 2:
      num_coeff = 3;
      indices_out[0] = 0;
      indices_out[1] = 1;
      indices_out[2] = 2;
      coeff_out[0] = 1.0/3.0;
      coeff_out[1] = 1.0/3.0;
      coeff_out[2] = 1.0/3.0;
      for (int i = 0; i < 3; ++i) { // for each mid-edge node
        if (nodeset.mid_edge_node(i)) { // if node is present
          indices_out[num_coeff] = i+3;
            // coeff for mid-edge node
          coeff_out[num_coeff] = 4.0/9.0;
            // adjust coeff for adj corner nodes
          coeff_out[i]       -= 2.0/9.0;
          coeff_out[(i+1)%3] -= 2.0/9.0;
            // update count
          ++num_coeff;
        }
      }
      break;
    default:
      MSQ_SETERR(err)(MsqError::UNSUPPORTED_ELEMENT,
                  "Request for dimension %d mapping function value"
                  "for a triangular element", loc.dimension);
  }
}

static inline void get_linear_derivatives( size_t* vertices,
                                           MsqVector<2>* derivs )
{
  vertices[0] = 0;
  vertices[1] = 1;
  vertices[2] = 2;
  derivs[0][0] = -1.0;
  derivs[0][1] = -1.0;
  derivs[1][0] =  1.0;
  derivs[1][1] =  0.0;
  derivs[2][0] =  0.0;
  derivs[2][1] =  1.0;
}  

static void derivatives_at_corner( unsigned corner,
                                   NodeSet nodeset,
                                   size_t* vertices,
                                   MsqVector<2>* derivs,
                                   size_t& num_vtx )
{
  num_vtx = 3;
  get_linear_derivatives( vertices, derivs );
  switch (corner) {
    case 0:
      if (nodeset.mid_edge_node(0)) {
        vertices[num_vtx] = 3;
        derivs[num_vtx][0] = 4.0;
        derivs[num_vtx][1] = 0.0;
        ++num_vtx;
        derivs[0][0] -= 2.0;
        derivs[1][0] -= 2.0;
      }
      if (nodeset.mid_edge_node(2)) {
        vertices[num_vtx] = 5;
        derivs[num_vtx][0] = 0.0;
        derivs[num_vtx][1] = 4.0;
        ++num_vtx;
        derivs[0][1] -= 2.0;
        derivs[2][1] -= 2.0;
      }
      break;
    
    case 1:
      if (nodeset.mid_edge_node(0)) {
        vertices[num_vtx] = 3;
        derivs[num_vtx][0] = -4.0;
        derivs[num_vtx][1] = -4.0;
        ++num_vtx;
        derivs[0][0] += 2.0;
        derivs[0][1] += 2.0;
        derivs[1][0] += 2.0;
        derivs[1][1] += 2.0;
      }
      if (nodeset.mid_edge_node(1)) {
        vertices[num_vtx] = 4;
        derivs[num_vtx][0] = 0.0;
        derivs[num_vtx][1] = 4.0;
        ++num_vtx;
        derivs[1][1] -= 2.0;
        derivs[2][1] -= 2.0;
      }
      break;
    
    case 2:
      if (nodeset.mid_edge_node(1)) {
        vertices[num_vtx] = 4;
        derivs[num_vtx][0] = 4.0;
        derivs[num_vtx][1] = 0.0;
        ++num_vtx;
        derivs[1][0] -= 2.0;
        derivs[2][0] -= 2.0;
      }
      if (nodeset.mid_edge_node(2)) {
        vertices[num_vtx] = 5;
        derivs[num_vtx][0] = -4.0;
        derivs[num_vtx][1] = -4.0;
        ++num_vtx;
        derivs[0][0] += 2.0;
        derivs[0][1] += 2.0;
        derivs[2][0] += 2.0;
        derivs[2][1] += 2.0;
      }
      break;
  }
}

static const double edr[3][3] = { { 0.0,-2.0, 2.0 },
                                  { 0.0, 2.0, 2.0 },
                                  { 0.0,-2.0,-2.0 } };
static const double eds[3][3] = { {-2.0,-2.0, 0.0 },
                                  { 2.0, 2.0, 0.0 },
                                  { 2.0,-2.0, 0.0 } };

static void derivatives_at_mid_edge( unsigned edge, 
                                     NodeSet nodeset,
                                     size_t* vertices,
                                     MsqVector<2>* derivs,
                                     size_t& num_vtx )
{
    // The mid-edge behavior is rather strange.
    // A corner vertex contributes to the jacobian
    // at the mid-point of the opposite edge unless
    // one, but *not* both of the adjacent mid-edge
    // nodes is present.
    
    // The value for each corner is incremented when 
    // a mid-side node is present.  If the value is
    // 0 when finished, the corner doesn't contribute.
    // Initialize values to 0 for corners adjacent to
    // edge so they are never zero.
  int corner_count[3] = { 1, 1, 1 };
  corner_count[(edge+2)%3] = -1;
  
    // begin with derivatives for linear tri
  double corner_derivs[3][2] = { {-1.0,-1.0 },
                                 { 1.0, 0.0 },
                                 { 0.0, 1.0 } };
  
    // do mid-side nodes first
  num_vtx = 0;
  for (unsigned i = 0; i < 3; ++i) {
    if (nodeset.mid_edge_node(i)) {
      vertices[num_vtx] =  i+3 ;
      derivs[num_vtx][0] =  edr[i][edge] ;
      derivs[num_vtx][1] =  eds[i][edge] ;
      ++num_vtx;

      int a = (i+1)%3;
      corner_derivs[i][0] -= 0.5 * edr[i][edge];
      corner_derivs[i][1] -= 0.5 * eds[i][edge];
      corner_derivs[a][0] -= 0.5 * edr[i][edge];
      corner_derivs[a][1] -= 0.5 * eds[i][edge];
      ++corner_count[i];
      ++corner_count[a];
    }
  }
  
    // now add corner nodes to list
  for (unsigned i = 0; i < 3; ++i) {
    if (corner_count[i]) {
      vertices[num_vtx] = i;
      derivs[num_vtx][0] = corner_derivs[i][0];
      derivs[num_vtx][1] = corner_derivs[i][1];
      ++num_vtx;
    }
  }
}

static const double fdr[] = { 0.0,     4.0/3.0, -4.0/3.0 };
static const double fds[] = {-4.0/3.0, 4.0/3.0,  0.0     };
       
static void derivatives_at_mid_elem( NodeSet nodeset,
                                     size_t* vertices,
                                     MsqVector<2>* derivs,
                                     size_t& num_vtx )
{
  get_linear_derivatives( vertices, derivs );
  num_vtx = 3;
  for (unsigned i = 0; i < 3; ++i) {
    if (nodeset.mid_edge_node(i)) {
      vertices[num_vtx] = i+3;
      derivs[num_vtx][0] = fdr[i];
      derivs[num_vtx][1] = fds[i];
      ++num_vtx;
       
      int a = (i+1)%3;
      derivs[i][0] -= 0.5 * fdr[i];
      derivs[i][1] -= 0.5 * fds[i];
      derivs[a][0] -= 0.5 * fdr[i];
      derivs[a][1] -= 0.5 * fds[i];
    }
  }
}

void TriLagrangeShape::derivatives( Sample loc,
                                    NodeSet nodeset,
                                    size_t* vertex_indices_out,
                                    MsqVector<2>* d_coeff_d_xi_out,
                                    size_t& num_vtx,
                                    MsqError& err ) const
{
  if (!nodeset.have_any_mid_node()) {
    num_vtx = 3;
    get_linear_derivatives( vertex_indices_out, d_coeff_d_xi_out );
    return;
  }

  if (nodeset.have_any_mid_face_node()) {
    MSQ_SETERR(err)("TriLagrangeShape does not support mid-element nodes",
                    MsqError::UNSUPPORTED_ELEMENT);
    return;
  }
  
  switch (loc.dimension) {
    case 0:
      derivatives_at_corner( loc.number, nodeset, vertex_indices_out, d_coeff_d_xi_out, num_vtx );
      break;
    case 1:
      derivatives_at_mid_edge( loc.number, nodeset, vertex_indices_out, d_coeff_d_xi_out, num_vtx );
      break;
    case 2:
      derivatives_at_mid_elem( nodeset, vertex_indices_out, d_coeff_d_xi_out, num_vtx );
      break;
    default:
      MSQ_SETERR(err)("Invalid/unsupported logical dimension",MsqError::INVALID_ARG);
  }
}


void TriLagrangeShape::ideal( Sample , 
                              MsqMatrix<3,2>& J,
                              MsqError&  ) const
{
//  const double g = 1.074569931823542; // sqrt(2/sqrt(3));
//  J(0,0) = g;    J(0,1) = 0.5*g;
//  J(1,0) = 0.0;  J(1,1) = 1.0/g;
//  J(2,0) = 0.0;  J(2,1) = 0.0;
    const double a = 0.70710678118654746; // 1/sqrt(2)
    const double b = 1.3160740129524924;  // 4th root of 3
    J(0,0) = -a/b; J(0,1) =  a/b;
    J(1,0) = -a*b; J(1,1) = -a*b;
    J(2,0) =  0.0; J(2,1) =  0.0;
}


} // namespace Mesquite
