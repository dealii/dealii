/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Lawrence Livermore National Laboratory.  Under 
    the terms of Contract B545069 with the University of Wisconsin -- 
    Madison, Lawrence Livermore National Laboratory retains certain
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

    kraftche@cae.wisc.edu    

  ***************************************************************** */
/** \file QuadLagrangeShape.cpp
 *  \author Jason Kraftcheck
 */
 
#include "QuadLagrangeShape.hpp"
#include "LinearQuadrilateral.hpp"
#include "MsqError.hpp"

namespace MESQUITE_NS {

EntityTopology QuadLagrangeShape::element_topology() const
  { return QUADRILATERAL; }
  
int QuadLagrangeShape::num_nodes() const
  { return 9; }

void QuadLagrangeShape::coefficients( Sample loc,
                                      NodeSet nodeset,
                                      double* coeff_out,
                                      size_t* indices_out,
                                      size_t& num_coeff,
                                      MsqError& err ) const
{
  if (!nodeset.have_any_mid_node()) {
    LinearQuadrilateral::coefficients( loc, nodeset, coeff_out, indices_out, num_coeff );
    return;
  }

  switch (loc.dimension) {
    case 0:
      num_coeff = 1;
      indices_out[0] = loc.number;
      coeff_out[0] = 1.0;
      break;
    case 1:
      coeff_out[0] = coeff_out[1] = coeff_out[2] =
      coeff_out[3] = coeff_out[4] = coeff_out[5] = 
      coeff_out[6] = coeff_out[7] = coeff_out[8] = 0.0;
      if (nodeset.mid_edge_node(loc.number)) {  
          // if mid-edge node is present
        num_coeff = 1;
        indices_out[0] = loc.number+4;
        coeff_out[0] = 1.0;
      }
      else {
          // If mid-edge node is not present, mapping function value
          // for linear edge is even weight of adjacent vertices.
        num_coeff = 2;
        indices_out[0] = loc.number;
        indices_out[1] = (loc.number+1)%4;
        coeff_out[0] = 0.5;
        coeff_out[1] = 0.5;
      }
      break;
    case 2:
      if (nodeset.mid_face_node(0)) { // if quad center node is present
        num_coeff = 1;
        indices_out[0] = 8;
        coeff_out[0] = 1.0;
      } 
      else {
          // for linear element, (no mid-edge nodes), all corners contribute 1/4.
        num_coeff = 4;
        indices_out[0] = 0;
        indices_out[1] = 1;
        indices_out[2] = 2;
        indices_out[3] = 3;
        coeff_out[0] = 0.25;
        coeff_out[1] = 0.25;
        coeff_out[2] = 0.25;
        coeff_out[3] = 0.25;
          // add in contribution for any mid-edge nodes present
        for (int i = 0; i < 4; ++i) { // for each edge
          if (nodeset.mid_edge_node(i))
          {
            indices_out[num_coeff] = i+4;
            coeff_out[num_coeff] = 0.5;
            coeff_out[ i     ] -= 0.25;
            coeff_out[(i+1)%4] -= 0.25;
            ++num_coeff;
          }
        }
      }
      break;
    default:
      MSQ_SETERR(err)(MsqError::UNSUPPORTED_ELEMENT,
                  "Request for dimension %d mapping function value"
                  "for a quadrilateral element", loc.dimension);
  }
}
     

static void derivatives_at_corner( unsigned corner, 
                                   NodeSet nodeset,
                                   size_t* vertices,
                                   MsqVector<2>* derivs,
                                   size_t& num_vtx )
{
  static const unsigned xi_adj_corners[]  = { 1, 0, 3, 2 };
  static const unsigned xi_adj_edges[]    = { 0, 0, 2, 2 };
  static const unsigned eta_adj_corners[] = { 3, 2, 1, 0 };
  static const unsigned eta_adj_edges[]   = { 3, 1, 1, 3 };
  
  static const double corner_xi[]  = { -1,  1,  1, -1 }; // xi values by corner
  static const double corner_eta[] = { -1, -1,  1,  1 }; // eta values by corner
  static const double other_xi[]   = {  1, -1, -1,  1 }; // xi values for adjacent corner in xi direction
  static const double other_eta[]  = {  1,  1, -1, -1 }; // eta values for adjcent corner in eta direction
  static const double mid_xi[]     = {  2, -2, -2,  2 }; // xi values for mid-node in xi direction
  static const double mid_eta[]    = {  2,  2, -2, -2 }; // eta values for mid-node in eta direction
  
  num_vtx = 3;
  vertices[0] = corner;
  vertices[1] = xi_adj_corners[corner];
  vertices[2] = eta_adj_corners[corner];

  derivs[0][0] = corner_xi [corner];
  derivs[0][1] = corner_eta[corner];
  derivs[1][0] = other_xi  [corner];
  derivs[1][1] = 0.0;
  derivs[2][0] = 0.0;
  derivs[2][1] = other_eta [corner];

  if (nodeset.mid_edge_node(xi_adj_edges[corner])) {
    vertices[num_vtx] = 4 + xi_adj_edges[corner];
    derivs[num_vtx][0] = 2.0*mid_xi[corner];
    derivs[num_vtx][1] = 0.0;
    derivs[0][0] -= mid_xi[corner];
    derivs[1][0] -= mid_xi[corner];
    ++num_vtx;
  }

  if (nodeset.mid_edge_node(eta_adj_edges[corner])) {
    vertices[num_vtx] = 4 + eta_adj_edges[corner];
    derivs[num_vtx][0] = 0.0;
    derivs[num_vtx][1] = 2.0*mid_eta[corner];
    derivs[0][1] -= mid_eta[corner];
    derivs[2][1] -= mid_eta[corner];
    ++num_vtx;
  }
}

static void derivatives_at_mid_edge( unsigned edge, 
                                     NodeSet nodeset,
                                     size_t* vertices,
                                     MsqVector<2>* derivs,
                                     size_t& num_vtx )
{
  static const double values[][9] = { {-0.5, -0.5, 0.5,  0.5, -1.0,  2.0,  1.0,  2.0,  4.0},
                                      {-0.5,  0.5, 0.5, -0.5, -2.0,  1.0, -2.0, -1.0, -4.0},
                                      {-0.5, -0.5, 0.5,  0.5, -1.0, -2.0,  1.0, -2.0, -4.0},
                                      {-0.5,  0.5, 0.5, -0.5,  2.0,  1.0,  2.0, -1.0,  4.0} };
  static const double edge_values[][2] = { {-1,  1},
                                           {-1,  1},
                                           { 1, -1},
                                           { 1, -1} }; 
  const unsigned prev_corner = edge;           // index of start vertex of edge
  const unsigned next_corner = (edge+1)%4;     // index of end vertex of edge
  const unsigned is_eta_edge = edge % 2;       // edge is xi = +/- 0
  const unsigned is_xi_edge  = 1 - is_eta_edge;// edge is eta = +/- 0
  //const unsigned mid_edge_node = edge + 4;     // mid-edge node index
  const unsigned prev_opposite = (prev_corner+3)%4; // index of corner adjacent to prev_corner
  const unsigned next_opposite = (next_corner+1)%4; // index of corner adjacent to next_corner;
 
    // First do derivatives along edge (e.g. wrt xi if edge is eta = +/-1)
  num_vtx = 2;
  vertices[0] = prev_corner;
  vertices[1] = next_corner;
  derivs[0][is_eta_edge] = edge_values[edge][0];
  derivs[0][is_xi_edge]  = 0.0;
  derivs[1][is_eta_edge] = edge_values[edge][1];
  derivs[1][is_xi_edge]  = 0.0;
    // That's it for the edge-direction derivatives.  No other vertices contribute.
    
    // Next handle the linear element case.  Handle this as a special case first,
    // so the generalized solution doesn't impact performance for linear elements
    // too much.
  if (!nodeset.have_any_mid_node()) {
    num_vtx = 4;
    vertices[2] = prev_opposite;
    vertices[3] = next_opposite;
    derivs[0][is_xi_edge] = values[edge][prev_corner];
    derivs[1][is_xi_edge] = values[edge][next_corner];
    derivs[2][is_xi_edge] = values[edge][prev_opposite];
    derivs[2][is_eta_edge] = 0.0;
    derivs[3][is_xi_edge] = values[edge][next_opposite];
    derivs[3][is_eta_edge] = 0.0;
    return;
  }
  
    // Initial (linear) contribution for each corner
  double v[4] = { values[edge][0], 
                  values[edge][1], 
                  values[edge][2], 
                  values[edge][3] };

    // If mid-face node is present
  double v8 = 0.0;
  if (nodeset.mid_face_node(0)) {
    v8 = values[edge][8];
    vertices[num_vtx] = 8;
    derivs[num_vtx][is_eta_edge] = 0.0;
    derivs[num_vtx][is_xi_edge] = v8;
    v[0] -= 0.25 * v8;
    v[1] -= 0.25 * v8;
    v[2] -= 0.25 * v8;
    v[3] -= 0.25 * v8;
    ++num_vtx;
  }

    // If mid-edge nodes are present
  for (unsigned i = 0; i < 4; ++i) {
    if (nodeset.mid_edge_node(i)) {
      const double value = values[edge][i+4] - 0.5 * v8;
      if (fabs(value) > 0.125) {
        v[ i     ] -= 0.5 * value;
        v[(i+1)%4] -= 0.5 * value;
        vertices[num_vtx] = i+4;
        derivs[num_vtx][is_eta_edge] = 0.0;
        derivs[num_vtx][is_xi_edge] = value;
        ++num_vtx;
      }
    }
  }

    // update values for adjacent corners
  derivs[0][is_xi_edge] = v[prev_corner];
  derivs[1][is_xi_edge] = v[next_corner];
    // do other two corners
  if (fabs(v[prev_opposite]) > 0.125) {
    vertices[num_vtx] = prev_opposite;
    derivs[num_vtx][is_eta_edge] = 0.0;
    derivs[num_vtx][is_xi_edge] = v[prev_opposite];
    ++num_vtx;
  }
  if (fabs(v[next_opposite]) > 0.125) {
    vertices[num_vtx] = next_opposite;
    derivs[num_vtx][is_eta_edge] = 0.0;
    derivs[num_vtx][is_xi_edge] = v[next_opposite];
    ++num_vtx;
  }
}


static void derivatives_at_mid_elem( NodeSet nodeset,
                                     size_t* vertices,
                                     MsqVector<2>* derivs,
                                     size_t& num_vtx )
{
    // fast linear case
    // This is provided as an optimization for linear elements.
    // If this block of code were removed, the general-case code
    // below should produce the same result.
  if (!nodeset.have_any_mid_node()) {
    num_vtx = 4;
    vertices[0] = 0; derivs[0][0] = -0.5; derivs[0][1] = -0.5;
    vertices[1] = 1; derivs[1][0] =  0.5; derivs[1][1] = -0.5;
    vertices[2] = 2; derivs[2][0] =  0.5; derivs[2][1] =  0.5;
    vertices[3] = 3; derivs[3][0] = -0.5; derivs[3][1] =  0.5;
    return;
  }
  
  num_vtx = 0;
  
    // N_0
  if (!nodeset.both_edge_nodes(0,3)) {  // if eiter adjacent mid-edge node is missing
    vertices[num_vtx] = 0;
    derivs[num_vtx][0] = nodeset.mid_edge_node(3) ? 0.0 : -0.5;
    derivs[num_vtx][1] = nodeset.mid_edge_node(0) ? 0.0 : -0.5;
    ++num_vtx;
  }
  
    // N_1
  if (!nodeset.both_edge_nodes(0,1)) {  // if eiter adjacent mid-edge node is missing
    vertices[num_vtx] = 1;
    derivs[num_vtx][0] = nodeset.mid_edge_node(1) ? 0.0 :  0.5;
    derivs[num_vtx][1] = nodeset.mid_edge_node(0) ? 0.0 : -0.5;
    ++num_vtx;
  }
  
    // N_2
  if (!nodeset.both_edge_nodes(1,2)) {  // if eiter adjacent mid-edge node is missing
    vertices[num_vtx] = 2;
    derivs[num_vtx][0] = nodeset.mid_edge_node(1) ? 0.0 :  0.5;
    derivs[num_vtx][1] = nodeset.mid_edge_node(2) ? 0.0 :  0.5;
    ++num_vtx;
  }
  
    // N_3
  if (!nodeset.both_edge_nodes(2,3)) {  // if eiter adjacent mid-edge node is missing
    vertices[num_vtx] = 3;
    derivs[num_vtx][0] = nodeset.mid_edge_node(3) ? 0.0 : -0.5;
    derivs[num_vtx][1] = nodeset.mid_edge_node(2) ? 0.0 :  0.5;
    ++num_vtx;
  }
  
    // N_4
  if (nodeset.mid_edge_node(0)) {
    vertices[num_vtx] = 4;
    derivs[num_vtx][0] =  0.0;
    derivs[num_vtx][1] = -1.0;
    ++num_vtx;
  }
  
    // N_5
  if (nodeset.mid_edge_node(1)) {
    vertices[num_vtx] = 5;
    derivs[num_vtx][0] =  1.0;
    derivs[num_vtx][1] =  0.0;
    ++num_vtx;
  }
  
    // N_6
  if (nodeset.mid_edge_node(2)) {
    vertices[num_vtx] = 6;
    derivs[num_vtx][0] =  0.0;
    derivs[num_vtx][1] =  1.0;
    ++num_vtx;
  }
  
    // N_7
  if (nodeset.mid_edge_node(3)) {
    vertices[num_vtx] = 7;
    derivs[num_vtx][0] = -1.0;
    derivs[num_vtx][1] =  0.0;
    ++num_vtx;
  }
  
    // N_8 (mid-quad node) never contributes to Jacobian at element center!!!
}

void QuadLagrangeShape::derivatives( Sample loc,
                                     NodeSet nodeset,
                                     size_t* vertex_indices_out,
                                     MsqVector<2>* d_coeff_d_xi_out,
                                     size_t& num_vtx,
                                     MsqError& err ) const
{
  if (!nodeset.have_any_mid_node()) {
    LinearQuadrilateral::derivatives( loc, nodeset, vertex_indices_out, d_coeff_d_xi_out, num_vtx );
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


void QuadLagrangeShape::ideal( Sample , 
                               MsqMatrix<3,2>& J,
                               MsqError&  ) const
{
  J(0,0) = J(1,1) = 1.0;
  J(0,1) = J(1,0) = 0.0;
  J(2,0) = J(2,1) = 0.0;
}

} // namespace Mesquite
