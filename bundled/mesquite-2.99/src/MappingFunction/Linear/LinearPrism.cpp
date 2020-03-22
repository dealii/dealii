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

    (2006) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file LinearPrism.cpp
 *  \brief mapping function for linear prism
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "MsqError.hpp"
#include "LinearPrism.hpp"

namespace MESQUITE_NS {

static const char* nonlinear_error 
 = "Attempt to use LinearPrism mapping function for a nonlinear element\n";


EntityTopology LinearPrism::element_topology() const
  { return PRISM; }
  
int LinearPrism::num_nodes() const
  { return 6; }
  
static const int edge_beg[] = { 0, 1, 2, 0, 1, 2, 3, 4, 5 };
static const int edge_end[] = { 1, 2, 0, 3, 4, 5, 4, 5, 3 };
static const int faces[5][5] = { { 4, 0, 1, 4, 3 },
                                 { 4, 1, 2, 5, 4 },
                                 { 4, 2, 0, 3, 5 },
                                 { 3, 0, 1, 2,-1 },
                                 { 3, 3, 4, 5,-1 } };

static void coefficients_at_corner( unsigned corner,
                                    double* coeff_out,
                                    size_t* indices_out,
                                    size_t& num_coeff )
{
  num_coeff = 1;
  indices_out[0] = corner;
  coeff_out[0] = 1.0;
}

static void coefficients_at_mid_edge( unsigned edge,
                                      double* coeff_out,
                                      size_t* indices_out,
                                      size_t& num_coeff )
{
  num_coeff = 2;
  indices_out[0] = edge_beg[edge];
  indices_out[1] = edge_end[edge];
  coeff_out[0] = 0.5;
  coeff_out[1] = 0.5;
}

static void coefficients_at_mid_face( unsigned face,
                                      double* coeff_out,
                                      size_t* indices_out,
                                      size_t& num_coeff )
{
  double f;
  if (faces[face][0] == 4) {
    num_coeff = 4;
    f = 0.25;
    indices_out[3] = faces[face][4];
    coeff_out[3] = f;
  }
  else {
    num_coeff = 3;
    f = MSQ_ONE_THIRD; 
  }
  
  coeff_out[0] = f;
  coeff_out[1] = f;
  coeff_out[2] = f;
  indices_out[0] = faces[face][1];
  indices_out[1] = faces[face][2];
  indices_out[2] = faces[face][3];
}

static void coefficients_at_mid_elem( double* coeff_out,
                                      size_t* indices_out,
                                      size_t& num_coeff )
{
  num_coeff = 6;
  const double sixth = 1.0/6.0;
  coeff_out[0] = sixth;
  coeff_out[1] = sixth;
  coeff_out[2] = sixth;
  coeff_out[3] = sixth;
  coeff_out[4] = sixth;
  coeff_out[5] = sixth;
  indices_out[0] = 0;
  indices_out[1] = 1;
  indices_out[2] = 2;
  indices_out[3] = 3;
  indices_out[4] = 4;
  indices_out[5] = 5;
}

void LinearPrism::coefficients( Sample loc,
                                NodeSet nodeset,
                                double* coeff_out,
                                size_t* indices_out,
                                size_t& num_coeff,
                                MsqError& err ) const
{
  if (nodeset.have_any_mid_node()) {
    MSQ_SETERR(err)(nonlinear_error, MsqError::UNSUPPORTED_ELEMENT );
    return;
  }
  
  switch (loc.dimension) {
    case 0:
      coefficients_at_corner( loc.number, coeff_out, indices_out, num_coeff );
      break;
    case 1:
      coefficients_at_mid_edge( loc.number, coeff_out, indices_out, num_coeff );
      break;
    case 2:
      coefficients_at_mid_face( loc.number, coeff_out, indices_out, num_coeff );
      break;
    case 3:
      coefficients_at_mid_elem( coeff_out, indices_out, num_coeff );
      break;
    default:
      MSQ_SETERR(err)("Invalid/unsupported logical dimension",MsqError::INVALID_ARG);
  }
}

                                

static void derivatives_at_corner( unsigned corner, 
                                   size_t* vertex_indices_out,
                                   MsqVector<3>* d_coeff_d_xi_out,
                                   size_t& num_vtx )
{
  int tri = (corner / 3); // 0 for xi=0, 1 for xi=1
  int tv = corner % 3;    // index of corner with xi=constant triangle

  num_vtx = 4;
    // three vertices within the xi=constant triangle
  vertex_indices_out[0] = 3*tri;
  vertex_indices_out[1] = 3*tri+1;
  vertex_indices_out[2] = 3*tri+2;
    // vertex adjacent to corner in other triangle
  vertex_indices_out[3] = 3 - 6*tri + corner;
  
    // three vertices within the xi=constant triangle
  d_coeff_d_xi_out[0][0] =  0.0;
  d_coeff_d_xi_out[0][1] = -1.0;
  d_coeff_d_xi_out[0][2] = -1.0;
  d_coeff_d_xi_out[1][0] =  0.0;
  d_coeff_d_xi_out[1][1] =  1.0;
  d_coeff_d_xi_out[1][2] =  0.0;
  d_coeff_d_xi_out[2][0] =  0.0;
  d_coeff_d_xi_out[2][1] =  0.0;
  d_coeff_d_xi_out[2][2] =  1.0;
    // fix dxi value for input corner
  d_coeff_d_xi_out[tv][0] = 2*tri - 1; 
    // vertex adjacent to corner in other triangle
  d_coeff_d_xi_out[3][0] =  1 - 2*tri;
  d_coeff_d_xi_out[3][1] =  0.0;
  d_coeff_d_xi_out[3][2] =  0.0;
}


static void derivatives_at_mid_edge( unsigned edge, 
                                     size_t* vertex_indices_out,
                                     MsqVector<3>* d_coeff_d_xi_out,
                                     size_t& num_vtx )
{
  int opp; // vertex opposite edge in same triagle
  
  switch (edge/3) {
    case 0:  // triangle at xi = 0
      opp = (edge+2)%3;

      num_vtx = 5;
        // vertices in this xi = 0 triagnle
      vertex_indices_out[0] = 0;
      vertex_indices_out[1] = 1;
      vertex_indices_out[2] = 2;
        // adjacent vertices in xi = 1 triangle
      vertex_indices_out[3] = 3 + edge;
      vertex_indices_out[4] = 3 + (edge+1)%3;

        // vertices in this xi = 0 triagnle
      d_coeff_d_xi_out[0][0] = -0.5;
      d_coeff_d_xi_out[0][1] = -1.0;
      d_coeff_d_xi_out[0][2] = -1.0;
      d_coeff_d_xi_out[1][0] = -0.5;
      d_coeff_d_xi_out[1][1] =  1.0;
      d_coeff_d_xi_out[1][2] =  0.0;
      d_coeff_d_xi_out[2][0] = -0.5;
      d_coeff_d_xi_out[2][1] =  0.0;
      d_coeff_d_xi_out[2][2] =  1.0;
        // clear dxi for vertex opposite edge in xi = 0 triangle
      d_coeff_d_xi_out[opp][0] = 0.0;
        // adjacent vertices in xi = 1 triangle
      d_coeff_d_xi_out[3][0] =  0.5;
      d_coeff_d_xi_out[3][1] =  0.0;
      d_coeff_d_xi_out[3][2] =  0.0;
      d_coeff_d_xi_out[4][0] =  0.5;
      d_coeff_d_xi_out[4][1] =  0.0;
      d_coeff_d_xi_out[4][2] =  0.0;
      break;

    case 1:  // lateral edges (not in either triangle)
      num_vtx = 6;
      vertex_indices_out[0] = 0;
      vertex_indices_out[1] = 1;
      vertex_indices_out[2] = 2;
      vertex_indices_out[3] = 3;
      vertex_indices_out[4] = 4;
      vertex_indices_out[5] = 5;
      
        // set all deta & dzeta values, zero all dxi values
      d_coeff_d_xi_out[0][0] =  0.0;
      d_coeff_d_xi_out[0][1] = -0.5;
      d_coeff_d_xi_out[0][2] = -0.5;
      d_coeff_d_xi_out[1][0] =  0.0;
      d_coeff_d_xi_out[1][1] =  0.5;
      d_coeff_d_xi_out[1][2] =  0.0;
      d_coeff_d_xi_out[2][0] =  0.0;
      d_coeff_d_xi_out[2][1] =  0.0;
      d_coeff_d_xi_out[2][2] =  0.5;
      d_coeff_d_xi_out[3][0] =  0.0;
      d_coeff_d_xi_out[3][1] = -0.5;
      d_coeff_d_xi_out[3][2] = -0.5;
      d_coeff_d_xi_out[4][0] =  0.0;
      d_coeff_d_xi_out[4][1] =  0.5;
      d_coeff_d_xi_out[4][2] =  0.0;
      d_coeff_d_xi_out[5][0] =  0.0;
      d_coeff_d_xi_out[5][1] =  0.0;
      d_coeff_d_xi_out[5][2] =  0.5;
      
        // set dxi values for end points of edge
      d_coeff_d_xi_out[(edge-3)][0] = -1;
      d_coeff_d_xi_out[ edge   ][0] =  1;
      break;
    
    case 2:  // triangle at xi = 1
      opp = (edge+2)%3;

      num_vtx = 5;
        // vertices in this xi = 1 triagnle
      vertex_indices_out[0] = 3;
      vertex_indices_out[1] = 4;
      vertex_indices_out[2] = 5;
        // adjacent vertices in xi = 1 triangle
      vertex_indices_out[3] = edge - 6;
      vertex_indices_out[4] = (edge-5)%3;

        // vertices in this xi = 1 triagnle
      d_coeff_d_xi_out[0][0] =  0.5;
      d_coeff_d_xi_out[0][1] = -1.0;
      d_coeff_d_xi_out[0][2] = -1.0;
      d_coeff_d_xi_out[1][0] =  0.5;
      d_coeff_d_xi_out[1][1] =  1.0;
      d_coeff_d_xi_out[1][2] =  0.0;
      d_coeff_d_xi_out[2][0] =  0.5;
      d_coeff_d_xi_out[2][1] =  0.0;
      d_coeff_d_xi_out[2][2] =  1.0;
        // clear dxi for vertex opposite edge in xi = 1 triangle
      d_coeff_d_xi_out[opp][0] = 0.0;
        // adjacent vertices in xi = 0 triangle
      d_coeff_d_xi_out[3][0] = -0.5;
      d_coeff_d_xi_out[3][1] =  0.0;
      d_coeff_d_xi_out[3][2] =  0.0;
      d_coeff_d_xi_out[4][0] = -0.5;
      d_coeff_d_xi_out[4][1] =  0.0;
      d_coeff_d_xi_out[4][2] =  0.0;
      break;
  }
}
static void derivatives_at_mid_face( unsigned face, 
                                     size_t* vertex_indices_out,
                                     MsqVector<3>* d_coeff_d_xi_out,
                                     size_t& num_vtx )
{
  num_vtx = 6;
  vertex_indices_out[0] = 0;
  vertex_indices_out[1] = 1;
  vertex_indices_out[2] = 2;
  vertex_indices_out[3] = 3;
  vertex_indices_out[4] = 4;
  vertex_indices_out[5] = 5;
  
  int opp; // start vtx of edge opposite from quad face
  int tri_offset; // offset in d_coeff_d_xi_out for triangle containing edge
  
  if (face < 3) { // quad face
      // set all values
    d_coeff_d_xi_out[0][0] = -0.5;
    d_coeff_d_xi_out[0][1] = -0.5;
    d_coeff_d_xi_out[0][2] = -0.5;
    d_coeff_d_xi_out[1][0] = -0.5;
    d_coeff_d_xi_out[1][1] =  0.5;
    d_coeff_d_xi_out[1][2] =  0.0;
    d_coeff_d_xi_out[2][0] = -0.5;
    d_coeff_d_xi_out[2][1] =  0.0;
    d_coeff_d_xi_out[2][2] =  0.5;
    d_coeff_d_xi_out[3][0] =  0.5;
    d_coeff_d_xi_out[3][1] = -0.5;
    d_coeff_d_xi_out[3][2] = -0.5;
    d_coeff_d_xi_out[4][0] =  0.5;
    d_coeff_d_xi_out[4][1] =  0.5;
    d_coeff_d_xi_out[4][2] =  0.0;
    d_coeff_d_xi_out[5][0] =  0.5;
    d_coeff_d_xi_out[5][1] =  0.0;
    d_coeff_d_xi_out[5][2] =  0.5;
      // clear dxi for ends of edge opposite from face
    opp = (face+2)%3;
    d_coeff_d_xi_out[opp][0] = 0.0;
    d_coeff_d_xi_out[(opp+3)][0] = 0.0;
  }
  else { // triangular faces
      // set all xi values, zero all other values
    const double third = 1./3;
    d_coeff_d_xi_out[0][0] = -third;
    d_coeff_d_xi_out[0][1] =  0;
    d_coeff_d_xi_out[0][2] =  0;
    d_coeff_d_xi_out[1][0] = -third;
    d_coeff_d_xi_out[1][1] =  0;
    d_coeff_d_xi_out[1][2] =  0;
    d_coeff_d_xi_out[2][0] = -third;
    d_coeff_d_xi_out[2][1] =  0;
    d_coeff_d_xi_out[2][2] =  0;
    d_coeff_d_xi_out[3][0] =  third;
    d_coeff_d_xi_out[3][1] =  0;
    d_coeff_d_xi_out[3][2] =  0;
    d_coeff_d_xi_out[4][0] =  third;
    d_coeff_d_xi_out[4][1] =  0;
    d_coeff_d_xi_out[4][2] =  0;
    d_coeff_d_xi_out[5][0] =  third;
    d_coeff_d_xi_out[5][1] =  0;
    d_coeff_d_xi_out[5][2] =  0;
      // set deta and dzeta values for vertices in same triangle as edge
    tri_offset = 3 * (face - 3);  // either 0 or 3
    d_coeff_d_xi_out[tri_offset][1] = -1.0;
    d_coeff_d_xi_out[tri_offset][2] = -1.0;
    d_coeff_d_xi_out[tri_offset+1][1] =  1.0;
    d_coeff_d_xi_out[tri_offset+2][2] =  1.0;
  }
}
static void derivatives_at_mid_elem( size_t* vertex_indices_out,
                                     MsqVector<3>* d_coeff_d_xi_out,
                                     size_t& num_vtx )
{
  const double third = 1./3;
  
  num_vtx = 6;;
  vertex_indices_out[0] = 0;
  vertex_indices_out[1] = 1;
  vertex_indices_out[2] = 2;
  vertex_indices_out[3] = 3;
  vertex_indices_out[4] = 4;
  vertex_indices_out[5] = 5;
  
  d_coeff_d_xi_out[0][0] = -third;
  d_coeff_d_xi_out[0][1] = -0.5;
  d_coeff_d_xi_out[0][2] = -0.5;
  d_coeff_d_xi_out[1][0] = -third;
  d_coeff_d_xi_out[1][1] =  0.5;
  d_coeff_d_xi_out[1][2] =  0.0;
  d_coeff_d_xi_out[2][0] = -third;
  d_coeff_d_xi_out[2][1] =  0.0;
  d_coeff_d_xi_out[2][2] =  0.5;
  d_coeff_d_xi_out[3][0] =  third;
  d_coeff_d_xi_out[3][1] = -0.5;
  d_coeff_d_xi_out[3][2] = -0.5;
  d_coeff_d_xi_out[4][0] =  third;
  d_coeff_d_xi_out[4][1] =  0.5;
  d_coeff_d_xi_out[4][2] =  0.0;
  d_coeff_d_xi_out[5][0] =  third;
  d_coeff_d_xi_out[5][1] =  0.0;
  d_coeff_d_xi_out[5][2] =  0.5;
}
  
void LinearPrism::derivatives( Sample loc,
                               NodeSet nodeset,
                               size_t* vertex_indices_out,
                               MsqVector<3>* d_coeff_d_xi_out,
                               size_t& num_vtx,
                               MsqError& err ) const
{
  if (nodeset.have_any_mid_node()) {
    MSQ_SETERR(err)(nonlinear_error, MsqError::UNSUPPORTED_ELEMENT );
    return;
  }
  
  switch (loc.dimension) {
    case 0:
      derivatives_at_corner( loc.number, vertex_indices_out, d_coeff_d_xi_out, num_vtx );
      break;
    case 1:
      derivatives_at_mid_edge( loc.number, vertex_indices_out, d_coeff_d_xi_out, num_vtx );
      break;
    case 2:
      derivatives_at_mid_face( loc.number, vertex_indices_out, d_coeff_d_xi_out, num_vtx );
      break;
    case 3:
      derivatives_at_mid_elem( vertex_indices_out, d_coeff_d_xi_out, num_vtx );
      break;
    default:
      MSQ_SETERR(err)("Invalid/unsupported logical dimension",MsqError::INVALID_ARG);
  }
}


void LinearPrism::ideal( Sample , 
                         MsqMatrix<3,3>& J,
                         MsqError&  ) const
{
  const double a = 0.52455753171082409;  // 2^(-2/3) * 3^(-1/6)
  const double b = 0.90856029641606983;  // a * sqrt(3) = 1/2 cbrt(6)

  J(0,0) = 2*a; J(0,1) = 0.0;  J(0,2) = 0.0;
  J(1,0) = 0.0; J(1,1) = 2*a;  J(1,2) = a;
  J(2,0) = 0.0; J(2,1) = 0.0;  J(2,2) = b;
}

} // namespace Mesquite
