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


/** \file TetLagrangeShape.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "TetLagrangeShape.hpp"
#include "MsqError.hpp"
#include <assert.h>

namespace MESQUITE_NS {

EntityTopology TetLagrangeShape::element_topology() const
  { return TETRAHEDRON; }
  
int TetLagrangeShape::num_nodes() const
  { return 10; }

NodeSet TetLagrangeShape::sample_points( NodeSet ns ) const
{
  if (ns.have_any_mid_node()) {
    ns.set_all_corner_nodes(TETRAHEDRON);
  }
  else {
    ns.clear();
    ns.set_mid_region_node();
  }
  return ns;
}

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
                                      NodeSet nodeset,
                                      double* coeff_out,
                                      size_t* indices_out,
                                      size_t& num_coeff )
{
  if (nodeset.mid_edge_node(edge)) { // if mid-edge node is present
    num_coeff = 1;
    indices_out[0] = 4+edge;
    coeff_out[0] = 1.0;
  }
  else { // no mid node on edge
    num_coeff = 2;
    coeff_out[0] = coeff_out[1] = 0.5;
    if (edge < 3) {
      indices_out[0] = edge;
      indices_out[1] = (edge+1) % 3;
    }
    else {
      indices_out[0] = edge-3;
      indices_out[1] = 3;
    }
  }
}

static void coefficients_at_mid_face( unsigned face, 
                                      NodeSet nodeset,
                                      double* coeff_out,
                                      size_t* indices_out,
                                      size_t& num_coeff )
{
  const double one_ninth = 1.0/9.0;
  const double two_ninth = 2.0/9.0;
  const double four_ninth = 4.0/9.0;
  
  if (face < 3) {
    const int next = (face+1)%3;
    indices_out[0] = face;
    indices_out[1] = next;
    indices_out[2] = 3;
    coeff_out[0] = -one_ninth;
    coeff_out[1] = -one_ninth;
    coeff_out[2] = -one_ninth;
    num_coeff = 3;
    if (nodeset.mid_edge_node(face)) {
      indices_out[num_coeff] = 4+face;
      coeff_out[num_coeff] = four_ninth;
      ++num_coeff;
    }
    else {
      coeff_out[0] += two_ninth;
      coeff_out[1] += two_ninth;
    }
    if (nodeset.mid_edge_node(3+next)) {
      indices_out[num_coeff] = 7+next;
      coeff_out[num_coeff] = four_ninth;
      ++num_coeff;
    }
    else {
      coeff_out[1] += two_ninth;
      coeff_out[2] += two_ninth;
    }
    if (nodeset.mid_edge_node(3+face)) {
      indices_out[num_coeff] = 7+face;
      coeff_out[num_coeff] = four_ninth;
      ++num_coeff;
    }
    else {
      coeff_out[0] += two_ninth;
      coeff_out[2] += two_ninth;
    }
  }
  else {
    assert( face == 3);
    indices_out[0] = 0;
    indices_out[1] = 1;
    indices_out[2] = 2;
    coeff_out[0] = -one_ninth;
    coeff_out[1] = -one_ninth;
    coeff_out[2] = -one_ninth;
    num_coeff = 3;
    if (nodeset.mid_edge_node(0)) {
      indices_out[num_coeff] = 4;
      coeff_out[num_coeff] = four_ninth;
      ++num_coeff;
    }
    else {
      coeff_out[0] += two_ninth;
      coeff_out[1] += two_ninth;
    }
    if (nodeset.mid_edge_node(1)) {
      indices_out[num_coeff] = 5;
      coeff_out[num_coeff] = four_ninth;
      ++num_coeff;
    }
    else {
      coeff_out[1] += two_ninth;
      coeff_out[2] += two_ninth;
    }
    if (nodeset.mid_edge_node(2)) {
      indices_out[num_coeff] = 6;
      coeff_out[num_coeff] = four_ninth;
      ++num_coeff;
    }
    else {
      coeff_out[2] += two_ninth;
      coeff_out[0] += two_ninth;
    }
  }
}

static void coefficients_at_mid_elem( NodeSet nodeset,
                                      double* coeff_out,
                                      size_t* indices_out,
                                      size_t& num_coeff )
{
  num_coeff = 4;
  indices_out[0] = 0;
  indices_out[1] = 1;
  indices_out[2] = 2;
  indices_out[3] = 3;
  coeff_out[0] = -0.125;
  coeff_out[1] = -0.125;
  coeff_out[2] = -0.125;
  coeff_out[3] = -0.125;
  if (nodeset.mid_edge_node(0)) {
    indices_out[num_coeff] = 4;
    coeff_out[num_coeff] = 0.25;
    ++num_coeff;
  }
  else {
    coeff_out[0] += 0.125;
    coeff_out[1] += 0.125;
  }
  if (nodeset.mid_edge_node(1)) {
    indices_out[num_coeff] = 5;
    coeff_out[num_coeff] = 0.25;
    ++num_coeff;
  }
  else {
    coeff_out[1] += 0.125;
    coeff_out[2] += 0.125;
  }
  if (nodeset.mid_edge_node(2)) {
    indices_out[num_coeff] = 6;
    coeff_out[num_coeff] = 0.25;
    ++num_coeff;
  }
  else {
    coeff_out[2] += 0.125;
    coeff_out[0] += 0.125;
  }
  if (nodeset.mid_edge_node(3)) {
    indices_out[num_coeff] = 7;
    coeff_out[num_coeff] = 0.25;
    ++num_coeff;
  }
  else {
    coeff_out[0] += 0.125;
    coeff_out[3] += 0.125;
  }
  if (nodeset.mid_edge_node(4)) {
    indices_out[num_coeff] = 8;
    coeff_out[num_coeff] = 0.25;
    ++num_coeff;
  }
  else {
    coeff_out[1] += 0.125;
    coeff_out[3] += 0.125;
  }
  if (nodeset.mid_edge_node(5)) {
    indices_out[num_coeff] = 9;
    coeff_out[num_coeff] = 0.25;
    ++num_coeff;
  }
  else {
    coeff_out[2] += 0.125;
    coeff_out[3] += 0.125;
  }
}

void TetLagrangeShape::coefficients( Sample loc,
                                     NodeSet nodeset,
                                     double* coeff_out,
                                     size_t* indices_out,
                                     size_t& num_coeff,
                                     MsqError& err ) const
{
  if (nodeset.have_any_mid_face_node() | nodeset.have_any_mid_region_node()) {
    MSQ_SETERR(err)("TetLagrangeShape does not support mid-face/mid-element nodes",
                    MsqError::UNSUPPORTED_ELEMENT);
    return;
  }
  
  switch (loc.dimension) {
    case 0:
      coefficients_at_corner( loc.number, coeff_out, indices_out, num_coeff );
      break;
    case 1:
      coefficients_at_mid_edge( loc.number, nodeset, coeff_out, indices_out, num_coeff );
      break;
    case 2:
      coefficients_at_mid_face( loc.number, nodeset, coeff_out, indices_out, num_coeff );
      break;
    case 3:
      coefficients_at_mid_elem( nodeset, coeff_out, indices_out, num_coeff );
      break;
    default:
      MSQ_SETERR(err)("Invalid/unsupported logical dimension",MsqError::INVALID_ARG);
  }
}

static void get_linear_derivatives( size_t* vertices,
                                    MsqVector<3>* derivs )
{
  vertices[0] = 0;
  derivs[0][0] = -1.0;
  derivs[0][1] = -1.0;
  derivs[0][2] = -1.0;
  
  vertices[1] = 1;
  derivs[1][0] = 1.0;
  derivs[1][1] = 0.0;
  derivs[1][2] = 0.0;
  
  vertices[2] = 2;
  derivs[2][0] = 0.0;
  derivs[2][1] = 1.0;
  derivs[2][2] = 0.0;
  
  vertices[3] = 3;
  derivs[3][0] = 0.0;
  derivs[3][1] = 0.0;
  derivs[3][2] = 1.0;
}

static const unsigned edges[][2] = { { 0, 1 },
                                     { 1, 2 },
                                     { 2, 0 },
                                     { 0, 3 },
                                     { 1, 3 },
                                     { 2, 3 } };

static void derivatives_at_corner( unsigned corner,
                                   NodeSet nodeset,
                                   size_t* vertices,
                                   MsqVector<3>* derivs,
                                   size_t& num_vtx )
{
    // begin with derivatives for linear tetrahedron
  num_vtx = 4;
  get_linear_derivatives( vertices, derivs );
  
    // adjust for the presence of mid-edge nodes
  switch (corner) {
    case 0:
      if (nodeset.mid_edge_node(0)) {
        vertices[num_vtx] = 4;
        derivs[num_vtx][0] = 4.0;
        derivs[num_vtx][1] = 0.0;
        derivs[num_vtx][2] = 0.0;
        derivs[0][0] -= 2.0;
        derivs[1][0] -= 2.0;
        ++num_vtx;
      }
      if (nodeset.mid_edge_node(2)) {
        vertices[num_vtx] = 6;
        derivs[num_vtx][0] = 0.0;
        derivs[num_vtx][1] = 4.0;
        derivs[num_vtx][2] = 0.0;
        derivs[0][1] -= 2.0;
        derivs[2][1] -= 2.0;
        ++num_vtx;
      }
      if (nodeset.mid_edge_node(3)) {
        vertices[num_vtx] = 7;
        derivs[num_vtx][0] = 0.0;
        derivs[num_vtx][1] = 0.0;
        derivs[num_vtx][2] = 4.0;
        derivs[0][2] -= 2.0;
        derivs[3][2] -= 2.0;
        ++num_vtx;
      }
      break;

    case 1:
      if (nodeset.mid_edge_node(0)) {
        vertices[num_vtx] = 4;
        derivs[num_vtx][0] = -4.0;
        derivs[num_vtx][1] = -4.0;
        derivs[num_vtx][2] = -4.0;
        derivs[0][0] += 2.0;
        derivs[0][1] += 2.0;
        derivs[0][2] += 2.0;
        derivs[1][0] += 2.0;
        derivs[1][1] += 2.0;
        derivs[1][2] += 2.0;
        ++num_vtx;
      }
      if (nodeset.mid_edge_node(1)) {
        vertices[num_vtx] = 5;
        derivs[num_vtx][0] = 0.0;
        derivs[num_vtx][1] = 4.0;
        derivs[num_vtx][2] = 0.0;
        derivs[1][1] -= 2.0;
        derivs[2][1] -= 2.0;
        ++num_vtx;
      }
      if (nodeset.mid_edge_node(4)) {
        vertices[num_vtx] = 8;
        derivs[num_vtx][0] = 0.0;
        derivs[num_vtx][1] = 0.0;
        derivs[num_vtx][2] = 4.0;
        derivs[1][2] -= 2.0;
        derivs[3][2] -= 2.0;
        ++num_vtx;
      }
      break;
  
    case 2:
      if (nodeset.mid_edge_node(1)) {
        vertices[num_vtx] = 5;
        derivs[num_vtx][0] = 4.0;
        derivs[num_vtx][1] = 0.0;
        derivs[num_vtx][2] = 0.0;
        derivs[1][0] -= 2.0;
        derivs[2][0] -= 2.0;
        ++num_vtx;
      }
      if (nodeset.mid_edge_node(2)) {
        vertices[num_vtx] = 6;
        derivs[num_vtx][0] = -4.0;
        derivs[num_vtx][1] = -4.0;
        derivs[num_vtx][2] = -4.0;
        derivs[0][0] += 2.0;
        derivs[0][1] += 2.0;
        derivs[0][2] += 2.0;
        derivs[2][0] += 2.0;
        derivs[2][1] += 2.0;
        derivs[2][2] += 2.0;
        ++num_vtx;
      }
      if (nodeset.mid_edge_node(5)) {
        vertices[num_vtx] = 9;
        derivs[num_vtx][0] = 0.0;
        derivs[num_vtx][1] = 0.0;
        derivs[num_vtx][2] = 4.0;
        derivs[2][2] -= 2.0;
        derivs[3][2] -= 2.0;
        ++num_vtx;
      }
      break;
  
    case 3:
      if (nodeset.mid_edge_node(3)) {
        vertices[num_vtx] = 7;
        derivs[num_vtx][0] = -4.0;
        derivs[num_vtx][1] = -4.0;
        derivs[num_vtx][2] = -4.0;
        derivs[0][0] += 2.0;
        derivs[0][1] += 2.0;
        derivs[0][2] += 2.0;
        derivs[3][0] += 2.0;
        derivs[3][1] += 2.0;
        derivs[3][2] += 2.0;
        ++num_vtx;
      }
      if (nodeset.mid_edge_node(4)) {
        vertices[num_vtx] = 8;
        derivs[num_vtx][0] = 4.0;
        derivs[num_vtx][1] = 0.0;
        derivs[num_vtx][2] = 0.0;
        derivs[1][0] -= 2.0;
        derivs[3][0] -= 2.0;
        ++num_vtx;
      }
      
      if (nodeset.mid_edge_node(5)) {
        vertices[num_vtx] = 9;
        derivs[num_vtx][0] = 0.0;
        derivs[num_vtx][1] = 4.0;
        derivs[num_vtx][2] = 0.0;
        derivs[2][1]-= 2.0;
        derivs[3][1]-= 2.0;
        ++num_vtx;
      }
      break;
  }
}
  
static void derivatives_at_mid_edge( unsigned edge,
                                     NodeSet nodeset,
                                     size_t* vertices,
                                     MsqVector<3>* derivs,
                                     size_t& num_vtx )
{
  int sign;
  num_vtx = 2;
  switch (edge) {
    case 0:
      vertices[0] = 0;
      derivs[0][0] = -1.0;
      derivs[0][1] = -1.0;
      derivs[0][2] = -1.0;
      
      vertices[1] = 1;
      derivs[1][0] = 1.0;
      derivs[1][1] = 0.0;
      derivs[1][2] = 0.0;
      
      if (nodeset.mid_edge_node(1) == nodeset.mid_edge_node(2)) {
        vertices[num_vtx] = 2;
        sign = 1 - 2*nodeset.mid_edge_node(1);
        derivs[num_vtx][0] =  0.0;
        derivs[num_vtx][1] = sign;
        derivs[num_vtx][2] =  0.0;
        ++num_vtx;
      }
      
      if (nodeset.mid_edge_node(3) == nodeset.mid_edge_node(4)) {
        vertices[num_vtx] = 3;
        sign = 1 - 2*nodeset.mid_edge_node(3);
        derivs[num_vtx][0] =  0.0;
        derivs[num_vtx][1] =  0.0;
        derivs[num_vtx][2] = sign;
        ++num_vtx;
      }
      
      if (nodeset.mid_edge_node(0)) {
        vertices[num_vtx] = 4;
        derivs[num_vtx][0] =  0.0;
        derivs[num_vtx][1] = -2.0;
        derivs[num_vtx][2] = -2.0;
        derivs[0][1] += 1.0;
        derivs[0][2] += 1.0;
        derivs[1][1] += 1.0;
        derivs[1][2] += 1.0;
        ++num_vtx;
      }
      
      if (nodeset.mid_edge_node(1)) {
        vertices[num_vtx] = 5;
        derivs[num_vtx][0] = 0.0;
        derivs[num_vtx][1] = 2.0;
        derivs[num_vtx][2] = 0.0;
        derivs[1][1] -= 1.0;
        ++num_vtx;
      }
      
      if (nodeset.mid_edge_node(2)) {
        vertices[num_vtx] = 6;
        derivs[num_vtx][0] = 0.0;
        derivs[num_vtx][1] = 2.0;
        derivs[num_vtx][2] = 0.0;
        derivs[0][1] -= 1.0;
        ++num_vtx;
      }
      
      if (nodeset.mid_edge_node(3)) {
        vertices[num_vtx] = 7;
        derivs[num_vtx][0] = 0.0;
        derivs[num_vtx][1] = 0.0;
        derivs[num_vtx][2] = 2.0;
        derivs[0][2] -= 1.0;
        ++num_vtx;
      }
      
      if (nodeset.mid_edge_node(4)) {
        vertices[num_vtx] = 8;
        derivs[num_vtx][0] = 0.0;
        derivs[num_vtx][1] = 0.0;
        derivs[num_vtx][2] = 2.0;
        derivs[1][2] -= 1.0;
        ++num_vtx;
      }
      break;
    
    case 1:
      vertices[0] = 1;
      derivs[0][0] = 1.0;
      derivs[0][1] = 0.0;
      derivs[0][2] = 0.0;
      
      vertices[1] = 2;
      derivs[1][0] = 0.0;
      derivs[1][1] = 1.0;
      derivs[1][2] = 0.0;
   
      if (nodeset.mid_edge_node(0) == nodeset.mid_edge_node(2)) {
        vertices[num_vtx] = 0;
        sign = 2*nodeset.mid_edge_node(0) - 1;
        derivs[num_vtx][0] = sign;
        derivs[num_vtx][1] = sign;
        derivs[num_vtx][2] = sign;
        ++num_vtx;
      }
      
      if (nodeset.mid_edge_node(4) == nodeset.mid_edge_node(5)) {
        vertices[num_vtx] = 3;
        sign = 1 - 2*nodeset.mid_edge_node(4);
        derivs[num_vtx][0] =  0.0;
        derivs[num_vtx][1] =  0.0;
        derivs[num_vtx][2] = sign;
        ++num_vtx;
      }
      
      if (nodeset.mid_edge_node(0)) {
        vertices[num_vtx] = 4;
        derivs[num_vtx][0] = -2.0;
        derivs[num_vtx][1] = -2.0;
        derivs[num_vtx][2] = -2.0;
        ++num_vtx;
        derivs[0][0] += 1.0;
        derivs[0][1] += 1.0;
        derivs[0][2] += 1.0;
     }
      
      if (nodeset.mid_edge_node(1)) {
        vertices[num_vtx] = 5;
        derivs[num_vtx][0] = 2.0;
        derivs[num_vtx][1] = 2.0;
        derivs[num_vtx][2] = 0.0;
        ++num_vtx;
        derivs[0][0] -= 1.0;
        derivs[0][1] -= 1.0;
        derivs[1][0] -= 1.0;
        derivs[1][1] -= 1.0;
      }
      
      if (nodeset.mid_edge_node(2)) {
        vertices[num_vtx] = 6;
        derivs[num_vtx][0] = -2.0;
        derivs[num_vtx][1] = -2.0;
        derivs[num_vtx][2] = -2.0;
        ++num_vtx;
        derivs[1][0] += 1.0;
        derivs[1][1] += 1.0;
        derivs[1][2] += 1.0;
      }
      
      if (nodeset.mid_edge_node(4)) {
        vertices[num_vtx] = 8;
        derivs[num_vtx][0] = 0.0;
        derivs[num_vtx][1] = 0.0;
        derivs[num_vtx][2] = 2.0;
        ++num_vtx;
        derivs[0][2] -= 1.0;
      }
      
      if (nodeset.mid_edge_node(5)) {
        vertices[num_vtx] = 9;
        derivs[num_vtx][0] = 0.0;
        derivs[num_vtx][1] = 0.0;
        derivs[num_vtx][2] = 2.0;
        ++num_vtx;
        derivs[1][2] -= 1.0;
      }
      break;
      
    case 2:
      vertices[0] = 0;
      derivs[0][0] = -1.0;
      derivs[0][1] = -1.0;
      derivs[0][2] = -1.0;
      
      vertices[1] = 2;
      derivs[1][0] = 0.0;
      derivs[1][1] = 1.0;
      derivs[1][2] = 0.0;
   
      if (nodeset.mid_edge_node(0) == nodeset.mid_edge_node(1)) {
        vertices[num_vtx] = 1;
        sign = 1 - 2*nodeset.mid_edge_node(0);
        derivs[num_vtx][0] = sign;
        derivs[num_vtx][1] =  0.0;
        derivs[num_vtx][2] =  0.0;
        ++num_vtx;
      }
      
      if (nodeset.mid_edge_node(3) == nodeset.mid_edge_node(5)) {
        vertices[num_vtx] = 3;
        sign = 1 - 2*nodeset.mid_edge_node(3);
        derivs[num_vtx][0] =  0.0;
        derivs[num_vtx][1] =  0.0;
        derivs[num_vtx][2] = sign;
        ++num_vtx;
      }
      
      if (nodeset.mid_edge_node(0)) {
        vertices[num_vtx] = 4;
        derivs[num_vtx][0] =  2.0;
        derivs[num_vtx][1] =  0.0;
        derivs[num_vtx][2] =  0.0;
        ++num_vtx;
        derivs[0][0] -= 1.0;
      }
      
      if (nodeset.mid_edge_node(1)) {
        vertices[num_vtx] = 5;
        derivs[num_vtx][0] =  2.0;
        derivs[num_vtx][1] =  0.0;
        derivs[num_vtx][2] =  0.0;
        ++num_vtx;
        derivs[1][0] -= 1.0;
      }
      
      if (nodeset.mid_edge_node(2)) {
        vertices[num_vtx] = 6;
        derivs[num_vtx][0] = -2.0;
        derivs[num_vtx][1] =  0.0;
        derivs[num_vtx][2] = -2.0;
        ++num_vtx;
        derivs[0][0] += 1.0;
        derivs[0][2] += 1.0;
        derivs[1][0] += 1.0;
        derivs[1][2] += 1.0;
      }
       
      if (nodeset.mid_edge_node(3)) {
        vertices[num_vtx] = 7;
        derivs[num_vtx][0] = 0.0;
        derivs[num_vtx][1] = 0.0;
        derivs[num_vtx][2] = 2.0;
        ++num_vtx;
        derivs[0][2] -= 1.0;
      }
      
      if (nodeset.mid_edge_node(5)) {
        vertices[num_vtx] = 9;
        derivs[num_vtx][0] = 0.0;
        derivs[num_vtx][1] = 0.0;
        derivs[num_vtx][2] = 2.0;
        ++num_vtx;
        derivs[1][2] -= 1.0;
      }
      break;
    
    case 3:
      vertices[0] = 0;
      derivs[0][0] = -1.0;
      derivs[0][1] = -1.0;
      derivs[0][2] = -1.0;
      
      vertices[1] = 3;
      derivs[1][0] = 0.0;
      derivs[1][1] = 0.0;
      derivs[1][2] = 1.0;
      
      if (nodeset.mid_edge_node(0) == nodeset.mid_edge_node(4)) {
        vertices[num_vtx] = 1;
        sign = 1 - 2*nodeset.mid_edge_node(0);
        derivs[num_vtx][0] = sign;
        derivs[num_vtx][1] =  0.0;
        derivs[num_vtx][2] =  0.0;
        ++num_vtx;
      }
      
      if (nodeset.mid_edge_node(2) == nodeset.mid_edge_node(5)) {
        vertices[num_vtx] = 2;
        sign = 1 - 2*nodeset.mid_edge_node(2);
        derivs[num_vtx][0] =  0.0;
        derivs[num_vtx][1] = sign;
        derivs[num_vtx][2] =  0.0;
        ++num_vtx;
      }

      if (nodeset.mid_edge_node(0)) {
        vertices[num_vtx] = 4;
        derivs[num_vtx][0] =  2.0;
        derivs[num_vtx][1] =  0.0;
        derivs[num_vtx][2] =  0.0;
        ++num_vtx;
        derivs[0][0] -= 1.0;
      }
      
      if (nodeset.mid_edge_node(2)) {
        vertices[num_vtx] = 6;
        derivs[num_vtx][0] =  0.0;
        derivs[num_vtx][1] =  2.0;
        derivs[num_vtx][2] =  0.0;
        ++num_vtx;
        derivs[0][1] -= 1.0;
      }
      
      if (nodeset.mid_edge_node(3)) {
        vertices[num_vtx] = 7;
        derivs[num_vtx][0] = -2.0;
        derivs[num_vtx][1] = -2.0;
        derivs[num_vtx][2] =  0.0;
        ++num_vtx;
        derivs[0][0] += 1.0;
        derivs[0][1] += 1.0;
        derivs[1][0] += 1.0;
        derivs[1][1] += 1.0;
      }
      
      if (nodeset.mid_edge_node(4)) {
        vertices[num_vtx] = 8;
        derivs[num_vtx][0] =  2.0;
        derivs[num_vtx][1] =  0.0;
        derivs[num_vtx][2] =  0.0;
        ++num_vtx;
        derivs[1][0] -= 1.0;
      }
      
      if (nodeset.mid_edge_node(5)) {
        vertices[num_vtx] = 9;
        derivs[num_vtx][0] =  0.0;
        derivs[num_vtx][1] =  2.0;
        derivs[num_vtx][2] =  0.0;
        ++num_vtx;
        derivs[1][1] -= 1.0;
      }
      break;
    
    case 4:
      vertices[0] = 1;
      derivs[0][0] = 1.0;
      derivs[0][1] = 0.0;
      derivs[0][2] = 0.0;
      
      vertices[1] = 3;
      derivs[1][0] = 0.0;
      derivs[1][1] = 0.0;
      derivs[1][2] = 1.0;

      if (nodeset.mid_edge_node(0) == nodeset.mid_edge_node(3)) {
        vertices[num_vtx] = 0;
        sign = 2*nodeset.mid_edge_node(0) - 1;
        derivs[num_vtx][0] = sign;
        derivs[num_vtx][1] = sign;
        derivs[num_vtx][2] = sign;
        ++num_vtx;
      }
      
      if (nodeset.mid_edge_node(1) == nodeset.mid_edge_node(5)) {
        vertices[num_vtx] = 2;
        sign = 1 - 2*nodeset.mid_edge_node(1);
        derivs[num_vtx][0] =  0.0;
        derivs[num_vtx][1] = sign;
        derivs[num_vtx][2] =  0.0;
        ++num_vtx;
      }

      if (nodeset.mid_edge_node(0)) {
        vertices[num_vtx] = 4;
        derivs[num_vtx][0] = -2.0;
        derivs[num_vtx][1] = -2.0;
        derivs[num_vtx][2] = -2.0;
        ++num_vtx;
        derivs[0][0] += 1.0;
        derivs[0][1] += 1.0;
        derivs[0][2] += 1.0;
      }
      
      if (nodeset.mid_edge_node(1)) {
        vertices[num_vtx] = 5;
        derivs[num_vtx][0] =  0.0;
        derivs[num_vtx][1] =  2.0;
        derivs[num_vtx][2] =  0.0;
        ++num_vtx;
        derivs[0][1] -= 1.0;
      }
 
      if (nodeset.mid_edge_node(3)) {
        vertices[num_vtx] = 7;
        derivs[num_vtx][0] = -2.0;
        derivs[num_vtx][1] = -2.0;
        derivs[num_vtx][2] = -2.0;
        ++num_vtx;
        derivs[1][0] += 1.0;
        derivs[1][1] += 1.0;
        derivs[1][2] += 1.0;
      }
      
      if (nodeset.mid_edge_node(4)) {
        vertices[num_vtx] = 8;
        derivs[num_vtx][0] = 2.0;
        derivs[num_vtx][1] =  0.0;
        derivs[num_vtx][2] = 2.0;
        ++num_vtx;
        derivs[0][0] -= 1.0;
        derivs[0][2] -= 1.0;
        derivs[1][0] -= 1.0;
        derivs[1][2] -= 1.0;
      }
      
      if (nodeset.mid_edge_node(5)) {
        vertices[num_vtx] = 9;
        derivs[num_vtx][0] =  0.0;
        derivs[num_vtx][1] =  2.0;
        derivs[num_vtx][2] =  0.0;
        ++num_vtx;
        derivs[1][1] -= 1.0;
      }
      break;
    
    case 5:
      vertices[0] = 2;
      derivs[0][0] = 0.0;
      derivs[0][1] = 1.0;
      derivs[0][2] = 0.0;
      
      vertices[1] = 3;
      derivs[1][0] = 0.0;
      derivs[1][1] = 0.0;
      derivs[1][2] = 1.0;

      if (nodeset.mid_edge_node(2) == nodeset.mid_edge_node(3)) {
        vertices[num_vtx] = 0;
        sign = 2*nodeset.mid_edge_node(2) - 1;
        derivs[num_vtx][0] = sign;
        derivs[num_vtx][1] = sign;
        derivs[num_vtx][2] = sign;
        ++num_vtx;
      }
      
      if (nodeset.mid_edge_node(1) == nodeset.mid_edge_node(4)) {
        vertices[num_vtx] = 1;
        sign = 1 - 2*nodeset.mid_edge_node(1);
        derivs[num_vtx][0] = sign;
        derivs[num_vtx][1] =  0.0;
        derivs[num_vtx][2] =  0.0;
        ++num_vtx;
      }
      
      if (nodeset.mid_edge_node(1)) {
        vertices[num_vtx] = 5;
        derivs[num_vtx][0] =  2.0;
        derivs[num_vtx][1] =  0.0;
        derivs[num_vtx][2] =  0.0;
        ++num_vtx;
        derivs[0][0] -= 1.0;
      }

      if (nodeset.mid_edge_node(2)) {
        vertices[num_vtx] = 6;
        derivs[num_vtx][0] = -2.0;
        derivs[num_vtx][1] = -2.0;
        derivs[num_vtx][2] = -2.0;
        ++num_vtx;
        derivs[0][0] += 1.0;
        derivs[0][1] += 1.0;
        derivs[0][2] += 1.0;
      }
      
      if (nodeset.mid_edge_node(3)) {
        vertices[num_vtx] = 7;
        derivs[num_vtx][0] = -2.0;
        derivs[num_vtx][1] = -2.0;
        derivs[num_vtx][2] = -2.0;
        ++num_vtx;
        derivs[1][0] += 1.0;
        derivs[1][1] += 1.0;
        derivs[1][2] += 1.0;
      }

      if (nodeset.mid_edge_node(4)) {
        vertices[num_vtx] = 8;
        derivs[num_vtx][0] =  2.0;
        derivs[num_vtx][1] =  0.0;
        derivs[num_vtx][2] =  0.0;
        ++num_vtx;
        derivs[1][0] -= 1.0;
      }
       
      if (nodeset.mid_edge_node(5)) {
        vertices[num_vtx] = 9;
        derivs[num_vtx][0] = 0.0;
        derivs[num_vtx][1] = 2.0;
        derivs[num_vtx][2] = 2.0;
        ++num_vtx;
        derivs[0][1] -= 1.0;
        derivs[0][2] -= 1.0;
        derivs[1][1] -= 1.0;
        derivs[1][2] -= 1.0;
      }
      break;
  }
}

// Derivatives of coefficients for mid-edge nodes

const double ft = 4.0/3.0;

const double ho_dr[6][4] = { { 0.,-ft, ft, 0. },
                             { 0., ft, ft, ft },
                             { 0.,-ft,-ft,-ft },
                             {-ft,-ft,-ft, 0. },
                             { ft, ft, ft, 0. },
                             { 0., 0., 0., 0. } };

const double ho_ds[6][4] = { {-ft,-ft, 0.,-ft },
                             { ft, ft, 0., ft },
                             { ft,-ft, 0., 0. },
                             {-ft,-ft,-ft, 0. },
                             { 0., 0., 0., 0. },
                             { ft, ft, ft, 0. } };

const double ho_dt[6][4] = { {-ft,-ft, 0.,-ft },
                             { 0., 0., 0., 0. },
                             { 0.,-ft,-ft,-ft },
                             { 0.,-ft, 0., ft },
                             { ft, ft, 0., ft },
                             { 0., ft, ft, ft } };


static void derivatives_at_mid_face( unsigned face,
                                     NodeSet nodeset,
                                     size_t* vertices,
                                     MsqVector<3>* derivs,
                                     size_t& num_vtx )
{
    // begin with derivatives for linear tetrahedron
  num_vtx = 4;
  get_linear_derivatives( vertices, derivs );
  
  for (unsigned i = 0; i < 6; ++i) 
    if (nodeset.mid_edge_node(i)) {
      vertices[num_vtx] = i+4;
      derivs[num_vtx][0] = ho_dr[i][face];
      derivs[num_vtx][1] = ho_ds[i][face];
      derivs[num_vtx][2] = ho_dt[i][face];
      ++num_vtx;
      int j = edges[i][0];
      derivs[j][0] -= 0.5*ho_dr[i][face];
      derivs[j][1] -= 0.5*ho_ds[i][face];
      derivs[j][2] -= 0.5*ho_dt[i][face];
      j = edges[i][1];
      derivs[j][0] -= 0.5*ho_dr[i][face];
      derivs[j][1] -= 0.5*ho_ds[i][face];
      derivs[j][2] -= 0.5*ho_dt[i][face];
    }
}

  
  // position (0->r, 1->s, 2->t) of zero-valued term for mid-edge node
static const int zeros[6] = { 0, 2, 1, 2, 1, 0 };
  // value of mid-edge terms
static const int signs[6] = { -1, 1, -1, -1, 1, 1 };

static void derivatives_at_mid_elem( NodeSet nodeset,
                                     size_t* vertices,
                                     MsqVector<3>* derivs,
                                     size_t& num_vtx )
{
                                    
  bool corners[4] = { false, false, false, false };
  double corner_vals[4][3] = { {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0} };
  
  num_vtx = 0;
  for (unsigned i = 4;  i < 10; ++i) {
    int sign = signs[i-4];
    int zero = zeros[i-4];
   
    if (nodeset.mid_edge_node(i-4)) {
      vertices[num_vtx] = i;
      derivs[num_vtx][0] = (double)sign;
      derivs[num_vtx][1] = (double)sign;
      derivs[num_vtx][2] = (double)sign;
      derivs[num_vtx][zero] = 0.0;
      ++num_vtx;
    }
    else {
      for (unsigned j = 0; j < 2; ++j) {
        int corner = edges[i-4][j];
        int v1 = (zero + 1) % 3;
        int v2 = (zero + 2) % 3;
        corners[corner] = true;
        corner_vals[corner][v1] += 0.5*sign;
        corner_vals[corner][v2] += 0.5*sign;
      }
    }
  }
  
  for (unsigned i = 0; i < 4; ++i)
    if (corners[i]) {
      vertices[num_vtx] = i;
      derivs[num_vtx][0] = corner_vals[i][0];
      derivs[num_vtx][1] = corner_vals[i][1];
      derivs[num_vtx][2] = corner_vals[i][2];
      ++num_vtx;
    }
}
    
void TetLagrangeShape::derivatives( Sample loc,
                                    NodeSet nodeset,
                                    size_t* vertex_indices_out,
                                    MsqVector<3>* d_coeff_d_xi_out,
                                    size_t& num_vtx,
                                    MsqError& err ) const
{
  if (!nodeset.have_any_mid_node()) {
    num_vtx = 4;
    get_linear_derivatives( vertex_indices_out, d_coeff_d_xi_out );
    return;
  }

  if (nodeset.have_any_mid_face_node() | nodeset.have_any_mid_region_node()) {
    MSQ_SETERR(err)("TetLagrangeShape does not support mid-face/mid-element nodes",
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
      derivatives_at_mid_face( loc.number, nodeset, vertex_indices_out, d_coeff_d_xi_out, num_vtx );
      break;
    case 3:
      derivatives_at_mid_elem( nodeset, vertex_indices_out, d_coeff_d_xi_out, num_vtx );
      break;
    default:
      MSQ_SETERR(err)("Invalid/unsupported logical dimension",MsqError::INVALID_ARG);
  }
}

void TetLagrangeShape::ideal( Sample , 
                              MsqMatrix<3,3>& J,
                              MsqError&  ) const
{
  const double a = 1.122462048309373;  // 6th root of 2
  const double b = 1.7320508075688772; // sqrt(3)
  const double c = 1.5874010519681994; // 2 to the 2/3
  J(0,0) = a;   J(0,1) = 0.5*a;   J(0,2) = 0.5*a;
  J(1,0) = 0.0; J(1,1) = 0.5*a*b; J(1,2) = 0.5*a/b;
  J(2,0) = 0.0; J(2,1) = 0.0;     J(2,2) = c/b;
}


} // namespace Mesquite
