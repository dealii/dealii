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


/** \file LinearPyramid.cpp
 *  \brief LinearPyramid implementation
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "LinearPyramid.hpp"
#include "MsqError.hpp"

namespace MESQUITE_NS {

static const char* nonlinear_error 
 = "Attempt to use LinearTriangle mapping function for a nonlinear element\n";

static inline void set_equal_derivatives( double value, 
                                          size_t* indices,
                                          MsqVector<3>* derivs,
                                          size_t& num_vtx )
{
  num_vtx = 5;
  indices[0] = 0;
  indices[1] = 1;
  indices[2] = 2;
  indices[3] = 3;
  indices[4] = 4;
    
  derivs[0][0] = -value;
  derivs[0][1] = -value;
  derivs[0][2] = -0.25;
  
  derivs[1][0] =  value;
  derivs[1][1] = -value;
  derivs[1][2] = -0.25;
  
  derivs[2][0] =  value;
  derivs[2][1] =  value;
  derivs[2][2] = -0.25;
  
  derivs[3][0] = -value;
  derivs[3][1] =  value;
  derivs[3][2] = -0.25;
  
  derivs[4][0] =  0.0;
  derivs[4][1] =  0.0;
  derivs[4][2] =  1.0;
}

static inline void set_edge_derivatives( unsigned base_corner,
                                         double value,
                                         size_t* indices,
                                         MsqVector<3>* derivs,
                                         size_t& num_vtx )
{
  const int direction = base_corner % 2;
  const int edge_beg =  base_corner;
  const int edge_end = (base_corner+1)%4;
  const int adj_end  = (base_corner+2)%4;
  const int adj_beg  = (base_corner+3)%4;
  const int dir_sign = 2*(edge_beg/2) - 1;
  const int oth_sign = 2*((edge_beg+1)/2%2) - 1;

  num_vtx = 5;
  indices[0] = edge_beg;
  indices[1] = edge_end;
  indices[2] = adj_end;
  indices[3] = adj_beg;
  indices[4] = 4;

  derivs[0][  direction] =  2 * dir_sign * value;
  derivs[0][1-direction] =      oth_sign * value;
  derivs[0][2] = -0.5;

  derivs[1][  direction] = -2 * dir_sign * value;
  derivs[1][1-direction] =      oth_sign * value;
  derivs[1][2] = -0.5;

  derivs[2][  direction] =  0.0;
  derivs[2][1-direction] = -oth_sign * value;
  derivs[2][2]           =  0.0;

  derivs[3][  direction] =  0.0;
  derivs[3][1-direction] = -oth_sign * value;
  derivs[3][2]           =  0.0;

  derivs[4][0] = 0.0;
  derivs[4][1] = 0.0;
  derivs[4][2] = 1.0;
}

static inline void set_corner_derivatives( unsigned corner,
                                           double value,
                                           size_t* indices,
                                           MsqVector<3>* derivs,
                                           size_t& num_vtx )
{
  const unsigned adj_in_xi = (5 - corner) % 4;
  const unsigned adj_in_eta = 3 - corner;

  const int dxi_sign  = 2*((corner+1)/2%2)-1;
  const int deta_sign = 2*(corner/2) - 1;
  const double dxi_value = dxi_sign * value;
  const double deta_value = deta_sign * value;

  num_vtx = 4;
  indices[0] = corner;
  indices[1] = adj_in_xi;
  indices[2] = adj_in_eta;
  indices[3] = 4;

  derivs[0][0] =  dxi_value;
  derivs[0][1] =  deta_value;
  derivs[0][2] = -1.0;

  derivs[1][0] = -dxi_value;
  derivs[1][1] =  0.0;
  derivs[1][2] =  0.0;

  derivs[2][0] =  0.0;
  derivs[2][1] = -deta_value;
  derivs[2][2] =  0.0;

  derivs[3][0] =  0.0;
  derivs[3][1] =  0.0;
  derivs[3][2] =  1.0;
}

EntityTopology LinearPyramid::element_topology() const
  { return PYRAMID; }
  
int LinearPyramid::num_nodes() const
  { return 5; }

NodeSet LinearPyramid::sample_points( NodeSet ) const
{
  NodeSet result;
  result.set_all_corner_nodes(PYRAMID);
  result.clear_corner_node(4);
  return result;
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
                                      double* coeff_out,
                                      size_t* indices_out,
                                      size_t& num_coeff )
{
  num_coeff = 2;
  coeff_out[0] = 0.5;
  coeff_out[1] = 0.5;
  
  if (edge < 4) {
    indices_out[0] = edge;
    indices_out[1] = (edge+1)%4;
  }
  else {
    indices_out[0] = edge-4;
    indices_out[1] = 4;
  }
}

static void coefficients_at_mid_face( unsigned face,
                                      double* coeff_out,
                                      size_t* indices_out,
                                      size_t& num_coeff )
{
  if (face == 4) {
    num_coeff = 4;
    coeff_out[0] = 0.25;
    coeff_out[1] = 0.25;
    coeff_out[2] = 0.25;
    coeff_out[3] = 0.25;
    indices_out[0] = 0;
    indices_out[1] = 1;
    indices_out[2] = 2;
    indices_out[3] = 3;
  }
  else {
    num_coeff = 3;
    indices_out[0] = face;
    indices_out[1] = (face+1)%4;
    indices_out[2] = 4;
    coeff_out[0] = 0.25;
    coeff_out[1] = 0.25;
    coeff_out[2] = 0.50;
  }
}

static void coefficients_at_mid_elem( double* coeff_out,
                                      size_t* indices_out,
                                      size_t& num_coeff )
{
  num_coeff = 5;
  coeff_out[0] = 0.125;
  coeff_out[1] = 0.125;
  coeff_out[2] = 0.125;
  coeff_out[3] = 0.125;
  coeff_out[4] = 0.500;
  indices_out[0] = 0;
  indices_out[1] = 1;
  indices_out[2] = 2;
  indices_out[3] = 3;
  indices_out[4] = 4;
}

void LinearPyramid::coefficients( Sample loc,
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

void LinearPyramid::derivatives( Sample loc,
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
      if (loc.number == 4) {
        set_equal_derivatives( 0.0, vertex_indices_out, d_coeff_d_xi_out, num_vtx );
      }
      else {
        set_corner_derivatives( loc.number, 1.0, vertex_indices_out, d_coeff_d_xi_out, num_vtx );
      }
      break;
    case 1:
      if (loc.number < 4) {
        set_edge_derivatives( loc.number, 0.50, vertex_indices_out, d_coeff_d_xi_out, num_vtx );
      }
      else {
        set_corner_derivatives( loc.number-4, 0.50, vertex_indices_out, d_coeff_d_xi_out, num_vtx );
      }    
      break;
    case 2:
      if (loc.number == 4) {
        set_equal_derivatives( 0.5, vertex_indices_out, d_coeff_d_xi_out, num_vtx );
      }
      else {
        set_edge_derivatives( loc.number, 0.25, vertex_indices_out, d_coeff_d_xi_out, num_vtx );
      }
      break;
    case 3:
      set_equal_derivatives( 0.25, vertex_indices_out, d_coeff_d_xi_out, num_vtx );
      break;
    default:
      MSQ_SETERR(err)("Invalid/unsupported logical dimension",MsqError::INVALID_ARG);
  }
}

void LinearPyramid::ideal( Sample location,
                           MsqMatrix<3,3>& J,
                           MsqError& err ) const
{
  // For an ideal element with unit edge length at the base and unit
  // height, the Jacobian matrix is:
  // | 1-zeta    0      1/2 - xi  |
  // |  0       1-zeta  1/2 - eta |
  // |  0        0       1        |
  //
  // The coefficient to produce a unit determinant
  // is therefore (1-zeta)^(-2/3).  
  // 
  // Thus the unit-determinate ideal element Jacobian
  // is, given alpha = (1-zeta)^(-1/3):
  //
  // | 1/alpha  0      alpha^2 (1/2 - xi) |
  // | 0       1/alpha alpha^2 (1/2 - eta)|
  // | 0        0      alpha^2            |
  //
  // There are only three zeta values of interest:
  //  zeta = 1 : the degenerate case
  //  zeta = 0 : both 1/alpha and alpha^2 are 1.0
  //  zeta = 1/2 : 1/alpha = 1/cbrt(2.0) and alpha^2 = 2*(1/alpha)
  
    // special case for apex
  if (location.dimension == 0 && location.number == 4) {
    J = MsqMatrix<3,3>(0.0);
    return;
  }
    
    // These are always zero
  J(0,1) = J(1,0) = J(2,0) = J(2,1) = 0.0;

    // Set diagonal terms and magnitude of terms in 3rd column based on zeta

    // All of the zeta=0 locations
  double f;
  if ( location.dimension == 0 || 
      (location.dimension == 1 && location.number < 4) ||
      (location.dimension == 2 && location.number == 4)) {
      J(0,0) = J(1,1) = J(2,2) = 1.0;
      f = 0.5;
  }
    // all of the zeta=1/2 locations
  else {
    f = J(0,0) = J(1,1) = 0.79370052598409979;
    J(2,2) = 2.0*f;
  }

    // Set terms in 3rd column based on xi,eta
  
    // The xi = eta = 0.5 locations (mid-element in xi and eta)
  if ( location.dimension == 3 ||
      (location.dimension == 2 && location.number == 4)) {
        J(0,2) = J(1,2) = 0.0;
  }
    // The corner locations
  else if ( location.dimension == 0 ||  
           (location.dimension == 1 && location.number >= 4)) {
    switch (location.number % 4) {
      case 0: J(0,2) =  f; J(1,2) =  f; break;
      case 1: J(0,2) = -f; J(1,2) =  f; break;
      case 2: J(0,2) = -f; J(1,2) = -f; break;
      case 3: J(0,2) =  f; J(1,2) = -f; break;
    }
  }
    // The mid-edge locations
  else {
    switch (location.number) {
      case 0: J(0,2) =  0; J(1,2) =  f; break;
      case 1: J(0,2) = -f; J(1,2) =  0; break;
      case 2: J(0,2) =  0; J(1,2) = -f; break;
      case 3: J(0,2) =  f; J(1,2) =  0; break;
    }
  }  
}

} // namespace Mesquite
