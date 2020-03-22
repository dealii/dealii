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

#include "Mesquite.hpp"
#include "MsqError.hpp"
#include "LinearTetrahedron.hpp"

namespace MESQUITE_NS {

static const char* nonlinear_error 
 = "Attempt to use LinearTetrahedron mapping function for a nonlinear element\n";
 
EntityTopology LinearTetrahedron::element_topology() const
  { return TETRAHEDRON; }
  
int LinearTetrahedron::num_nodes() const
  { return 4; }

NodeSet LinearTetrahedron::sample_points( NodeSet ) const
{
  NodeSet result;
  result.set_mid_region_node();
  return result;
}

static const unsigned faces[][3] = { { 1, 0, 3 },
                                     { 3, 2, 1 },
                                     { 0, 2, 3 },
                                     { 0, 1, 2 } };

void LinearTetrahedron::coefficients( Sample location,
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
  
  switch (location.dimension) {
    case 0:
      num_coeff = 1;
      indices_out[0] = location.number;
      coeff_out[0] = 1.0;
      break;
    case 1:
      num_coeff = 2;
      coeff_out[0] = 0.5;
      coeff_out[1] = 0.5;
      if (location.number < 3) {
        indices_out[0] = location.number;
        indices_out[1] = (location.number+1)%3;
      }
      else {
        indices_out[0] = location.number - 3;
        indices_out[1] = 3;
      }
      break;
    case 2:
      num_coeff = 3;
      indices_out[0] = faces[location.number][0];
      indices_out[1] = faces[location.number][1];
      indices_out[2] = faces[location.number][2];
      coeff_out[0] = coeff_out[1] = coeff_out[2] = coeff_out[3] = MSQ_ONE_THIRD;
      break;
    case 3:
      num_coeff = 4;
      indices_out[0] = 0;
      indices_out[1] = 1;
      indices_out[2] = 2;
      indices_out[3] = 3;
      coeff_out[0] = coeff_out[1] = coeff_out[2] = coeff_out[3] = 0.25;
      break;
    default:
      MSQ_SETERR(err)("Invalid/unsupported logical dimension",MsqError::INVALID_ARG);
  }
}

void LinearTetrahedron::derivatives( Sample ,
                                     NodeSet nodeset,
                                     size_t* vertices,
                                     MsqVector<3>* coeff_derivs,
                                     size_t& num_vtx,
                                     MsqError& err ) const
{
  if (nodeset.have_any_mid_node()) {
    MSQ_SETERR(err)(nonlinear_error, MsqError::UNSUPPORTED_ELEMENT );
    return;
  }
  else {
    num_vtx = 4;
    vertices[0] = 0;
    vertices[1] = 1;
    vertices[2] = 2;
    vertices[3] = 3;
    
    coeff_derivs[0][0] = -1.0;
    coeff_derivs[0][1] = -1.0;
    coeff_derivs[0][2] = -1.0;

    coeff_derivs[1][0] =  1.0;
    coeff_derivs[1][1] =  0.0;
    coeff_derivs[1][2] =  0.0;

    coeff_derivs[2][0] =  0.0;
    coeff_derivs[2][1] =  1.0;
    coeff_derivs[2][2] =  0.0;

    coeff_derivs[3][0] =  0.0;
    coeff_derivs[3][1] =  0.0;
    coeff_derivs[3][2] =  1.0;
  }
}

void LinearTetrahedron::ideal( Sample , 
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
