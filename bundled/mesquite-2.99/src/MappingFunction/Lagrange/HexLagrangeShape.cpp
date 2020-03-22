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
/** \file HexLagrangeShape.cpp
 *  \author Nicholas Voshell
 */
 
#include "HexLagrangeShape.hpp"
#include "MsqError.hpp"

namespace MESQUITE_NS {

EntityTopology HexLagrangeShape::element_topology() const
{ return HEXAHEDRON; }
  
int HexLagrangeShape::num_nodes() const
{ return 27; }

void HexLagrangeShape::coefficients( Sample loc,
                                     NodeSet nodeset,
                                     double* coeff_out,
                                     size_t* indices_out,
                                     size_t& num_coeff,
                                     MsqError& err ) const
{
  if (nodeset.num_nodes() != 27) {
    MSQ_SETERR(err)("Mapping function supports only 27-node hexahedra with no slaved nodes.",
                    MsqError::UNSUPPORTED_ELEMENT);
    return;
  }

  switch (loc.dimension) {
  case 0: //Corner sample point - assume that it is there always
      num_coeff = 1;
      indices_out[0] = loc.number;
      coeff_out[0] = 1.0;
      break;
    case 1: //Line sample point - check if it is there,
      num_coeff = 1;
      indices_out[0] = loc.number+8;
      coeff_out[0] = 1.0;
      break;
    case 2: //Face sample point - check if it is there,
      num_coeff = 1;
      indices_out[0] = loc.number+20;
      coeff_out[0] = 1.0;
      break;
    case 3: //Center sample point - check if it is there, 
      num_coeff = 1;
      indices_out[0] = 26;
      coeff_out[0] = 1.0;
      break;
    default:
      MSQ_SETERR(err)(MsqError::UNSUPPORTED_ELEMENT,
                  "Request for dimension %d mapping function value"
                  "for a quadratic hexahedral element", loc.dimension);
  }

}
     

void HexLagrangeShape::derivatives( Sample loc,
                                    NodeSet nodeset,
				    size_t* vertices,
                                    MsqVector<3>* derivs,
                                    size_t& num_vtx,
                                    MsqError& err ) const
{
  if (nodeset.num_nodes() != 27) {
    MSQ_SETERR(err)("Mapping function supports only 27-node hexahedra with no slaved nodes.",
                    MsqError::UNSUPPORTED_ELEMENT);
    return;
  }

  //r coordinate
  const int HLS1[]= {1, 3, 3, 1, 1, 3, 3, 1, 2, 3, 2, 1, 1, 3, 3, 1, 2, 3, 2, 1, 2, 3, 2, 1, 2, 2, 2};
  //s coordinate
  const int HLS2[]= {1, 1, 3, 3, 1, 1, 3, 3, 1, 2, 3, 2, 1, 1, 3, 3, 1, 2, 3, 2, 1, 2, 3, 2, 2, 2, 2};
  //t coordinate
  const int HLS3[]= {1, 1, 1, 1, 3, 3, 3, 3, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 2, 2, 2, 2, 1, 3, 2};

  const int HLSup1[]= {12, 13, 14, 15, -1, -1, -1, -1, 20, 21, 22, 23, 4, 5, 6, 7, -1, -1, -1, -1, 16, 17, 18, 19, 26, -1, 25};
  const int HLSup2[]= {4, 5, 6, 7, -1, -1, -1, -1, 16, 17, 18, 19, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 25, -1, -1};
  const int HLSdn1[]= {-1, -1, -1, -1, 12, 13, 14, 15, -1, -1, -1, -1, 0, 1, 2, 3, 20, 21, 22, 23, 8, 9, 10, 11, -1, 26, 24};
  const int HLSdn2[]= {-1, -1, -1, -1, 0, 1, 2, 3, -1, -1, -1, -1, -1, -1, -1, -1, 8, 9, 10, 11, -1, -1, -1, -1, -1, 24, -1};

  const int HLSlt1[]= {-1, 8, 10, -1, -1, 16, 18, -1, 0, 24, 3, -1, -1, 20, 22, -1, 4, 25, 7, -1, 12, 26, 15, -1, 11, 19, 23};
  const int HLSlt2[]= {-1, 0, 3, -1, -1, 4, 7, -1, -1, 11, -1, -1, -1, 12, 15, -1, -1, 19, -1, -1, -1, 23, -1, -1, -1, -1, -1};
  const int HLSrt1[]= {8, -1, -1, 10, 16, -1, -1, 18, 1, -1, 2, 24, 20, -1, -1, 22, 5, -1, 6, 25, 13, -1, 14, 26, 9, 17, 21};
  const int HLSrt2[]= {1, -1, -1, 2, 5, -1, -1, 6, -1, -1, -1, 9, 13, -1, -1, 14, -1, -1, -1, 17, -1, -1, -1, 21, -1, -1, -1};

  const int HLSft1[]= {-1, -1, 9, 11, -1, -1, 17, 19, -1, 1, 24, 0, -1, -1, 21, 23, -1, 5, 25, 4, -1, 13, 26, 12, 8, 16, 20};
  const int HLSft2[]= {-1, -1, 1, 0, -1, -1, 5, 4, -1, -1, 8, -1, -1, -1, 13, 12, -1, -1, 16, -1, -1, -1, 20, -1, -1, -1, -1};
  const int HLSbk1[]= {11, 9, -1, -1, 19, 17, -1, -1, 24, 2, -1, 3, 23, 21, -1, -1, 25, 6, -1, 7, 26, 14, -1, 15, 10, 18, 22};
  const int HLSbk2[]= {3, 2, -1, -1, 7, 6, -1, -1, 10, -1, -1, -1, 15, 14, -1, -1, 18, -1, -1, -1, 22, -1, -1, -1, -1, -1, -1};

  double entries[] = {0, -3, 0, 3};

  int location=loc.number;

  switch (loc.dimension) {
    case 1:
      location+=8;
      break;
    case 2:
      location+=20;
      break;
    case 3:
      location+=26;
      break;
  }

  num_vtx=0;
  
  if (location != 26) {
    vertices[num_vtx]=location;
    derivs[num_vtx][0]=entries[HLS1[location]];
    derivs[num_vtx][1]=entries[HLS2[location]]; 
    derivs[num_vtx][2]=entries[HLS3[location]];  
    num_vtx++;
  }
  
  if(HLSup1[location]!=-1){
    vertices[num_vtx]=HLSup1[location];
    if(HLS3[location]==2) {
      derivs[num_vtx][0]=0; derivs[num_vtx][1]=0; derivs[num_vtx][2]=1;
    } else {
      derivs[num_vtx][0]=0; derivs[num_vtx][1]=0; derivs[num_vtx][2]=4;
    }
    num_vtx++;
  }
  
  if(HLSdn1[location]!=-1){
    vertices[num_vtx]=HLSdn1[location];
    if(HLS3[location]==2) {
      derivs[num_vtx][0]=0; derivs[num_vtx][1]=0; derivs[num_vtx][2]=-1;
    } else {
      derivs[num_vtx][0]=0; derivs[num_vtx][1]=0; derivs[num_vtx][2]=-4;
    }
    num_vtx++;
  }
  
  if(HLSup2[location]!=-1){
  vertices[num_vtx]=HLSup2[location];
  derivs[num_vtx][0]=0; derivs[num_vtx][1]=0; derivs[num_vtx][2]=-1;
  num_vtx++;
  }

  if(HLSdn2[location]!=-1){
  vertices[num_vtx]=HLSdn2[location];
  derivs[num_vtx][0]=0; derivs[num_vtx][1]=0; derivs[num_vtx][2]=1;
  num_vtx++;
  }
  
  if(HLSlt1[location]!=-1){
    vertices[num_vtx]=HLSlt1[location];
    if(HLS1[location]==2) {
      derivs[num_vtx][0]=-1; derivs[num_vtx][1]=0; derivs[num_vtx][2]=0;
    } else {
      derivs[num_vtx][0]=-4; derivs[num_vtx][1]=0; derivs[num_vtx][2]=0;
    }
    num_vtx++;
  }
  
  if(HLSrt1[location]!=-1){
    vertices[num_vtx]=HLSrt1[location];
    if(HLS1[location]==2) {
      derivs[num_vtx][0]=1; derivs[num_vtx][1]=0; derivs[num_vtx][2]=0;
    } else {
      derivs[num_vtx][0]=4; derivs[num_vtx][1]=0; derivs[num_vtx][2]=0;
    }
    num_vtx++;
  }
  
  if(HLSlt2[location]!=-1){
  vertices[num_vtx]=HLSlt2[location];
  derivs[num_vtx][0]=1; derivs[num_vtx][1]=0; derivs[num_vtx][2]=0;
  num_vtx++;
  }

  if(HLSrt2[location]!=-1){
  vertices[num_vtx]=HLSrt2[location];
  derivs[num_vtx][0]=-1; derivs[num_vtx][1]=0; derivs[num_vtx][2]=0;
  num_vtx++;
  }

  if(HLSft1[location]!=-1){
    vertices[num_vtx]=HLSft1[location];
    if(HLS2[location]==2) {
      derivs[num_vtx][0]=0; derivs[num_vtx][1]=-1; derivs[num_vtx][2]=0;
    } else {
      derivs[num_vtx][0]=0; derivs[num_vtx][1]=-4; derivs[num_vtx][2]=0;
    }
    num_vtx++;
  }
  
  if(HLSbk1[location]!=-1){
    vertices[num_vtx]=HLSbk1[location];
    if(HLS2[location]==2) {
      derivs[num_vtx][0]=0; derivs[num_vtx][1]=1; derivs[num_vtx][2]=0;
    } else {
      derivs[num_vtx][0]=0; derivs[num_vtx][1]=4; derivs[num_vtx][2]=0;
    }
    num_vtx++;
  }

  if(HLSft2[location]!=-1){
  vertices[num_vtx]=HLSft2[location];
  derivs[num_vtx][0]=0; derivs[num_vtx][1]=1; derivs[num_vtx][2]=0;
  num_vtx++;
  }

  if(HLSbk2[location]!=-1){
  vertices[num_vtx]=HLSbk2[location];
  derivs[num_vtx][0]=0; derivs[num_vtx][1]=-1; derivs[num_vtx][2]=0;
  num_vtx++;
  }

}


void HexLagrangeShape::ideal( Sample , 
                              MsqMatrix<3,3>& J,
                              MsqError&  ) const
{
  J = MsqMatrix<3,3>(1.0);
}

} // namespace Mesquite
