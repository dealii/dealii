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
  \file   AspectRatioGammaQualityMetric.cpp
   
  \brief Evaluates the Aspect Ratio Gamma metric for two- and
  three-diminsional simplicial elements.
  \author Michael Brewer
  \date   2002-05-09
*/

#include "AspectRatioGammaQualityMetric.hpp"
#include <math.h>
#include "Vector3D.hpp"
#include "MsqMeshEntity.hpp"
#include "PatchData.hpp"

#include <vector>
using std::vector;

using namespace Mesquite;

const double fourDivRootThree = 4.0/sqrt(3.0);
const double twelveDivRootTwo = 12.0/sqrt(2.0);


std::string AspectRatioGammaQualityMetric::get_name() const
  { return "AspectRatioGamma"; }
  
int AspectRatioGammaQualityMetric::get_negate_flag() const
  { return 1; }

//note that we can define this metric for other element types?
//!Evaluate aspect ratio gamma on ``element''
bool AspectRatioGammaQualityMetric::evaluate(PatchData &pd,
                                             size_t elem_index,
                                             double &fval,
                                             MsqError &err)
{
  MsqMeshEntity& element = pd.element_by_index( elem_index );
  EntityTopology entity = element.get_element_type();
  double vol;
  Vector3D cross, normal(0,0,0);
  fval = MSQ_MAX_CAP;
  
    //get element's nodes
  vert.clear();
  pd.get_element_vertex_coordinates(elem_index, vert, err);  MSQ_ERRZERO(err);
  
  switch(entity)
  {
    case TRIANGLE:
        //area
      cross = (vert[1] - vert[0]) * (vert[2] - vert[0]);
      vol= cross.length() / 2.0;
      if (vol < MSQ_MIN)
        return false;
      
      if (pd.domain_set()) { // need domain to check for inverted elements
        pd.get_domain_normal_at_corner( elem_index, 0, normal, err ); MSQ_ERRZERO(err);
        if ((cross % normal) < -MSQ_MIN)
          return false;
      }
   
        // sum squares of edge lengths
      fval = ((vert[1] - vert[0]).length_squared()
            + (vert[2] - vert[0]).length_squared()
            + (vert[1] - vert[2]).length_squared());
        // average
      fval /= 3.0;
      
         //normalize to equil. and div by area
      fval /= vol * fourDivRootThree;
      
      break;
      
    case TETRAHEDRON:
      vol = (vert[1] - vert[0]) % ((vert[2] - vert[0]) * (vert[3] - vert[0])) / 6.0;
      if (vol < MSQ_MIN)  // zero for degenerate and negative for inverted
        return false;
      

        // sum squares of edge lengths
      fval = (vert[1] - vert[0]).length_squared()
           + (vert[2] - vert[0]).length_squared()
           + (vert[3] - vert[0]).length_squared()
           + (vert[2] - vert[1]).length_squared()
           + (vert[3] - vert[1]).length_squared()
           + (vert[3] - vert[2]).length_squared();
        // average
      fval /= 6.0;
      
      fval *= sqrt(fval);
        //normalize to equil. and div by volume
      fval /= vol * twelveDivRootTwo;

      break;
    default:
      MSQ_SETERR(err)(MsqError::UNSUPPORTED_ELEMENT, 
                     "Entity type %d is not valid for Aspect Ratio Gamma\n", 
                     (int)entity);
      return false;
  };
  
  return true;
}
