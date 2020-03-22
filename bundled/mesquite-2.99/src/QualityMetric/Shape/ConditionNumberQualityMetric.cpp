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
  \file   ConditionNumberQualityMetric.cpp
  \brief  

  \author Michael Brewer
  \date   2002-06-9
*/
#include <vector>
#include "ConditionNumberQualityMetric.hpp"
#include <math.h>
#include "Vector3D.hpp"
#include "ConditionNumberFunctions.hpp"

using namespace Mesquite;

ConditionNumberQualityMetric::ConditionNumberQualityMetric()
  : AveragingQM(QualityMetric::LINEAR)
  { }

std::string ConditionNumberQualityMetric::get_name() const
  { return "Condition Number"; }

int ConditionNumberQualityMetric::get_negate_flag() const
  { return 1; }

bool ConditionNumberQualityMetric::evaluate( PatchData& pd, 
                                     size_t handle, 
                                     double& fval, 
                                     MsqError& err )
{
  const MsqMeshEntity *const element = &pd.element_by_index(handle);
  bool return_flag;
  double met_vals[MSQ_MAX_NUM_VERT_PER_ENT];
  fval=MSQ_MAX_CAP;
  const size_t* v_i = element->get_vertex_index_array();
    //only 3 temp_vec will be sent to cond-num calculator, but the
    //additional vector3Ds may be needed during the calculations
  Vector3D temp_vec[6];
  const MsqVertex *vertices=pd.get_vertex_array(err);
  EntityTopology type = element->get_element_type();
  switch(type){
    case TRIANGLE:
      temp_vec[0]=vertices[v_i[1]]-vertices[v_i[0]];
      temp_vec[2]=vertices[v_i[2]]-vertices[v_i[0]];
        //make relative to equilateral
      temp_vec[1]=((2*temp_vec[2])-temp_vec[0])*MSQ_SQRT_THREE_INV;
      return_flag=condition_number_2d(temp_vec,handle,pd,fval,err); MSQ_ERRZERO(err);
      return return_flag;
    case QUADRILATERAL:
      temp_vec[0]=vertices[v_i[1]]-vertices[v_i[0]];
      temp_vec[1]=vertices[v_i[3]]-vertices[v_i[0]];
      return_flag=condition_number_2d(temp_vec,handle,pd,met_vals[0],err);  MSQ_ERRZERO(err);
      if(!return_flag)
        return return_flag;
      temp_vec[0]=vertices[v_i[2]]-vertices[v_i[1]];
      temp_vec[1]=vertices[v_i[0]]-vertices[v_i[1]];
      return_flag=condition_number_2d(temp_vec,handle,pd,met_vals[1],err);  MSQ_ERRZERO(err);
      if(!return_flag)
        return return_flag;
      temp_vec[0]=vertices[v_i[3]]-vertices[v_i[2]];
      temp_vec[1]=vertices[v_i[1]]-vertices[v_i[2]];
      return_flag=condition_number_2d(temp_vec,handle,pd,met_vals[2],err);  MSQ_ERRZERO(err);
      if(!return_flag)
        return return_flag;
      temp_vec[0]=vertices[v_i[0]]-vertices[v_i[3]];
      temp_vec[1]=vertices[v_i[2]]-vertices[v_i[3]];
      return_flag=condition_number_2d(temp_vec,handle,pd,met_vals[3],err);  MSQ_ERRZERO(err);
      fval = average_metrics(met_vals, 4, err);
      return return_flag;
    case TETRAHEDRON:
      temp_vec[0]=vertices[v_i[1]]-vertices[v_i[0]];
      temp_vec[3]=vertices[v_i[2]]-vertices[v_i[0]];
      temp_vec[4]=vertices[v_i[3]]-vertices[v_i[0]];
        //transform to equilateral tet
      temp_vec[1]=((2*temp_vec[3])-temp_vec[0])/MSQ_SQRT_THREE;
      temp_vec[2]=((3*temp_vec[4])-temp_vec[0]-temp_vec[3])/
        (MSQ_SQRT_THREE*MSQ_SQRT_TWO);
      return_flag=condition_number_3d(temp_vec,pd,fval,err);  MSQ_ERRZERO(err);
      return return_flag;

    case HEXAHEDRON:
        //transform to v_i[0]
      temp_vec[0]=vertices[v_i[1]]-vertices[v_i[0]];
      temp_vec[1]=vertices[v_i[3]]-vertices[v_i[0]];
      temp_vec[2]=vertices[v_i[4]]-vertices[v_i[0]];
      return_flag=condition_number_3d(temp_vec,pd,met_vals[0],err);  MSQ_ERRZERO(err);
      if(!return_flag)
        return return_flag;
      temp_vec[0]=vertices[v_i[2]]-vertices[v_i[1]];
      temp_vec[1]=vertices[v_i[0]]-vertices[v_i[1]];
      temp_vec[2]=vertices[v_i[5]]-vertices[v_i[1]];
      return_flag=condition_number_3d(temp_vec,pd,met_vals[1],err);  MSQ_ERRZERO(err);
      if(!return_flag)
        return return_flag;
      temp_vec[0]=vertices[v_i[3]]-vertices[v_i[2]];
      temp_vec[1]=vertices[v_i[1]]-vertices[v_i[2]];
      temp_vec[2]=vertices[v_i[6]]-vertices[v_i[2]];
      return_flag=condition_number_3d(temp_vec,pd,met_vals[2],err);  MSQ_ERRZERO(err);
      if(!return_flag)
        return return_flag;
      temp_vec[0]=vertices[v_i[0]]-vertices[v_i[3]];
      temp_vec[1]=vertices[v_i[2]]-vertices[v_i[3]];
      temp_vec[2]=vertices[v_i[7]]-vertices[v_i[3]];
      return_flag=condition_number_3d(temp_vec,pd,met_vals[3],err);  MSQ_ERRZERO(err);
      if(!return_flag)
        return return_flag;
      temp_vec[0]=vertices[v_i[7]]-vertices[v_i[4]];
      temp_vec[1]=vertices[v_i[5]]-vertices[v_i[4]];
      temp_vec[2]=vertices[v_i[0]]-vertices[v_i[4]];
      return_flag=condition_number_3d(temp_vec,pd,met_vals[4],err);  MSQ_ERRZERO(err);
      if(!return_flag)
        return return_flag;
      temp_vec[0]=vertices[v_i[4]]-vertices[v_i[5]];
      temp_vec[1]=vertices[v_i[6]]-vertices[v_i[5]];
      temp_vec[2]=vertices[v_i[1]]-vertices[v_i[5]];
      return_flag=condition_number_3d(temp_vec,pd,met_vals[5],err);  MSQ_ERRZERO(err);
      if(!return_flag)
        return return_flag;
      temp_vec[0]=vertices[v_i[5]]-vertices[v_i[6]];
      temp_vec[1]=vertices[v_i[7]]-vertices[v_i[6]];
      temp_vec[2]=vertices[v_i[2]]-vertices[v_i[6]];
      return_flag=condition_number_3d(temp_vec,pd,met_vals[6],err);  MSQ_ERRZERO(err);
      if(!return_flag)
        return return_flag;
      temp_vec[0]=vertices[v_i[6]]-vertices[v_i[7]];
      temp_vec[1]=vertices[v_i[4]]-vertices[v_i[7]];
      temp_vec[2]=vertices[v_i[3]]-vertices[v_i[7]];
      return_flag=condition_number_3d(temp_vec,pd,met_vals[7],err);  MSQ_ERRZERO(err);
      fval=average_metrics(met_vals, 8, err);  MSQ_ERRZERO(err);
      return return_flag;

    case PYRAMID:
    {
      unsigned num_adj;
      const unsigned* adj_idx;
      return_flag = true;
      for (size_t j = 0; j < 4; ++j) // skip fifth vertex (apex)
      {
        adj_idx = TopologyInfo::adjacent_vertices( type, j, num_adj );
        assert( num_adj == 3 );
        
        temp_vec[0] = vertices[v_i[adj_idx[0]]] - vertices[v_i[j]];
        temp_vec[1] = vertices[v_i[adj_idx[1]]] - vertices[v_i[j]];
          // calculate last vect map to right tetrahedron
        temp_vec[3] = vertices[v_i[adj_idx[2]]] - vertices[v_i[adj_idx[0]]];
        temp_vec[4] = vertices[v_i[adj_idx[2]]] - vertices[v_i[adj_idx[1]]];
        temp_vec[2] = 0.5 * (temp_vec[3] + temp_vec[4]);
        
        return_flag = return_flag && condition_number_3d( temp_vec, pd, met_vals[j], err );
      }
      fval = average_metrics( met_vals, 4, err ); MSQ_ERRZERO(err);
      return return_flag;
    }
    
    case PRISM:
    {
      unsigned num_adj;
      const unsigned* adj_idx;
      return_flag = true;
      for (size_t j = 0; j < 6; ++j) 
      {
        adj_idx = TopologyInfo::adjacent_vertices( type, j, num_adj );
        assert( num_adj == 3 );
        
        temp_vec[0] = vertices[v_i[adj_idx[0]]] - vertices[v_i[j]];
        temp_vec[1] = vertices[v_i[adj_idx[1]]] - vertices[v_i[j]];
        temp_vec[2] = vertices[v_i[adj_idx[2]]] - vertices[v_i[j]];
          // map to right tetrahedron
        temp_vec[1] += vertices[v_i[adj_idx[1]]];
        temp_vec[1] -= vertices[v_i[adj_idx[0]]];
        temp_vec[1] *= MSQ_SQRT_THREE_INV;
        
        return_flag = return_flag && condition_number_3d( temp_vec, pd, met_vals[j], err );
      }
      fval = average_metrics( met_vals, 6, err ); MSQ_ERRZERO(err);
      return return_flag;
    }
        
    default:
       MSQ_SETERR(err)( MsqError::UNSUPPORTED_ELEMENT,
                      "Unsupported cell type (%ld) for Condition Number quality metric.",
                       type);

     fval=MSQ_MAX_CAP;
     return false;
  }// end switch over element type
  return false;
}


