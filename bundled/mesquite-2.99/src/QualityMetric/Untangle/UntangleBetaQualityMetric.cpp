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
/*! \file   UntangleBetaQualityMetric.cpp
  UntangleBeta is an untangle quality metric which can be used to evaluate
  the quality of two- or three-dimensional elements.

  \author Michael Brewer
  \date   2002-10-10
*/

#include "UntangleBetaQualityMetric.hpp"
#include "Vector3D.hpp"
#include "PatchData.hpp"
#include "MsqError.hpp"

using namespace Mesquite;

static inline void untangle_function_2d( double beta,
                                         const Vector3D temp_vec[],
                                         size_t e_ind,
                                         PatchData &pd, 
                                         double &fval, 
                                         MsqError &err)
{
  Vector3D surface_normal;
  pd.get_domain_normal_at_element(e_ind,surface_normal,err); MSQ_ERRRTN(err);
  Vector3D cross_vec=temp_vec[0]*temp_vec[1];
    //cout<<"\nsurface_normal "<<surface_normal;
    //cout<<"\cross_vec "<<cross_vec;
  double temp_var=cross_vec.length();
  if(cross_vec%surface_normal<0.0){
    temp_var*=-1;
  }
  temp_var -= beta;
    //cout<<"temp_var == "<<temp_var;
  fval=0.0;
  if(temp_var<0.0){
     fval=fabs(temp_var)-temp_var;
  }
    //  cout<<"\nfval == "<<fval<<"  e_ind "<<e_ind;
}

static inline void untangle_function_3d( double beta,
                                         const Vector3D temp_vec[],
                                         double &fval )
{
  double temp_var=temp_vec[0]%(temp_vec[1]*temp_vec[2]);
  temp_var-=beta;
  fval=0.0;
  if(temp_var<0.0){
    fval=fabs(temp_var)-temp_var;
  }
}

/*! \fn UntangleBetaQualityMetric::UntangleBetaQualityMetric(double bet)
  \brief For untangle beta, the constructor defaults to the SUM
  averaging method, and to the ELEMENT_VERTICES evaluation mode.
*/
UntangleBetaQualityMetric::UntangleBetaQualityMetric(double bet)
  : AveragingQM( RMS ),
    mBeta( bet )
  { }

std::string UntangleBetaQualityMetric::get_name() const
  { return "Untangle Beta"; }

int UntangleBetaQualityMetric::get_negate_flag() const
  { return 1; }

bool UntangleBetaQualityMetric::evaluate( PatchData& pd, 
                                          size_t e_ind, 
                                          double& fval, 
                                          MsqError& err )
{
  
  double met_vals[MSQ_MAX_NUM_VERT_PER_ENT];
  fval=MSQ_MAX_CAP;
  const MsqMeshEntity* element = &pd.element_by_index(e_ind);
  const size_t* v_i = element->get_vertex_index_array();
    //only 3 temp_vec will be sent to untangle calculator, but the
    //additional vector3Ds may be needed during the calculations
  Vector3D temp_vec[5];
  const MsqVertex *vertices=pd.get_vertex_array(err);  MSQ_ERRZERO(err);
  EntityTopology type = element->get_element_type();
  switch(type){
    case TRIANGLE:
      temp_vec[0]=vertices[v_i[1]]-vertices[v_i[0]];
      temp_vec[2]=vertices[v_i[2]]-vertices[v_i[0]];
        //make relative to equilateral
      temp_vec[1]=((2*temp_vec[2])-temp_vec[0])*MSQ_SQRT_THREE_INV;
      untangle_function_2d( mBeta, temp_vec, e_ind, pd, fval, err );
      break;
    case QUADRILATERAL:
      temp_vec[0]=vertices[v_i[1]]-vertices[v_i[0]];
      temp_vec[1]=vertices[v_i[3]]-vertices[v_i[0]];
      untangle_function_2d( mBeta, temp_vec, e_ind, pd, met_vals[0], err);  MSQ_ERRZERO(err);
      
      temp_vec[0]=vertices[v_i[2]]-vertices[v_i[1]];
      temp_vec[1]=vertices[v_i[0]]-vertices[v_i[1]];
      untangle_function_2d( mBeta, temp_vec, e_ind, pd, met_vals[1], err);  MSQ_ERRZERO(err);
      
      temp_vec[0]=vertices[v_i[3]]-vertices[v_i[2]];
      temp_vec[1]=vertices[v_i[1]]-vertices[v_i[2]];
      untangle_function_2d( mBeta, temp_vec, e_ind, pd, met_vals[2], err);  MSQ_ERRZERO(err);
      
      temp_vec[0]=vertices[v_i[0]]-vertices[v_i[3]];
      temp_vec[1]=vertices[v_i[2]]-vertices[v_i[3]];
      untangle_function_2d( mBeta, temp_vec, e_ind, pd, met_vals[3], err);  MSQ_ERRZERO(err);
      fval=average_metrics(met_vals, 4, err);  MSQ_ERRZERO(err);
      break;
    case TETRAHEDRON:
      temp_vec[0]=vertices[v_i[1]]-vertices[v_i[0]];
      temp_vec[3]=vertices[v_i[2]]-vertices[v_i[0]];
      temp_vec[4]=vertices[v_i[3]]-vertices[v_i[0]];
        //transform to equilateral tet
      temp_vec[1]=((2*temp_vec[3])-temp_vec[0])/MSQ_SQRT_THREE;
      temp_vec[2]=((3*temp_vec[4])-temp_vec[0]-temp_vec[3])/
        (MSQ_SQRT_THREE*MSQ_SQRT_TWO);
      untangle_function_3d( mBeta, temp_vec,fval );
      break;
    case HEXAHEDRON:
        //transform to v_i[0]
      temp_vec[0]=vertices[v_i[1]]-vertices[v_i[0]];
      temp_vec[1]=vertices[v_i[3]]-vertices[v_i[0]];
      temp_vec[2]=vertices[v_i[4]]-vertices[v_i[0]];
      untangle_function_3d( mBeta, temp_vec, met_vals[0] );
      
      temp_vec[0]=vertices[v_i[2]]-vertices[v_i[1]];
      temp_vec[1]=vertices[v_i[0]]-vertices[v_i[1]];
      temp_vec[2]=vertices[v_i[5]]-vertices[v_i[1]];
      untangle_function_3d( mBeta, temp_vec, met_vals[1] );
      
      temp_vec[0]=vertices[v_i[3]]-vertices[v_i[2]];
      temp_vec[1]=vertices[v_i[1]]-vertices[v_i[2]];
      temp_vec[2]=vertices[v_i[6]]-vertices[v_i[2]];
      untangle_function_3d( mBeta, temp_vec, met_vals[2] );
      
      temp_vec[0]=vertices[v_i[0]]-vertices[v_i[3]];
      temp_vec[1]=vertices[v_i[2]]-vertices[v_i[3]];
      temp_vec[2]=vertices[v_i[7]]-vertices[v_i[3]];
      untangle_function_3d( mBeta, temp_vec, met_vals[3] );
      
      temp_vec[0]=vertices[v_i[7]]-vertices[v_i[4]];
      temp_vec[1]=vertices[v_i[5]]-vertices[v_i[4]];
      temp_vec[2]=vertices[v_i[0]]-vertices[v_i[4]];
      untangle_function_3d( mBeta, temp_vec,met_vals[4] );
      
      temp_vec[0]=vertices[v_i[4]]-vertices[v_i[5]];
      temp_vec[1]=vertices[v_i[6]]-vertices[v_i[5]];
      temp_vec[2]=vertices[v_i[1]]-vertices[v_i[5]];
      untangle_function_3d( mBeta, temp_vec, met_vals[5] );
      
      temp_vec[0]=vertices[v_i[5]]-vertices[v_i[6]];
      temp_vec[1]=vertices[v_i[7]]-vertices[v_i[6]];
      temp_vec[2]=vertices[v_i[2]]-vertices[v_i[6]];
      untangle_function_3d( mBeta, temp_vec, met_vals[6] );
      
      temp_vec[0]=vertices[v_i[6]]-vertices[v_i[7]];
      temp_vec[1]=vertices[v_i[4]]-vertices[v_i[7]];
      temp_vec[2]=vertices[v_i[3]]-vertices[v_i[7]];
      untangle_function_3d( mBeta, temp_vec, met_vals[7] );
      fval=average_metrics(met_vals, 8, err);  MSQ_ERRZERO(err);
      break;
    case PYRAMID:
      for (unsigned i = 0; i < 4; ++i) {
        temp_vec[0] = vertices[v_i[(i+1)%4]] - vertices[v_i[i]];
        temp_vec[1] = vertices[v_i[(i+3)%4]] - vertices[v_i[i]];
        temp_vec[3] = vertices[v_i[4]] - vertices[v_i[i]];
        temp_vec[2] = (4*temp_vec[3] - temp_vec[0] - temp_vec[1])/
                      (2.0 - MSQ_SQRT_TWO);
        untangle_function_3d( mBeta, temp_vec, met_vals[i] );
      }
      fval = average_metrics( met_vals, 4, err ); MSQ_ERRZERO(err);
      break;
    case PRISM:
      temp_vec[0] = vertices[v_i[1]] - vertices[v_i[0]];
      temp_vec[3] = vertices[v_i[2]] - vertices[v_i[0]];
      temp_vec[1] =((2*temp_vec[3])-temp_vec[0])*MSQ_SQRT_THREE_INV;
      temp_vec[2] = vertices[v_i[3]] - vertices[v_i[0]];
      untangle_function_3d( mBeta, temp_vec, met_vals[0] );
    
      temp_vec[0] = vertices[v_i[2]] - vertices[v_i[1]];
      temp_vec[3] = vertices[v_i[0]] - vertices[v_i[1]];
      temp_vec[1] =((2*temp_vec[3])-temp_vec[0])*MSQ_SQRT_THREE_INV;
      temp_vec[2] = vertices[v_i[4]] - vertices[v_i[1]];
      untangle_function_3d( mBeta, temp_vec, met_vals[1] );
    
      temp_vec[0] = vertices[v_i[0]] - vertices[v_i[2]];
      temp_vec[3] = vertices[v_i[1]] - vertices[v_i[2]];
      temp_vec[1] =((2*temp_vec[3])-temp_vec[0])*MSQ_SQRT_THREE_INV;
      temp_vec[2] = vertices[v_i[5]] - vertices[v_i[2]];
      untangle_function_3d( mBeta, temp_vec, met_vals[2] );
    
      temp_vec[0] = vertices[v_i[5]] - vertices[v_i[3]];
      temp_vec[3] = vertices[v_i[4]] - vertices[v_i[3]];
      temp_vec[1] =((2*temp_vec[3])-temp_vec[0])*MSQ_SQRT_THREE_INV;
      temp_vec[2] = vertices[v_i[0]] - vertices[v_i[3]];
      untangle_function_3d( mBeta, temp_vec, met_vals[3] );
    
      temp_vec[0] = vertices[v_i[3]] - vertices[v_i[4]];
      temp_vec[3] = vertices[v_i[5]] - vertices[v_i[4]];
      temp_vec[1] =((2*temp_vec[3])-temp_vec[0])*MSQ_SQRT_THREE_INV;
      temp_vec[2] = vertices[v_i[1]] - vertices[v_i[4]];
      untangle_function_3d( mBeta, temp_vec, met_vals[4] );
    
      temp_vec[0] = vertices[v_i[4]] - vertices[v_i[5]];
      temp_vec[3] = vertices[v_i[3]] - vertices[v_i[5]];
      temp_vec[1] =((2*temp_vec[3])-temp_vec[0])*MSQ_SQRT_THREE_INV;
      temp_vec[2] = vertices[v_i[2]] - vertices[v_i[5]];
      untangle_function_3d( mBeta, temp_vec, met_vals[5] );
    
      fval=average_metrics(met_vals, 6, err);  MSQ_ERRZERO(err);
      break;
    default:
      MSQ_SETERR(err)( MsqError::UNSUPPORTED_ELEMENT,
                       "Unsupported cell type (%ld) for Untangle quality metric.",
                       type);


      fval=MSQ_MAX_CAP;
      return false;
  }// end switch over element type
  return true;
}





     
