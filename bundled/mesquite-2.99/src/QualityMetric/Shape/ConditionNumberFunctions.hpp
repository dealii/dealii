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
// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file ConditionNumberFunctions.hpp

Header file for the Mesquite::ShapeQualityMetric class

  \author Thomas Leurent
  \date   2002-09-01
 */


#ifndef ConditionNumberFunctions_hpp
#define ConditionNumberFunctions_hpp

#include "Mesquite.hpp"
#include "MsqError.hpp"
#include "QualityMetric.hpp"
#include "PatchData.hpp"

namespace MESQUITE_NS
{
   static inline bool condition_number_2d( const Vector3D temp_vec[],
                                           size_t e_ind,
                                           PatchData &pd,
                                           double &fval,
                                           MsqError &err)
   {   
       //norm squared of J
     double term1=temp_vec[0]%temp_vec[0]+temp_vec[1]%temp_vec[1];

     Vector3D unit_surf_norm;
     pd.get_domain_normal_at_element(e_ind,unit_surf_norm,err);MSQ_ERRZERO(err);
     unit_surf_norm.normalize();

     // det J
     double temp_var=unit_surf_norm%(temp_vec[0]*temp_vec[1]);
     
     // revert to old, non-barrier form
     if (temp_var <= 0.0)
       return false;
     fval = term1 / (2.0*temp_var);
     return true;
     
     /*
     double h;
     double delta=pd.get_barrier_delta(err); MSQ_ERRZERO(err);

     // Note: technically, we want delta=eta*tau-max
     //       whereas the function above gives delta=eta*alpha-max
     //      
     //       Because the only requirement on eta is eta << 1,
     //       and because tau-max = alpha-max/0.707 we can
     //       ignore the discrepancy

     if (delta==0) { 
        if (temp_var < MSQ_DBL_MIN ) {
           return false;
        }
        else {
           h=temp_var;
        }

     // Note: when delta=0, the vertex_barrier_function
     //       formally gives h=temp_var as well.
     //       We just do it this way to avoid any 
     //       roundoff issues.
     // Also: when delta=0, this metric is identical
     //       to the original condition number with
     //       the barrier at temp_var=0
     }
     else {
        h = QualityMetric::vertex_barrier_function(temp_var,delta);

        if (h<MSQ_DBL_MIN && fabs(temp_var) > MSQ_DBL_MIN ) { 
          h = delta*delta/fabs(temp_var); }
        // Note: Analytically, h is strictly positive, but
        //       it can be zero numerically if temp_var
        //       is a large negative number 
        //       In the case h=0, we use a different analytic
        //       approximation to compute h.
     }

     if (h<MSQ_DBL_MIN) {
       MSQ_SETERR(err)( "Barrier function is zero due to excessively large "
                        "negative area compared to delta. /n Try to untangle "
                        "mesh another way. ", MsqError::INVALID_MESH);
       return false;
     }

     fval=term1/(2*h);

     if (fval>MSQ_MAX_CAP) {
        fval=MSQ_MAX_CAP;
     }
     return true;
     */
   }
   

   //} //namespace

   static inline bool condition_number_3d( const Vector3D temp_vec[],
                                           PatchData &pd,
                                           double &fval,
                                           MsqError &err)
   {   
       //norm squared of J
     double term1=temp_vec[0]%temp_vec[0]+
        temp_vec[1]%temp_vec[1]+
        temp_vec[2]%temp_vec[2];
       //norm squared of adjoint of J
     double term2=(temp_vec[0]*temp_vec[1])%
        (temp_vec[0]*temp_vec[1])+
        (temp_vec[1]*temp_vec[2])%
        (temp_vec[1]*temp_vec[2])+
        (temp_vec[2]*temp_vec[0])%
        (temp_vec[2]*temp_vec[0]);
       //det of J
     double temp_var=temp_vec[0]%(temp_vec[1]*temp_vec[2]);
     
      // revert to old, non-barrier formulation
     if (temp_var <= 0.0)
       return false;
     fval = sqrt( term1 * term2 ) / (3.0*temp_var);
     return true;
     
     /*
     double h;
     double delta=pd.get_barrier_delta(err); MSQ_ERRZERO(err);

     // Note: technically, we want delta=eta*tau-max
     //       whereas the function above gives delta=eta*alpha-max
     //      
     //       Because the only requirement on eta is eta << 1,
     //       and because tau-max = alpha-max/0.707 we can
     //       ignore the discrepancy

     if (delta==0) { 
        if (temp_var < MSQ_DBL_MIN ) {
           return false;
        }
        else {
           h=temp_var;
        }

     // Note: when delta=0, the vertex_barrier_function
     //       formally gives h=temp_var as well.
     //       We just do it this way to avoid any 
     //       roundoff issues.
     // Also: when delta=0, this metric is identical
     //       to the original condition number with
     //       the barrier at temp_var=0

     }
     else {
        h = QualityMetric::vertex_barrier_function(temp_var,delta);
 
        if (h<MSQ_DBL_MIN && fabs(temp_var) > MSQ_DBL_MIN ) { 
          h = delta*delta/fabs(temp_var); }
 
        // Note: Analytically, h is strictly positive, but
        //       it can be zero numerically if temp_var
        //       is a large negative number 
        //       In the h=0, we use a different analytic
        //       approximation to compute h.
     }
     if (h<MSQ_DBL_MIN) {
       MSQ_SETERR(err)("Barrier function is zero due to excessively large "
                       "negative area compared to delta. /n Try to untangle "
                       "mesh another way. ", MsqError::INVALID_MESH);
       return false;
     }

     fval=sqrt(term1*term2)/(3*h);

     if (fval>MSQ_MAX_CAP) {
        fval=MSQ_MAX_CAP;
     }
     return true;
     */
   }
   

} //namespace


#endif // ConditionNumberFunctions_hpp
