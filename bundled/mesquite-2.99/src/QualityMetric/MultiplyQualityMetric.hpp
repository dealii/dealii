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

/*! \file MultiplyQualityMetric.hpp
\brief
Header file for the Mesquite::MultiplyQualityMetric class

  \author Michael Brewer
  \date   2002-09-05
 */


#ifndef MultiplyQualityMetric_hpp
#define MultiplyQualityMetric_hpp

#include "Mesquite.hpp"
#include "QualityMetric.hpp"

namespace MESQUITE_NS
{
   class Vector3D;
   class MsqError;
   class PatchData;
   class MsqMeshEntity;
   class MsqVertex;
   
     /*! \class MultiplyQualityMetric
       \brief Combines two quality metrics (qMetric1 and qMetric2 defined
       in the parent class CompositeQualityMetric) by multiplication for two-
       and three-diminsional elements.  Note:  This function should not
       be used to combine a node-based metric with an element-based
       metric.  
     */
   class MultiplyQualityMetric : public QualityMetric
   {
  public:
       /*! Ensures that qm1 and qm2 are not NULL.  If either qm1 or qm2
         are valid only on a feasible region, then the composite
         metric's feasibility flag is set to one.  If qm1 and qm2 have
         different negateFlags, then a warning is printed, and the composite
         metric's negate flag is set to one.  Otherwise, the composite
         metric's negateFlag is set to qm1's negateFlag (and, thus, qm2's
         negateFlag).  
       */
     MultiplyQualityMetric(QualityMetric* qm1, QualityMetric* qm2,
                           MsqError &err);
     
       // virtual destructor ensures use of polymorphism during destruction
     virtual ~MultiplyQualityMetric();
     
     MetricType get_metric_type() const;

     std::string get_name() const;

     int get_negate_flag() const;

     QualityMetric* get_first_metric() const { return &metric1; }
     QualityMetric* get_second_metric() const { return &metric2; }

     virtual
     void get_evaluations( PatchData& pd, 
                           std::vector<size_t>& handles, 
                           bool free_vertices_only,
                           MsqError& err );

     virtual
     bool evaluate( PatchData& pd, 
                    size_t handle, 
                    double& value, 
                    MsqError& err );


     virtual
     bool evaluate_with_indices( PatchData& pd,
                    size_t handle,
                    double& value,
                    std::vector<size_t>& indices,
                    MsqError& err );

     virtual
     bool evaluate_with_gradient( PatchData& pd,
                    size_t handle,
                    double& value,
                    std::vector<size_t>& indices,
                    std::vector<Vector3D>& gradient,
                    MsqError& err );


     virtual
     bool evaluate_with_Hessian( PatchData& pd,
                    size_t handle,
                    double& value,
                    std::vector<size_t>& indices,
                    std::vector<Vector3D>& gradient,
                    std::vector<Matrix3D>& Hessian,
                    MsqError& err );
     
  private:

     QualityMetric& metric1;
     QualityMetric& metric2;
     mutable std::vector<size_t> mHandles;
     mutable std::vector<size_t> indices1, indices2;
     mutable std::vector<Vector3D> grad1, grad2;
     mutable std::vector<Matrix3D> Hess1, Hess2;
  };
   

} //namespace


#endif // MultiplyQualityMetric_hpp







