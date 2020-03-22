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

/*! \file PowerQualityMetric.hpp
\brief
Header file for the Mesquite::PowerQualityMetric class

  \author Michael Brewer
  \date   2002-09-05
 */


#ifndef PowerQualityMetric_hpp
#define PowerQualityMetric_hpp

#include "Mesquite.hpp"
#include "QualityMetric.hpp"
#include "Exponent.hpp"

namespace MESQUITE_NS
{
  class MsqError;
  class Vector3D;

     /*! \class PowerQualityMetric
       \brief Raises a single quality metrics (qMetric1) to an arbitrary
       power (a double value, scaleAlpha) for two- and three-diminsional
       elements.  
     */
   class PowerQualityMetric : public QualityMetric
   {
  public:
       /*! Ensures that qm1 is not NULL.  If qm1 is only valid
         on a certain feasible, then the composite metric has the same
         constraint.  The composite metric also has the same negate flag
         as qm1.
       */
     PowerQualityMetric(QualityMetric* qm1, double pow_double );
     
     
       // virtual destructor ensures use of polymorphism during destruction
     virtual ~PowerQualityMetric();
     
    MetricType get_metric_type() const
      { return mMetric.get_metric_type(); }

     std::string get_name() const;

     int get_negate_flag() const;

     QualityMetric* get_metric() const { return &mMetric; }

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
  
    QualityMetric& mMetric;
    Mesquite::Exponent mPower;
   };
   

} //namespace


#endif // PowerQualityMetric_hpp







