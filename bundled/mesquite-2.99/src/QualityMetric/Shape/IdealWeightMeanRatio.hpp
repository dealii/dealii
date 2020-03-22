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

/*! \file IdealWeightMeanRatio.hpp

Header file for the Mesquite::IdealWeightMeanRatio class

\author Michael Brewer
\author Todd Munson
\date   2002-11-11
 */


#ifndef IdealWeightMeanRatio_hpp
#define IdealWeightMeanRatio_hpp

#include "Mesquite.hpp"
#include "MsqError.hpp"
#include "ElementQM.hpp"
#include "AveragingQM.hpp"
#include "Vector3D.hpp"
#include "Matrix3D.hpp"
#include "Exponent.hpp"

namespace MESQUITE_NS
{
   /*! \class IdealWeightMeanRatio
     \brief Computes the mean ratio quality metric
     of given element.
     
     The metric does not use the sample point functionality or the
     compute_weighted_jacobian.  It evaluates the metric at
     the element vertices, and uses the isotropic ideal element.
     It does require a feasible region, and the metric needs
     to be maximized.
   */
   class IdealWeightMeanRatio : public ElementQM, public AveragingQM
   {
   public:
 
     MESQUITE_EXPORT IdealWeightMeanRatio()
      : AveragingQM( QualityMetric::LINEAR ),
        a2Con(2.0), 
        b2Con(-1.0), 
        c2Con(1.0),
        a3Con(3.0),
        b3Con(-1.0),
        c3Con(2.0/3.0)
      { }
     
     //! virtual destructor ensures use of polymorphism during destruction
     MESQUITE_EXPORT virtual ~IdealWeightMeanRatio() {
     }
     
     
     virtual std::string get_name() const;

      //! 1 if metric should be minimized, -1 if metric should be maximized.
     virtual int get_negate_flag() const;

     virtual
     bool evaluate( PatchData& pd, 
                    size_t handle, 
                    double& value, 
                    MsqError& err );
      
     virtual
     bool evaluate_with_gradient( PatchData& pd,
                    size_t handle,
                    double& value,
                    std::vector<size_t>& indices,
                    std::vector<Vector3D>& gradient,
                    MsqError& err );

     virtual
     bool evaluate_with_Hessian_diagonal( PatchData& pd,
                    size_t handle,
                    double& value,
                    std::vector<size_t>& indices,
                    std::vector<Vector3D>& gradient,
                    std::vector<SymMatrix3D>& Hessian,
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
      // arrays used in Hessian computations 
      // We allocate them here, so that one allocation only is done.
      // This gives a big computation speed increase.
      Vector3D mCoords[4]; // Vertex coordinates for the (decomposed) elements
      Vector3D mGradients[32]; // Gradient of metric with respect to the coords
      Matrix3D mHessians[80]; // Hessian of metric with respect to the coords
      double   mMetrics[8]; // Metric values for the (decomposed) elements
      
      const double a2Con;
      const Exponent b2Con;
      const Exponent c2Con;
      
      const double a3Con;
      const Exponent b3Con;
      const Exponent c3Con;
   };
} //namespace


#endif // IdealWeightMeanRatio_hpp


