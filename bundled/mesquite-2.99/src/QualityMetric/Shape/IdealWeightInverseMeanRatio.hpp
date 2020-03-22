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

/*! \file IdealWeightInverseMeanRatio.hpp

Header file for the Mesquite::IdealWeightInverseMeanRatio class

\author Michael Brewer
\author Thomas Leurent
\date   2002-06-19
 */


#ifndef IdealWeightInverseMeanRatio_hpp
#define IdealWeightInverseMeanRatio_hpp

#include "Mesquite.hpp"
#include "ElementQM.hpp"
#include "AveragingQM.hpp"
#include "Vector3D.hpp"
#include "Matrix3D.hpp"
#include "Exponent.hpp"

namespace MESQUITE_NS
{
    class MsqMeshEntity;
    class PatchData;
    class MsqError;

   /*! \class IdealWeightInverseMeanRatio
     \brief Computes the inverse mean ratio of given element.

     The metric does not use the sample point functionality or the
     compute_weighted_jacobian.  It evaluates the metric at
     the element vertices, and uses the isotropic ideal element.
     Optionally, the metric computation can be raised to the
     'pow_dbl' power.  This does not necessarily raise the metric
     value to the 'pow_dbl' power but instead raises each local
     metric.  For example, if the corner inverse mean ratios of a quadraliteral
     element were m1,m2,m3, and m4 and we set pow_dbl=2 and
     used linear averaging, the metric value would then be
     m = .25(m1*m1 + m2*m2 + m3*m3 + m4*m4).  The metric does
     require a feasible region, and the metric needs to be minimized
     if pow_dbl is greater than zero and maximized if pow_dbl
     is less than zero.  pow_dbl being equal to zero is invalid.
   */
   class IdealWeightInverseMeanRatio : public ElementQM, public AveragingQM
   {
   public:
      MESQUITE_EXPORT IdealWeightInverseMeanRatio(MsqError& err, double power = 1.0);
      MESQUITE_EXPORT IdealWeightInverseMeanRatio();

      //! virtual destructor ensures use of polymorphism during destruction
      MESQUITE_EXPORT virtual ~IdealWeightInverseMeanRatio() {
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
       //! Sets the power value in the metric computation.
     void set_metric_power(double pow_dbl, MsqError& err);
     
      // arrays used in Hessian computations 
      // We allocate them here, so that one allocation only is done.
      // This gives a big computation speed increase.
      Vector3D mCoords[4]; // Vertex coordinates for the (decomposed) elements
      Vector3D mGradients[32]; // Gradient of metric with respect to the coords
      Matrix3D mHessians[80]; // Hessian of metric with respect to the coords
      double   mMetrics[8]; // Metric values for the (decomposed) elements
       //variables used in the definition of the metric (2d and 3d)
     double a2Con;
     Exponent b2Con;
     Exponent c2Con;
     
     double a3Con;
     Exponent b3Con;
     Exponent c3Con;
   };
} //namespace


#endif // IdealWeightInverseMeanRatio_hpp
