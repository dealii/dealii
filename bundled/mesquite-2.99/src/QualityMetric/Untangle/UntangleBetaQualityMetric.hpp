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

/*! \file UntangleBetaQualityMetric.hpp

Header file for the Mesquite::UntangleBetaQualityMetric class

  \author Michael Brewer
  \date   2002-09-10
 */

#ifndef UNTANGLE_BETA_QUALITY_METRIC_HPP
#define UNTANGLE_BETA_QUALITY_METRIC_HPP

#include "Mesquite.hpp"
#include "ElementQM.hpp"
#include "AveragingQM.hpp"
namespace MESQUITE_NS
{
     /*! \class UntangleBetaQualityMetric
       \brief The untangle beta quality metric.
       
       Given a scalar value beta and local signed element volume alpha_i,
       define delta_i to be alpha_i minus beta.  The Untangle beta value
       is then defined as square root of the sum over sample points
       of the absolute value of delta_i minus delta_i, difference squared.
       That is, the root mean square of the difference, abs(delta_i) minus
       delta_i.

       The constructor defaults to RMS AveragingMethod and
       ELEMENT_VERTICES evaluationMode.  The default beta value is
       .05.
     */
   
   class UntangleBetaQualityMetric : public ElementQM, public AveragingQM
   {
   public:
     
     MESQUITE_EXPORT UntangleBetaQualityMetric(double bet=0.05);

       // virtual destructor ensures use of polymorphism during destruction
     MESQUITE_EXPORT virtual ~UntangleBetaQualityMetric()
        {}
     
       /*!Function to allow users to set the beta value after the
         metric has already been created. */
     MESQUITE_EXPORT void set_beta(double beta_in)
       {mBeta = beta_in;}
       /*!Function to allow the user to check the value of beta.*/
     MESQUITE_EXPORT double get_beta()
       {return mBeta;}
       
     MESQUITE_EXPORT virtual
     std::string get_name() const;
     
     MESQUITE_EXPORT virtual
     int get_negate_flag() const;
     
     MESQUITE_EXPORT virtual
     bool evaluate( PatchData& pd, 
                    size_t handle, 
                    double& value, 
                    MsqError& err );

   private:
     double mBeta;
   };
   
} //namespace


#endif // UntangleBetaQualityMetric_hpp


