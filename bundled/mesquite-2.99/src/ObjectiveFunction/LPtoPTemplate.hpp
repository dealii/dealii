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

/*! \file LPtoPTemplate.hpp
  \brief Header file for the Mesquite::LPtoPTemplate class
 \author Michael Brewer
 \author Thomas Leurent
  \date   2002-05-23
 */


#ifndef LPtoPTemplate_hpp
#define LPtoPTemplate_hpp

#include "Mesquite.hpp"
#include "ObjectiveFunctionTemplate.hpp"

namespace MESQUITE_NS
{
  class Matrix3D;

     /*! \class LPtoPTemplate
       \brief Calculates the L_p objective function raised to the pth
       power.  That is, sums the p_th powers of (the absolute value of)
       the quality metric values.

       \todo MB. Suggestions made by Todd Munson:
       a)  There is an inconsistent use of fabs.  The hessian evaluation
       when using the one norm does not take the absolute value, while the
       gradient does.
       b)  The analytic gradient and hessian evaluations are incorrect when
       the quality metric changes sign due to taking the absolute value.
       The negative of the element gradient and hessian also needs to be
       taken.
       c)  Done.  The analytic gradient and hessian evaluations are
       incorrect when the negate flag is set to -1.  The negative
       of the element gradient and hessian also needs to be taken
       in this case.
       d)  The malloc in the concrete_eval routine should be removed.

     */
  class LPtoPTemplate :public ObjectiveFunctionTemplate
  {
  public:
	MESQUITE_EXPORT
    LPtoPTemplate(QualityMetric *, short, MsqError &);
	MESQUITE_EXPORT
    LPtoPTemplate( short, QualityMetric* );
     
	MESQUITE_EXPORT
    virtual ~LPtoPTemplate();
    
	MESQUITE_EXPORT
    virtual void clear();
    
	MESQUITE_EXPORT
    virtual bool evaluate( EvalType type, 
                           PatchData& pd,
                           double& value_out,
                           bool free,
                           MsqError& err ); 
     
	MESQUITE_EXPORT
    virtual bool evaluate_with_gradient( EvalType type, 
                                         PatchData& pd,
                                         double& value_out,
                                         std::vector<Vector3D>& grad_out,
                                         MsqError& err ); 

	MESQUITE_EXPORT
    virtual bool evaluate_with_Hessian_diagonal( EvalType type, 
                                        PatchData& pd,
                                        double& value_out,
                                        std::vector<Vector3D>& grad_out,
                                        std::vector<SymMatrix3D>& hess_diag_out,
                                        MsqError& err ); 

	MESQUITE_EXPORT
    virtual bool evaluate_with_Hessian( EvalType type, 
                                        PatchData& pd,
                                        double& value_out,
                                        std::vector<Vector3D>& grad_out,
                                        MsqHessian& Hessian_out,
                                        MsqError& err ); 

	MESQUITE_EXPORT
    virtual ObjectiveFunction* clone() const;

       /*!Use set_dividing_by_n to control whether this objective
         function divides it's final value by the number of
         metric values used to compute the objective function
         value.  That is, if the associated metric is element
         based, the obejctive function value is divided by
         the number of elements.  If it is vertex based, the
         objective function is divided by the number of vertices.
         If this function is passed 'true', the function value
         will be scale.  If it is passed false, the function
         value will not be scaled.*/
	MESQUITE_EXPORT
    void set_dividing_by_n(bool d_bool){dividingByN=d_bool;}
     
   private:
     double get_value( double power_sum, size_t count, EvalType type, 
		       size_t& global_count, MsqError& err );
     
       //! The metric value entries are raised to the pVal power
     short pVal;
       //! dividingByN is true if we are dividing the objective function
       //! by the number of metric values.
     bool dividingByN;
     
     
     size_t mCount;    /**< The number of accumulated entires */
     double mPowSum;   /**< The accumulated sum of values */
     size_t saveCount; /**< Saved count from previous patch */
     double savePowSum;/**< Saved sum from previous patch */
     
     /** Temporary storage for qm sample handles */
     mutable std::vector<size_t> qmHandles;
     /** Temporary storage for qm vertex indices */
     mutable std::vector<size_t> mIndices;
     /** Temporary storage for qm gradient */
     mutable std::vector<Vector3D> mGradient;
     /** Temporary storage for qm Hessian diagonal blocks */
     mutable std::vector<SymMatrix3D> mDiag;
      /** Temporary storage for qm Hessian */
     mutable std::vector<Matrix3D> mHessian;
   };
   
}//namespace

#endif // LPtoPTemplate_hpp

