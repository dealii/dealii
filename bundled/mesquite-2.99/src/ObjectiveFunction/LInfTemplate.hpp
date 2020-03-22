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

/*! \file LInfTemplate.hpp

Header file for the Mesquite::LInfTemplate class

  \author Michael Brewer
  \date   2002-07-3
 */


#ifndef LInfTemplate_hpp
#define LInfTemplate_hpp

#include "Mesquite.hpp"
#include "ObjectiveFunctionTemplate.hpp"

namespace MESQUITE_NS
{
  /*! \class LInfTemplate
    \brief Computes the L_infinity objective function for a given patch,
    i.e., LInfTemplate::concrete_evaluate returns the maximum absolute value of
    the quality metric values  on 'patch'.
  */
   class LInfTemplate :public ObjectiveFunctionTemplate
   {
   public:
     MESQUITE_EXPORT LInfTemplate(QualityMetric *);
     MESQUITE_EXPORT virtual ~LInfTemplate();
     MESQUITE_EXPORT virtual bool evaluate( EvalType type, 
                                            PatchData& pd,
                                            double& value_out,
                                            bool free,
                                            MsqError& err ); 
     MESQUITE_EXPORT virtual ObjectiveFunction* clone() const;
     MESQUITE_EXPORT virtual void clear();
   private:
     /** Temporary storage for qm sample handles */
     mutable std::vector<size_t> qmHandles;
   };
}//namespace
#endif // LInfTemplate_hpp
