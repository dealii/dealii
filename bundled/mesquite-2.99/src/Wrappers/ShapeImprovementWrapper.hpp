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
//
//   SUMMARY: 
//     USAGE:
//
//    AUTHOR: Thomas Leurent <tleurent@mcs.anl.gov>
//       ORG: Argonne National Laboratory
//    E-MAIL: tleurent@mcs.anl.gov
//
// ORIG-DATE: 14-Nov-02 at 16:51:36
//  LAST-MOD: 22-May-03 at 09:04:07 by Michael Brewer


/*! \file ShapeImprovementWrapper.hpp
  ShapeImprovementWrapper header file.

*/
// DESCRIP-END.
//


#ifndef ShapeImprovementWrapper_hpp
#define ShapeImprovementWrapper_hpp

#include "Mesquite.hpp"
#include "Wrapper.hpp"

namespace MESQUITE_NS { 
  /*! \class ShapeImprovementWrapper
       \brief Wrapper which performs a Feasible Newton solve using
       an \f$\ell_2^2 \f$ objective function template with inverse mean
       ratio.
       
     */
  class ShapeImprovementWrapper : public Wrapper {
     
  public:  
      //Constructor sets the instructions in the queue.
    MESQUITE_EXPORT
    ShapeImprovementWrapper(MsqError& err,
                            double cpu_time = 0.0, 
                            double grad_norm =1.e-6,
                            int parallel_iterations = 10);
      //Constructor sets the instructions in the queue.
    MESQUITE_EXPORT
    ShapeImprovementWrapper(double cpu_time = 0.0, 
                            double grad_norm =1.e-6,
                            int parallel_iterations = 10);

  protected:
  
    MESQUITE_EXPORT
    void run_wrapper( MeshDomainAssoc* mesh_and_domain,
                      ParallelMesh* pmesh,
                      Settings* settings,
                      QualityAssessor* qa,
                      MsqError& err );
    
      
  private:

    double maxTime, gradNorm;
      // constants
    const double untBeta;
    const double successiveEps;
    int parallelIterations;
  };
  
  
} // namespace

#endif // ShapeImprovementWrapper_hpp
