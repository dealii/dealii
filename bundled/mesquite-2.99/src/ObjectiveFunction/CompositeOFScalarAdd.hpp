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

/*! \file CompositeOFScalarAdd.hpp

Header file for the Mesquite:: CompositeOFScalarAdd class

  \author Michael Brewer
  \date   2002-06-24
 */


#ifndef CompositeOFScalarAdd_hpp
#define CompositeOFScalarAdd_hpp

#include "Mesquite.hpp"
#include "ObjectiveFunction.hpp"

namespace MESQUITE_NS
{
   /*! \class CompositeOFScalarAdd.
       \brief Adds a scalar to a given ObjectiveFunction.
     */
   class MsqMeshEntity;
   class PatchData;
   class MsqError;
   
   class MESQUITE_EXPORT CompositeOFScalarAdd : public ObjectiveFunction
   {
   public:
     CompositeOFScalarAdd(double, ObjectiveFunction*, bool delete_OF = false);
     virtual ~CompositeOFScalarAdd();
     
      //!\brief Called at start of instruction queue processing
      //!
      //! Do any preliminary global initialization, consistency checking,
      //! etc. 
     virtual void initialize_queue( MeshDomainAssoc* mesh_and_domain,
                                    const Settings* settings,
                                    MsqError& err );

     virtual bool initialize_block_coordinate_descent( MeshDomainAssoc* mesh_and_domain,
                                                       const Settings* settings,
                                                       PatchSet* user_set,
                                                       MsqError& err );

     virtual bool evaluate( EvalType type, 
                            PatchData& pd,
                            double& value_out,
                            bool free,
                            MsqError& err ); 

     virtual bool evaluate_with_gradient( EvalType type, 
                                          PatchData& pd,
                                          double& value_out,
                                          std::vector<Vector3D>& grad_out,
                                          MsqError& err ); 

     virtual bool evaluate_with_Hessian_diagonal( EvalType type, 
                                        PatchData& pd,
                                        double& value_out,
                                        std::vector<Vector3D>& grad_out,
                                        std::vector<SymMatrix3D>& hess_diag_out,
                                        MsqError& err ); 
    
     virtual bool evaluate_with_Hessian( EvalType type, 
                                         PatchData& pd,
                                         double& value_out,
                                         std::vector<Vector3D>& grad_out,
                                         MsqHessian& Hessian_out,
                                         MsqError& err ); 

     virtual ObjectiveFunction* clone() const;

     virtual void clear();
     
     virtual int min_patch_layers() const;
     
   private:
     double mAlpha;
     ObjectiveFunction* objFunc;
     bool deleteObjFunc;
   };
}//namespace
#endif //  CompositeOFScalarAdd_hpp
