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
//    AUTHOR: Thomas Leurent <tleurent@mcs.anl.gov>
//       ORG: Argonne National Laboratory
//    E-MAIL: tleurent@mcs.anl.gov
//
// ORIG-DATE: 15-Jan-03 at 08:05:56
//  LAST-MOD: 23-May-03 at 11:20:14 by Thomas Leurent
//
// DESCRIPTION:
// ============
/*!
  \file   FeasibleNewton.hpp
  \brief  

  The FeasibleNewton Class implements the newton non-linear programming algorythm
  in order to move a free vertex to an optimal position given an
  ObjectiveFunction object and a QualityMetric object.

  \author Thomas Leurent
  \author Todd Munson
  \date   2003-01-15
*/
// DESCRIP-END.
//

#ifndef MSQ_FeasibleNewton_hpp 
#define MSQ_FeasibleNewton_hpp

#include "Mesquite.hpp"
#include "VertexMover.hpp"
#include "MsqHessian.hpp"
#include "PatchSetUser.hpp"

namespace MESQUITE_NS
{
  class ObjectiveFunction;


  /*! \class FeasibleNewton

      \brief High Performance implementation of the Feasible Newton algorythm.

      Consider our non-linear objective function
      \f$ f: I\!\!R^{3N} \rightarrow I\!\!R \f$ where \f$ N \f$
      is the number of vertices of the mesh, and \f$ 3N \f$ is therefore the number
      of degrees of freedom of the mesh.
      The Taylor expansion of \f$ f \f$ around the point \f$ x_0 \f$ is 
      \f[ f(x_0+d) = f(x_0) + \nabla f(x_0)d + \frac{1}{2} d^T\nabla^2 f(x_0)d
          + ...  \;\;\; .\f]

      Each iteration of the Newton algorithm tries to find a descent vector that
      minimizes the above quadratic approximation, i.e. it looks for
      \f[ \min_{d} q(d;x_0) = f(x_0) + \nabla f(x_0)d + \frac{1}{2} d^T\nabla^2 f(x_0)d
          \;\; . \f]
      We know that if a quadratic function has a finite minimum, it is reached at the
      point where the function gradient is null and that the function Hessian
      is then positive definite. 
      Therefore we are looking for \f$ d \f$ such that \f$ \nabla q(d;x_0) =0 \f$. We have
      \f[ \nabla q(d;x_0) = \nabla f(x_0) + \nabla^2 f(x_0)d \;\;, \f]
      therefore we must solve for \f$ d \f$ the system
      \f[ \nabla^2 f(x_0)d = -\nabla f(x_0) \;\; . \f]

      We assume that the Hessian is positive definite and we use the conjugate gradient
      algebraic solver to solve the above system. If the conjugate gradient solver finds
      a direction of negative curvature, the Hessian was not positive definite and we take
      a step in that direction of negative curvature, which is a descent direction. 
  */ 
  class FeasibleNewton : public VertexMover, public PatchSetUser
  {
  public:
    MESQUITE_EXPORT FeasibleNewton(ObjectiveFunction* of);

    MESQUITE_EXPORT virtual ~FeasibleNewton()
    { delete coordsMem; }

    /*! Sets a minimum value for the gradient. If the gradient is below that value,
      we stop iterating. */  
    MESQUITE_EXPORT void set_lower_gradient_bound(double gradc){
        convTol=gradc;}
    
    PatchSet* get_patch_set();
    
    MESQUITE_EXPORT std::string get_name() const;
    
  protected:
    virtual void initialize(PatchData &pd, MsqError &err);
    virtual void optimize_vertex_positions(PatchData &pd,
                                           MsqError &err);
    virtual void initialize_mesh_iteration(PatchData &pd, MsqError &err);
    virtual void terminate_mesh_iteration(PatchData &pd, MsqError &err);
    virtual void cleanup();

  private:
    double convTol;
    MsqHessian mHessian;
    PatchDataVerticesMemento* coordsMem;
    bool havePrintedDirectionMessage;
  };
  
}

#endif // MSQ_FeasibleNewton_hpp 
