/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2007 Sandia National Laboratories.  Developed at the
    University of Wisconsin--Madison under SNL contract number
    624796.  The U.S. Government and the University of Wisconsin
    retain certain rights to this software.

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

    (2008) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file TrustRegion.hpp
 *  \brief Port Todd Munson's trust region solver to Mesquite
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_TRUST_REGION_HPP
#define MSQ_TRUST_REGION_HPP

#include "Mesquite.hpp"
#include "VertexMover.hpp"
#include "PatchSetUser.hpp"
#include "MsqHessian.hpp"

namespace MESQUITE_NS {

class PatchDataVerticesMemento;

class TrustRegion : public VertexMover, public PatchSetUser
{
  public:
  
    MESQUITE_EXPORT TrustRegion( ObjectiveFunction* of );
    
    MESQUITE_EXPORT virtual ~TrustRegion();
    
    PatchSet* get_patch_set();
    
    MESQUITE_EXPORT std::string get_name() const;
    
  protected:
    
    virtual void initialize( PatchData& pd, MsqError& err );
    virtual void optimize_vertex_positions( PatchData& pd, MsqError& err );
    virtual void initialize_mesh_iteration( PatchData& pd, MsqError& err );
    virtual void terminate_mesh_iteration( PatchData& pd, MsqError& err );
    virtual void cleanup();
    
  private:
  
    void compute_preconditioner( MsqError& err );
    void apply_preconditioner( Vector3D* z, Vector3D* r, MsqError& err );
  
    PatchDataVerticesMemento* mMemento;
    MsqHessian mHess;
    std::vector<Vector3D> mGrad;
    std::vector<Vector3D> wVect, zVect, dVect, pVect, rVect;
    std::vector<double> preCond;
};

} // namespace Mesquite

#endif
