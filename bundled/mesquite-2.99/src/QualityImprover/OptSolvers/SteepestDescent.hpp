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
/*!
  \file   SteepestDescent.hpp
  \brief  

  The SteepestDescent Class implements the steepest descent algorithm in
  order to move a free vertex to an optimal position given an
  ObjectiveFunction object and a QualityMetric object.

  \author Thomas Leurent
  \date   2002-06-13
*/

#ifndef Mesquite_SteepestDescent_hpp 
#define Mesquite_SteepestDescent_hpp

#include "Mesquite.hpp"
#include "VertexMover.hpp"
#include "PatchSetUser.hpp"

namespace MESQUITE_NS
{
  class ObjectiveFunction;

  /*! \class SteepestDescent

      This is a very basic implementation of the steepest descent optimization algorythm.
   */ 
  class SteepestDescent : public VertexMover, public PatchSetUser
  {
  public:
    MESQUITE_EXPORT 
    SteepestDescent(ObjectiveFunction* of);

    MESQUITE_EXPORT virtual
    ~SteepestDescent() { }
    
    MESQUITE_EXPORT virtual
    std::string get_name() const;
    
    MESQUITE_EXPORT virtual
    PatchSet* get_patch_set();
    
    MESQUITE_EXPORT
    bool project_gradient() const 
      { return projectGradient; }
    
    MESQUITE_EXPORT
    void project_gradient( bool yesno ) 
      { projectGradient = yesno; }
    
    //bool cosine_projection_step() const 
    //  { return cosineStep; }
    //
    //void cosine_projection_step( bool yesno ) 
    //  { cosineStep = yesno; }
    
  protected:
    MESQUITE_EXPORT virtual
    void initialize( PatchData &pd, MsqError &err );
    MESQUITE_EXPORT virtual
    void optimize_vertex_positions( PatchData &pd, MsqError &err );
    MESQUITE_EXPORT virtual
    void initialize_mesh_iteration( PatchData &pd, MsqError &err );
    MESQUITE_EXPORT virtual
    void terminate_mesh_iteration( PatchData &pd, MsqError &err );
    MESQUITE_EXPORT virtual
    void cleanup();
  
  private:
    bool projectGradient;
    //bool cosineStep;
    SteepestDescent(const SteepestDescent &pd); //disable copying
    SteepestDescent& operator=(const SteepestDescent &pd);  //disable assignment

  };
  
}

#endif
