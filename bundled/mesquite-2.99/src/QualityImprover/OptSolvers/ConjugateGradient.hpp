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
  \file   ConjugateGradient.hpp
  \brief 

  Conjugate Gradient minimization method ...

  \author Michael Brewer
  \date   2002-06/19
*/

#ifndef Mesquite_ConjugateGradient_hpp 
#define Mesquite_ConjugateGradient_hpp
#include "Mesquite.hpp"
#include "VertexMover.hpp"
#include "PatchData.hpp"
#include "PatchSetUser.hpp"


namespace MESQUITE_NS
{
  class ObjectiveFunction;


  /*! \class ConjugateGradient
    \brief Optimizes the objective function using the Polack-Ribiere scheme.
   */ 
  class ConjugateGradient : public VertexMover, public PatchSetUser
  {
  public:
    MESQUITE_EXPORT ConjugateGradient(ObjectiveFunction* objective);
    MESQUITE_EXPORT ConjugateGradient(ObjectiveFunction* objective, 
                                      MsqError &err);

    MESQUITE_EXPORT virtual ~ConjugateGradient();
    
    MESQUITE_EXPORT virtual std::string get_name() const;
    
    virtual PatchSet* get_patch_set();
    
      //!Just for debugging purposes or for obtaining more data
      //! during the optimization process.
    MESQUITE_EXPORT void set_debugging_level(int new_lev)
      {
        conjGradDebug=new_lev;
      }
    
  protected:
      
      //!Initialize data for smoothing process
    virtual void initialize(PatchData &pd, MsqError &err);
 
    virtual void optimize_vertex_positions(PatchData &pd, MsqError &err);

    virtual void initialize_mesh_iteration(PatchData &pd, MsqError &err);
    
    virtual void terminate_mesh_iteration(PatchData &pd, MsqError &err);
    
      //!Delete arrays initially created in initialize().
    virtual void cleanup();
   
      //!Returns the step distance to take in the search direction.
    double get_step(PatchData &pd, double f0,int &j, MsqError &err);
      
      //!Culls the vertex list free_vertex_list.     
      //void cull_list(PatchData &pd, double beta, MsqError &err);
    
     
private:
    std::vector<Vector3D> fGrad, pGrad, fNewGrad;
    PatchDataVerticesMemento* pMemento;
      //just for debugging
    int conjGradDebug;
  };

  

}
 
#endif
