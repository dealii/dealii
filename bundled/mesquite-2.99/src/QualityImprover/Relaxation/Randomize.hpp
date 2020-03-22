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
  \file   Randomize.hpp
  \brief  

  The Randomize Class implements the Randomize Vertex Mover
  for a patch with one free vertex. 

  \author Michael Brewer      
  \date   2002-10-27
*/

#ifndef Mesquite_Randomize_hpp 
#define Mesquite_Randomize_hpp

#include "Mesquite.hpp"
#include "VertexMover.hpp"
#include "VertexPatches.hpp"
#include <vector>

namespace MESQUITE_NS
{

  /*! \class Randomize
   \brief Randomly perftubs the (un-culled) vertices.
  */ 
  class Randomize : public VertexMover 
  {
  public:
      //!Constructor defaulting mPercent to .05.
    MESQUITE_EXPORT Randomize();
      //!Constructor allowing user to set mPercent
    MESQUITE_EXPORT Randomize(double percent);

    MESQUITE_EXPORT virtual ~Randomize();
    
    MESQUITE_EXPORT virtual std::string get_name() const;

    MESQUITE_EXPORT virtual PatchSet* get_patch_set();
    
  protected:
    virtual void initialize(PatchData &pd, MsqError &err);
    virtual void optimize_vertex_positions(PatchData &pd,
                                         MsqError &err);
    virtual void initialize_mesh_iteration(PatchData &pd, MsqError &err);
    virtual void terminate_mesh_iteration(PatchData &pd, MsqError &err);
    virtual void cleanup();
  private:
      //! \param The percentage of the scale factor each vertex will be moved.
    double mPercent;
    std::vector<size_t> adjVtxList;
    VertexPatches patchSet;
  };


  
}

#endif
