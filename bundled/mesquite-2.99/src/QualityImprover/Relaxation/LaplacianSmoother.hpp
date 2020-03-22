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
  \file   LaplacianSmoother.hpp
  \brief  

  The LaplacianSmoother Class implements the Laplacian smoothing
  for a patch with one free vertex. 

  \author Thomas Leurent
  \date   2002-01-17
*/

#ifndef Mesquite_LaplacianSmoother_hpp 
#define Mesquite_LaplacianSmoother_hpp

#include "Mesquite.hpp"
#include "RelaxationSmoother.hpp"

namespace MESQUITE_NS
{

  /*! \class LaplacianSmoother
    Moves free center vertex to the average of the neighboring vertices.
   */  
  class LaplacianSmoother : public RelaxationSmoother 
  {
  public:
    /**
     *\param OF ObjectiveFunction used by some termination criteria
     */
    MESQUITE_EXPORT LaplacianSmoother( ObjectiveFunction* OF = NULL ) 
      : RelaxationSmoother(OF) {}
    
    MESQUITE_EXPORT virtual ~LaplacianSmoother();
    MESQUITE_EXPORT virtual std::string get_name() const;
  protected:
    MESQUITE_EXPORT virtual 
    void optimize_vertex_positions( PatchData &pd, MsqError &err );
  private:
    std::vector<size_t> adjVtxList;
  };
}

#endif
