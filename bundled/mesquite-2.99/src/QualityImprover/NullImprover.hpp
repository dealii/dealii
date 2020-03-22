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
#ifndef MESQUITE_NULL_IMPROVER_HPP
#define MESQUITE_NULL_IMPROVER_HPP
/*!
  \file   NullImprover.hpp
  \brief  The NullImprover Class is a do-nothing VertexMover.  It just
          loops over the mesh without doing any real work.
          It is used to test functions
          found in VertexMover, such as loop_over_mesh().
          
  \author Darryl Melander
  \date   2002-12-10
*/

#include "VertexMover.hpp"

namespace MESQUITE_NS
{
  class NullImprover : public VertexMover
  {
  protected:
    virtual void initialize(PatchData &, MsqError &)
      {}
    virtual void cleanup()
      {}
    virtual void optimize_vertex_positions(PatchData &, 
                                           MsqError &)
      {}
    virtual void initialize_mesh_iteration(PatchData &, 
                                           MsqError &)
      {}
    virtual void terminate_mesh_iteration(PatchData &, 
                                          MsqError &)
      {}
  };
}

#endif
