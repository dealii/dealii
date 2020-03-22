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
  \file   LaplacianSmoother.cpp
  \brief  

  The LaplacianSmoother Class is the concrete class
  that performs Laplacian Smoothing

  \author Thomas Leurent
  \date   2002-01-17
*/

#include "LaplacianSmoother.hpp"
#include "PatchData.hpp"
#include "MsqVertex.hpp"

namespace MESQUITE_NS {

std::string LaplacianSmoother::get_name() const { return "LaplacianSmoother"; }

LaplacianSmoother::~LaplacianSmoother() 
{ }    

void LaplacianSmoother::optimize_vertex_positions( PatchData &pd, 
                                                   MsqError &err )
{
  assert(pd.num_free_vertices() == 1);
  const size_t center_vtx_index = 0;
  
  adjVtxList.clear();
  pd.get_adjacent_vertex_indices( center_vtx_index, adjVtxList, err );
  MSQ_ERRRTN(err);
  
  if (adjVtxList.empty())
    return;
  
  const MsqVertex* verts = pd.get_vertex_array(err);
  const size_t n = adjVtxList.size();
  
  Vector3D new_pos = verts[ adjVtxList[0] ];
  for (size_t i = 1; i < n; ++i)
    new_pos += verts[ adjVtxList[i] ];
  new_pos *= 1.0/n;
  pd.set_vertex_coordinates( new_pos, center_vtx_index, err );
  pd.snap_vertex_to_domain( center_vtx_index, err );  MSQ_ERRRTN(err);
}

} // namespace Mesquite

