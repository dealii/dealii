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
  \file   Randomize.cpp
  \brief  

  The Randomize Class is the concrete class that randomizes
  the vertex positions.          

  \author Michael Brewer
  \date   2002-10-27
*/


#include "Randomize.hpp"
#include "MsqFreeVertexIndexIterator.hpp"
#include "MsqDebug.hpp"
#include <math.h>

using namespace Mesquite;

Randomize::~Randomize() {}

std::string Randomize::get_name() const { return "Randomize"; }

PatchSet* Randomize::get_patch_set() { return &patchSet; }

Randomize::Randomize() : mPercent(0.5), patchSet( 1, true )
  { }

Randomize::Randomize(double percent) : mPercent(percent), patchSet( 1, true )
  { }
  
void Randomize::initialize(PatchData &/*pd*/, MsqError &)
{ }

void Randomize::initialize_mesh_iteration(PatchData &/*pd*/, MsqError &/*err*/)
{
  //  cout << "- Executing Randomize::iteration_complete()\n";
}

  
      /*!Function calculates a scale factor for the patch, then moves
        the incident vertex randomly in each of the three coordinate
        directions (relative to the scale factor multiplied by mPercent).
      */
static inline void randomize_vertex(PatchData &pd,
                                    size_t free_ind,
                                    double percent,
                                    MsqError &err)
{
  size_t i;
  short j;
  const MsqVertex* verts = pd.get_vertex_array(err);  MSQ_ERRRTN(err);
  const size_t num_vtx = pd.num_nodes();
    //a scale w.r.t. the patch size
  double scale_factor=0.0;
    //a "random" number between -1 and 1
  double rand_double=0.0;
    //a "random" int
  int rand_int=0;
  if (num_vtx<=1){
    MSQ_PRINT(1)("WARNING: Number of incident vertex is zero.  Returning.\n");
    return;
  }

  for (i=0;i<num_vtx;++i){
    if(i != free_ind)
      scale_factor+=(verts[i]-verts[free_ind]).length();
  }
  scale_factor/=( (double) num_vtx - 1.0 );    
  Vector3D delta;
  for (j=0;j<3;++j){
    rand_int = rand();
      //number between 0 and 1000
    rand_int = rand_int%1000;
      //number between -1 and 1
    rand_double = (((double) rand_int)/500.0)-1.0;
    delta[j] = scale_factor * rand_double * percent;
  }
  pd.move_vertex( delta, free_ind, err );

  return;
}


void Randomize::optimize_vertex_positions(PatchData &pd, 
                                                MsqError &err)
{
    //cout << "- Executing Randomize::optimize_vertex_position()\n";

  // gets the array of coordinates for the patch and print it 
  // does the randomize smooth
  MsqFreeVertexIndexIterator free_iter(pd, err); MSQ_ERRRTN(err);
  free_iter.reset();
  free_iter.next();
    //find the free vertex.
  int m=free_iter.value();
  randomize_vertex(pd, m, mPercent, err);  MSQ_ERRRTN(err);
  pd.snap_vertex_to_domain(m,err); MSQ_ERRRTN(err);
}
  
void Randomize::terminate_mesh_iteration(PatchData &/*pd*/, MsqError &/*err*/)
{
  //  cout << "- Executing Randomize::iteration_complete()\n";
}
  
void Randomize::cleanup()
{
  //  cout << "- Executing Randomize::iteration_end()\n";
}
  

