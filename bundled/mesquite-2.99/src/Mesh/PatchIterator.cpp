/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Lawrence Livermore National Laboratory.  Under 
    the terms of Contract B545069 with the University of Wisconsin -- 
    Madison, Lawrence Livermore National Laboratory retains certain
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

    (2006) kraftche@cae.wisc.edu    

  ***************************************************************** */

#ifndef MSQ_PATCH_ITERATOR_CPP
#define MSQ_PATCH_ITERATOR_CPP

#include "PatchIterator.hpp"
#include "MsqError.hpp"
#include "PatchData.hpp"

namespace MESQUITE_NS {

bool PatchIterator::get_next_patch( PatchData& pd, MsqError& err )
{
  if (!patches.size()) {
    if (pd.get_mesh() != patchSet->get_mesh()) {
      if (pd.get_mesh() == 0)
        pd.set_mesh( patchSet->get_mesh() );
      else if(patchSet->get_mesh() == 0)
        patchSet->set_mesh( pd.get_mesh() );
      else {
        MSQ_SETERR(err)("PatchSet and PatchData do not share the same Mesh instance.",
                        MsqError::INVALID_STATE);
        return false;
      }
    }
    patchSet->get_patch_handles( patches, err ); MSQ_ERRZERO(err);
    current = patches.begin();
  } 
  
  if (current == patches.end())
    return false;
  
  patchSet->get_patch( *current, elems, verts, err ); MSQ_ERRZERO(err);
  pd.set_mesh_entities( elems, verts, err ); MSQ_ERRZERO(err);
  ++current;
  return true;
}

void PatchIterator::reset()
{
  current = patches.begin();
}


} // namespace Mesquite

#endif
