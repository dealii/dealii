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


/** \file ExtraData.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_EXTRA_DATA_HPP
#define MSQ_EXTRA_DATA_HPP

#include "Mesquite.hpp"
#include <stdlib.h>

namespace MESQUITE_NS {

class PatchData;
class MsqError;

/**\brief Object used to attach auxiliary data to PatchData */
class ExtraData {
  public:
    ExtraData( PatchData& patch );

    virtual ~ExtraData();

    PatchData* get_patch_data() const { return patchPtr; }

      /**\brief Notify that the owning PatchData is being destroyed.
       *
       * Notify an ExtraData object that the patch it is attached to
       * is being destroyed.  The ExtraData will have been removed
       * from the PatchData before being notified.  Therefore 
       * this->get_patch_data() will return NULL.
       *
       * ExtraData instances will also be notified via this method
       * and removed and removed from the PatchData if the attached
       * Mesquite::Mesh or Mesquite::MeshDomain instance is changed.
       */
    virtual void notify_patch_destroyed( ) = 0;

      /**\brief Notify that the patch (mesh) in the PatchData is changing
       *
       * Notify attached ExtraData that the PatchData is being changed
       * to contain an new patch (set of mesh entities.)
       */
    virtual void notify_new_patch( ) = 0;

      /**\brief Nofity that a subpatch is being created from this patch.
       *\param sub_patch The new, populated subpatch
       *\param vertex_index_map The indices in the original patch for each
       *          vertex in the subpatch.
       *\param element_index_map The indices in the original patch for each
       *          element in the subpatch.
       */
    virtual void notify_sub_patch( PatchData& sub_patch,
                                   const size_t* vertex_index_map,
                                   const size_t* element_index_map,
                                   MsqError& err ) = 0;

  private:
    friend class PatchData;
    ExtraData* patchNext;
    PatchData* patchPtr;
};



} // namespace Mesquite

#endif
