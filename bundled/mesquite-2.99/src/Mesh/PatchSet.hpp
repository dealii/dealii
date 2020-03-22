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

#ifndef MSQ_PATCH_SET_HPP
#define MSQ_PATCH_SET_HPP

/** \file PatchSet.hpp
 *  \brief 
 *  \author Jason Kraftcheck
 */

#include <vector>

#include "Mesquite.hpp"
#include "MeshInterface.hpp"

namespace MESQUITE_NS {

class MsqError;

/**\brief Specify a division of the Mesh into working patches
 *
 * This class provides an interface for specifying how Mesquite
 * will divide the active mesh into working patches.
 */
class MESQUITE_EXPORT PatchSet 
{
  public:
  
    typedef void* PatchHandle;
    
    inline PatchSet() {}
    
    /**\brief Declare destructor virtual */
    virtual ~PatchSet();
  
    /**\brief Specify the working Mesh */
    inline void set_mesh( Mesh* mesh ) { myMesh = mesh; }
    
    /**\brief Get a list of handles, one for each patch */
    virtual void get_patch_handles( std::vector<PatchHandle>& patch_handles_out,
                                    MsqError& err ) = 0;
    
    /**\brief Get the mesh entities in a patch
     *
     * Given one of the handles returned by get_patch_handles(),
     * return the mesh entities in the corresponding patch.
     *\param patch_handle one of the handles returned by get_patch_handles()
     *\param elem_handles_out the list of elements in the mesh
     *\param free_vertices_out The list of vertices interior to the patch
     *                         If this list is empty, it is assumed that all
     *                         vertices in the closure of the elements are 
     *                         free.  
     */
    virtual void get_patch( PatchHandle patch_handle,
                            std::vector<Mesh::ElementHandle>& elem_handles_out,
                            std::vector<Mesh::VertexHandle>& free_vertices_out,
                            MsqError& err ) = 0;
    
      /**\brief get the Mesh object passed to set_mesh() */                        
    inline Mesh* get_mesh() const { return myMesh; }
    
  private:
  
      /**\brief disallow copying*/
    PatchSet( const PatchSet& );
      /**\brief disallow assignment*/
    PatchSet& operator=(const PatchSet& );
    
    Mesh* myMesh;
};

} // namespace Mesquite

#endif
