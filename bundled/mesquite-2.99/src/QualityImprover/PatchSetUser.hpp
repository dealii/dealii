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

#ifndef MSQ_PATCH_SET_USER_HPP
#define MSQ_PATCH_SET_USER_HPP


#include "Mesquite.hpp"
#include "VertexPatches.hpp"
#include "GlobalPatch.hpp"

namespace MESQUITE_NS {

/**\brief Utility class for handling variable patch types
 *
 * Common implementation for classes supporting variable patch types.
 */
class PatchSetUser
{
public:

    PatchSetUser( bool defaultGlobal ) :
      myVertexPatches( 1, true )
      { 
        if(defaultGlobal)
          activePatchSet = &myGlobalPatch;
        else
          activePatchSet = &myVertexPatches; 
      }
    
    PatchSetUser( PatchSet* my_patch_set ) :
      myVertexPatches( 1, true ),
      activePatchSet( my_patch_set )
      {}
    
    virtual ~PatchSetUser();
      
      /**\brief Use a single patch representing the entire mesh */
    void use_global_patch() 
      { activePatchSet = &myGlobalPatch; }
    
      /**\brief Using a single patch representing the entire mesh */
    bool using_global_patch() const 
      { return activePatchSet == &myGlobalPatch; }
    
      /**\brief Construct a patch for each free vertex in the mesh 
       *\param num_layers Number of layers of adjacent elements to
       *                  include in the patch.  If unsure, use 1.
       */
    void use_element_on_vertex_patch( unsigned num_layers = 1 )
      { 
        activePatchSet = &myVertexPatches;
        myVertexPatches.set_num_layers( num_layers );
      }
    
      /**\brief True if a patch will be constructed for each free vertex */
    bool using_element_on_vertex_patch() const
      { return activePatchSet == &myVertexPatches; }
    
      /**\brief Get the number of layers of elements that will be included
       *        in an element-on-vertex patch
       */
    unsigned num_element_on_vertex_layers() const
      { return myVertexPatches.get_num_layers(); }
      
    virtual PatchSet* get_patch_set()
      { return activePatchSet; }
      
    void use_patch_set( PatchSet* patch_set )
      { activePatchSet = patch_set; }

private:

    VertexPatches myVertexPatches;
    GlobalPatch   myGlobalPatch;
    PatchSet* activePatchSet;
};



} // namespace Mesquite

#endif
