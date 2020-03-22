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

/*! \file MsqVertex.hpp
  \brief Mesquite's vertex object.

  \author Darryl Melander
  \author Thomas Leurent
*/
#ifndef MSQVERTEX_HPP
#define MSQVERTEX_HPP

#include "Mesquite.hpp"
#include "Vector3D.hpp"

namespace MESQUITE_NS
{
    /*!
      \class MsqVertex
      \brief MsqVertex is the Mesquite object that stores information about
      the vertices in the mesh.

      This class has no virtual destructor for performance reasons.
      !!! Make sure NOT to delete a MsqVertex object from a pointer
          to Vector3D !!!
    */
  class MESQUITE_EXPORT MsqVertex : public Vector3D
   {
   public:
       //!Construct vertex using three doubles.
     MsqVertex(double x, double y, double z) 
       : Vector3D(x, y, z), vertexBitFlags(0)
       {}
     
       //!Construct vertex using Vector3D.
     MsqVertex(const Vector3D &vec) 
       : Vector3D(vec), vertexBitFlags(0)
       {}
     
       //!Construct default vertex with coordinates (0.0,0.0,0.0)
     MsqVertex() 
       : Vector3D(0,0,0), vertexBitFlags(0)
       {}

       //!Construct default vertex with coordinates (0.0,0.0,0.0)
     MsqVertex(const MsqVertex& rhs) 
       : Vector3D(rhs), vertexBitFlags(rhs.vertexBitFlags)
       {}

       //! Initializes with coordinates. Sets tag data/pointer to 0.
     MsqVertex& operator=(const Vector3D& rhs)
       { Vector3D::operator=(rhs);
         return *this; }
     
       // This allows for 8 flag bits.
       // I don't think we'll want more than that (yet).
     typedef unsigned char FlagMask;
     
       //! \enum FlagMaskID
       //!   Those are the available flags... currently only return
       //!   is_free.
       //!   Developers: The values used in that enum are used by a bitset,
       //!               so they have to be 2-based (2,4,8,16,32, ...)
     enum FlagMaskID
     {
       MSQ_HARD_FIXED  = 1<<0, //!< vertex is always fixed. This can only be set on and never off.
       MSQ_DEPENDENT   = 1<<1, //!< higher-order node w/ position determined by mapping function
       MSQ_CULLED      = 1<<2, //!< vertex is fixed. This flag can be set on and off. 
       MSQ_PATCH_FIXED = 1<<3, //!< vertex is fixed only because it is on patch boundary (not by app request)
       MSQ_MARK        = 1<<4, //!< arbitrary mark for use by code - clear before using
       MSQ_FIXED = (MSQ_HARD_FIXED|MSQ_CULLED|MSQ_PATCH_FIXED)
     };
       //!Returns true if vertex is ``free''.
     bool is_free_vertex() const
       { return (vertexBitFlags & (MSQ_HARD_FIXED|MSQ_PATCH_FIXED)) == 0; }
     
     void set_soft_fixed_flag()
       { vertexBitFlags|=MSQ_CULLED; }
     
     void remove_soft_fixed_flag()
       { vertexBitFlags &= (~MSQ_CULLED); }
     
     void set_hard_fixed_flag()
       { vertexBitFlags|=MSQ_HARD_FIXED; }
     
     void set_vertex_flag(FlagMaskID flag)
       { vertexBitFlags|=flag; }
     
     void remove_vertex_flag(FlagMaskID flag)
       { vertexBitFlags &= (~flag); }
     
     bool is_flag_set(FlagMaskID flag) const
       { return (vertexBitFlags & flag) != 0; }
    
     FlagMask get_flags() const
      { return vertexBitFlags; }
    
     FlagMask& flags() 
      { return vertexBitFlags; }
    
     void set_flags( FlagMask flags )
      { vertexBitFlags = flags; }
     
   private:
     FlagMask vertexBitFlags;
   };

} //namespace
  

#endif // MsqVertex_hpp
