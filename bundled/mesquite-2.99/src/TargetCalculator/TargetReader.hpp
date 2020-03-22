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


/** \file TargetReader.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_TARGET_READER_HPP
#define MSQ_TARGET_READER_HPP

#include "Mesquite.hpp"
#include "TargetCalculator.hpp"
#include "ExtraDataUser.hpp"
#include "MeshInterface.hpp"

#include <vector>
#include <string>

namespace MESQUITE_NS {

class PatchData;
class MsqError;

/**\brief Internal structure used by TargetReader
 *
 * Store per-PatchData information.
 */
struct TargetReaderData {
  std::vector<TagHandle> handles2D, handles3D; //< tag handles, indexed by #tags/elem
  std::vector< MsqMatrix<3,3> > targets3D; //< cached values for last element
  std::vector< MsqMatrix<2,2> > targets2D; //< cached values for last element
  std::vector< MsqMatrix<3,2> > targetsSurface; //< cached values for last element
  size_t elementIndex;                      //< element for which values are cached.
};
  

/**\brief Read targets from tag data */
class TargetReader 
 : public TargetCalculator, 
   private ExtraDataUser<TargetReaderData>
{
  public:
  
    MESQUITE_EXPORT 
    TargetReader( bool oriented_2D_targets,
                  std::string tag_base_name = "MSQ_TARGET_MATRIX" );
    
    MESQUITE_EXPORT virtual 
    ~TargetReader();
    
    MESQUITE_EXPORT virtual 
    bool get_3D_target( PatchData &pd,
                        size_t element,
                        Sample sample,
                        MsqMatrix<3,3>& W_out,
                        MsqError& err );
    
    MESQUITE_EXPORT virtual 
    bool get_2D_target( PatchData &pd,
                        size_t element,
                        Sample sample,
                        MsqMatrix<2,2>& W_out,
                        MsqError& err );
    
    MESQUITE_EXPORT virtual 
    bool get_surface_target( PatchData &pd,
                             size_t element,
                             Sample sample,
                             MsqMatrix<3,2>& W_out,
                             MsqError& err );
    
    MESQUITE_EXPORT virtual
    bool have_surface_orient() const
      { return orient2D; }

  private:
  
    virtual void notify_patch_destroyed( TargetReaderData& data );
    virtual void notify_new_patch( PatchData& pd, TargetReaderData& data );
    virtual void notify_sub_patch( PatchData& pd, TargetReaderData& data,
                                   PatchData& subpatch, const size_t* vert_map,
                                   const size_t* elem_map, MsqError& err );

    std::string tagBaseName; //!< Base name for tags used to store targets
    bool orient2D;           //!< 2D targets included orientation (3x2 rather than 2x2)
};


} // namespace Mesquite

#endif
