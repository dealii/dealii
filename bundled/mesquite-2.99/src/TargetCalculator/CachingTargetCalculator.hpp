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


/** \file CachingTargetCalculator.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_CACHING_TARGET_CALCULATOR_HPP
#define MSQ_CACHING_TARGET_CALCULATOR_HPP

#include "Mesquite.hpp"
#include "TargetCalculator.hpp"
#include "ExtraDataUser.hpp"
#include "MsqMatrix.hpp"

#include <vector>

namespace MESQUITE_NS {

struct CachedTargetData {
  std::vector<size_t> elementOffsets;
  std::vector< MsqMatrix<3,3> > targets3D;
  std::vector< MsqMatrix<3,2> > targetsSurface;
  std::vector< MsqMatrix<2,2> > targets2D;
  bool has_data() const { return !elementOffsets.empty(); }
  void clear() 
    { 
      elementOffsets.clear();
      targets3D.clear(); 
      targets2D.clear(); 
      targetsSurface.clear(); 
    }
};

/**\brief Cache target matrices on PatchData
 *
 * This class is a decorator for concrete TargetCalculator
 * that caches the previously calculated targets.  This class
 * should not be used if the target matrices produced by the
 * concrete target calculator are a function of the vertex
 * positions in the active mesh.
 */
class CachingTargetCalculator : public TargetCalculator, 
                                private ExtraDataUser<CachedTargetData>
{
public:
  
  CachingTargetCalculator( TargetCalculator* cached )
    : cachedCalculator( cached )
      {}
  
  virtual ~CachingTargetCalculator();
  
  TargetCalculator* get_cached_calculator() const
    { return cachedCalculator; }
  
  virtual bool get_3D_target( PatchData& pd, 
                              size_t element,
                              Sample sample,
                              MsqMatrix<3,3>& W_out,
                              MsqError& err );
  
  virtual bool get_surface_target( PatchData& pd,
                              size_t element,
                              Sample sample,
                              MsqMatrix<3,2>& W_out,
                              MsqError& err );
  
  virtual bool get_2D_target( PatchData& pd,
                              size_t element,
                              Sample sample,
                              MsqMatrix<2,2>& W_out,
                              MsqError& err );

  virtual bool have_surface_orient() const
    { return cachedCalculator->have_surface_orient(); }

protected:
                            
  void notify_patch_destroyed( CachedTargetData& d );

  void notify_sub_patch( PatchData& orig_patch,
                         CachedTargetData& data,
                         PatchData& sub_patch,
                         const size_t* vertex_index_map,
                         const size_t* element_index_map,
                         MsqError& err );
  
  void notify_new_patch( PatchData& patch, CachedTargetData& data );

private:
  
  TargetCalculator *const cachedCalculator;
};

} // namespace Mesquite

#endif
