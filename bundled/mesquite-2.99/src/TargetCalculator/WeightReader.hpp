/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Sandia National Laboratories.  Developed at the
    University of Wisconsin--Madison under SNL contract number
    624796.  The U.S. Government and the University of Wisconsin
    retain certain rights to this software.

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


/** \file WeightReader.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_WEIGHT_READER_HPP
#define MSQ_WEIGHT_READER_HPP

#include "Mesquite.hpp"
#include "WeightCalculator.hpp"
#include "ExtraDataUser.hpp"
#include "MeshInterface.hpp"

#include <vector>
#include <string>

namespace MESQUITE_NS {

class PatchData;
class MsqError;

/**\brief Internal structure used by WeightReader
 *
 * Store per-PatchData information.
 */
struct WeightReaderData {
  std::vector<TagHandle> handles; //< tag handles, indexed by #tags/elem
  std::vector< double > weights;  //< cached values for last element
  size_t elementIndex;                //< element for which values are cached.
};
  

/**\brief Read targets from tag data */
class WeightReader 
 : public WeightCalculator, 
   private ExtraDataUser<WeightReaderData>
{
  public:
  
    MESQUITE_EXPORT
    WeightReader( std::string tag_base_name = "MSQ_TARGET_WEIGHT"  );
    
    MESQUITE_EXPORT virtual 
    ~WeightReader();
    
    MESQUITE_EXPORT virtual 
    double get_weight( PatchData &pd,
                       size_t element,
                       Sample Sample,
                       MsqError& err );
    
  private:
  
    virtual void notify_patch_destroyed( WeightReaderData& data );
    virtual void notify_new_patch( PatchData& pd, WeightReaderData& data );
    virtual void notify_sub_patch( PatchData& pd, 
                                   WeightReaderData& data,
                                   PatchData& subpatch, 
                                   const size_t* vertex_map,
                                   const size_t* element_map,
                                   MsqError& err );

    std::string tagBaseName;
};


} // namespace Mesquite

#endif
