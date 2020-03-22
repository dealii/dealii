/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2005 Lawrence Livermore National Laboratory.  Under 
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

    kraftche@cae.wisc.edu    

  ***************************************************************** */

#ifndef MSQ_VTK_TYPE_INFO_HPP
#define MSQ_VTK_TYPE_INFO_HPP

#include "Mesquite.hpp"
#include <sys/types.h>
#include <vector>

namespace MESQUITE_NS {

class MsqError;

/** A struct describing a mapping between a Mesquite type/node_count
 *  combination an a VTK element type.*/
struct MESQUITE_EXPORT VtkTypeInfo {
  
  const char* name;
  
  unsigned vtkType;        //!< The VTK type number
  EntityTopology msqType;  //!< The Mesquite element topology for the VTK type
  unsigned numNodes;       //!< The number of nodes in the VTK type.
  const unsigned* vtkConnOrder;  /**< NULL if VTK node ordering is the same as 
                            *   Mesquite's internal ordering.  If non-null,
                            *   an array of length VtkTypeInfo::numNodes, indexed
                            *   with the Mesquite connectivity offset and
                            *   containing the corresponding VTK connectivity
                            *   offset
                            */
  
    /** Get VtkTypeInfo from VTK type number */
  static const VtkTypeInfo* find_type( unsigned vtk_type, MsqError& err );
    /** Get VtkTypeInfo from Mesquite type and number of nodes */
  static const VtkTypeInfo* find_type( EntityTopology msq_type,
                                       unsigned num_nodes,
                                       MsqError& err );
   
    /** Reorder element connectivty list for writing to a VTK file
     *
     * If the VTK node ordering is the same as Mesquite's node ordering
     * for the type, the input list is not changed.  If the canonical
     * ordering for the elment differs, the passed list will be reordered
     * for writing to a VTK file.
     */
  void mesquiteToVtkOrder( std::vector<size_t>& connectivity_list ) const; 
};


} // namespace Mesquite

#endif
