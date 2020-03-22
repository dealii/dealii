/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Lawrence Livermore National Laboratory.  Under 
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

/*!
  \file   MsqIBase.cpp
  \brief  


  \author Jason Kraftcheck
  \date   2007-08-14
*/

#include "iBase.h"
#include "MsqIBase.hpp"
#include "MsqDebug.hpp"

namespace MESQUITE_NS {

std::string process_itaps_error( int ierr )
{
  std::string result( "ITAPS ERROR: " );
  switch (ierr) {
    case iBase_MESH_ALREADY_LOADED:      result += "File Already Loaded";   break;
    case iBase_FILE_NOT_FOUND:           result += "File Not Found";        break;
    case iBase_FILE_WRITE_ERROR:         result += "File Write Error";      break;
    case iBase_NIL_ARRAY:                result += "NULL Array";            break;
    case iBase_BAD_ARRAY_SIZE:           result += "Bad Array Size";        break;
    case iBase_BAD_ARRAY_DIMENSION:      result += "Bad Array Dimension";   break;
    case iBase_INVALID_ENTITY_HANDLE:    result += "Invalid Handle";        break;
    case iBase_INVALID_ENTITY_COUNT:     result += "Invalid Count";         break;
    case iBase_INVALID_ENTITY_TYPE:      result += "Invalid Type";          break;
    case iBase_INVALID_ENTITY_TOPOLOGY:  result += "Invalid Topology";      break;
    case iBase_BAD_TYPE_AND_TOPO:        result += "Invalid Type";          break;
    case iBase_ENTITY_CREATION_ERROR:    result += "Creation Failed";       break;
    case iBase_INVALID_TAG_HANDLE:       result += "Invalid Tag";           break;
    case iBase_TAG_NOT_FOUND:            result += "Tag Not Found";         break;
    case iBase_TAG_ALREADY_EXISTS:       result += "Tag Exists";            break;
    case iBase_TAG_IN_USE:               result += "Tag In Use";            break;
    case iBase_INVALID_ENTITYSET_HANDLE: result += "Invalid Handle";        break;
    case iBase_INVALID_ITERATOR_HANDLE:  result += "Invalid Iterator";      break;
    case iBase_INVALID_ARGUMENT:         result += "Invalid Argument";      break;
    case iBase_MEMORY_ALLOCATION_FAILED: result += "Out of Memory";         break;
    case iBase_NOT_SUPPORTED:            result += "Not Supported";         break;
    default:                             result += "Uknown/Internal Error"; break;
  }
  MSQ_DBGOUT(1) << result << std::endl;
  return result;
}

} // namespace Mesquite

