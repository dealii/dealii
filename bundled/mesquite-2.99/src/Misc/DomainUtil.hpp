/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2010 Sandia National Laboratories.  Developed at the
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

    (2010) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file DomainUtil.hpp
 *  \brief Utility functions for use in build-in geometric domains
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_DOMAIN_UTIL_HPP
#define MSQ_DOMAIN_UTIL_HPP

#include "Mesquite.hpp"
#include "MeshInterface.hpp"
#include <vector>
#include <stdlib.h> // size_t

namespace MESQUITE_NS {

class MsqVertex;
class Vector3D;
class Mesh;

namespace DomainUtil {

void bounding_box( const MsqVertex* vertex_array,
                   size_t num_vertices,
                   Vector3D& box_min,
                   Vector3D& box_max );

double max_box_extent( const MsqVertex* vertex_array, size_t num_vertices );

inline double 
default_tolerance( const MsqVertex* vertex_array, size_t num_vertices )
  { return 1e-3 * max_box_extent( vertex_array, num_vertices ); }

void get_fixed_vertices( Mesh* mesh_instance,
                         const Mesh::VertexHandle* vertex_array,
                         size_t num_vertices,
                         std::vector<Mesh::VertexHandle>& fixed_handles_out,
                         MsqError& err );

bool non_colinear_vertices( const MsqVertex* vertex_array,
                            size_t num_vertices,
                            Vector3D coords_out[3],
                            double epsilon );

bool non_coplanar_vertices( const MsqVertex* vertex_array,
                            size_t num_vertices,
                            Vector3D coords_out[4],
                            double epsilon );


} // namespace DomainUtil

} // namespace MESQUITE_NS

#endif
