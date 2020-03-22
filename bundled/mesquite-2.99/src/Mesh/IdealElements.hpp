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


/** \file IdealElements.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_IDEAL_ELEMENTS_HPP
#define MSQ_IDEAL_ELEMENTS_HPP

#include "Mesquite.hpp"

namespace MESQUITE_NS {

class Vector3D;

/**\brief Get ideal element with unit edge length 
 *
 * Get list of vertex coordinates for an ideal element with it's
 * centroid at the origin and all edges of unit length.  Surface
 * elements lie in the XY plane.
 *
 *\param type the type of the element to obtain.
 *\param unit_height_pyramid If true, ideal pyramid has it's height equal
 *          to the length of an edge of the base, rather than the default
 *          of equilateral triangular faces.
 *\return corner vertex coordinates in canonical order.
 */
const Vector3D* unit_edge_element( EntityTopology type, bool unit_height_pyramid = false);

/**\brief Get ideal element with unit area or volume
 *
 * Get list of vertex coordinates for an ideal element with it's
 * centroid at the origin and unit area/volume.  Surface
 * elements lie in the XY plane.
 *
 *\param type the type of the element to obtain.
 *\param unit_height_pyramid If true, ideal pyramid has it's height equal
 *          to the length of an edge of the base, rather than the default
 *          of equilateral triangular faces.
 *\return corner vertex coordinates in canonical order.
 */
const Vector3D* unit_element( EntityTopology type, bool unit_height_pyramid = false );

} // namespace Mesquite

#endif
