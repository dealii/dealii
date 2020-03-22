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


/** \file JacobianCalculator.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_JACOBIAN_CALCULATOR_HPP
#define MSQ_JACOBIAN_CALCULATOR_HPP

#include "Mesquite.hpp"
#include "MsqMatrix.hpp"
#include "Vector3D.hpp"
#include "NodeSet.hpp"

#include <vector>

namespace MESQUITE_NS {

class MappingFunction2D;
class MappingFunction3D;

/**\brief Calculate Jacobian matrices given vertex coordinates and MappingFunction
*/
class JacobianCalculator
{
public:    
  /**\brief Calculate Jacobian for surface element
   *
   * Calculate the Jacobian matrix at a specified location in a surface element.
   *\param mf The mapping function
   *\param ho_bits bit mask indicating which higher-order nodes are present in the element
   *               (zero for a linear element)
   *\param location Logical position within element at which to evaluate Jacobian
   *\param J_out The resulting Jacobian matrix.
   */
  void get_Jacobian_2D( const MappingFunction2D* mf,
                        NodeSet ho_bits,
                        Sample location,
                        const Vector3D* vertex_coords,
                        size_t num_vertex,
                        MsqMatrix<3,2>& J_out,
                        MsqError& err );
  
  /**\brief Calculate Jacobian for volume element
   *
   * Calculate the Jacobian matrix at a specified location in a volume element.
   *\param mf The mapping function
   *\param ho_bits bit mask indicating which higher-order nodes are present in the element
   *               (zero for a linear element)
   *\param location Logical position within element at which to evaluate Jacobian
   *\param J_out The resulting Jacobian matrix.
   */
  void get_Jacobian_3D( const MappingFunction3D* mf,
                        NodeSet ho_bits,
                        Sample location,
                        const Vector3D* vertex_coords,
                        size_t num_vertex,
                        MsqMatrix<3,3>& J_out,
                        MsqError& err );
private:

  enum { MAX_ELEM_NODES = 27 };
  size_t mIndices[MAX_ELEM_NODES];
  MsqVector<3> mDerivs3D[MAX_ELEM_NODES];
  MsqVector<2> mDerivs2D[MAX_ELEM_NODES];
};




} // namespace Mesquite

#endif
