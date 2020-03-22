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

#ifndef MSQ_LINEAR_HEXAHEDRON_HPP
#define MSQ_LINEAR_HEXAHEDRON_HPP

#include "Mesquite.hpp"
#include "MappingFunction.hpp"

namespace MESQUITE_NS {

/**\brief Linear mapping function for a hexahedral element
 *
 * This class implements the MappingFunction interface, providing
 * a Linear shape function for hexahedral elements.
 *
 * \f$\vec{x}(\xi,\eta,\zeta) = \sum_{i=0}^{7} N_i(\xi,\eta,\zeta) \vec{x_i}\f$
 * 
 * \f$N_0(\xi,\eta,\zeta) = (1-\xi)(1-\eta)(1-\zeta)\f$
 * 
 * \f$N_1(\xi,\eta,\zeta) =    \xi (1-\eta)(1-\zeta)\f$
 * 
 * \f$N_2(\xi,\eta,\zeta) =    \xi    \eta (1-\zeta)\f$
 * 
 * \f$N_3(\xi,\eta,\zeta) = (1-\xi)   \eta (1-\zeta)\f$
 * 
 * \f$N_4(\xi,\eta,\zeta) = (1-\xi)(1-\eta)   \zeta \f$
 * 
 * \f$N_5(\xi,\eta,\zeta) =    \xi (1-\eta)   \zeta \f$
 * 
 * \f$N_6(\xi,\eta,\zeta) =    \xi    \eta    \zeta \f$
 * 
 * \f$N_7(\xi,\eta,\zeta) = (1-\xi)   \eta    \zeta \f$
 * 
 */ 
class MESQUITE_EXPORT LinearHexahedron : public MappingFunction3D
{
public:

  virtual
  EntityTopology element_topology() const;
  
  virtual
  int num_nodes() const;

  virtual 
  void coefficients( Sample location,
                     NodeSet nodeset,
                     double* coeff_out,
                     size_t* indices_out,
                     size_t& num_coeff_out,
                     MsqError& err ) const;
  
  virtual 
  void derivatives( Sample location,
                    NodeSet nodeset,
                    size_t* vertex_indices_out,
                    MsqVector<3>* d_coeff_d_xi_out,
                    size_t& num_vtx,
                    MsqError& err ) const;

  virtual
  void ideal( Sample location, 
              MsqMatrix<3,3>& jacobian_out,
              MsqError& err ) const;
};



} // namespace Mesquite

#endif
