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

#ifndef MSQ_LINEAR_TETRAHEDRON_HPP
#define MSQ_LINEAR_TETRAHEDRON_HPP

#include "Mesquite.hpp"
#include "MappingFunction.hpp"

namespace MESQUITE_NS {

/**\brief Linear mapping function for a tetrahedral element
 *
 * \f$\vec{x}(\vec{\xi})=\sum_{i=1}^4 N_i(\vec{\xi})\vec{x_i}\f$
 * - \f$N_1(\vec{\xi})=1-\xi_1-\xi_2-\xi_3\f$
 * - \f$N_2(\vec{\xi})=\xi_1\f$
 * - \f$N_3(\vec{\xi})=\xi_2\f$
 * - \f$N_4(\vec{\xi})=\xi_3\f$
 */ 
class MESQUITE_EXPORT LinearTetrahedron : public MappingFunction3D
{
public:

  virtual
  EntityTopology element_topology() const;
  
  virtual
  int num_nodes() const;
  
  virtual
  NodeSet sample_points( NodeSet higher_order_nodes ) const;
  

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
