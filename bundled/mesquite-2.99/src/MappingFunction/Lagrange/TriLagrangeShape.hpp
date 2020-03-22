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


/** \file TriLagrangeShape.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_TRI_LAGRANGE_SHAPE_HPP
#define MSQ_TRI_LAGRANGE_SHAPE_HPP

#include "Mesquite.hpp"
#include "MappingFunction.hpp"

namespace MESQUITE_NS {

/**\brief Lagrange shape function for triangle elements
 *
 * This class implements the MappingFunction interface, providing
 * a Lagrange shape function for 6-node triangle.
 *
 * \f$\vec{x}(r,s) = sum_{i=0}^{n-1} N_i(r,s) \vec{x_i}\f$
 * 
 * \f$N_1 = t(2t - 1)\f$
 * 
 * \f$N_2 = r(2r - 1)\f$
 * 
 * \f$N_3 = s(2s - 1)\f$
 * 
 * \f$N_4 = 4rt\f$
 * 
 * \f$N_5 = 4rs\f$
 * 
 * \f$N_6 = 4st\f$
 *
 * \f$t = 1 - r - s\f$
 */
class MESQUITE_EXPORT TriLagrangeShape : public MappingFunction2D
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
                    MsqVector<2>* d_coeff_d_xi_out,
                    size_t& num_vtx,
                    MsqError& err ) const;

  virtual
  void ideal( Sample location, 
              MsqMatrix<3,2>& jacobian_out,
              MsqError& err ) const;
};



} // namespace Mesquite

#endif
