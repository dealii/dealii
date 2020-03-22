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

#ifndef MSQ_QUAD_LAGRANGE_SHAPE_HPP
#define MSQ_QUAD_LAGRANGE_SHAPE_HPP

/** \file QuadLagrangeShape.hpp
 *  \brief Lagrange mapping funtion for a 9-node quad.
 *  \author Jason Kraftcheck
 */
 
#include "MappingFunction.hpp"

namespace MESQUITE_NS {

/**\brief Lagrange shape function for 9-node quadrilateral elements
 *
 * This class implements the MappingFunction interface, providing
 * a Lagrange shape function for quadrilateral elements.
 *
 * \f$\vec{x}(\xi,\eta) = \sum_{i=0}^{n-1} N_i(\xi,\eta) \vec{x_i}\f$
 * 
 * \f$N_a = l^2_b(\xi) l^2_c(\eta)\f$
 *
 * \f$l^2_1(\xi) = (\xi - 1) (2 \xi - 1)\f$
 *
 * \f$l^2_2(\xi) = 4 \xi (1 - \xi)\f$
 *
 * \f$l^2_3(\xi) = \xi (2 \xi - 1)\f$
 *
 * \f$\begin{array}{ccc}
 *    a & b & c \\ \hline
 *    0 & 1 & 1 \\
 *    1 & 3 & 1 \\
 *    2 & 3 & 3 \\
 *    3 & 1 & 3 \\
 *    4 & 2 & 1 \\
 *    5 & 3 & 2 \\
 *    6 & 2 & 3 \\
 *    7 & 1 & 2 \\
 *    8 & 2 & 2 \end{array}\f$
 *    
 */
class MESQUITE_EXPORT QuadLagrangeShape : public MappingFunction2D
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
