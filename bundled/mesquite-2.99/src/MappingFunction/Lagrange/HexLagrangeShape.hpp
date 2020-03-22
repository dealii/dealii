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

#ifndef MSQ_HEX_LAGRANGE_SHAPE_HPP
#define MSQ_HEX_LAGRANGE_SHAPE_HPP

/** \file HexLagrangeShape.hpp
 *  \brief Lagrange mapping funtion for a 27-node hex.
 *  \author Nicholas Voshell
 */
 
#include "MappingFunction.hpp"

namespace MESQUITE_NS {

/**\brief Lagrange shape function for 27-node hexahedral elements
 *
 * This class implements the MappingFunction interface, providing
 * a Lagrange shape function for hexahedral elements.
 *
 * \f$\vec{x}(\xi,\eta) = \sum_{i=0}^{n-1} N_i(r, s, t) \vec{x_i}\f$
 * 
 * \f$N_a = l^2_b(r) l^2_c(s) l^2_d(t)\f$
 *
 * \f$l^2_1(\xi) = (\xi - 1) (2 \xi - 1)\f$
 *
 * \f$l^2_2(\xi) = 4 \xi (1 - \xi)\f$
 *
 * \f$l^2_3(\xi) = \xi (2 \xi - 1)\f$
 *
 * Mesquite / iTAPS / MBCN order
 *
 * \f$\begin{array}{cccc}
 *    a & b & c & d \\ \hline
 *    0 & 1 & 1 & 1 \\
 *    1 & 3 & 1 & 1 \\
 *    2 & 3 & 3 & 1 \\
 *    3 & 1 & 3 & 1 \\
 *    4 & 1 & 1 & 3 \\
 *    5 & 3 & 1 & 3 \\
 *    6 & 3 & 3 & 3 \\
 *    7 & 1 & 3 & 3 \\
 *    8 & 2 & 1 & 1 \\
 *    9 & 3 & 2 & 1 \\
 *    10 & 2 & 3 & 1 \\
 *    11 & 1 & 2 & 1 \\
 *    12 & 1 & 1 & 2 \\
 *    13 & 3 & 1 & 2 \\
 *    14 & 3 & 3 & 2 \\
 *    15 & 1 & 3 & 2 \\
 *    16 & 2 & 1 & 3 \\
 *    17 & 3 & 2 & 3 \\
 *    18 & 2 & 3 & 3 \\
 *    19 & 1 & 2 & 3 \\
 *    20 & 2 & 1 & 2 \\
 *    21 & 3 & 2 & 2 \\
 *    22 & 2 & 3 & 2 \\
 *    23 & 1 & 2 & 2 \\
 *    24 & 2 & 2 & 1 \\
 *    25 & 2 & 2 & 3 \\
 *    26 & 2 & 2 & 2 \end{array}\f$
 *
 *  MBCN ORDER
 *
 *     7 - 18- 6
 *    /       /
 *   19  25  17
 *  /       /
 * 4 - 16- 6    
 *
 *     15- 22- 14
 *    /       /
 *   23  26  21
 *  /       /
 * 12- 20- 13    
 *
 *     3 - 10- 2
 *    /       /
 *   11  24  9    t s
 *  /       /     |/
 * 0 - 8 - 1      +-r
 */
class MESQUITE_EXPORT HexLagrangeShape : public MappingFunction3D
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
