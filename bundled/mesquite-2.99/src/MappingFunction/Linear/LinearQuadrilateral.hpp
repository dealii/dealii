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

#ifndef MSQ_LINEAR_QUADRILATERAL_HPP
#define MSQ_LINEAR_QUADRILATERAL_HPP

/** \file LinearQuadrilateral.hpp
 *  \brief Linear mapping funtion for Quadrilateral elements.
 *  \author Jason Kraftcheck
 */
 
#include "Mesquite.hpp"
#include "MappingFunction.hpp"

namespace MESQUITE_NS {

/**\brief Linear shape function for quadrilateral elements
 *
 * This class implements the MappingFunction interface, providing
 * a Linear shape function for quadrilateral elements.
 *
 * \f$\vec{x}(\xi,\eta) = \sum_{i=0}^{3} N_i(\xi,\eta) \vec{x_i}\f$
 * 
 * \f$N_0(\xi,\eta) = (1-\xi)(1-\eta)\f$
 * 
 * \f$N_1(\xi,\eta) =    \xi (1-\eta)\f$
 * 
 * \f$N_2(\xi,\eta) =    \xi    \eta \f$
 * 
 * \f$N_3(\xi,\eta) = (1-\xi)   \eta\f$
 */
class MESQUITE_EXPORT LinearQuadrilateral : public MappingFunction2D
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
  
  static 
  void coefficients( Sample location,
                     NodeSet nodeset,
                     double* coeff_out,
                     size_t* indices_out,
                     size_t& num_coeff_out );
  
  virtual 
  void derivatives( Sample location,
                    NodeSet nodeset,
                    size_t* vertex_indices_out,
                    MsqVector<2>* d_coeff_d_xi_out,
                    size_t& num_vtx,
                    MsqError& err ) const;

  static
  void derivatives( Sample location,
                    NodeSet nodeset,
                    size_t* vertex_indices_out,
                    MsqVector<2>* d_coeff_d_xi_out,
                    size_t& num_vtx );

  virtual
  void ideal( Sample location, 
              MsqMatrix<3,2>& jacobian_out,
              MsqError& err ) const;
};

} // namespace Mesquite

#endif
