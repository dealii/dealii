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

/** \file LinearPyramid.hpp
 *  \brief mapping function for linear prism
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_LINEAR_PYRAMID_HPP
#define MSQ_LINEAR_PYRAMID_HPP

#include "Mesquite.hpp"
#include "MappingFunction.hpp"

namespace MESQUITE_NS {

/**\brief Linear mapping function for a pyramid element
 *
 * \f$\vec{x}(\vec{\xi})=\sum_{i=0}^5 N_i(\vec{\xi})\vec{x_i}\f$
 * - \f$ N_0(\vec{\xi})=(1-\xi)(1-\eta)(1-\zeta) \f$
 * - \f$ N_1(\vec{\xi})=   \xi (1-\eta)(1-\zeta) \f$
 * - \f$ N_2(\vec{\xi})=   \xi    \eta (1-\zeta) \f$
 * - \f$ N_3(\vec{\xi})=(1-\xi)   \eta (1-\zeta) \f$
 * - \f$ N_4(\vec{\xi})=\zeta \f$
 *
 * Where the mid-edge coordinates are:
 * - \f$ \vec{\xi}_{e0}=( \frac{1}{2}, 0,           0) \f$
 * - \f$ \vec{\xi}_{e1}=( 1,           \frac{1}{2}, 0) \f$
 * - \f$ \vec{\xi}_{e2}=( \frac{1}{2}, 1,           0) \f$
 * - \f$ \vec{\xi}_{e3}=( 0,           \frac{1}{2}, 0) \f$
 * - \f$ \vec{\xi}_{e4}=( 0,           0,           \frac{1}{2}) \f$
 * - \f$ \vec{\xi}_{e5}=( 1,           0,           \frac{1}{2}) \f$
 * - \f$ \vec{\xi}_{e6}=( 1,           1,           \frac{1}{2}) \f$
 * - \f$ \vec{\xi}_{e7}=( 0,           1,           \frac{1}{2}) \f$
 *
 * and mid-face the coordinates are:
 * - \f$ \vec{\xi}_{f0}=( \frac{1}{2}, 0,           \frac{1}{2}) \f$
 * - \f$ \vec{\xi}_{f1}=( 1,           \frac{1}{2}, \frac{1}{2}) \f$
 * - \f$ \vec{\xi}_{f2}=( \frac{1}{2}, 1,           \frac{1}{2}) \f$
 * - \f$ \vec{\xi}_{f3}=( 0,           \frac{1}{2}, \frac{1}{2}) \f$
 * - \f$ \vec{\xi}_{f4}=( \frac{1}{2}, \frac{1}{2}, 0)           \f$
 *
 * and the mid-element coorindates are
 * - \f$ \vec{\xi}_m=( \frac{1}{2}, \frac{1}{2}, \frac{1}{2}) \f$
 */ 
class MESQUITE_EXPORT LinearPyramid : public MappingFunction3D
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
              MsqMatrix<3,3>& J,
              MsqError& err ) const;
};



} // namespace Mesquite

#endif
