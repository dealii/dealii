/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2009 Sandia National Laboratories.  Developed at the
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

    (2009) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file TetDihedralWeight.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_TET_DIHEDRAL_WEIGHT_HPP
#define MSQ_TET_DIHEDRAL_WEIGHT_HPP

#include "Mesquite.hpp"
#include "WeightCalculator.hpp"

namespace MESQUITE_NS {

class ReferenceMesh;

/**\brief Calculate weight as a function of maximum dihedral angle 
 *
 * \f$ c_k = 1/(1 + \exp(-0.4394*(d_max - mCutoff))) \f*
 */
class TetDihedralWeight : public WeightCalculator
{
public:

  TetDihedralWeight( double cutoff = 135, double a = 0.4395 )
    : refMesh(0), mCutoff(cutoff), mA(a)
    {}

  TetDihedralWeight( ReferenceMesh* mesh, double cutoff = 135, double a = 0.4395 )
    : refMesh(mesh), mCutoff(cutoff), mA(a)
    {}

  /**\brief Get target metric weight
   *
   *\param pd      The current PatchData
   *\param element The index an element within the patch data.
   *\param sample  The sample point in the element.
   */
  virtual double get_weight( PatchData& pd, 
                             size_t element,
                             Sample sample,
                             MsqError& err );

  private:
    ReferenceMesh* refMesh;
    double mCutoff, mA;
};

} // namespace MESQUITE_NS

#endif
