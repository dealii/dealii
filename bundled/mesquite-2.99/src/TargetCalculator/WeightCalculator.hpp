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


/** \file WeightCalculator.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_WEIGHT_CALCULATOR_HPP
#define MSQ_WEIGHT_CALCULATOR_HPP

#include "Mesquite.hpp"
#include "Sample.hpp"
#include <stddef.h>

namespace MESQUITE_NS {

class PatchData;
class MsqError;
class Mesh;
class MeshDomain;
class MeshDomainAssoc;
class Settings;

class MESQUITE_EXPORT WeightCalculator
{
public:

  virtual ~WeightCalculator();
  
   //!\brief Called at start of instruction queue processing
   //!
   //! Do any preliminary global initialization, consistency checking,
   //! etc.  Default implementation does nothing.
  virtual void initialize_queue( MeshDomainAssoc* mesh_and_domain,
                                 const Settings* settings,
                                 MsqError& err );


  /**\brief Get target metric weight
   *
   *\param pd      The current PatchData
   *\param element The index an element within the patch data.
   *\param sample  The sample point in the element.
   */
  virtual double get_weight( PatchData& pd, 
                             size_t element,
                             Sample sample,
                             MsqError& err ) = 0;
};



} // namespace Mesquite

#endif
