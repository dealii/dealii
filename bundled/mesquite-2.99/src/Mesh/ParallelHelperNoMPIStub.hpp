/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Sandia Corporation and Argonne National
    Laboratory.  Under the terms of Contract DE-AC04-94AL85000 
    with Sandia Corporation, the U.S. Government retains certain 
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
 
    diachin2@llnl.gov, djmelan@sandia.gov, mbrewer@sandia.gov, 
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov      
   
  ***************************************************************** */
/*!
  \file   ParallelHelperNoMPIStub.hpp
  \brief  

  Standin for a couple stubbed-out functions for MPI rank and size.

  \author Stephen Kennon
  \date   2012-01-11
*/

#ifndef Mesquite_ParallelHelperNoMPIStub_hpp 
#define Mesquite_ParallelHelperNoMPIStub_hpp

#include "Mesquite.hpp"

namespace MESQUITE_NS
{

  int get_parallel_rank();
  int get_parallel_size();
  double reduce_parallel_max(double value);
  void parallel_barrier();
  
} // namespace

#endif // Mesquite_ParallelHelperNoMPIStub_hpp
