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
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov,
    kraftche@cae.wisc.edu
   
  ***************************************************************** */
/*!
  \file   ParallelMeshInterface.hpp
  \brief  This file contains the parallel Mesquite mesh interface.
          Many users will want to implement a concrete class derived from
          the ParallelMeshInterface class to access their mesh.

*/
#ifndef MESQUITE_PARALLEL_MESH_INTERFACE_HPP
#define MESQUITE_PARALLEL_MESH_INTERFACE_HPP

#include "MeshInterface.hpp"

namespace MESQUITE_NS
{
  class ParallelHelper;
    
  /*! \class ParallelMesh \brief Mesquite::ParallelMesh is an abstract class
   *  which defines required methods required for using Mesquite in parallel.
   *  It derives from the Mesquite::Mesh interface so the user must provide
   *  implementations of the pure virtual methods in both Mesquite::Mesh as
   *  well as those defined here.
   */
  class MESQUITE_EXPORT ParallelMesh : virtual public Mesh
  {
  public:

    /*! Get global ids for given vertices.
     */
    virtual void vertices_get_global_id ( const VertexHandle vert_array[],
                                          size_t global_id[],
					  size_t num_vtx,
					  MsqError& err) = 0;     
     
    /*! Get processor ids for given vertices.
     */
    virtual void vertices_get_processor_id ( const VertexHandle vert_array[],
                                             int proc_id[],                                          
					     size_t num_vtx,
					     MsqError& err) = 0;

    /*! Set parallel helper
     */
    virtual void set_parallel_helper(ParallelHelper* helper) {
      this->helper = helper;
    }

    /*! Get parallel helper
     */
    virtual ParallelHelper* get_parallel_helper() {
      return helper;
    }

  protected:
    ParallelHelper* helper;
      //! Don't allow a ParallelMesh to be deleted directly.
    virtual ~ParallelMesh()
      {}
  };
} // namespace

#endif
