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
  \file   ParallelHelperInterface.hpp
  \brief  

  \author Martin Isenburg
  \date   2008-03-25
*/

#ifndef Mesquite_ParallelHelperInterface_hpp 
#define Mesquite_ParallelHelperInterface_hpp

#include "ParallelMeshInterface.hpp"

namespace MESQUITE_NS
{
  /*! \class ParallelHelper \brief Mesquite::ParallelHelper is an abstract class
   *  which defines methods used by VertexMover and QualityAssessor to parallelize
   *  operations so the details and includes for MPI communication do not need to
   *  be exposed when running Mesquite in serial mode.
   */
  class MESQUITE_EXPORT ParallelHelper
  {
  public:
    // function called by application during set-up
    virtual void set_parallel_mesh(ParallelMesh* mesh) = 0;
    virtual void set_communicator(size_t comm) = 0;
    virtual void set_communication_model(int model, MsqError&) = 0;
    virtual void set_generate_random_numbers(int grn, MsqError&) = 0;

    virtual ~ParallelHelper() {}

  protected:
    friend class VertexMover;
    // functions called by VertexMover::loop_over_mesh()
    virtual void smoothing_init(MsqError&) = 0;
    virtual void compute_first_independent_set(std::vector<Mesh::VertexHandle>& fixed_vertices) = 0;
    virtual void communicate_first_independent_set(MsqError&) = 0;
    virtual bool compute_next_independent_set() = 0;
    virtual bool get_next_partition_boundary_vertex(Mesquite::Mesh::VertexHandle& vertex_handle) = 0;
    virtual void communicate_next_independent_set(MsqError&) = 0;
    virtual void smoothing_close(MsqError&) = 0;

  protected:
    friend class QualityAssessor;
    // functions called by QualityAssessor::loop_over_mesh()
    virtual int get_rank() const = 0;
    virtual int get_nprocs() const = 0;
    virtual bool is_our_element(Mesh::ElementHandle element_handle, MsqError&) const = 0;
    virtual bool is_our_vertex(Mesh::VertexHandle vertex_handle, MsqError&) const = 0;
    virtual void communicate_min_max_to_all(double* minimum, double* maximum, MsqError&) const = 0;
    virtual void communicate_min_max_to_zero(double* minimum, double* maximum, MsqError&) const = 0;
    virtual void communicate_sums_to_zero(size_t* freeElementCount, int* invertedElementCount, size_t* elementCount, int* invertedSampleCount, size_t* sampleCount, long unsigned int* count, long unsigned int* invalid, double* sum, double *sqrSum, MsqError&) const = 0;
    virtual void communicate_power_sum_to_zero(double* pMean, MsqError&) const = 0;
    virtual void communicate_histogram_to_zero(std::vector<int> &histogram, MsqError&) const = 0;
    virtual void communicate_all_true( bool& value, MsqError& err ) const = 0;
    virtual void communicate_any_true( bool& value, MsqError& err ) const = 0;
  };
} // namespace

#endif
