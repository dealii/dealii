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
  \file   ParallelHelper.hpp
  \brief  

  Implements ParallelHelper Class 

  \author Martin Isenburg
  \date   2008-03-04
*/

#ifndef Mesquite_ParallelHelper_hpp 
#define Mesquite_ParallelHelper_hpp

#include "ParallelHelperInterface.hpp"

#include <vector>
#include <map>

namespace MESQUITE_NS
{
  typedef struct VertexIdMapKey {
    size_t glob_id;
    int proc_id;
  } VertexIdMapKey;
  
  struct VertexIdLessFunc {
    bool operator()( const VertexIdMapKey &that1, const VertexIdMapKey& that2 ) const
    {
      return ( (that1.proc_id < that2.proc_id) || ((that1.proc_id==that2.proc_id)&&(that1.glob_id<that2.glob_id)) );
    }
  };
  
  typedef std::map<VertexIdMapKey,int,VertexIdLessFunc> VertexIdMap;

  int get_parallel_rank();
  int get_parallel_size();
  double reduce_parallel_max(double value);
  void parallel_barrier();

  class ParallelHelperImpl : public ParallelHelper
  {
  public:

    enum CommunicationModel{ TrulyNonBlocking = 0, TrulyNonBlockingAvoidAllReduce,
			     NonBlocking, NonBlockingAvoidAllReduce,
			     Blocking, BlockingAvoidAllReduce };

    ParallelHelperImpl();
    ~ParallelHelperImpl();

    // function called by application during set-up
    void set_parallel_mesh(ParallelMesh* mesh);
    void set_communicator(size_t comm);
    void set_communicator(const void* comm) 
      { set_communicator( reinterpret_cast<size_t>(comm) ); }
    void set_communication_model(int model, MsqError& err);
    void set_generate_random_numbers(int grn, MsqError& err);

  protected:
    friend class VertexMover;
    // functions called by VertexMover::loop_over_mesh()
    void smoothing_init( MsqError& err );
    void compute_first_independent_set(std::vector<Mesh::VertexHandle>& fixed_vertices);
    void communicate_first_independent_set(MsqError& err);
    bool compute_next_independent_set();
    bool get_next_partition_boundary_vertex(Mesquite::Mesh::VertexHandle& vertex_handle);
    void communicate_next_independent_set(MsqError& err);
    void smoothing_close(MsqError& err);

  protected:
    friend class QualityAssessor;
    // functions called by QualityAssessor::loop_over_mesh()
    int get_rank() const;
    int get_nprocs() const;
    bool is_our_element(Mesh::ElementHandle element_handle, MsqError& err) const;
    bool is_our_vertex(Mesh::VertexHandle vertex_handle, MsqError& err) const;
    void communicate_min_max_to_all(double* minimum, double* maximum, MsqError& ) const;
    void communicate_min_max_to_zero(double* minimum, double* maximum, MsqError&) const;
    void communicate_sums_to_zero(size_t* freeElementCount, int* invertedElementCount, size_t* elementCount, int* invertedSampleCount, size_t* sampleCount, long unsigned int* count, long unsigned int* invalid, double* sum, double *sqrSum, MsqError&) const;
    void communicate_power_sum_to_zero(double* pMean, MsqError&) const;
    void communicate_histogram_to_zero(std::vector<int> &histogram, MsqError&) const;
    void communicate_any_true( bool& value, MsqError& err ) const;
    void communicate_all_true( bool& value, MsqError& err ) const;

  private:
    ParallelMesh* mesh;

    size_t communicator;
    int communication_model;

    int rank;
    int nprocs;

    // variables for VertexMover::loop_over_mesh()
    int generate_random_numbers;
    std::vector<Mesquite::Mesh::VertexHandle> vertices;
    int num_vertex;
    std::vector<char> vtx_in_partition_boundary;
    int num_vtx_partition_boundary;
    int num_vtx_partition_boundary_local;
    int num_vtx_partition_boundary_remote;
    std::vector<Mesquite::Mesh::VertexHandle> part_vertices;
    std::vector<int> part_proc_owner;
    std::vector<size_t> part_gid;
    std::vector<int> part_smoothed_flag;
    std::vector<double> part_rand_number;
    int num_exportVtx;
    std::vector<size_t> exportVtxGIDs;
    std::vector<int> exportVtxLIDs;
    std::vector<int> exportProc;
    std::vector<bool> in_independent_set;
    VertexIdMap vid_map;
    int total_num_vertices_to_smooth;
    int total_num_vertices_to_recv;
    std::vector<int> neighbourProcSend;
    std::vector<int> neighbourProcRecv;
    std::vector<int> neighbourProcSendRemain;
    std::vector<int> neighbourProcRecvRemain;
    int num_already_smoothed_vertices;
    int num_already_recv_vertices;
    std::vector< std::vector<int> > vtx_off_proc_list;
    std::vector<int> neighbourProc;
    int iteration;
    int global_work_remains;
    int next_vtx_partition_boundary;
    /* for exchanging unused ghost node information */
    int unghost_num_vtx;
    std::vector<Mesquite::Mesh::VertexHandle> unghost_vertices;
    int unghost_num_procs;
    std::vector<int> unghost_procs;
    std::vector<int> unghost_procs_num_vtx;
    std::vector<int> unghost_procs_offset;
    int update_num_vtx;
    std::vector<size_t> update_gid;
    int update_num_procs;
    std::vector<int> update_procs;
    std::vector<int> update_procs_num_vtx;
    std::vector<int> update_procs_offset;

    // functions for VertexMover::loop_over_mesh()
    void compute_independent_set();
    int comm_smoothed_vtx_b(MsqError& err);
    int comm_smoothed_vtx_b_no_all(MsqError& err);
    int comm_smoothed_vtx_nb(MsqError& err);
    int comm_smoothed_vtx_nb_no_all(MsqError& err);
    int comm_smoothed_vtx_tnb(MsqError& err);
    int comm_smoothed_vtx_tnb_no_all(MsqError& err);
  };
  
} // namespace
#endif // Mesquite_ParallelHelper_hpp
