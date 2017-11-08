// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#include <deal.II/base/point.h>
#include <deal.II/base/bounding_box.h>
#include <deal.II/base/mpi.h>
#include <deal.II/distributed/grid_tools.h>

DEAL_II_NAMESPACE_OPEN

#ifdef DEAL_II_WITH_MPI

namespace parallel
{
  namespace GridTools
  {
    template<int spacedim>
    std::vector< std::vector< BoundingBox<spacedim> > >
    exchange_local_bounding_boxes(const std::vector< BoundingBox<spacedim> > &local_bboxes,
                                  MPI_Comm                                    mpi_communicator)
    {
#ifndef DEAL_II_WITH_MPI
      (void)local_bboxes;
      (void)mpi_communicator;
      Assert(false, ExcMessage("parallel::GridTools::exchange_local_bounding_boxes() requires MPI."));
#else
      // Step 1: preparing data to be sent
      unsigned int n_bboxes = local_bboxes.size();
      // Dimension of the array to be exchanged (number of double)
      int n_local_data = 2*spacedim*n_bboxes;
      // data array stores each entry of each point describing the bounding boxes
      std::vector<double> loc_data_array(n_local_data);
      for (unsigned int i=0; i<n_bboxes; ++i)
        for (unsigned int d=0; d < spacedim; ++d)
          {
            // Extracting the coordinates of each boundary point
            loc_data_array[2*i*spacedim + d] = local_bboxes[i].get_boundary_points().first[d];
            loc_data_array[2*i*spacedim + spacedim + d] = local_bboxes[i].get_boundary_points().second[d];
          }

      // Step 2: exchanging the size of local data
      unsigned int n_procs = Utilities::MPI::n_mpi_processes(mpi_communicator);

      // Vector to store the size of loc_data_array for every process
      std::vector<int> size_all_data(n_procs);

      // Exchanging the number of bboxes
      MPI_Allgather(&n_local_data, 1, MPI_INT,
                    &(size_all_data[0]), 1, MPI_INT,
                    mpi_communicator);

      // Now computing the the displacement, relative to recvbuf,
      // at which to store the incoming data
      std::vector<int> rdispls(n_procs);
      rdispls[0] = 0;
      for (unsigned int i=1; i < n_procs; ++i)
        rdispls[i] = rdispls[i-1] + size_all_data[i-1];

      // Step 3: exchange the data and bounding boxes:
      // Allocating a vector to contain all the received data
      std::vector<double> data_array(rdispls.back() + size_all_data.back());

      MPI_Allgatherv(&(loc_data_array[0]), n_local_data, MPI_DOUBLE,
                     &(data_array[0]), &(size_all_data[0]),
                     &(rdispls[0]), MPI_DOUBLE, mpi_communicator);

      // Step 4: create the array of bboxes for output
      std::vector< std::vector< BoundingBox<spacedim> > > global_bboxes(n_procs);
      unsigned int begin_idx = 0;
      for (unsigned int i=0; i < n_procs; ++i)
        {
          // Number of local bounding boxes
          unsigned int n_bbox_i = size_all_data[i]/(spacedim*2);
          global_bboxes[i].resize(n_bbox_i);
          for (unsigned int bbox=0; bbox<n_bbox_i; ++bbox)
            {
              Point<spacedim> p1,p2; // boundary points for bbox
              for (unsigned int d=0; d<spacedim; ++d)
                {
                  p1[d] = data_array[ begin_idx + 2*bbox*spacedim + d];
                  p2[d] = data_array[ begin_idx + 2*bbox*spacedim + spacedim + d];
                }
              BoundingBox<spacedim> loc_bbox(std::make_pair(p1,p2));
              global_bboxes[i][bbox] = loc_bbox;
            }
          // Shifting the first index to the start of the next vector
          begin_idx += size_all_data[i];
        }
      return global_bboxes;
#endif // DEAL_II_WITH_MPI
    }
  }
}

// explicit instantiations
#include "grid_tools.inst"

#endif // DEAL_II_WITH_MPI
DEAL_II_NAMESPACE_CLOSE
