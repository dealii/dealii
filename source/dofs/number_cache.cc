// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/mpi.h>

#include <deal.II/dofs/number_cache.h>

#include <numeric>


DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace DoFHandlerImplementation
  {
    NumberCache::NumberCache()
      : n_global_dofs(0)
      , n_locally_owned_dofs(0)
    {}



    NumberCache::NumberCache(const types::global_dof_index n_global_dofs)
      : n_global_dofs(n_global_dofs)
      , n_locally_owned_dofs(n_global_dofs)
      , locally_owned_dofs(complete_index_set(n_global_dofs))
      , n_locally_owned_dofs_per_processor(1, n_global_dofs)
      , locally_owned_dofs_per_processor(1, complete_index_set(n_global_dofs))
    {}



    NumberCache::NumberCache(
      const std::vector<IndexSet> &locally_owned_dofs_per_processor,
      const unsigned int           my_rank)
      : locally_owned_dofs_per_processor(locally_owned_dofs_per_processor)
    {
      const unsigned int n_procs = locally_owned_dofs_per_processor.size();

      // compress IndexSet representation before using it for anything else
      for (unsigned int p = 0; p < n_procs; ++p)
        locally_owned_dofs_per_processor[p].compress();

      n_locally_owned_dofs_per_processor.resize(n_procs);
      for (unsigned int p = 0; p < n_procs; ++p)
        n_locally_owned_dofs_per_processor[p] =
          locally_owned_dofs_per_processor[p].n_elements();
      n_locally_owned_dofs = n_locally_owned_dofs_per_processor[my_rank];
      locally_owned_dofs   = locally_owned_dofs_per_processor[my_rank];


      n_global_dofs =
        std::accumulate(n_locally_owned_dofs_per_processor.begin(),
                        n_locally_owned_dofs_per_processor.end(),
                        types::global_dof_index(0));
    }



    void
    NumberCache::clear()
    {
      n_global_dofs        = 0;
      n_locally_owned_dofs = 0;
      locally_owned_dofs.clear();
      n_locally_owned_dofs_per_processor.clear();
      locally_owned_dofs_per_processor.clear();
    }



    std::vector<types::global_dof_index>
    NumberCache::get_n_locally_owned_dofs_per_processor(
      const MPI_Comm mpi_communicator) const
    {
      const unsigned int n_procs =
        Utilities::MPI::job_supports_mpi() ?
          Utilities::MPI::n_mpi_processes(mpi_communicator) :
          1;
      if (n_global_dofs == 0)
        return std::vector<types::global_dof_index>();
      else if (n_locally_owned_dofs_per_processor.empty() == false)
        {
          AssertDimension(n_locally_owned_dofs_per_processor.size(), n_procs);
          return n_locally_owned_dofs_per_processor;
        }
      else
        {
          std::vector<types::global_dof_index> result(n_procs,
                                                      n_locally_owned_dofs);
#ifdef DEAL_II_WITH_MPI
          if (n_procs > 1)
            MPI_Allgather(DEAL_II_MPI_CONST_CAST(&n_locally_owned_dofs),
                          1,
                          DEAL_II_DOF_INDEX_MPI_TYPE,
                          result.data(),
                          1,
                          DEAL_II_DOF_INDEX_MPI_TYPE,
                          mpi_communicator);
#endif
          return result;
        }
    }



    std::vector<IndexSet>
    NumberCache::get_locally_owned_dofs_per_processor(
      const MPI_Comm mpi_communicator) const
    {
      AssertDimension(locally_owned_dofs.size(), n_global_dofs);
      const unsigned int n_procs =
        Utilities::MPI::job_supports_mpi() ?
          Utilities::MPI::n_mpi_processes(mpi_communicator) :
          1;
      if (n_global_dofs == 0)
        return std::vector<IndexSet>();
      else if (locally_owned_dofs_per_processor.empty() == false)
        {
          AssertDimension(locally_owned_dofs_per_processor.size(), n_procs);
          return locally_owned_dofs_per_processor;
        }
      else
        {
          std::vector<IndexSet> locally_owned_dofs_per_processor(
            n_procs, locally_owned_dofs);

#ifdef DEAL_II_WITH_MPI
          if (n_procs > 1)
            {
              // this step is substantially more complicated because indices
              // might be distributed arbitrarily among the processors. Here we
              // have to serialize the IndexSet objects and shop them across the
              // network.
              std::vector<char> my_data;
              {
#  ifdef DEAL_II_WITH_ZLIB

                boost::iostreams::filtering_ostream out;
                out.push(boost::iostreams::gzip_compressor(
                  boost::iostreams::gzip_params(
                    boost::iostreams::gzip::best_compression)));
                out.push(boost::iostreams::back_inserter(my_data));

                boost::archive::binary_oarchive archive(out);

                archive << locally_owned_dofs;
                out.flush();
#  else
                std::ostringstream              out;
                boost::archive::binary_oarchive archive(out);
                archive << locally_owned_dofs;
                const std::string &s = out.str();
                my_data.reserve(s.size());
                my_data.assign(s.begin(), s.end());
#  endif
              }

              // determine maximum size of IndexSet
              const unsigned int max_size =
                Utilities::MPI::max(my_data.size(), mpi_communicator);

              // as the MPI_Allgather call will be reading max_size elements,
              // and as this may be past the end of my_data, we need to increase
              // the size of the local buffer. This is filled with zeros.
              my_data.resize(max_size);

              std::vector<char> buffer(max_size * n_procs);
              const int         ierr = MPI_Allgather(my_data.data(),
                                             max_size,
                                             MPI_BYTE,
                                             buffer.data(),
                                             max_size,
                                             MPI_BYTE,
                                             mpi_communicator);
              AssertThrowMPI(ierr);

              for (unsigned int i = 0; i < n_procs; ++i)
                if (i == Utilities::MPI::this_mpi_process(mpi_communicator))
                  locally_owned_dofs_per_processor[i] = locally_owned_dofs;
                else
                  {
                    // copy the data previously received into a stringstream
                    // object and then read the IndexSet from it
                    std::string decompressed_buffer;

                    // first decompress the buffer
                    {
#  ifdef DEAL_II_WITH_ZLIB

                      boost::iostreams::filtering_ostream decompressing_stream;
                      decompressing_stream.push(
                        boost::iostreams::gzip_decompressor());
                      decompressing_stream.push(
                        boost::iostreams::back_inserter(decompressed_buffer));

                      decompressing_stream.write(&buffer[i * max_size],
                                                 max_size);
#  else
                      decompressed_buffer.assign(&buffer[i * max_size],
                                                 max_size);
#  endif
                    }

                    // then restore the object from the buffer
                    std::istringstream              in(decompressed_buffer);
                    boost::archive::binary_iarchive archive(in);

                    archive >> locally_owned_dofs_per_processor[i];
                  }
            }
#endif
          return locally_owned_dofs_per_processor;
        }
    }


    std::size_t
    NumberCache::memory_consumption() const
    {
      return MemoryConsumption::memory_consumption(n_global_dofs) +
             MemoryConsumption::memory_consumption(n_locally_owned_dofs) +
             MemoryConsumption::memory_consumption(locally_owned_dofs) +
             MemoryConsumption::memory_consumption(
               n_locally_owned_dofs_per_processor) +
             MemoryConsumption::memory_consumption(
               locally_owned_dofs_per_processor);
    }
  } // namespace DoFHandlerImplementation
} // namespace internal

DEAL_II_NAMESPACE_CLOSE
