// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2019 by the deal.II authors
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

#include <deal.II/base/partitioner.h>
#include <deal.II/base/partitioner.templates.h>

DEAL_II_NAMESPACE_OPEN



namespace Utilities
{
  namespace MPI
  {
    void
    Partitioner::initialize_import_indices_plain_dev() const
    {
      const unsigned int n_import_targets = import_targets_data.size();
      import_indices_plain_dev.reserve(n_import_targets);
      for (unsigned int i = 0; i < n_import_targets; i++)
        {
          // Expand the indices on the host
          std::vector<std::pair<unsigned int, unsigned int>>::const_iterator
            my_imports = import_indices_data.begin() +
                         import_indices_chunks_by_rank_data[i],
            end_my_imports = import_indices_data.begin() +
                             import_indices_chunks_by_rank_data[i + 1];
          std::vector<unsigned int> import_indices_plain_host;
          for (; my_imports != end_my_imports; ++my_imports)
            {
              const unsigned int chunk_size =
                my_imports->second - my_imports->first;
              for (unsigned int j = 0; j < chunk_size; ++j)
                import_indices_plain_host.push_back(my_imports->first + j);
            }

          // Move the indices to the device
          import_indices_plain_dev.emplace_back(std::make_pair(
            std::unique_ptr<unsigned int[], void (*)(unsigned int *)>(
              nullptr, Utilities::CUDA::delete_device_data<unsigned int>),
            import_indices_plain_host.size()));

          import_indices_plain_dev[i].first.reset(
            Utilities::CUDA::allocate_device_data<unsigned int>(
              import_indices_plain_dev[i].second));
          Utilities::CUDA::copy_to_dev(import_indices_plain_host,
                                       import_indices_plain_dev[i].first.get());
        }
    }
  } // namespace MPI
} // namespace Utilities

// explicit instantiations from .templates.h file
#include "partitioner.cuda.inst"

DEAL_II_NAMESPACE_CLOSE
