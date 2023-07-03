// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2023 by the deal.II authors
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

#include <deal.II/base/config.h>

#include <deal.II/base/geometry_info.h>

#include <deal.II/grid/data_serializer.h>
#include <deal.II/grid/tria_accessor.h>

#include <utility>
#include <vector>

DEAL_II_NAMESPACE_OPEN

template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
CellAttachedDataSerializer<dim, spacedim>::CellAttachedDataSerializer()
  : variable_size_data_stored(false)
{}


template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
void CellAttachedDataSerializer<dim, spacedim>::pack_data(
  const std::vector<cell_relation_t> &cell_relations,
  const std::vector<
    typename internal::CellAttachedData<dim, spacedim>::pack_callback_t>
    &pack_callbacks_fixed,
  const std::vector<
    typename internal::CellAttachedData<dim, spacedim>::pack_callback_t>
    &             pack_callbacks_variable,
  const MPI_Comm &mpi_communicator)
{
  Assert(src_data_fixed.size() == 0,
         ExcMessage("Previously packed data has not been released yet!"));
  Assert(src_sizes_variable.size() == 0, ExcInternalError());

  const unsigned int n_callbacks_fixed    = pack_callbacks_fixed.size();
  const unsigned int n_callbacks_variable = pack_callbacks_variable.size();

  // Store information that we packed variable size data in
  // a member variable for later.
  variable_size_data_stored = (n_callbacks_variable > 0);

  // If variable transfer is scheduled, we will store the data size that
  // each variable size callback function writes in this auxiliary
  // container. The information will be stored by each cell in this vector
  // temporarily.
  std::vector<unsigned int> cell_sizes_variable_cumulative(
    n_callbacks_variable);

  // Prepare the buffer structure, in which each callback function will
  // store its data for each active cell.
  // The outmost shell in this container construct corresponds to the
  // data packed per cell. The next layer resembles the data that
  // each callback function packs on the corresponding cell. These
  // buffers are chains of chars stored in an std::vector<char>.
  // A visualisation of the data structure:
  /* clang-format off */
      // |             cell_1                | |             cell_2                | ...
      // ||  callback_1  ||  callback_2  |...| ||  callback_1  ||  callback_2  |...| ...
      // |||char|char|...|||char|char|...|...| |||char|char|...|||char|char|...|...| ...
  /* clang-format on */
  std::vector<std::vector<std::vector<char>>> packed_fixed_size_data(
    cell_relations.size());
  std::vector<std::vector<std::vector<char>>> packed_variable_size_data(
    variable_size_data_stored ? cell_relations.size() : 0);

  //
  // --------- Pack data for fixed and variable size transfer ---------
  //
  // Iterate over all cells, call all callback functions on each cell,
  // and store their data in the corresponding buffer scope.
  {
    auto cell_rel_it           = cell_relations.cbegin();
    auto data_cell_fixed_it    = packed_fixed_size_data.begin();
    auto data_cell_variable_it = packed_variable_size_data.begin();
    for (; cell_rel_it != cell_relations.cend(); ++cell_rel_it)
      {
        const auto &dealii_cell = cell_rel_it->first;
        const auto &cell_status = cell_rel_it->second;

        // Assertions about the tree structure.
        switch (cell_status)
          {
            case CellStatus::cell_will_persist:
            case CellStatus::cell_will_be_refined:
              // double check the condition that we will only ever attach
              // data to active cells when we get here
              Assert(dealii_cell->is_active(), ExcInternalError());
              break;

            case CellStatus::children_will_be_coarsened:
              // double check the condition that we will only ever attach
              // data to cells with children when we get here. however, we
              // can only tolerate one level of coarsening at a time, so
              // check that the children are all active
              Assert(dealii_cell->is_active() == false, ExcInternalError());
              for (unsigned int c = 0;
                   c < GeometryInfo<dim>::max_children_per_cell;
                   ++c)
                Assert(dealii_cell->child(c)->is_active(), ExcInternalError());
              break;

            case CellStatus::cell_invalid:
              // do nothing on invalid cells
              break;

            default:
              Assert(false, ExcInternalError());
              break;
          }

        // Reserve memory corresponding to the number of callback
        // functions that will be called.
        // If variable size transfer is scheduled, we need to leave
        // room for an array that holds information about how many
        // bytes each of the variable size callback functions will
        // write.
        // On cells flagged with CellStatus::cell_invalid, only its CellStatus
        // will be stored.
        const unsigned int n_fixed_size_data_sets_on_cell =
          1 + ((cell_status == CellStatus::cell_invalid) ?
                 0 :
                 ((variable_size_data_stored ? 1 : 0) + n_callbacks_fixed));
        data_cell_fixed_it->resize(n_fixed_size_data_sets_on_cell);

        // We continue with packing all data on this specific cell.
        auto data_fixed_it = data_cell_fixed_it->begin();

        // First, we pack the CellStatus information.
        // to get consistent data sizes on each cell for the fixed size
        // transfer, we won't allow compression
        *data_fixed_it =
          Utilities::pack(cell_status, /*allow_compression=*/false);
        ++data_fixed_it;

        // Proceed with all registered callback functions.
        // Skip cells with the CellStatus::cell_invalid flag.
        if (cell_status != CellStatus::cell_invalid)
          {
            // Pack fixed size data.
            for (auto callback_it = pack_callbacks_fixed.cbegin();
                 callback_it != pack_callbacks_fixed.cend();
                 ++callback_it, ++data_fixed_it)
              {
                *data_fixed_it = (*callback_it)(dealii_cell, cell_status);
              }

            // Pack variable size data.
            // If we store variable size data, we need to transfer
            // the sizes of each corresponding callback function
            // via fixed size transfer as well.
            if (variable_size_data_stored)
              {
                const unsigned int n_variable_size_data_sets_on_cell =
                  ((cell_status == CellStatus::cell_invalid) ?
                     0 :
                     n_callbacks_variable);
                data_cell_variable_it->resize(
                  n_variable_size_data_sets_on_cell);

                auto callback_it       = pack_callbacks_variable.cbegin();
                auto data_variable_it  = data_cell_variable_it->begin();
                auto sizes_variable_it = cell_sizes_variable_cumulative.begin();
                for (; callback_it != pack_callbacks_variable.cend();
                     ++callback_it, ++data_variable_it, ++sizes_variable_it)
                  {
                    *data_variable_it =
                      (*callback_it)(dealii_cell, cell_status);

                    // Store data sizes for each callback function first.
                    // Make it cumulative below.
                    *sizes_variable_it = data_variable_it->size();
                  }

                // Turn size vector into its cumulative representation.
                std::partial_sum(cell_sizes_variable_cumulative.begin(),
                                 cell_sizes_variable_cumulative.end(),
                                 cell_sizes_variable_cumulative.begin());

                // Serialize cumulative variable size vector
                // value-by-value. This way we can circumvent the overhead
                // of storing the container object as a whole, since we
                // know its size by the number of registered callback
                // functions.
                data_fixed_it->resize(n_callbacks_variable *
                                      sizeof(unsigned int));
                for (unsigned int i = 0; i < n_callbacks_variable; ++i)
                  std::memcpy(&(data_fixed_it->at(i * sizeof(unsigned int))),
                              &(cell_sizes_variable_cumulative.at(i)),
                              sizeof(unsigned int));

                ++data_fixed_it;
              }

            // Double check that we packed everything we wanted
            // in the fixed size buffers.
            Assert(data_fixed_it == data_cell_fixed_it->end(),
                   ExcInternalError());
          }

        ++data_cell_fixed_it;

        // Increment the variable size data iterator
        // only if we actually pack this kind of data
        // to avoid getting out of bounds.
        if (variable_size_data_stored)
          ++data_cell_variable_it;
      } // loop over cell_relations
  }

  //
  // ----------- Gather data sizes for fixed size transfer ------------
  //
  // Generate a vector which stores the sizes of each callback function,
  // including the packed CellStatus transfer.
  // Find the very first cell that we wrote to with all callback
  // functions (i.e. a cell that was not flagged with CellStatus::cell_invalid)
  // and store the sizes of each buffer.
  //
  // To deal with the case that at least one of the processors does not
  // own any cell at all, we will exchange the information about the data
  // sizes among them later. The code in between is still well-defined,
  // since the following loops will be skipped.
  std::vector<unsigned int> local_sizes_fixed(
    1 + n_callbacks_fixed + (variable_size_data_stored ? 1 : 0));
  for (const auto &data_cell : packed_fixed_size_data)
    {
      if (data_cell.size() == local_sizes_fixed.size())
        {
          auto sizes_fixed_it = local_sizes_fixed.begin();
          auto data_fixed_it  = data_cell.cbegin();
          for (; data_fixed_it != data_cell.cend();
               ++data_fixed_it, ++sizes_fixed_it)
            {
              *sizes_fixed_it = data_fixed_it->size();
            }

          break;
        }
    }

  // Check if all cells have valid sizes.
  for (auto data_cell_fixed_it = packed_fixed_size_data.cbegin();
       data_cell_fixed_it != packed_fixed_size_data.cend();
       ++data_cell_fixed_it)
    {
      Assert((data_cell_fixed_it->size() == 1) ||
               (data_cell_fixed_it->size() == local_sizes_fixed.size()),
             ExcInternalError());
    }

  // Share information about the packed data sizes
  // of all callback functions across all processors, in case one
  // of them does not own any cells at all.
  std::vector<unsigned int> global_sizes_fixed(local_sizes_fixed.size());
  Utilities::MPI::max(local_sizes_fixed, mpi_communicator, global_sizes_fixed);

  // Construct cumulative sizes, since this is the only information
  // we need from now on.
  sizes_fixed_cumulative.resize(global_sizes_fixed.size());
  std::partial_sum(global_sizes_fixed.begin(),
                   global_sizes_fixed.end(),
                   sizes_fixed_cumulative.begin());

  //
  // ---------- Gather data sizes for variable size transfer ----------
  //
  if (variable_size_data_stored)
    {
      src_sizes_variable.reserve(packed_variable_size_data.size());
      for (const auto &data_cell : packed_variable_size_data)
        {
          int variable_data_size_on_cell = 0;

          for (const auto &data : data_cell)
            variable_data_size_on_cell += data.size();

          src_sizes_variable.push_back(variable_data_size_on_cell);
        }
    }

  //
  // ------------------------ Build buffers ---------------------------
  //
  const unsigned int expected_size_fixed =
    cell_relations.size() * sizes_fixed_cumulative.back();
  const unsigned int expected_size_variable =
    std::accumulate(src_sizes_variable.begin(),
                    src_sizes_variable.end(),
                    std::vector<int>::size_type(0));

  // Move every piece of packed fixed size data into the consecutive
  // buffer.
  src_data_fixed.reserve(expected_size_fixed);
  for (const auto &data_cell_fixed : packed_fixed_size_data)
    {
      // Move every fraction of packed data into the buffer
      // reserved for this particular cell.
      for (const auto &data_fixed : data_cell_fixed)
        std::move(data_fixed.begin(),
                  data_fixed.end(),
                  std::back_inserter(src_data_fixed));

      // If we only packed the CellStatus information
      // (i.e. encountered a cell flagged CellStatus::cell_invalid),
      // fill the remaining space with invalid entries.
      // We can skip this if there is nothing else to pack.
      if ((data_cell_fixed.size() == 1) && (sizes_fixed_cumulative.size() > 1))
        {
          const std::size_t bytes_skipped =
            sizes_fixed_cumulative.back() - sizes_fixed_cumulative.front();

          src_data_fixed.insert(src_data_fixed.end(),
                                bytes_skipped,
                                static_cast<char>(-1)); // invalid_char
        }
    }

  // Move every piece of packed variable size data into the consecutive
  // buffer.
  if (variable_size_data_stored)
    {
      src_data_variable.reserve(expected_size_variable);
      for (const auto &data_cell : packed_variable_size_data)
        {
          // Move every fraction of packed data into the buffer
          // reserved for this particular cell.
          for (const auto &data : data_cell)
            std::move(data.begin(),
                      data.end(),
                      std::back_inserter(src_data_variable));
        }
    }

  // Double check that we packed everything correctly.
  Assert(src_data_fixed.size() == expected_size_fixed, ExcInternalError());
  Assert(src_data_variable.size() == expected_size_variable,
         ExcInternalError());
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
void CellAttachedDataSerializer<dim, spacedim>::unpack_cell_status(
  std::vector<
    typename CellAttachedDataSerializer<dim, spacedim>::cell_relation_t>
    &cell_relations) const
{
  Assert(sizes_fixed_cumulative.size() > 0,
         ExcMessage("No data has been packed!"));
  if (cell_relations.size() > 0)
    {
      Assert(dest_data_fixed.size() > 0,
             ExcMessage("No data has been received!"));
    }

  // Size of CellStatus object that will be unpacked on each cell.
  const unsigned int size = sizes_fixed_cumulative.front();

  // Iterate over all cells and overwrite the CellStatus
  // information from the transferred data.
  // Proceed buffer iterator position to next cell after
  // each iteration.
  auto cell_rel_it   = cell_relations.begin();
  auto dest_fixed_it = dest_data_fixed.cbegin();
  for (; cell_rel_it != cell_relations.end();
       ++cell_rel_it, dest_fixed_it += sizes_fixed_cumulative.back())
    {
      cell_rel_it->second = // cell_status
        Utilities::unpack<CellStatus>(dest_fixed_it,
                                      dest_fixed_it + size,
                                      /*allow_compression=*/false);
    }
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
void CellAttachedDataSerializer<dim, spacedim>::unpack_data(
  const std::vector<
    typename CellAttachedDataSerializer<dim, spacedim>::cell_relation_t>
    &                cell_relations,
  const unsigned int handle,
  const std::function<
    void(const cell_iterator &,
         const CellStatus &,
         const boost::iterator_range<std::vector<char>::const_iterator> &)>
    &unpack_callback) const
{
  // We decode the handle returned by register_data_attach() back into
  // a format we can use. All even handles belong to those callback
  // functions which write/read variable size data, all odd handles
  // interact with fixed size buffers.
  const bool         callback_variable_transfer = (handle % 2 == 0);
  const unsigned int callback_index             = handle / 2;

  // Cells will always receive fixed size data (i.e., CellStatus
  // information), but not necessarily variable size data (e.g., with a
  // ParticleHandler a cell might not contain any particle at all).
  // Thus it is sufficient to check if fixed size data has been received.
  Assert(sizes_fixed_cumulative.size() > 0,
         ExcMessage("No data has been packed!"));
  if (cell_relations.size() > 0)
    {
      Assert(dest_data_fixed.size() > 0,
             ExcMessage("No data has been received!"));
    }

  std::vector<char>::const_iterator dest_data_it;
  std::vector<char>::const_iterator dest_sizes_cell_it;

  // Depending on whether our callback function unpacks fixed or
  // variable size data, we have to pursue different approaches
  // to localize the correct fraction of the buffer from which
  // we are allowed to read.
  unsigned int offset         = numbers::invalid_unsigned_int;
  unsigned int size           = numbers::invalid_unsigned_int;
  unsigned int data_increment = numbers::invalid_unsigned_int;

  if (callback_variable_transfer)
    {
      // For the variable size data, we need to extract the
      // data size from the fixed size buffer on each cell.
      //
      // We packed this information last, so the last packed
      // object in the fixed size buffer corresponds to the
      // variable data sizes.
      //
      // The last entry of sizes_fixed_cumulative corresponds
      // to the size of all fixed size data packed on the cell.
      // To get the offset for the last packed object, we need
      // to get the next-to-last entry.
      const unsigned int offset_variable_data_sizes =
        sizes_fixed_cumulative[sizes_fixed_cumulative.size() - 2];

      // This iterator points to the data size that the
      // callback_function packed for each specific cell.
      // Adjust buffer iterator to the offset of the callback
      // function so that we only have to advance its position
      // to the next cell after each iteration.
      dest_sizes_cell_it = dest_data_fixed.cbegin() +
                           offset_variable_data_sizes +
                           callback_index * sizeof(unsigned int);

      // Let the data iterator point to the correct buffer.
      dest_data_it = dest_data_variable.cbegin();
    }
  else
    {
      // For the fixed size data, we can get the information about
      // the buffer location on each cell directly from the
      // sizes_fixed_cumulative vector.
      offset         = sizes_fixed_cumulative[callback_index];
      size           = sizes_fixed_cumulative[callback_index + 1] - offset;
      data_increment = sizes_fixed_cumulative.back();

      // Let the data iterator point to the correct buffer.
      // Adjust buffer iterator to the offset of the callback
      // function so that we only have to advance its position
      // to the next cell after each iteration.
      if (cell_relations.begin() != cell_relations.end())
        dest_data_it = dest_data_fixed.cbegin() + offset;
    }

  // Iterate over all cells and unpack the transferred data.
  auto cell_rel_it   = cell_relations.begin();
  auto dest_sizes_it = dest_sizes_variable.cbegin();
  for (; cell_rel_it != cell_relations.end(); ++cell_rel_it)
    {
      const auto &dealii_cell = cell_rel_it->first;
      const auto &cell_status = cell_rel_it->second;

      if (callback_variable_transfer)
        {
          // Update the increment according to the whole data size
          // of the current cell.
          data_increment = *dest_sizes_it;

          if (cell_status != CellStatus::cell_invalid)
            {
              // Extract the corresponding values for offset and size from
              // the cumulative sizes array stored in the fixed size
              // buffer.
              if (callback_index == 0)
                offset = 0;
              else
                std::memcpy(&offset,
                            &(*(dest_sizes_cell_it - sizeof(unsigned int))),
                            sizeof(unsigned int));

              std::memcpy(&size, &(*dest_sizes_cell_it), sizeof(unsigned int));

              size -= offset;

              // Move the data iterator to the corresponding position
              // of the callback function and adjust the increment
              // accordingly.
              dest_data_it += offset;
              data_increment -= offset;
            }

          // Advance data size iterators to the next cell, avoid iterating
          // past the end of dest_sizes_cell_it
          if (cell_rel_it != cell_relations.end() - 1)
            dest_sizes_cell_it += sizes_fixed_cumulative.back();
          ++dest_sizes_it;
        }

      switch (cell_status)
        {
          case CellStatus::cell_will_persist:
          case CellStatus::children_will_be_coarsened:
            unpack_callback(dealii_cell,
                            cell_status,
                            boost::make_iterator_range(dest_data_it,
                                                       dest_data_it + size));
            break;

          case CellStatus::cell_will_be_refined:
            unpack_callback(dealii_cell->parent(),
                            cell_status,
                            boost::make_iterator_range(dest_data_it,
                                                       dest_data_it + size));
            break;

          case CellStatus::cell_invalid:
            // Skip this cell.
            break;

          default:
            Assert(false, ExcInternalError());
            break;
        }

      if (cell_rel_it != cell_relations.end() - 1)
        dest_data_it += data_increment;
    }
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
void CellAttachedDataSerializer<dim, spacedim>::save(
  const unsigned int global_first_cell,
  const unsigned int global_num_cells,
  const std::string &filename,
  const MPI_Comm &   mpi_communicator) const
{
  Assert(sizes_fixed_cumulative.size() > 0,
         ExcMessage("No data has been packed!"));

#ifdef DEAL_II_WITH_MPI
  // Large fractions of this function have been copied from
  // DataOutInterface::write_vtu_in_parallel.
  // TODO: Write general MPIIO interface.

  const int myrank  = Utilities::MPI::this_mpi_process(mpi_communicator);
  const int mpisize = Utilities::MPI::n_mpi_processes(mpi_communicator);

  if (mpisize > 1)
    {
      const unsigned int bytes_per_cell = sizes_fixed_cumulative.back();

      //
      // ---------- Fixed size data ----------
      //
      {
        const std::string fname_fixed = std::string(filename) + "_fixed.data";

        MPI_Info info;
        int      ierr = MPI_Info_create(&info);
        AssertThrowMPI(ierr);

        MPI_File fh;
        ierr = MPI_File_open(mpi_communicator,
                             fname_fixed.c_str(),
                             MPI_MODE_CREATE | MPI_MODE_WRONLY,
                             info,
                             &fh);
        AssertThrowMPI(ierr);

        ierr = MPI_File_set_size(fh, 0); // delete the file contents
        AssertThrowMPI(ierr);
        // this barrier is necessary, because otherwise others might already
        // write while one core is still setting the size to zero.
        ierr = MPI_Barrier(mpi_communicator);
        AssertThrowMPI(ierr);
        ierr = MPI_Info_free(&info);
        AssertThrowMPI(ierr);
        // ------------------

        // Write cumulative sizes to file.
        // Since each processor owns the same information about the data
        // sizes, it is sufficient to let only the first processor perform
        // this task.
        if (myrank == 0)
          {
            ierr = Utilities::MPI::LargeCount::File_write_at_c(
              fh,
              0,
              sizes_fixed_cumulative.data(),
              sizes_fixed_cumulative.size(),
              MPI_UNSIGNED,
              MPI_STATUS_IGNORE);
            AssertThrowMPI(ierr);
          }

        // Write packed data to file simultaneously.
        const MPI_Offset size_header =
          sizes_fixed_cumulative.size() * sizeof(unsigned int);

        // Make sure we do the following computation in 64bit integers to be
        // able to handle 4GB+ files:
        const MPI_Offset my_global_file_position =
          size_header +
          static_cast<MPI_Offset>(global_first_cell) * bytes_per_cell;

        ierr =
          Utilities::MPI::LargeCount::File_write_at_c(fh,
                                                      my_global_file_position,
                                                      src_data_fixed.data(),
                                                      src_data_fixed.size(),
                                                      MPI_BYTE,
                                                      MPI_STATUS_IGNORE);
        AssertThrowMPI(ierr);

        ierr = MPI_File_close(&fh);
        AssertThrowMPI(ierr);
      }



      //
      // ---------- Variable size data ----------
      //
      if (variable_size_data_stored)
        {
          const std::string fname_variable =
            std::string(filename) + "_variable.data";

          MPI_Info info;
          int      ierr = MPI_Info_create(&info);
          AssertThrowMPI(ierr);

          MPI_File fh;
          ierr = MPI_File_open(mpi_communicator,
                               fname_variable.c_str(),
                               MPI_MODE_CREATE | MPI_MODE_WRONLY,
                               info,
                               &fh);
          AssertThrowMPI(ierr);

          ierr = MPI_File_set_size(fh, 0); // delete the file contents
          AssertThrowMPI(ierr);
          // this barrier is necessary, because otherwise others might already
          // write while one core is still setting the size to zero.
          ierr = MPI_Barrier(mpi_communicator);
          AssertThrowMPI(ierr);
          ierr = MPI_Info_free(&info);
          AssertThrowMPI(ierr);

          // Write sizes of each cell into file simultaneously.
          {
            const MPI_Offset my_global_file_position =
              static_cast<MPI_Offset>(global_first_cell) * sizeof(unsigned int);

            // It is very unlikely that a single process has more than
            // 2 billion cells, but we might as well check.
            AssertThrow(src_sizes_variable.size() <
                          static_cast<std::size_t>(
                            std::numeric_limits<int>::max()),
                        ExcNotImplemented());

            ierr = Utilities::MPI::LargeCount::File_write_at_c(
              fh,
              my_global_file_position,
              src_sizes_variable.data(),
              src_sizes_variable.size(),
              MPI_INT,
              MPI_STATUS_IGNORE);
            AssertThrowMPI(ierr);
          }

          // Gather size of data in bytes we want to store from this
          // processor and compute the prefix sum. We do this in 64 bit
          // to avoid overflow for files larger than 4GB:
          const std::uint64_t size_on_proc = src_data_variable.size();
          std::uint64_t       prefix_sum   = 0;
          ierr                             = MPI_Exscan(&size_on_proc,
                            &prefix_sum,
                            1,
                            MPI_UINT64_T,
                            MPI_SUM,
                            mpi_communicator);
          AssertThrowMPI(ierr);

          const MPI_Offset my_global_file_position =
            static_cast<MPI_Offset>(global_num_cells) * sizeof(unsigned int) +
            prefix_sum;

          // Write data consecutively into file.
          ierr = Utilities::MPI::LargeCount::File_write_at_c(
            fh,
            my_global_file_position,
            src_data_variable.data(),
            src_data_variable.size(),
            MPI_BYTE,
            MPI_STATUS_IGNORE);
          AssertThrowMPI(ierr);


          ierr = MPI_File_close(&fh);
          AssertThrowMPI(ierr);
        }
    } // if (mpisize > 1)
  else
#endif
    {
      (void)global_first_cell;
      (void)global_num_cells;

      //
      // ---------- Fixed size data ----------
      //
      {
        const std::string fname_fixed = std::string(filename) + "_fixed.data";

        std::ofstream file(fname_fixed, std::ios::binary | std::ios::out);

        // Write header data.
        file.write(reinterpret_cast<const char *>(
                     sizes_fixed_cumulative.data()),
                   sizes_fixed_cumulative.size() * sizeof(unsigned int));

        // Write packed data.
        file.write(reinterpret_cast<const char *>(src_data_fixed.data()),
                   src_data_fixed.size() * sizeof(char));

        file.close();
      }

      //
      // ---------- Variable size data ----------
      //
      if (variable_size_data_stored)
        {
          const std::string fname_variable =
            std::string(filename) + "_variable.data";

          std::ofstream file(fname_variable, std::ios::binary | std::ios::out);

          // Write header data.
          file.write(reinterpret_cast<const char *>(src_sizes_variable.data()),
                     src_sizes_variable.size() * sizeof(int));

          // Write packed data.
          file.write(reinterpret_cast<const char *>(src_data_variable.data()),
                     src_data_variable.size() * sizeof(char));

          file.close();
        }
    }
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
void CellAttachedDataSerializer<dim, spacedim>::load(
  const unsigned int global_first_cell,
  const unsigned int global_num_cells,
  const unsigned int local_num_cells,
  const std::string &filename,
  const unsigned int n_attached_deserialize_fixed,
  const unsigned int n_attached_deserialize_variable,
  const MPI_Comm &   mpi_communicator)
{
  Assert(dest_data_fixed.size() == 0,
         ExcMessage("Previously loaded data has not been released yet!"));

  variable_size_data_stored = (n_attached_deserialize_variable > 0);

#ifdef DEAL_II_WITH_MPI
  // Large fractions of this function have been copied from
  // DataOutInterface::write_vtu_in_parallel.
  // TODO: Write general MPIIO interface.

  const int mpisize = Utilities::MPI::n_mpi_processes(mpi_communicator);

  if (mpisize > 1)
    {
      //
      // ---------- Fixed size data ----------
      //
      {
        const std::string fname_fixed = std::string(filename) + "_fixed.data";

        MPI_Info info;
        int      ierr = MPI_Info_create(&info);
        AssertThrowMPI(ierr);

        MPI_File fh;
        ierr = MPI_File_open(
          mpi_communicator, fname_fixed.c_str(), MPI_MODE_RDONLY, info, &fh);
        AssertThrowMPI(ierr);

        ierr = MPI_Info_free(&info);
        AssertThrowMPI(ierr);

        // Read cumulative sizes from file.
        // Since all processors need the same information about the data
        // sizes, let each of them retrieve it by reading from the same
        // location in the file.
        sizes_fixed_cumulative.resize(1 + n_attached_deserialize_fixed +
                                      (variable_size_data_stored ? 1 : 0));
        ierr = Utilities::MPI::LargeCount::File_read_at_c(
          fh,
          0,
          sizes_fixed_cumulative.data(),
          sizes_fixed_cumulative.size(),
          MPI_UNSIGNED,
          MPI_STATUS_IGNORE);
        AssertThrowMPI(ierr);

        // Allocate sufficient memory.
        const unsigned int bytes_per_cell = sizes_fixed_cumulative.back();
        dest_data_fixed.resize(static_cast<size_t>(local_num_cells) *
                               bytes_per_cell);

        // Read packed data from file simultaneously.
        const MPI_Offset size_header =
          sizes_fixed_cumulative.size() * sizeof(unsigned int);

        // Make sure we do the following computation in 64bit integers to be
        // able to handle 4GB+ files:
        const MPI_Offset my_global_file_position =
          size_header +
          static_cast<MPI_Offset>(global_first_cell) * bytes_per_cell;

        ierr =
          Utilities::MPI::LargeCount::File_read_at_c(fh,
                                                     my_global_file_position,
                                                     dest_data_fixed.data(),
                                                     dest_data_fixed.size(),
                                                     MPI_BYTE,
                                                     MPI_STATUS_IGNORE);
        AssertThrowMPI(ierr);


        ierr = MPI_File_close(&fh);
        AssertThrowMPI(ierr);
      }

      //
      // ---------- Variable size data ----------
      //
      if (variable_size_data_stored)
        {
          const std::string fname_variable =
            std::string(filename) + "_variable.data";

          MPI_Info info;
          int      ierr = MPI_Info_create(&info);
          AssertThrowMPI(ierr);

          MPI_File fh;
          ierr = MPI_File_open(mpi_communicator,
                               fname_variable.c_str(),
                               MPI_MODE_RDONLY,
                               info,
                               &fh);
          AssertThrowMPI(ierr);

          ierr = MPI_Info_free(&info);
          AssertThrowMPI(ierr);

          // Read sizes of all locally owned cells.
          dest_sizes_variable.resize(local_num_cells);

          const MPI_Offset my_global_file_position_sizes =
            static_cast<MPI_Offset>(global_first_cell) * sizeof(unsigned int);

          ierr = Utilities::MPI::LargeCount::File_read_at_c(
            fh,
            my_global_file_position_sizes,
            dest_sizes_variable.data(),
            dest_sizes_variable.size(),
            MPI_INT,
            MPI_STATUS_IGNORE);
          AssertThrowMPI(ierr);


          // Compute my data size in bytes and compute prefix sum. We do this
          // in 64 bit to avoid overflow for files larger than 4 GB:
          const std::uint64_t size_on_proc =
            std::accumulate(dest_sizes_variable.begin(),
                            dest_sizes_variable.end(),
                            0ULL);

          std::uint64_t prefix_sum = 0;
          ierr                     = MPI_Exscan(&size_on_proc,
                            &prefix_sum,
                            1,
                            MPI_UINT64_T,
                            MPI_SUM,
                            mpi_communicator);
          AssertThrowMPI(ierr);

          const MPI_Offset my_global_file_position =
            static_cast<MPI_Offset>(global_num_cells) * sizeof(unsigned int) +
            prefix_sum;

          dest_data_variable.resize(size_on_proc);

          ierr = Utilities::MPI::LargeCount::File_read_at_c(
            fh,
            my_global_file_position,
            dest_data_variable.data(),
            dest_data_variable.size(),
            MPI_BYTE,
            MPI_STATUS_IGNORE);
          AssertThrowMPI(ierr);

          ierr = MPI_File_close(&fh);
          AssertThrowMPI(ierr);
        }
    }
  else // if (mpisize > 1)
#endif
    {
      (void)global_first_cell;
      (void)global_num_cells;
      (void)local_num_cells;

      //
      // ---------- Fixed size data ----------
      //
      {
        const std::string fname_fixed = std::string(filename) + "_fixed.data";

        std::ifstream file(fname_fixed, std::ios::binary | std::ios::in);
        sizes_fixed_cumulative.resize(1 + n_attached_deserialize_fixed +
                                      (variable_size_data_stored ? 1 : 0));

        // Read header data.
        file.read(reinterpret_cast<char *>(sizes_fixed_cumulative.data()),
                  sizes_fixed_cumulative.size() * sizeof(unsigned int));

        const unsigned int bytes_per_cell = sizes_fixed_cumulative.back();
        dest_data_fixed.resize(static_cast<size_t>(local_num_cells) *
                               bytes_per_cell);

        // Read packed data.
        file.read(reinterpret_cast<char *>(dest_data_fixed.data()),
                  dest_data_fixed.size() * sizeof(char));

        file.close();
      }

      //
      // ---------- Variable size data ----------
      //
      if (variable_size_data_stored)
        {
          const std::string fname_variable =
            std::string(filename) + "_variable.data";

          std::ifstream file(fname_variable, std::ios::binary | std::ios::in);

          // Read header data.
          file.read(reinterpret_cast<char *>(dest_sizes_variable.data()),
                    dest_sizes_variable.size() * sizeof(int));

          // Read packed data.
          file.read(reinterpret_cast<char *>(dest_data_variable.data()),
                    dest_data_variable.size() * sizeof(char));

          file.close();
        }
    }
}


template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
void CellAttachedDataSerializer<dim, spacedim>::clear()
{
  variable_size_data_stored = false;

  // free information about data sizes
  sizes_fixed_cumulative.clear();
  sizes_fixed_cumulative.shrink_to_fit();

  // free fixed size transfer data
  src_data_fixed.clear();
  src_data_fixed.shrink_to_fit();

  dest_data_fixed.clear();
  dest_data_fixed.shrink_to_fit();

  // free variable size transfer data
  src_sizes_variable.clear();
  src_sizes_variable.shrink_to_fit();

  src_data_variable.clear();
  src_data_variable.shrink_to_fit();

  dest_sizes_variable.clear();
  dest_sizes_variable.shrink_to_fit();

  dest_data_variable.clear();
  dest_data_variable.shrink_to_fit();
}

// explicit instantiations
#include "data_serializer.inst"

DEAL_II_NAMESPACE_CLOSE
