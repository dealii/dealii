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

#ifndef dealii_tria_data_serializer_h
#define dealii_tria_data_serializer_h

#include <deal.II/base/config.h>

#include <deal.II/base/mpi.h>
#include <deal.II/base/mpi_large_count.h>
#include <deal.II/base/utilities.h>

#include <deal.II/grid/cell_status.h>
#include <deal.II/grid/tria_iterator_selector.h>

#include <boost/range/iterator_range.hpp>

#include <algorithm>
#include <fstream>
#include <functional>
#include <utility>
#include <vector>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  /**
   * A structure that binds information about data attached to cells.
   */
  template <int dim, int spacedim = dim>
  DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
  struct CellAttachedData
  {
    using cell_iterator = TriaIterator<CellAccessor<dim, spacedim>>;

    /**
     * Number of functions that get attached to the Triangulation through
     * register_data_attach() for example SolutionTransfer.
     */
    unsigned int n_attached_data_sets;

    /**
     * Number of functions that need to unpack their data after a call from
     * load().
     */
    unsigned int n_attached_deserialize;

    using pack_callback_t =
      std::function<std::vector<char>(cell_iterator, CellStatus)>;

    /**
     * These callback functions will be stored in the order in which they
     * have been registered with the register_data_attach() function.
     */
    std::vector<pack_callback_t> pack_callbacks_fixed;
    std::vector<pack_callback_t> pack_callbacks_variable;
  };
} // namespace internal

/**
 * A structure that stores information about the data that has been, or
 * will be, attached to cells via the register_data_attach() function
 * and later retrieved via notify_ready_to_unpack().
 *
 * This internalclass is dedicated to the data serialization and transfer across
 * repartitioned meshes and to/from the file system.
 *
 * It is designed to store all data buffers intended for serialization.
 */
template <int dim, int spacedim = dim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
class CellAttachedDataSerializer
{
public:
  using cell_iterator = TriaIterator<CellAccessor<dim, spacedim>>;

  /**
   * Auxiliary data structure for assigning a CellStatus to a deal.II cell
   * iterator. For an extensive description of the former, see the
   * documentation for the member function register_data_attach().
   */
  using cell_relation_t = typename std::pair<cell_iterator, CellStatus>;

  CellAttachedDataSerializer();

  /**
   * Prepare data serialization by calling the pack callback functions on each
   * cell in @p cell_relations.
   *
   * All registered callback functions in @p pack_callbacks_fixed will write
   * into the fixed size buffer, whereas each entry of @p pack_callbacks_variable
   * will write its data into the variable size buffer.
   */
  void
  pack_data(
    const std::vector<cell_relation_t> &cell_relations,
    const std::vector<
      typename internal::CellAttachedData<dim, spacedim>::pack_callback_t>
      &pack_callbacks_fixed,
    const std::vector<
      typename internal::CellAttachedData<dim, spacedim>::pack_callback_t>
      &             pack_callbacks_variable,
    const MPI_Comm &mpi_communicator);

  /**
   * Unpack the CellStatus information on each entry of
   * @p cell_relations.
   *
   * Data has to be previously transferred with execute_transfer()
   * or deserialized from the file system via load().
   */
  void
  unpack_cell_status(std::vector<cell_relation_t> &cell_relations) const;

  /**
   * Unpack previously serialized data on each cell registered in
   * @p cell_relations with the provided @p unpack_callback function.
   *
   * The parameter @p handle corresponds to the position where the
   * @p unpack_callback function is allowed to read from the memory. Its
   * value needs to be in accordance with the corresponding pack_callback
   * function that has been registered previously.
   *
   * Data has to be previously transferred with execute_transfer()
   * or deserialized from the file system via load().
   */
  void
  unpack_data(
    const std::vector<cell_relation_t> &cell_relations,
    const unsigned int                  handle,
    const std::function<
      void(const cell_iterator &,
           const CellStatus &,
           const boost::iterator_range<std::vector<char>::const_iterator> &)>
      &unpack_callback) const;

  /**
   * Serialize data to file system.
   *
   * The data will be written in a separate file, whose name
   * consists of the stem @p filename and an attached identifier
   * <tt>_fixed.data</tt> for fixed size data and <tt>_variable.data</tt>
   * for variable size data.
   *
   * If MPI support is enabled, all processors write into these files
   * simultaneously via MPIIO. Each processor's position to write to will be
   * determined from the provided input parameters.
   *
   * Data has to be previously packed with pack_data().
   */
  void
  save(const unsigned int global_first_cell,
       const unsigned int global_num_cells,
       const std::string &filename,
       const MPI_Comm &   mpi_communicator) const;

  /**
   * Deserialize data from file system.
   *
   * The data will be read from separate file, whose name
   * consists of the stem @p filename and an attached identifier
   * <tt>_fixed.data</tt> for fixed size data and <tt>_variable.data</tt>
   * for variable size data.
   * The @p n_attached_deserialize_fixed and @p n_attached_deserialize_variable
   * parameters are required to gather the memory offsets for each
   * callback.
   *
   * If MPI support is enabled, all processors read from these files
   * simultaneously via MPIIO. Each processor's position to read from will be
   * determined from the provided input arguments.
   *
   * After loading, unpack_data() needs to be called to finally
   * distribute data across the associated triangulation.
   */
  void
  load(const unsigned int global_first_cell,
       const unsigned int global_num_cells,
       const unsigned int local_num_cells,
       const std::string &filename,
       const unsigned int n_attached_deserialize_fixed,
       const unsigned int n_attached_deserialize_variable,
       const MPI_Comm &   mpi_communicator);

  /**
   * Clears all containers and associated data, and resets member
   * values to their default state.
   *
   * Frees memory completely.
   */
  void
  clear();

  /**
   * Flag that denotes if variable size data has been packed.
   */
  bool variable_size_data_stored;

  /**
   * Cumulative size in bytes that those functions that have called
   * register_data_attach() want to attach to each cell. This number
   * only pertains to fixed-sized buffers where the data attached to
   * each cell has exactly the same size.
   *
   * The last entry of this container corresponds to the data size
   * packed per cell in the fixed size buffer (which can be accessed
   * calling <tt>sizes_fixed_cumulative.back()</tt>).
   */
  std::vector<unsigned int> sizes_fixed_cumulative;

  /**
   * Consecutive buffers designed for the fixed size serialization
   * functions.
   */
  std::vector<char> src_data_fixed;
  std::vector<char> dest_data_fixed;

  /**
   * Consecutive buffers designed for the variable size serialization
   * functions.
   */
  std::vector<int>  src_sizes_variable;
  std::vector<int>  dest_sizes_variable;
  std::vector<char> src_data_variable;
  std::vector<char> dest_data_variable;
};

DEAL_II_NAMESPACE_CLOSE

#endif
