// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#ifndef dealii__patch_storage_h
#define dealii__patch_storage_h



#include <deal.II/base/data_out_base.h>
#include <deal.II/base/enable_observer_pointer.h>
#include <deal.II/base/iterator_range.h>
#include <deal.II/base/partitioner.h>
#include <deal.II/base/types.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/vector.h>

#include <deal.II/matrix_free/matrix_free.h>



DEAL_II_NAMESPACE_OPEN

#ifndef DOXYGEN

namespace internal
{
  /*
   * apply permutation to the given vector.
   */
  template <class T>
  void
  reorder(std::vector<T> &vA, std::vector<size_t> vOrder)
  {
    AssertDimension(vA.size(), vOrder.size());

    // for all elements to put in place
    for (size_t i = 0; i < vA.size(); ++i)
      {
        // while vOrder[i] is not yet in place
        // every swap places at least one element in it's proper place
        while (vOrder[i] != vOrder[vOrder[i]])
          {
            std::swap(vA[vOrder[i]], vA[vOrder[vOrder[i]]]);
            std::swap(vOrder[i], vOrder[vOrder[i]]);
          }
      }
  }



  // Computes index of given vertex within a cell
  template <int dim>
  unsigned int
  compute_vertex_index(const typename Triangulation<dim>::cell_iterator &cell,
                       const types::global_vertex_index                 &vertex)
  {
    for (const auto &i : GeometryInfo<dim>::vertex_indices())
      if (cell->vertex_index(i) == vertex)
        return i;
    Assert(false, ExcMessage("Vertex not found on the given cell!"));
    return numbers::invalid_unsigned_int;
  }

  // helper functions for RegularVertexPatch constructor
  template <int dim>
  std::vector<typename Triangulation<dim>::cell_iterator>
  order_patch(
    const std::set<typename Triangulation<dim>::cell_iterator> &patch_cells,
    const types::global_vertex_index                           &vertex_index);

  // rotate patch to minimize cell rotations
  void
  orient_patch2D(
    std::vector<typename Triangulation<2>::cell_iterator> &patch_cells,
    const types::global_vertex_index                      &vertex);

  template <int dim>
  void
  rotate_patch(
    std::vector<typename Triangulation<dim>::cell_iterator> &patch_cells);


  // Explicit specialization declaration BEFORE any use
  template <>
  void
  rotate_patch<2>(
    std::vector<typename Triangulation<2>::cell_iterator> &patch_cells);


  namespace GaussSeidel
  {
    class TaskInfoDummy
    {
    public:
      using cell_index_type_scalar = unsigned int;
      const static constexpr unsigned int invalid_category =
        std::numeric_limits<unsigned int>::max();
      using cell_index_type_vector = std::array<unsigned int, 1>;
      const static constexpr unsigned int n_categories = 1;

      TaskInfoDummy()
        : current_communicate_index(0)
      {}
      unsigned int current_communicate_index;

      template <typename number>
      void
      communicate(const unsigned int,
                  const int,
                  ::dealii::LinearAlgebra::distributed::Vector<number> &) const
      {}

      template <typename... Args>
      void
      initialize(const Args &...args)
      {
        (void)sizeof...(args);
      }


      unsigned int
      determine_minimum_category_of_patch(
        const std::set<types::global_dof_index> &,
        const std::vector<unsigned int> &) const
      {
        return 0;
      }

      static unsigned int
      determine_cell_partition(const unsigned int &)
      {
        return 0;
      }
    };

    unsigned int
    min_category(const unsigned int a, const unsigned int b);

  } // namespace GaussSeidel
} // namespace internal


#endif // DOXYGEN

/**
 * Represents a regular patch, i.e., a patch centered at a vertex
 * with exactly 2^dim cells.
 *
 * Stores the cell indices and potentially orientation information. Provides
 * methods to check constructibility and access cell data.
 *
 * @tparam dim The spatial dimension.
 */
template <int dim>
struct RegularVertexPatch
{
  /**
   * Type used to represent a unique index for each cell.
   */
  using CellIndex = unsigned int;

  /**
   * Alias for the cell iterator type of the underlying Triangulation.
   */
  using CellIterator = typename Triangulation<dim>::cell_iterator;

  /**
   * Type used to represent the orientation of a cell within a regular patch.
   */
  using CellOrientation = unsigned int;

  const static constexpr int dimension        = dim;
  const static bool          is_constant_size = true;
  const static unsigned int  n_cells          = 1 << dim;


  /**
   *  Constructor for a RegularVertexPatch.
   *
   * Orders the cells within the patch based on the central vertex index and
   * potentially generates orientation information. Asserts that the input
   * `patch` set contains exactly `n_cells`.
   *
   * @param patch A set of `CellIndex` objects representing the cells in the
   * patch.
   * @param vertex_index The global index of the central vertex.
   * @param index2cell A function object that converts a `CellIndex` to a
   * `CellIterator`.
   */
  RegularVertexPatch(
    const std::set<CellIndex>                            &patch,
    const types::global_vertex_index                     &vertex_index,
    const std::function<CellIterator(const CellIndex &)> &index2cell);



  /**
   *  Returns the number of cells in the patch (always `n_cells`).
   */
  constexpr unsigned int
  size() const
  {
    return n_cells;
  }
  /**
   *  Checks if this patch conflicts with another patch for parallel
   * processing (e.g., based on shared DoFs).
   * @param other The other RegularVertexPatch to check against.
   * @return `true` if there is no conflict, `false` otherwise.
   */
  bool
  has_conflict_with(const RegularVertexPatch &other) const;


  /**
   *  Provides read-only access to the array of cell indices forming
   * the patch.
   * @return A constant reference to the array of `CellIndex`.
   */
  const auto &
  get_cells() const
  {
    return cells;
  }


  /**
   *  Gets the orientation information for a specific cell within the
   * patch.
   * @param cell_index The local index of the cell within the patch (0 to
   * n_cells-1).
   * @return A constant reference to the `CellOrientation`.
   */
  const auto &
  get_orientation(const unsigned int &cell_index) const
  {
    AssertIndexRange(cell_index, size());
    return orientations[cell_index];
  }


  /**
   *  Static method to check if a given set of cells can form a
   * regular patch.
   *
   * Currently, this simply checks if the number of cells in the input `patch`
   * is equal to `n_cells`.
   *
   * @param patch A set of `CellIndex` objects.
   * @param vertex_index The global index of the potential central vertex
   * (unused).
   * @param index2cell A function object to convert `CellIndex` to
   * `CellIterator` (unused).
   * @return `true` if the patch can be constructed as a RegularVertexPatch,
   * `false` otherwise.
   */
  static bool
  is_constructible(
    const std::set<CellIndex>                            &patch,
    const types::global_vertex_index                     &vertex_index,
    const std::function<CellIterator(const CellIndex &)> &index2cell)
  {
    (void)vertex_index;
    (void)index2cell;
    if (patch.size() == n_cells)
      return true;
    return false;
  }


  /**
   *  Checks if any cell within this patch is a ghost cell on the
   * current MPI process.
   * @return `true` if the patch contains at least one ghost cell, `false`
   * otherwise.
   */
  bool
  is_partially_ghosted() const
  {
    return partially_ghosted;
  }

private:
  std::array<CellIndex, n_cells>       cells;
  std::array<CellOrientation, n_cells> orientations;


  bool partially_ghosted;
};

/**
 *  Represents a general patch, i.e., a patch that is not regular
 * (typically near boundaries or excluded vertices).
 *
 * Stores a variable number of cell indices.
 *
 * @tparam dim The spatial dimension.
 */
template <int dim>
struct GeneralVertexPatch
{
  /**
   * Type used to represent a unique index for each cell.
   */
  using CellIndex = unsigned int;

  /**
   * Alias for the cell iterator type of the underlying Triangulation.
   */
  using CellIterator = typename Triangulation<dim>::cell_iterator;

  const static constexpr int dimension        = dim;
  const static bool          is_constant_size = false;

  /**
   *  Constructor for a GeneralVertexPatch.
   *
   * @param patch A set of `CellIndex` objects representing the cells in the
   * patch.
   * @param vertex_index The global index of the central vertex (unused in
   * current implementation).
   * @param index2cell A function object that converts a `CellIndex` to a
   * `CellIterator` (unused in current implementation).
   */
  GeneralVertexPatch(
    const std::set<CellIndex>                            &patch,
    const types::global_vertex_index                     &vertex_index,
    const std::function<CellIterator(const CellIndex &)> &index2cell);

  /**
   *  Returns the number of cells in the patch.
   */
  unsigned int
  size() const
  {
    return cells.size();
  }

  /**
   *  Provides read-only access to the vector of cell indices forming
   * the patch.
   * @return A constant reference to the vector of `CellIndex`.
   */
  const auto &
  get_cells() const
  {
    return cells;
  }

  // Add other necessary members and methods for GeneralVertexPatch
  // For example:
  // bool has_conflict_with(const GeneralVertexPatch &other) const;
  // bool is_partially_ghosted() const;

private:
  std::vector<CellIndex> cells;
  // Add other necessary private members
};

/**
 * Manages the storage and categorization of cell patches centered around
 * vertices.
 *
 * This class identifies patches of cells surrounding each vertex in a mesh
 * level managed by a MatrixFree object. It distinguishes between "regular"
 * patches (those with the expected number of cells, 2^dim) and "general"
 * patches (e.g., near boundaries). Patches can be categorized for parallel
 * processing, typically for patch-based smoothers like Gauss-Seidel.
 *
 * The class provides mechanisms to iterate over patches, access patch data,
 * and handle communication for parallel operations.
 *
 * @tparam MFType The type of the MatrixFree object used to access cell and DoF
 * information. Expected to provide types like `dimension`, `value_type`,
 * `vectorized_value_type`, and methods to access DoF handlers and cell
 * iterators.
 */
template <class MFType>
class PatchStorage : EnableObserverPointer
{
public:
  /**
   * Alias for the template parameter `MFType`.
   */
  using MatrixFreeType = MFType;


  /**
   * The spatial dimension.
   */
  const static constexpr unsigned int dim = MatrixFreeType::dimension;

  /**
   * The scalar value type used by the MatrixFree object (e.g., float,
   * double).
   */
  using value_type = typename MatrixFreeType::value_type;


  /**
   * The vectorized value type used by the MatrixFree object.
   */
  using vectorized_value_type = typename MatrixFreeType::vectorized_value_type;


  /**
   * The number of lanes in the vectorized type.
   */
  const static constexpr unsigned int n_lanes = vectorized_value_type::size();


  /**
   * The expected number of cells in a regular patch (2^dim).
   */
  const static constexpr unsigned int n_patch_cells = 1 << dim;


  /**
   * Type used to represent a unique index for each cell within the
   * context of the associated MatrixFree object. This is typically calculated
   * as `batch_index * n_lanes + lane_index`.
   */
  using CellIndex = unsigned int;


  /**
   * Alias for the cell iterator type of the underlying Triangulation.
   */
  using CellIterator = typename Triangulation<dim>::cell_iterator;

  /**
   * Type used to represent the orientation of a cell within a regular
   * patch (currently unused).
   */
  using CellOrientation = unsigned int;

  /**
   * Structure to hold additional data for patch generation, such as
   * vertices to exclude.
   */
  struct AdditionalData
  {
    /**
     * Default constructor. Initializes the excluded vertices set with
     * an invalid index.
     */
    AdditionalData()
      : excluded_vertices{numbers::invalid_unsigned_int}
    {}

    /**
     * A set of global vertex indices that should be excluded during
     * patch generation.
     */
    std::set<types::global_vertex_index> excluded_vertices;
  };


  /**
   * Alias for RegularVertexPatch with the current dimension.
   */
  using RegularPatch = RegularVertexPatch<dim>;

  /**
   * Alias for GeneralVertexPatch with the current dimension.
   */
  using GeneralPatch = GeneralVertexPatch<dim>;


  /**
   *  Alias for the task information type used for parallel scheduling
   * (e.g., Gauss-Seidel coloring). At the moment MPI is disabled.
   */
  using TaskInfoType = internal::GaussSeidel::TaskInfoDummy;
  // LinearAlgebra::GaussSeidel::TaskInfo<value_type, vectorized_value_type>;

  /**
   *  Type representing a range of patch indices [begin, end).
   */
  using PatchRange = const std::pair<std::size_t, std::size_t>;
  /**
   *  Type representing a category assigned to a patch, often used for
   * grouping or identification.
   */
  using PatchCategory = unsigned int;



  /**
   *  Constructor.
   * @param mf A shared pointer to the constant MatrixFree object.
   */
  PatchStorage(const std::shared_ptr<const MatrixFreeType> &mf);


  /**
   *  Initializes the PatchStorage by generating and storing patches.
   *
   * Identifies all vertex-centered patches on the specified level of the
   * MatrixFree object. Filters patches based on ownership and MPI rank.
   * Initializes internal data structures like partitioners and task info.
   *
   * @param data Additional data, e.g., vertices to exclude.
   */
  void
  initialize(const AdditionalData &data = AdditionalData());


  /**
   * Assigns a category to each regular patch based on a user-provided
   * function. @warning This function may change indices of patches.
   *
   * @param category_function A function that takes a `const RegularPatch &`
   * and returns a `PatchCategory`.
   */
  void
  categorize_patches(
    std::function<PatchCategory(const RegularPatch &)> category_function);


  /**
   *  Clears all stored patch data and resets the state to
   * uninitialized.
   */
  void
  clear();


  /**
   *  Executes a given function (`patch_worker`) over ranges of patches,
   * handling parallel communication.
   *
   * This function iterates through the parallel categories determined by the
   * `TaskInfoType`. For each category, it performs necessary communication
   * (e.g., updating ghost values) and then calls the `patch_worker` with the
   * range of regular patches belonging to that category.
   *
   * @tparam OutVector The type of the output/solution vector.
   * @tparam InVector The type of the input/rhs vector.
   * @param patch_worker The function to execute for each patch range. It
   * takes the `PatchStorage` instance, output vector, input vector, and the
   * `PatchRange` as arguments.
   * @param solution The output/solution vector. Its ghost values will be
   * updated during the loop.
   * @param rhs The input/rhs vector. Its ghost values must be up-to-date
   * before calling this function.
   * @param do_forward Flag indicating the direction of sweep (currently only
   * `true` is supported).
   */
  template <class OutVector, class InVector>
  void
  patch_loop(const std::function<void(const PatchStorage<MFType> &,
                                      OutVector &,
                                      const InVector &,
                                      const PatchRange &)> &patch_worker,
             OutVector                                     &solution,
             const InVector                                &rhs,
             const bool &do_forward = true) const;


  /**
   *  Helper function to adjusts the partitioner of a
   * distributed vector if it doesn't match the internal partitioner for the
   * given component.
   *
   * This ensures that ghost value communication works correctly within the
   * `patch_loop`. If the vector's partitioner is different, the vector is
   * reinitialized with the correct partitioner, preserving its locally owned
   * data.
   *
   * @tparam number The value type of the vector.
   * @param component The component index for which to check the partitioner.
   * @param vec The distributed vector to potentially adjust.
   */
  template <typename number>
  void
  adjust_ghost_range_if_necessary(
    const unsigned int                                component,
    const LinearAlgebra::distributed::Vector<number> &vec) const;


  /**
   * Gets a constant reference to the regular patch at the given global
   * index `i`. Used by FEPatchEvaluation
   *
   * The index `i` ranges from 0 to `n_patches() - 1`.
   * @param i The global index of the regular patch.
   * @return A constant reference to the `RegularPatch`.
   */
  const RegularPatch &
  get_regular_patch(const std::size_t &i) const;


  /**
   *  Gets the category assigned to the regular patch at the given global
   * index `i`.
   *
   * Requires `categorize_patches()` to have been called first.
   * @param i The global index of the regular patch.
   * @return A constant reference to the `PatchCategory`.
   */
  const PatchCategory &
  get_regular_patch_category(const std::size_t &i) const;


  /**
   *  Returns the total number of regular patches stored on the current
   * MPI process.
   */
  std::size_t
  n_patches() const;


  /**
   *  Gets a constant reference to the shared pointer holding the
   * MatrixFree object.
   */
  const std::shared_ptr<const MatrixFreeType> &
  get_matrix_free() const;


  /**
   *  Converts a `CellIndex` back to a `CellIterator`.
   * @param index The `CellIndex` to convert.
   * @return The corresponding `CellIterator`.
   */
  inline CellIterator
  index2cell(const CellIndex &index) const;


  /**
   *  Converts a MatrixFree batch index and lane index into a unique
   * `CellIndex`.
   * @param cell_batch_index The index of the cell batch in the MatrixFree
   * object.
   * @param lane_index The index of the lane within the batch.
   * @return The calculated `CellIndex`.
   */
  inline CellIndex
  batch2index(const unsigned int &cell_batch_index,
              const unsigned int &lane_index) const;


  /**
   *  Converts a `CellIndex` back into its corresponding MatrixFree batch
   * index and lane index.
   * @param cell_index The `CellIndex` to convert.
   * @return A pair containing the batch index (first) and lane index (second).
   */
  inline std::pair<unsigned int, unsigned int>
  index2batch(const CellIndex &cell_index) const;


  /**
   *  Outputs the geometry of all stored regular patches to VTU/PVTU
   * files.
   *
   * Each cell within each patch is written as a separate element in the VTU
   * file. Data associated with each cell includes patch index, parallel
   * category, local cell index within the patch, global cell index, MPI rank,
   * and user-defined category (if available).
   *
   * @param filename_without_extension The base name for the output files (e.g.,
   * "patches"). `.procXXXX.vtu` and `.pvtu` will be appended.
   */
  void
  output_patches(const std::string &filename_without_extension) const;


  /**
   *  Outputs the center points (approximated by one vertex) of all
   * stored regular patches to VTU/PVTU files.
   *
   * Each patch is represented by a single point (vertex) in the output file.
   * Data associated with each point includes patch index, parallel category,
   * MPI rank, and user-defined category (if available).
   *
   * @param filename_without_extension The base name for the output files (e.g.,
   * "patch_centers"). `.procXXXX.vtu` and `.pvtu` will be appended.
   */
  void
  output_centerpoints(const std::string &filename_without_extension) const;

private:
  std::shared_ptr<const MatrixFreeType> matrix_free;
  const Triangulation<dim>             &triangulation;

  const MPI_Comm     mpi_communicator;
  const unsigned int this_mpi_process;
  const unsigned int n_mpi_process;

  const unsigned int level;
  const unsigned int n_components;



  std::map<types::global_vertex_index, std::set<CellIndex>>
  generate_patches();


  void
  push_back_patch(const std::set<CellIndex>        &patch,
                  const types::global_vertex_index &vertex_index);


  inline std::set<types::global_dof_index>
  collect_patch_dof_indices(const std::set<CellIndex> &patch_cells,
                            const unsigned int        &component) const;


  // patches[parallel cat][patch_index]
  std::array<std::vector<RegularPatch>, TaskInfoType::n_categories>

    regular_patches;
  std::array<std::vector<PatchCategory>, TaskInfoType::n_categories>
    regular_patch_categories;
  std::array<std::vector<GeneralPatch>, TaskInfoType::n_categories>
    other_patches;



  std::vector<std::shared_ptr<const Utilities::MPI::Partitioner>> partitioners;
  std::vector<unsigned int>         process_colors;
  mutable std::vector<TaskInfoType> task_infos;


  AdditionalData additional_data;

  bool is_initialized;
  bool are_patches_categorized;
};


#ifndef DOXYGEN
//---------------------------------------------------------------------------


// ============================================================================
// RegularVertexPatch constructor definitions
// ============================================================================
// dim==2
template <int dim>
RegularVertexPatch<dim>::RegularVertexPatch(
  const std::set<CellIndex>                            &patch,
  const types::global_vertex_index                     &vertex_index,
  const std::function<CellIterator(const CellIndex &)> &index2cell)
{
  if constexpr (dim == 2)
    {
      std::map<CellIterator, CellIndex> iterator2index;
      std::set<CellIterator>            cells_iterators;

      for (auto &cell_index : patch)
        {
          iterator2index[index2cell(cell_index)] = cell_index;
          cells_iterators.insert(index2cell(cell_index));
        }

      std::vector<CellIterator> ordered_patch_cells =
        internal::order_patch<dim>(cells_iterators, vertex_index);

      internal::reorder(ordered_patch_cells, {3, 2, 1, 0});

      std::vector<CellIndex> ordered_patch;
      for (auto &cell : ordered_patch_cells)
        ordered_patch.push_back(iterator2index.at(cell));


      AssertDimension(ordered_patch.size(), ordered_patch_cells.size());

      for (unsigned int i = 0; i < n_cells; ++i)
        cells[i] = ordered_patch[i];


      partially_ghosted = false;
      for (auto &cell : ordered_patch_cells)
        if (!cell->is_locally_owned_on_level())
          partially_ghosted = true;
    }

  if constexpr (dim == 3)
    {
      const static std::array<std::size_t, n_cells> vindex2position = {
        {7, 6, 5, 4, 3, 2, 1}}; // TODO: Check this order {7, 6, 5, 4, 3, 2,
                                // 1, 0}?

      std::vector<bool> used_cell(n_cells, false);

      for (const auto cell_index : patch)
        {
          const auto cell = index2cell(cell_index);
          const auto v_index =
            internal::compute_vertex_index<dim>(cell, vertex_index);

          const auto &cell__inpatch_index = vindex2position[v_index];
          Assert(
            used_cell[cell__inpatch_index] == false,
            ExcMessage(
              "You have tried to construct patch from rotated cells, that is currently not implemented"));
          used_cell[cell__inpatch_index] = true;
          cells[cell__inpatch_index]     = cell_index;
        }

      partially_ghosted = false;
      for (auto &cell_index : cells)
        if (!index2cell(cell_index)->is_locally_owned_on_level())
          partially_ghosted = true;
    }
}

// ============================================================================
// GeneralVertexPatch constructor definition
// ============================================================================
template <int dim>
GeneralVertexPatch<dim>::GeneralVertexPatch(
  const std::set<CellIndex> &patch,
  const types::global_vertex_index & /*vertex_index*/,
  const std::function<CellIterator(const CellIndex &)> & /*index2cell*/)
{
  // Example implementation: just copy the cell indices
  cells.assign(patch.begin(), patch.end());
  // TODO: Add logic to determine if partially ghosted, etc.
}


// ============================================================================
// PatchStorage member function definitions
// ============================================================================

template <class MFType>
PatchStorage<MFType>::PatchStorage(
  const std::shared_ptr<const MatrixFreeType> &mf)
  : EnableObserverPointer()
  , matrix_free(mf)
  , triangulation(matrix_free->get_dof_handler().get_triangulation())
  , mpi_communicator(
      matrix_free->get_vector_partitioner()->get_mpi_communicator())
  , this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator))
  , n_mpi_process(Utilities::MPI::n_mpi_processes(mpi_communicator))
  , level(matrix_free->get_mg_level())
  , n_components(matrix_free->n_components())
  , task_infos(n_components)
  , is_initialized(false)
  , are_patches_categorized(false)
{
  Assert(matrix_free->get_dof_handler().get_triangulation().n_global_levels() >
           0,
         ExcInternalError());
  Assert(level != numbers::invalid_unsigned_int, ExcInternalError());
}



template <class MFType>
inline typename PatchStorage<MFType>::CellIterator
PatchStorage<MFType>::index2cell(const CellIndex &index) const
{
  const unsigned int cell_batch_index = index / n_lanes;
  const unsigned int lane_index       = index % n_lanes;
  return matrix_free->get_cell_iterator(cell_batch_index, lane_index);
}



template <class MFType>
inline typename PatchStorage<MFType>::CellIndex
PatchStorage<MFType>::batch2index(const unsigned int &cell_batch_index,
                                  const unsigned int &lane_index) const
{
  return cell_batch_index * n_lanes + lane_index;
}



template <class MFType>
inline std::pair<unsigned int, unsigned int>
PatchStorage<MFType>::index2batch(const CellIndex &cell_index) const
{
  const unsigned int cell_batch_index = cell_index / n_lanes;
  const unsigned int lane_index       = cell_index % n_lanes;
  return std::make_pair(cell_batch_index, lane_index);
}



template <class MFType>
inline void
PatchStorage<MFType>::initialize(const AdditionalData &data)
{
  this->additional_data = data;
  // TODO: check if all DoFHandlers have the same triangulation
  // Assert....

  Assert(static_cast<unsigned int>(level) < triangulation.n_global_levels(),
         ExcInternalError());



  std::map<types::global_vertex_index, std::set<CellIndex>> vertex_to_cell_map =
    generate_patches();


  // remove patches with no owned cells
  // and patches that other processes takes care of:
  // If patch is shared among several processes,
  // the processor with lowest ID takes whole patch.
  // That is not optimal TODO: improve that
  // Count  patches.
  {
    unsigned int n_patches = 0;
    for (auto iterator = vertex_to_cell_map.begin();
         iterator != vertex_to_cell_map.end();
         ++iterator)
      {
        std::set<CellIndex> &patch = iterator->second;

        std::set<CellIterator> patch_cells;
        for (const CellIndex &cells_indices : patch)
          patch_cells.insert(index2cell(cells_indices));


        bool         has_owned = false;
        unsigned int owned_by  = numbers::invalid_unsigned_int;
        for (const CellIterator &cell : patch_cells)
          {
            if (cell->is_locally_owned_on_level())
              {
                has_owned = true;
              }
            if (cell->level_subdomain_id() < owned_by &&
                cell->level_subdomain_id() <= n_mpi_process)
              {
                owned_by = cell->level_subdomain_id();
              }
          }
        if (false == has_owned ||
            Utilities::MPI::this_mpi_process(mpi_communicator) != owned_by)
          patch.clear();
        else
          ++n_patches;
      }
    Assert(n_patches <= triangulation.n_vertices(), ExcInternalError());
  }


  // compute influence and dependency regions
  for (unsigned int component = 0; component < n_components; ++component)
    {
      const auto &dof_handler = matrix_free->get_dof_handler(component);

      IndexSet locally_relevant_dofs;
      DoFTools::extract_locally_relevant_level_dofs(dof_handler,
                                                    level,
                                                    locally_relevant_dofs);



      std::set<types::global_dof_index> locally_relevant_dofs_set;
      for (auto [vertex_index, patch] : vertex_to_cell_map)
        {
          const auto local_dofs = collect_patch_dof_indices(patch, component);
          locally_relevant_dofs_set.insert(local_dofs.begin(),
                                           local_dofs.end());
        }

      IndexSet dependency_region(dof_handler.n_dofs());
      dependency_region.add_indices(locally_relevant_dofs_set.begin(),
                                    locally_relevant_dofs_set.end());

      // TODO: is that true?
      IndexSet influence_region = dependency_region;

      partitioners.push_back(std::make_shared<Utilities::MPI::Partitioner>(
        dof_handler.locally_owned_mg_dofs(level),
        locally_relevant_dofs,
        mpi_communicator));
      // TODO: direction of sweep: Do_forward etc
      task_infos[component].initialize(*matrix_free,
                                       dependency_region,
                                       influence_region,
                                       process_colors,
                                       true);
    }


  for (auto iterator = vertex_to_cell_map.begin();
       iterator != vertex_to_cell_map.end();
       ++iterator)
    {
      std::set<CellIndex>              &patch        = iterator->second;
      const types::global_vertex_index &vertex_index = iterator->first;
      const unsigned int                n_cells      = patch.size();

      // Drop empty patches
      if (n_cells == 0)
        continue;


      push_back_patch(patch, vertex_index);
    }
  // sort following the order of MatrixFree object to increase cache
  // efficiency
  for (unsigned int i = 0; i < TaskInfoType::n_categories; ++i)
    std::sort(regular_patches[i].begin(),
              regular_patches[i].end(),
              [&](const auto a, const auto b) {
                return a.get_cells() < b.get_cells();
              });

  is_initialized = true;
}



template <class MFType>
std::map<types::global_vertex_index,
         std::set<typename PatchStorage<MFType>::CellIndex>>
PatchStorage<MFType>::generate_patches()
{
  //    unsigned int n_vertex_patches = 0;
  std::map<types::global_vertex_index, std::set<CellIndex>> vertex_to_cell_map;

  for (unsigned int i = 0;
       i < matrix_free->n_cell_batches() + matrix_free->n_ghost_cell_batches();
       ++i)
    for (unsigned int j = 0;
         j < matrix_free->n_active_entries_per_cell_batch(i);
         ++j)
      {
        const auto &cell = matrix_free->get_cell_iterator(i, j);

        for (const auto v : cell->vertex_indices())
          vertex_to_cell_map[cell->vertex_index(v)].insert(batch2index(i, j));
      }

  return vertex_to_cell_map;
}



template <class MFType>
inline void
PatchStorage<MFType>::push_back_patch(
  const std::set<typename PatchStorage<MFType>::CellIndex> &patch,
  const types::global_vertex_index                         &vertex_index)
{
  const std::function<CellIterator(const CellIndex &)> index2cell_local =
    [&](const CellIndex &index) -> CellIterator {
    return this->index2cell(index);
  };

  unsigned int category = TaskInfoType::invalid_category;
  for (unsigned int component = 0; component < n_components; ++component)
    {
      const auto local_dofs = collect_patch_dof_indices(patch, component);

      category = internal::GaussSeidel::min_category(
        category,
        task_infos[component].determine_minimum_category_of_patch(
          local_dofs, process_colors));
    }
  const unsigned int partition =
    TaskInfoType::determine_cell_partition(category);

  if (RegularPatch::is_constructible(patch, vertex_index, index2cell_local))
    {
      RegularPatch ordered_patch(patch, vertex_index, index2cell_local);
      regular_patches[partition].push_back(ordered_patch);
    }
  else
    {
      // Construct and push back a GeneralPatch
      GeneralPatch general_patch(patch, vertex_index, index2cell_local);
      other_patches[partition].push_back(general_patch);
    }
}



template <class MFType>
inline std::set<types::global_dof_index>
PatchStorage<MFType>::collect_patch_dof_indices(
  const std::set<typename PatchStorage<MFType>::CellIndex> &patch_cells,
  const unsigned int                                       &component) const
{
  const auto &dof_handler = matrix_free->get_dof_handler(component);
  std::set<types::global_dof_index> local_dofs;

  for (auto cell_index : patch_cells)
    {
      const auto &cell = index2cell(cell_index);
      typename DoFHandler<dim>::level_cell_iterator dof_cell(&triangulation,
                                                             cell->level(),
                                                             cell->index(),
                                                             &dof_handler);

      std::vector<types::global_dof_index> local_dof_indices(
        dof_cell->get_fe().n_dofs_per_cell());
      dof_cell->get_mg_dof_indices(local_dof_indices);
      local_dofs.insert(local_dof_indices.begin(), local_dof_indices.end());
    }
  return local_dofs;
}



template <class MFType>
void
PatchStorage<MFType>::categorize_patches(
  std::function<PatchCategory(const RegularPatch &)> category_function)
{
  Assert(is_initialized, ExcNotInitialized());

  for (unsigned int cat = 0; cat < TaskInfoType::n_categories; ++cat)
    {
      regular_patch_categories[cat].reserve(regular_patches[cat].size());
      for (const auto &patch : regular_patches[cat])
        {
          regular_patch_categories[cat].push_back(category_function(patch));
        }
    }
  are_patches_categorized = true;
}



template <class MFType>
template <typename number>
void
PatchStorage<MFType>::adjust_ghost_range_if_necessary(
  const unsigned int                                component,
  const LinearAlgebra::distributed::Vector<number> &vec) const
{
  if (vec.get_partitioner().get() == partitioners[component].get())
    return;
  LinearAlgebra::distributed::Vector<number> copy_vec(vec);
  const_cast<LinearAlgebra::distributed::Vector<number> &>(vec).reinit(
    partitioners[component]);
  const_cast<LinearAlgebra::distributed::Vector<number> &>(vec)
    .copy_locally_owned_data_from(copy_vec);
}



template <class MFType>
template <class OutVector, class InVector>
inline void
PatchStorage<MFType>::patch_loop(
  const std::function<void(const PatchStorage<MFType> &,
                           OutVector &,
                           const InVector &,
                           const PatchRange &)> &patch_worker,
  OutVector                                     &solution,
  const InVector                                &rhs,
  const bool                                    &do_forward) const
{
  Assert(is_initialized, ExcNotInitialized());
  Assert(do_forward == true, ExcNotImplemented());
  (void)do_forward;

  // fixme: if block, loop over blocks:
  // fixme: remove from here, only assert that solution and rhs are compatible
  adjust_ghost_range_if_necessary(0, rhs);
  adjust_ghost_range_if_necessary(0, solution);

  rhs.update_ghost_values();
  solution.set_ghost_state(true);

  solution.update_ghost_values();

  std::size_t current_begin = 0;


  for (unsigned int cat = 0; cat < TaskInfoType::n_categories; ++cat)
    {
      for (unsigned int component = 0; component < n_components; ++component)
        {
          const int communication_channel = component; // TODO
          task_infos[component].communicate(cat,
                                            communication_channel,
                                            solution);
        }
      if (regular_patches[cat].size() != 0)
        {
          PatchRange patch_range(current_begin,
                                 current_begin + regular_patches[cat].size());
          patch_worker(*this, solution, rhs, patch_range);

          current_begin += regular_patches[cat].size();
        }
    }

  for (unsigned int component = 0; component < n_components; ++component)
    {
      const int communication_channel = component; // TODO
      task_infos[component].communicate(TaskInfoType::n_categories,
                                        communication_channel,
                                        solution);
    }


  solution.zero_out_ghost_values();
  //    solution.compress(VectorOperation::insert);
}



template <class MFType>
const typename PatchStorage<MFType>::RegularPatch &
PatchStorage<MFType>::get_regular_patch(const std::size_t &i) const
{
  std::size_t ii = i;
  for (unsigned int cat_index = 0; cat_index < regular_patches.size();
       ++cat_index)
    {
      const auto &patch_cat = regular_patches[cat_index];
      if (ii < patch_cat.size())
        return patch_cat[ii];
      else
        ii -= patch_cat.size();
    }

  Assert(false, ExcInternalError());
  // Need to return something, but this path should not be reached.
  // Returning the first patch of the first category as a fallback.
  return regular_patches[0][0];
}

template <class MFType>
std::size_t
PatchStorage<MFType>::n_patches() const
{
  std::size_t n_patches = 0;
  for (const auto &patch_cat : regular_patches)
    n_patches += patch_cat.size();
  return n_patches;
}

template <class MFType>
const std::shared_ptr<const typename PatchStorage<MFType>::MatrixFreeType> &
PatchStorage<MFType>::get_matrix_free() const
{
  return matrix_free;
}

template <class MFType>
const typename PatchStorage<MFType>::PatchCategory &
PatchStorage<MFType>::get_regular_patch_category(const std::size_t &i) const
{
  Assert(are_patches_categorized, ExcNotInitialized());
  AssertDimension(regular_patches.size(), regular_patch_categories.size());

  std::size_t ii = i;
  for (unsigned int cat_index = 0; cat_index < regular_patches.size();
       ++cat_index)
    {
      AssertDimension(regular_patches[cat_index].size(),
                      regular_patch_categories[cat_index].size());
      const auto &patch_cat = regular_patch_categories[cat_index];
      if (ii < patch_cat.size())
        return patch_cat[ii];
      else
        ii -= patch_cat.size();
    }

  AssertThrow(false, ExcInternalError());
  return regular_patch_categories[0][0];
}

template <class MFType>
void
PatchStorage<MFType>::output_patches(
  const std::string &filename_without_extension) const
{
  using CellOutData = DataOutBase::Patch<dim, dim>;
  std::vector<CellOutData> patches_out;

  std::vector<std::string> data_names;
  data_names.emplace_back("patch_index");
  data_names.emplace_back("parallel category");
  data_names.emplace_back("local cell index");
  data_names.emplace_back("global cell index");
  data_names.emplace_back("MPIRank");
  data_names.emplace_back("Category");

  const unsigned n_datasets          = 6;
  unsigned int   patch_counter       = 0;
  unsigned int   patch_index_counter = 0;
  for (unsigned int cat_index = 0; cat_index < regular_patches.size();
       ++cat_index)
    {
      const auto &patch_cat = regular_patches[cat_index];
      for (unsigned int patch_index = 0; patch_index < patch_cat.size();
           ++patch_index)
        {
          const auto &patch = patch_cat[patch_index];
          for (unsigned int cell_index_within_patch = 0;
               cell_index_within_patch < patch.get_cells().size();
               ++cell_index_within_patch)
            {
              const auto &cell_index =
                patch.get_cells()[cell_index_within_patch];
              const auto &cell = index2cell(cell_index);
              CellOutData cell_out;
              for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_cell;
                   ++i)
                cell_out.vertices[i] = cell->vertex(i);
              cell_out.patch_index    = patch_counter;
              cell_out.reference_cell = cell->reference_cell();
              cell_out.data.reinit(n_datasets, 1 << dim);
              for (unsigned int i = 0; i < 1 << dim; ++i)
                {
                  cell_out.data(0, i) = patch_index_counter;
                  cell_out.data(1, i) = cat_index;
                  cell_out.data(2, i) = cell_index_within_patch;
                  cell_out.data(3, i) = cell_index;
                  cell_out.data(4, i) = this_mpi_process;
                  if (are_patches_categorized)
                    cell_out.data(5, i) =
                      get_regular_patch_category(patch_index_counter);
                  else
                    cell_out.data(5, i) =
                      std::numeric_limits<double>::quiet_NaN();
                }

              cell_out.n_subdivisions = 1;
              patches_out.push_back(cell_out);
              patch_counter++;
            }
          ++patch_index_counter;
        }
    }

  std::vector<std::string> piece_names(n_mpi_process);
  for (unsigned int i = 0; i < n_mpi_process; ++i)
    piece_names[i] = filename_without_extension + ".proc" +
                     Utilities::int_to_string(i, 4) + ".vtu";
  std::string new_file = piece_names[this_mpi_process];

  std::string out_pvtu = filename_without_extension + ".pvtu";

  std::ofstream out(new_file);
  std::vector<
    std::tuple<unsigned int,
               unsigned int,
               std::string,
               DataComponentInterpretation::DataComponentInterpretation>>
    vector_data_ranges;

  DataOutBase::VtkFlags vtu_flags;

  DataOutBase::write_vtu(
    patches_out, data_names, vector_data_ranges, vtu_flags, out);

  if (this_mpi_process == 0)
    {
      std::ofstream pvtu_output(out_pvtu);
      std::ostream &pvtu_out_steam = pvtu_output;
      DataOutBase::write_pvtu_record(
        pvtu_out_steam, piece_names, data_names, vector_data_ranges, vtu_flags);
    }
}



template <class MFType>
void
PatchStorage<MFType>::output_centerpoints(
  const std::string &filename_without_extension) const
{
  using PointOutData = DataOutBase::Patch<0, dim>;
  std::vector<PointOutData> centerpoints_out;

  std::vector<std::string> data_names;
  data_names.emplace_back("patch_index");
  data_names.emplace_back("parallel category");
  data_names.emplace_back("MPIRank");
  data_names.emplace_back("Category");

  const unsigned n_datasets    = 4;
  unsigned int   patch_counter = 0;

  for (unsigned int cat_index = 0; cat_index < regular_patches.size();
       ++cat_index)
    {
      const auto &patch_cat = regular_patches[cat_index];
      for (unsigned int patch_index = 0; patch_index < patch_cat.size();
           ++patch_index)
        {
          const auto        &patch                   = patch_cat[patch_index];
          const unsigned int cell_index_within_patch = 0;
          const auto  &cell_index = patch.get_cells()[cell_index_within_patch];
          const auto  &cell       = index2cell(cell_index);
          PointOutData vertex_out;

          unsigned int last_vertex_index =
            GeometryInfo<dim>::vertices_per_cell - 1;

          // fixme: replace with finding the common vertex
          vertex_out.vertices[0] = cell->vertex(last_vertex_index);

          vertex_out.patch_index = patch_counter;
          // vertex_out.reference_cell = ReferenceCells::Vertex.get_type();

          vertex_out.data.reinit(n_datasets, 1);

          vertex_out.data(0, 0) = patch_counter;
          vertex_out.data(1, 0) = cat_index;
          vertex_out.data(2, 0) = this_mpi_process;
          if (are_patches_categorized)
            vertex_out.data(3, 0) = get_regular_patch_category(patch_counter);
          else
            vertex_out.data(3, 0) = std::numeric_limits<double>::quiet_NaN();

          // vertex_out.n_subdivisions = 1;
          centerpoints_out.push_back(vertex_out);
          ++patch_counter;
        }
    }

  std::vector<std::string> piece_names(n_mpi_process);
  for (unsigned int i = 0; i < n_mpi_process; ++i)
    piece_names[i] = filename_without_extension + ".proc" +
                     Utilities::int_to_string(i, 4) + ".vtu";
  std::string new_file = piece_names[this_mpi_process];

  std::string out_pvtu = filename_without_extension + ".pvtu";


  std::ofstream out(new_file);
  std::vector<
    std::tuple<unsigned int,
               unsigned int,
               std::string,
               DataComponentInterpretation::DataComponentInterpretation>>
    vector_data_ranges;

  DataOutBase::VtkFlags vtu_flags;

  DataOutBase::write_vtu(
    centerpoints_out, data_names, vector_data_ranges, vtu_flags, out);
  if (this_mpi_process == 0)
    {
      std::ofstream pvtu_output(out_pvtu);
      std::ostream &pvtu_out_steam = pvtu_output;
      DataOutBase::write_pvtu_record(
        pvtu_out_steam, piece_names, data_names, vector_data_ranges, vtu_flags);
    }
}

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE



#endif // dealii__patch_storage_h
