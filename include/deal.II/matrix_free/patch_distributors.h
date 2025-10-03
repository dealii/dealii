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


#ifndef dealii__patch_distributors
#define dealii__patch_distributors

#include <deal.II/base/config.h>

#include <array>
#include <utility> // Required for std::pair
#include <vector>



DEAL_II_NAMESPACE_OPEN


/**
 * @brief This namespace provides tools for distributing and gathering data
 * between patch vectors and local cell storage in the context of matrix-free
 * vertex-patch smoothers
 */
namespace PatchDistributors
{

  /**
   * @brief Dynamic patch distributor for matrix-free patch operations.
   *
   * This class provides tools that
   * describe the layout of local patches (collections of 2^dim neighbor cells)
   * used by FEPatchEvaluation. The implementation aims to avoid
   * large precomputed index tables by evaluating global/local indices on the
   * fly and to support arbitrary polynomial degrees >= 1.
   *
   *
   * All sizes are compile-time constants where feasible. The distributor
   * does not precompute or store large index tables: global/local indices
   * for patch entries are evaluated on the fly, keeping memory usage low.
   * It supports arbitrary polynomial degree values (subject to the
   * degree >= 1 assertion).
   *
   */
  template <int dim, int degree>
  class Dynamic
  {
  public:
    const static int dimension = dim;

    static_assert(dim == 1 || dim == 2 || dim == 3,
                  "This distributor only supports dim = 1, 2, or 3");
    static_assert(degree >= 1, "This distributor only supports degree >= 1");

    const constexpr static std::size_t n_cells_1d          = 2;
    const constexpr static std::size_t n_dofs_per_patch_1d = 2 * degree - 1;
    const constexpr static std::size_t n_dofs_per_line     = degree + 1;

    const constexpr static std::size_t skip_size = n_dofs_per_line - 2;

    /**
     * Iterate over all patches managed by this distributor and invoke the
     * appropriate user-provided callback for each patch.
     *
     * The function visits every patch known to the distributor and dispatches
     * it to one of the three callables depending on the patch state:
     * - regular_operation: called for patches that should be processed normally
     *   (typically locally owned patches),
     * - duplicates_operation: called for patches that represent
     * duplicated/shared data and require duplicate-specific handling,
     * - skipped_operation: called for patches that must be ignored/skipped.
     *
     * Each argument must be a callable (function, functor or lambda) compatible
     * with the distributor's patch descriptor. The exact parameter types and
     * return values are determined by those callables. The iteration order is
     * unspecified. This method is const and does not modify the distributor.
     */
    template <typename RegularOperation,
              typename DuplicatesOperation,
              typename SkippedOperation>
    void
    loop(const RegularOperation    &regular_operation,
         const DuplicatesOperation &duplicates_operation,
         const SkippedOperation    &skipped_operation) const;



    /**
     * Return the number of degrees of freedom per patch.
     *
     * This value is a compile-time constant.
     */
    constexpr std::size_t
    n_patch_dofs() const
    {
      return n_patch_dofs_static;
    }

    const static unsigned int n_patch_dofs_static =
      Utilities::fixed_power<dim>(n_dofs_per_patch_1d);


  private:
    /**
     *
     * The base case for the recursion is a 1D entry (a line). For each degree
     * of freedom on the line, it determines if the DoF is regular, a
     * duplicate, or should be skipped, and then invokes the corresponding
     * operation (`regular_operation`, `duplicates_operation`, or
     * `skipped_operation`).
     *
     * For higher-dimensional entries (dim > 1), the function recursively calls
     * itself for each of the sub-entities that form the boundary of the
     * current entry. For a plane it will loop over all lines in the plane.
     *
     * The function updates the `dof_cell` and `dof_patch` counters as it
     * processes the degrees of freedom.
     *
     * @param[in] regular_operation Functor to be called for regular DoFs.
     * @param[in] duplicates_operation Functor to be called for duplicate DoFs.
     * @param[in] skipped_operation Functor to be called for skipped DoFs.
     * @param[in] cell_index The index of the cell being processed.
     * @param[in,out] dof_cell A running counter for the DoF index within the
     *                         cell.
     * @param[in,out] dof_patch A running counter for the DoF index within the
     *                          patch.
     * @param[in] skip_only A flag that, if true,indicates that all DoF on the
     *                     current entry should be treated as skipped. This is
     * used to optimize the recursion for entries that are fully skipped.
     */
    template <int entry_dim,
              typename RegularOperation,
              typename DuplicatesOperation,
              typename SkippedOperation>
    void
    do_entry(const DuplicatesOperation &regular_operation,
             const SkippedOperation    &duplicates_operation,
             const RegularOperation    &skipped_operation,
             const std::size_t         &cell_index,
             std::size_t               &dof_cell,
             std::size_t               &dof_patch,
             const bool                &skip_only) const;
  };



  namespace internal
  {
    template <int dim, int degree>
    struct CellPatchLookup
    {};
  } // namespace internal


  /**
   * @brief Provides lookup tables and methods for distributing and gathering
   * data between a patch vector and local cell storage.
   *
   * This struct relies on a specialization of `CellPatchLookup` for the given
   * `dim` and `degree` to provide the necessary mapping information between
   * patch DoFs and cell-local DoFs. It offers methods to loop over these
   * mappings and perform common operations like distributing patch data to
   * local cells or gathering local cell data into a patch vector.
   *
   * @tparam dim The spatial dimension.
   * @tparam degree The polynomial degree of the finite element.
   */
  template <int dim, int degree>
  struct Lookup
  {
    using Cell2Patch = internal::CellPatchLookup<dim, degree>;

    static const constexpr unsigned int n_cells = Cell2Patch::n_cells;
    static const constexpr unsigned int n_patch_dofs_static =
      Cell2Patch::n_patch_dofs;


    /**
     * Loops over the cell-to-patch mappings, applying different
     * operations based on whether a DoF is non-overlapping, overlapping, or
     * skipped.
     *
     * This is the core loop used by `distribute_patch_to_local` and
     * `gather_local_to_patch` in FEPatchEvaluation. It iterates through all
     * cells in the patch. For each cell, it first calls `regular_operation` for
     * DoFs that are unique to this cell within the patch context
     * (non-overlapping). Then, it calls `overlapping_operation` for DoFs that
     * are shared with other cells in the patch. Finally, it iterates through
     * all cells again and calls `skipped_operation` for DoFs that exist in the
     * cell's local storage but are not part of the patch vector (e.g., interior
     * DoFs when the patch represents a vertex or edge).
     *
     * @tparam RegularOperation Functor type for non-overlapping DoFs. Must be
     * callable with `(unsigned int patch_index, unsigned int cell, unsigned int
     * cell_index)`.
     * @tparam OverlappingOperation Functor type for overlapping DoFs. Must be
     * callable with `(unsigned int patch_index, unsigned int cell, unsigned int
     * cell_index)`.
     * @tparam SkippedOperation Functor type for skipped DoFs. Must be callable
     * with `(unsigned int cell, unsigned int cell_index)`.
     * @param regular_operation Operation applied to non-overlapping DoFs and first occurrences
     * of shared ones.
     * @param overlapping_operation Operation applied to subsequent occurrences of DoFs that
     * are shared between several cells.
     * @param skipped_operation Operation applied to DoFs that are not included in the patch.
     */
    template <typename RegularOperation,
              typename OverlappingOperation,
              typename SkippedOperation>
    inline void
    loop(const RegularOperation     &regular_operation,
         const OverlappingOperation &overlapping_operation,
         const SkippedOperation     &skipped_operation) const
    {
      for (unsigned int c = 0; c < Cell2Patch::n_cells; ++c)
        {
          unsigned int i = 0;
          for (; i < Cell2Patch::n_non_overlapping[c]; ++i)
            {
              const auto c2p = Cell2Patch::cell_to_patch[c][i];
              regular_operation(c2p.second, c, c2p.first);
            }
          for (; i < Cell2Patch::n_cell_to_patch[c]; ++i)
            {
              const auto c2p = Cell2Patch::cell_to_patch[c][i];
              overlapping_operation(c2p.second, c, c2p.first);
            }
        }
      for (unsigned int c = 0; c < Cell2Patch::n_cells; ++c)
        for (unsigned int i = 0; i < Cell2Patch::n_skipped_dofs[c]; ++i)
          {
            const auto &cell_dof = Cell2Patch::skipped_dofs[c][i];
            skipped_operation(c, cell_dof);
          }
    }


    /**
     * @brief Returns the total number of degrees of freedom in the patch vector.
     *
     * This corresponds to the size of the `patch_storage` expected by the
     * `distribute_patch_to_local` and `gather_local_to_patch` methods.
     */
    unsigned int
    n_patch_dofs() const
    {
      return n_patch_dofs_static;
    }
  };


//---------------------------------------------------------------------------
#ifndef DOXYGEN



  template <int dim, int degree>
  template <int entry_dim,
            typename RegularOperation,
            typename DuplicatesOperation,
            typename SkippedOperation>
  void
  Dynamic<dim, degree>::do_entry(const DuplicatesOperation &regular_operation,
                                 const SkippedOperation &duplicates_operation,
                                 const RegularOperation &skipped_operation,
                                 const std::size_t      &cell_index,
                                 std::size_t            &dof_cell,
                                 std::size_t            &dof_patch,
                                 const bool             &skip_only) const
  {
    const bool is_boundary_left = (cell_index & (1 << (entry_dim - 1))) == 0;



    const constexpr std::size_t patch_skip = []() constexpr -> std::size_t {
      static_assert(entry_dim >= 1 && entry_dim <= dim,
                    "Implemented for entry_dim in [1, dim]");

      if constexpr (entry_dim == 1)
        return static_cast<std::size_t>(1);
      else if constexpr (entry_dim == 2)
        return skip_size;

      return skip_size * n_dofs_per_patch_1d;
    }();



    auto boundary_operation = [&]() {
      if constexpr (entry_dim == 1)
        {
          skipped_operation(cell_index, dof_cell);
          ++dof_cell;
        }
      else
        {
          do_entry<entry_dim - 1>(regular_operation,
                                  duplicates_operation,
                                  skipped_operation,
                                  cell_index,
                                  dof_cell,
                                  dof_patch,
                                  true);
        }
    };

    auto middle_operation = [&]() {
      if constexpr (entry_dim == 1)
        {
          if (skip_only)
            skipped_operation(cell_index, dof_cell);
          else
            {
              regular_operation(dof_patch, cell_index, dof_cell);
              dof_patch += patch_skip;
            }
          ++dof_cell;
        }
      else
        {
          do_entry<entry_dim - 1>(regular_operation,
                                  duplicates_operation,
                                  skipped_operation,
                                  cell_index,
                                  dof_cell,
                                  dof_patch,
                                  skip_only);
          if (!skip_only)
            dof_patch += patch_skip;
        }
    };

    auto overlap_operation = [&]() {
      if constexpr (entry_dim == 1)
        {
          if (skip_only)
            skipped_operation(cell_index, dof_cell);
          else
            {
              duplicates_operation(dof_patch, cell_index, dof_cell);
              dof_patch += patch_skip;
            }
          ++dof_cell;
        }
      else
        {
          do_entry<entry_dim - 1>(duplicates_operation,
                                  duplicates_operation,
                                  skipped_operation,
                                  cell_index,
                                  dof_cell,
                                  dof_patch,
                                  skip_only);
          if (!skip_only)
            dof_patch += patch_skip;
        }
    };


    if (is_boundary_left)
      boundary_operation();
    else
      overlap_operation();

    for (std::size_t d = 1; d < n_dofs_per_line - 1; ++d)
      middle_operation();

    if (is_boundary_left)
      middle_operation();
    else
      boundary_operation();
  }



  template <int dim, int degree>
  template <typename RegularOperation,
            typename DuplicatesOperation,
            typename SkippedOperation>
  void
  Dynamic<dim, degree>::loop(const RegularOperation    &regular_operation,
                             const DuplicatesOperation &duplicates_operation,
                             const SkippedOperation    &skipped_operation) const
  {
    // Warning: @mwichro spent more then a week on optimizing this loop, but it
    // still remains singifficantly slower that the version with precomputed
    // tables. In fact, lookup tables are almost as fast as fast as direct
    // array access.
    //  So if you want to improve its performance, be aware that it is not easy.

    for (std::size_t k = 0; k < (dim > 2 ? n_cells_1d : 1); ++k)
      for (std::size_t j = 0; j < (dim > 1 ? n_cells_1d : 1); ++j)
        for (std::size_t i = 0; i < n_cells_1d; ++i)
          {
            const std::size_t cell_index =
              i + n_cells_1d * (j + n_cells_1d * k);

            // TODO: make it a compile-time lookup table
            std::size_t patch_index_i = i * skip_size;
            std::size_t patch_index_j = j * skip_size;
            std::size_t patch_index_k = k * skip_size;
            std::size_t dof_patch =
              patch_index_i + n_dofs_per_patch_1d * patch_index_j +
              n_dofs_per_patch_1d * n_dofs_per_patch_1d * patch_index_k;

            std::size_t dof_cell = 0;

            do_entry<dim>(regular_operation,
                          duplicates_operation,
                          skipped_operation,
                          cell_index,
                          dof_cell,
                          dof_patch,
                          false);
          }
  }
  //---------------------------------------------------------------------------


  namespace internal
  {
    template <>
    struct CellPatchLookup<1, 1>
    {
      using Index = unsigned int;
      using Pair  = std::pair<Index, Index>;

      static constexpr Index                n_cells        = 2;
      static constexpr Index                n_patch_dofs   = 1;
      static constexpr std::array<Index, 2> n_skipped_dofs = {{1, 1}};

      static constexpr std::array<std::array<Index, 1>, 2> skipped_dofs = {
        {{{0}}, {{1}}}};


      static constexpr std::array<Index, 2> n_cell_to_patch = {{1, 1}};

      static constexpr std::array<std::array<Pair, 1>, 2> cell_to_patch = {
        {{{Pair{1, 0}}}, {{Pair{0, 0}}}}};

      static constexpr std::array<Index, 2> n_non_overlapping = {{1, 0}};
    };


    template <>
    struct CellPatchLookup<1, 2>
    {
      using Index = unsigned int;
      using Pair  = std::pair<Index, Index>;

      static constexpr Index                n_cells        = 2;
      static constexpr Index                n_patch_dofs   = 3;
      static constexpr std::array<Index, 2> n_skipped_dofs = {{1, 1}};

      static constexpr std::array<std::array<Index, 1>, 2> skipped_dofs = {
        {{{0}}, {{2}}}};


      static constexpr std::array<Index, 2> n_cell_to_patch = {{2, 2}};

      static constexpr std::array<std::array<Pair, 2>, 2> cell_to_patch = {
        {{{Pair{1, 0}, Pair{2, 1}}}, {{Pair{1, 2}, Pair{0, 1}}}}};

      static constexpr std::array<Index, 2> n_non_overlapping = {{2, 1}};
    };


    template <>
    struct CellPatchLookup<1, 3>
    {
      using Index = unsigned int;
      using Pair  = std::pair<Index, Index>;

      static constexpr Index                n_cells        = 2;
      static constexpr Index                n_patch_dofs   = 5;
      static constexpr std::array<Index, 2> n_skipped_dofs = {{1, 1}};

      static constexpr std::array<std::array<Index, 1>, 2> skipped_dofs = {
        {{{0}}, {{3}}}};


      static constexpr std::array<Index, 2> n_cell_to_patch = {{3, 3}};

      static constexpr std::array<std::array<Pair, 3>, 2> cell_to_patch = {
        {{{Pair{1, 0}, Pair{2, 1}, Pair{3, 2}}},
         {{Pair{1, 3}, Pair{2, 4}, Pair{0, 2}}}}};

      static constexpr std::array<Index, 2> n_non_overlapping = {{3, 2}};
    };


    template <>
    struct CellPatchLookup<2, 1>
    {
      using Index = unsigned int;
      using Pair  = std::pair<Index, Index>;

      static constexpr Index                n_cells        = 4;
      static constexpr Index                n_patch_dofs   = 1;
      static constexpr std::array<Index, 4> n_skipped_dofs = {{3, 3, 3, 3}};

      static constexpr std::array<std::array<Index, 3>, 4> skipped_dofs = {
        {{{0, 1, 2}}, {{0, 1, 3}}, {{0, 2, 3}}, {{1, 2, 3}}}};


      static constexpr std::array<Index, 4> n_cell_to_patch = {{1, 1, 1, 1}};

      static constexpr std::array<std::array<Pair, 1>, 4> cell_to_patch = {
        {{{Pair{3, 0}}}, {{Pair{2, 0}}}, {{Pair{1, 0}}}, {{Pair{0, 0}}}}};

      static constexpr std::array<Index, 4> n_non_overlapping = {{1, 0, 0, 0}};
    };



    template <>
    struct CellPatchLookup<2, 2>
    {
      using Index                         = unsigned int;
      using Pair                          = std::pair<Index, Index>;
      static constexpr Index n_cells      = 4;
      static constexpr Index n_patch_dofs = 9;

      static constexpr std::array<Index, 4> n_skipped_dofs = {{5, 5, 5, 5}};

      static constexpr std::array<std::array<Index, 5>, 4> skipped_dofs = {
        {{{0, 1, 2, 3, 6}},
         {{0, 1, 2, 5, 8}},
         {{0, 3, 6, 7, 8}},
         {{2, 5, 6, 7, 8}}}};


      static constexpr std::array<Index, 4> n_cell_to_patch = {{4, 4, 4, 4}};

      static constexpr std::array<std::array<Pair, 4>, 4> cell_to_patch = {
        {{{Pair{4, 0}, Pair{5, 1}, Pair{7, 3}, Pair{8, 4}}},
         {{Pair{4, 2}, Pair{7, 5}, Pair{3, 1}, Pair{6, 4}}},
         {{Pair{4, 6}, Pair{5, 7}, Pair{1, 3}, Pair{2, 4}}},
         {{Pair{4, 8}, Pair{0, 4}, Pair{1, 5}, Pair{3, 7}}}}};

      static constexpr std::array<Index, 4> n_non_overlapping = {{4, 2, 2, 1}};
    };



    template <>
    struct CellPatchLookup<2, 3>
    {
      using Index = unsigned int;
      using Pair  = std::pair<Index, Index>;

      static constexpr Index n_cells      = 4;
      static constexpr Index n_patch_dofs = 25;

      static constexpr std::array<Index, 4> n_skipped_dofs = {{7, 7, 7, 7}};

      static constexpr std::array<std::array<Index, 7>, 4> skipped_dofs = {
        {{{0, 1, 2, 3, 4, 8, 12}},
         {{0, 1, 2, 3, 7, 11, 15}},
         {{0, 4, 8, 12, 13, 14, 15}},
         {{3, 7, 11, 12, 13, 14, 15}}}};


      static constexpr std::array<Index, 4> n_cell_to_patch = {{9, 9, 9, 9}};

      static constexpr std::array<std::array<Pair, 9>, 4> cell_to_patch = {
        {{{Pair{5, 0},
           Pair{6, 1},
           Pair{7, 2},
           Pair{9, 5},
           Pair{10, 6},
           Pair{11, 7},
           Pair{13, 10},
           Pair{14, 11},
           Pair{15, 12}}},
         {{Pair{5, 3},
           Pair{6, 4},
           Pair{9, 8},
           Pair{10, 9},
           Pair{13, 13},
           Pair{14, 14},
           Pair{4, 2},
           Pair{8, 7},
           Pair{12, 12}}},
         {{Pair{5, 15},
           Pair{6, 16},
           Pair{7, 17},
           Pair{9, 20},
           Pair{10, 21},
           Pair{11, 22},
           Pair{1, 10},
           Pair{2, 11},
           Pair{3, 12}}},
         {{Pair{5, 18},
           Pair{6, 19},
           Pair{9, 23},
           Pair{10, 24},
           Pair{0, 12},
           Pair{1, 13},
           Pair{2, 14},
           Pair{4, 17},
           Pair{8, 22}}}}};

      static constexpr std::array<Index, 4> n_non_overlapping = {{9, 6, 6, 4}};
    };

    template <>
    struct CellPatchLookup<2, 4>
    {
      using Index                    = unsigned int;
      using Pair                     = std::pair<Index, Index>;
      static constexpr Index n_cells = 4;

      static constexpr Index n_patch_dofs = 49;

      static constexpr std::array<Index, 4> n_skipped_dofs = {{9, 9, 9, 9}};

      static constexpr std::array<std::array<Index, 9>, 4> skipped_dofs = {
        {{{0, 1, 2, 3, 4, 5, 10, 15, 20}},
         {{0, 1, 2, 3, 4, 9, 14, 19, 24}},
         {{0, 5, 10, 15, 20, 21, 22, 23, 24}},
         {{4, 9, 14, 19, 20, 21, 22, 23, 24}}}};


      static constexpr std::array<Index, 4> n_cell_to_patch = {
        {16, 16, 16, 16}};

      static constexpr std::array<std::array<Pair, 16>, 4> cell_to_patch = {
        {{{Pair{6, 0},
           Pair{7, 1},
           Pair{8, 2},
           Pair{9, 3},
           Pair{11, 7},
           Pair{12, 8},
           Pair{13, 9},
           Pair{14, 10},
           Pair{16, 14},
           Pair{17, 15},
           Pair{18, 16},
           Pair{19, 17},
           Pair{21, 21},
           Pair{22, 22},
           Pair{23, 23},
           Pair{24, 24}}},
         {{Pair{6, 4},
           Pair{7, 5},
           Pair{8, 6},
           Pair{11, 11},
           Pair{12, 12},
           Pair{13, 13},
           Pair{16, 18},
           Pair{17, 19},
           Pair{18, 20},
           Pair{21, 25},
           Pair{22, 26},
           Pair{23, 27},
           Pair{5, 3},
           Pair{10, 10},
           Pair{15, 17},
           Pair{20, 24}}},
         {{Pair{6, 28},
           Pair{7, 29},
           Pair{8, 30},
           Pair{9, 31},
           Pair{11, 35},
           Pair{12, 36},
           Pair{13, 37},
           Pair{14, 38},
           Pair{16, 42},
           Pair{17, 43},
           Pair{18, 44},
           Pair{19, 45},
           Pair{1, 21},
           Pair{2, 22},
           Pair{3, 23},
           Pair{4, 24}}},
         {{Pair{6, 32},
           Pair{7, 33},
           Pair{8, 34},
           Pair{11, 39},
           Pair{12, 40},
           Pair{13, 41},
           Pair{16, 46},
           Pair{17, 47},
           Pair{18, 48},
           Pair{0, 24},
           Pair{1, 25},
           Pair{2, 26},
           Pair{3, 27},
           Pair{5, 31},
           Pair{10, 38},
           Pair{15, 45}}}}};

      static constexpr std::array<Index, 4> n_non_overlapping = {
        {16, 12, 12, 9}};
    };



    template <>
    struct CellPatchLookup<3, 1>
    {
      using Index = unsigned int;
      using Pair  = std::pair<Index, Index>;

      static constexpr Index                n_cells        = 8;
      static constexpr Index                n_patch_dofs   = 1;
      static constexpr std::array<Index, 8> n_skipped_dofs = {
        {7, 7, 7, 7, 7, 7, 7, 7}};

      static constexpr std::array<std::array<Index, 7>, 8> skipped_dofs = {
        {{{0, 1, 2, 3, 4, 5, 6}},
         {{0, 1, 2, 3, 4, 5, 7}},
         {{0, 1, 2, 3, 4, 6, 7}},
         {{0, 1, 2, 3, 5, 6, 7}},
         {{0, 1, 2, 4, 5, 6, 7}},
         {{0, 1, 3, 4, 5, 6, 7}},
         {{0, 2, 3, 4, 5, 6, 7}},
         {{1, 2, 3, 4, 5, 6, 7}}}};


      static constexpr std::array<Index, 8> n_cell_to_patch = {
        {1, 1, 1, 1, 1, 1, 1, 1}};

      static constexpr std::array<std::array<Pair, 1>, 8> cell_to_patch = {
        {{{Pair{7, 0}}},
         {{Pair{6, 0}}},
         {{Pair{5, 0}}},
         {{Pair{4, 0}}},
         {{Pair{3, 0}}},
         {{Pair{2, 0}}},
         {{Pair{1, 0}}},
         {{Pair{0, 0}}}}};

      static constexpr std::array<Index, 8> n_non_overlapping = {
        {1, 0, 0, 0, 0, 0, 0, 0}};
    };



    template <>
    struct CellPatchLookup<3, 2>
    {
      using Index = unsigned int;
      using Pair  = std::pair<Index, Index>;

      static constexpr Index                n_cells        = 8;
      static constexpr Index                n_patch_dofs   = 27;
      static constexpr std::array<Index, 8> n_skipped_dofs = {
        {19, 19, 19, 19, 19, 19, 19, 19}};

      static constexpr std::array<std::array<Index, 19>, 8> skipped_dofs = {
        {{{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 15, 18, 19, 20, 21, 24}},
         {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 14, 17, 18, 19, 20, 23, 26}},
         {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 15, 16, 17, 18, 21, 24, 25, 26}},
         {{0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 14, 15, 16, 17, 20, 23, 24, 25, 26}},
         {{0,
           1,
           2,
           3,
           6,
           9,
           10,
           11,
           12,
           15,
           18,
           19,
           20,
           21,
           22,
           23,
           24,
           25,
           26}},
         {{0,
           1,
           2,
           5,
           8,
           9,
           10,
           11,
           14,
           17,
           18,
           19,
           20,
           21,
           22,
           23,
           24,
           25,
           26}},
         {{0,
           3,
           6,
           7,
           8,
           9,
           12,
           15,
           16,
           17,
           18,
           19,
           20,
           21,
           22,
           23,
           24,
           25,
           26}},
         {{2,
           5,
           6,
           7,
           8,
           11,
           14,
           15,
           16,
           17,
           18,
           19,
           20,
           21,
           22,
           23,
           24,
           25,
           26}}}};


      static constexpr std::array<Index, 8> n_cell_to_patch = {
        {8, 8, 8, 8, 8, 8, 8, 8}};

      static constexpr std::array<std::array<Pair, 8>, 8> cell_to_patch = {
        {{{Pair{13, 0},
           Pair{14, 1},
           Pair{16, 3},
           Pair{17, 4},
           Pair{22, 9},
           Pair{23, 10},
           Pair{25, 12},
           Pair{26, 13}}},
         {{Pair{13, 2},
           Pair{16, 5},
           Pair{22, 11},
           Pair{25, 14},
           Pair{12, 1},
           Pair{15, 4},
           Pair{21, 10},
           Pair{24, 13}}},
         {{Pair{13, 6},
           Pair{14, 7},
           Pair{22, 15},
           Pair{23, 16},
           Pair{10, 3},
           Pair{11, 4},
           Pair{19, 12},
           Pair{20, 13}}},
         {{Pair{13, 8},
           Pair{22, 17},
           Pair{9, 4},
           Pair{10, 5},
           Pair{12, 7},
           Pair{18, 13},
           Pair{19, 14},
           Pair{21, 16}}},
         {{Pair{13, 18},
           Pair{14, 19},
           Pair{16, 21},
           Pair{17, 22},
           Pair{4, 9},
           Pair{5, 10},
           Pair{7, 12},
           Pair{8, 13}}},
         {{Pair{13, 20},
           Pair{16, 23},
           Pair{3, 10},
           Pair{4, 11},
           Pair{6, 13},
           Pair{7, 14},
           Pair{12, 19},
           Pair{15, 22}}},
         {{Pair{13, 24},
           Pair{14, 25},
           Pair{1, 12},
           Pair{2, 13},
           Pair{4, 15},
           Pair{5, 16},
           Pair{10, 21},
           Pair{11, 22}}},
         {{Pair{13, 26},
           Pair{0, 13},
           Pair{1, 14},
           Pair{3, 16},
           Pair{4, 17},
           Pair{9, 22},
           Pair{10, 23},
           Pair{12, 25}}}}};

      static constexpr std::array<Index, 8> n_non_overlapping = {
        {8, 4, 4, 2, 4, 2, 2, 1}};
    };



    template <>
    struct CellPatchLookup<3, 3>
    {
      using Index = unsigned int;
      using Pair  = std::pair<Index, Index>;

      static constexpr Index                n_cells        = 8;
      static constexpr Index                n_patch_dofs   = 125;
      static constexpr std::array<Index, 8> n_skipped_dofs = {
        {37, 37, 37, 37, 37, 37, 37, 37}};

      static constexpr std::array<std::array<Index, 37>, 8> skipped_dofs = {
        {{{0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12,
           13, 14, 15, 16, 17, 18, 19, 20, 24, 28, 32, 33, 34,
           35, 36, 40, 44, 48, 49, 50, 51, 52, 56, 60}},
         {{0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12,
           13, 14, 15, 16, 17, 18, 19, 23, 27, 31, 32, 33, 34,
           35, 39, 43, 47, 48, 49, 50, 51, 55, 59, 63}},
         {{0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12,
           13, 14, 15, 16, 20, 24, 28, 29, 30, 31, 32, 36, 40,
           44, 45, 46, 47, 48, 52, 56, 60, 61, 62, 63}},
         {{0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12,
           13, 14, 15, 19, 23, 27, 28, 29, 30, 31, 35, 39, 43,
           44, 45, 46, 47, 51, 55, 59, 60, 61, 62, 63}},
         {{0,  1,  2,  3,  4,  8,  12, 16, 17, 18, 19, 20, 24,
           28, 32, 33, 34, 35, 36, 40, 44, 48, 49, 50, 51, 52,
           53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63}},
         {{0,  1,  2,  3,  7,  11, 15, 16, 17, 18, 19, 23, 27,
           31, 32, 33, 34, 35, 39, 43, 47, 48, 49, 50, 51, 52,
           53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63}},
         {{0,  4,  8,  12, 13, 14, 15, 16, 20, 24, 28, 29, 30,
           31, 32, 36, 40, 44, 45, 46, 47, 48, 49, 50, 51, 52,
           53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63}},
         {{3,  7,  11, 12, 13, 14, 15, 19, 23, 27, 28, 29, 30,
           31, 35, 39, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52,
           53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63}}}};


      static constexpr std::array<Index, 8> n_cell_to_patch = {
        {27, 27, 27, 27, 27, 27, 27, 27}};

      static constexpr std::array<std::array<Pair, 27>, 8> cell_to_patch = {
        {{{Pair{21, 0},  Pair{22, 1},  Pair{23, 2},  Pair{25, 5},  Pair{26, 6},
           Pair{27, 7},  Pair{29, 10}, Pair{30, 11}, Pair{31, 12}, Pair{37, 25},
           Pair{38, 26}, Pair{39, 27}, Pair{41, 30}, Pair{42, 31}, Pair{43, 32},
           Pair{45, 35}, Pair{46, 36}, Pair{47, 37}, Pair{53, 50}, Pair{54, 51},
           Pair{55, 52}, Pair{57, 55}, Pair{58, 56}, Pair{59, 57}, Pair{61, 60},
           Pair{62, 61}, Pair{63, 62}}},
         {{Pair{21, 3},  Pair{22, 4},  Pair{25, 8},  Pair{26, 9},  Pair{29, 13},
           Pair{30, 14}, Pair{37, 28}, Pair{38, 29}, Pair{41, 33}, Pair{42, 34},
           Pair{45, 38}, Pair{46, 39}, Pair{53, 53}, Pair{54, 54}, Pair{57, 58},
           Pair{58, 59}, Pair{61, 63}, Pair{62, 64}, Pair{20, 2},  Pair{24, 7},
           Pair{28, 12}, Pair{36, 27}, Pair{40, 32}, Pair{44, 37}, Pair{52, 52},
           Pair{56, 57}, Pair{60, 62}}},
         {{Pair{21, 15}, Pair{22, 16}, Pair{23, 17}, Pair{25, 20}, Pair{26, 21},
           Pair{27, 22}, Pair{37, 40}, Pair{38, 41}, Pair{39, 42}, Pair{41, 45},
           Pair{42, 46}, Pair{43, 47}, Pair{53, 65}, Pair{54, 66}, Pair{55, 67},
           Pair{57, 70}, Pair{58, 71}, Pair{59, 72}, Pair{17, 10}, Pair{18, 11},
           Pair{19, 12}, Pair{33, 35}, Pair{34, 36}, Pair{35, 37}, Pair{49, 60},
           Pair{50, 61}, Pair{51, 62}}},
         {{Pair{21, 18}, Pair{22, 19}, Pair{25, 23}, Pair{26, 24}, Pair{37, 43},
           Pair{38, 44}, Pair{41, 48}, Pair{42, 49}, Pair{53, 68}, Pair{54, 69},
           Pair{57, 73}, Pair{58, 74}, Pair{16, 12}, Pair{17, 13}, Pair{18, 14},
           Pair{20, 17}, Pair{24, 22}, Pair{32, 37}, Pair{33, 38}, Pair{34, 39},
           Pair{36, 42}, Pair{40, 47}, Pair{48, 62}, Pair{49, 63}, Pair{50, 64},
           Pair{52, 67}, Pair{56, 72}}},
         {{Pair{21, 75},  Pair{22, 76},  Pair{23, 77},  Pair{25, 80},
           Pair{26, 81},  Pair{27, 82},  Pair{29, 85},  Pair{30, 86},
           Pair{31, 87},  Pair{37, 100}, Pair{38, 101}, Pair{39, 102},
           Pair{41, 105}, Pair{42, 106}, Pair{43, 107}, Pair{45, 110},
           Pair{46, 111}, Pair{47, 112}, Pair{5, 50},   Pair{6, 51},
           Pair{7, 52},   Pair{9, 55},   Pair{10, 56},  Pair{11, 57},
           Pair{13, 60},  Pair{14, 61},  Pair{15, 62}}},
         {{Pair{21, 78},  Pair{22, 79},  Pair{25, 83},  Pair{26, 84},
           Pair{29, 88},  Pair{30, 89},  Pair{37, 103}, Pair{38, 104},
           Pair{41, 108}, Pair{42, 109}, Pair{45, 113}, Pair{46, 114},
           Pair{4, 52},   Pair{5, 53},   Pair{6, 54},   Pair{8, 57},
           Pair{9, 58},   Pair{10, 59},  Pair{12, 62},  Pair{13, 63},
           Pair{14, 64},  Pair{20, 77},  Pair{24, 82},  Pair{28, 87},
           Pair{36, 102}, Pair{40, 107}, Pair{44, 112}}},
         {{Pair{21, 90},  Pair{22, 91},  Pair{23, 92},  Pair{25, 95},
           Pair{26, 96},  Pair{27, 97},  Pair{37, 115}, Pair{38, 116},
           Pair{39, 117}, Pair{41, 120}, Pair{42, 121}, Pair{43, 122},
           Pair{1, 60},   Pair{2, 61},   Pair{3, 62},   Pair{5, 65},
           Pair{6, 66},   Pair{7, 67},   Pair{9, 70},   Pair{10, 71},
           Pair{11, 72},  Pair{17, 85},  Pair{18, 86},  Pair{19, 87},
           Pair{33, 110}, Pair{34, 111}, Pair{35, 112}}},
         {{Pair{21, 93},  Pair{22, 94},  Pair{25, 98},  Pair{26, 99},
           Pair{37, 118}, Pair{38, 119}, Pair{41, 123}, Pair{42, 124},
           Pair{0, 62},   Pair{1, 63},   Pair{2, 64},   Pair{4, 67},
           Pair{5, 68},   Pair{6, 69},   Pair{8, 72},   Pair{9, 73},
           Pair{10, 74},  Pair{16, 87},  Pair{17, 88},  Pair{18, 89},
           Pair{20, 92},  Pair{24, 97},  Pair{32, 112}, Pair{33, 113},
           Pair{34, 114}, Pair{36, 117}, Pair{40, 122}}}}};

      static constexpr std::array<Index, 8> n_non_overlapping = {
        {27, 18, 18, 12, 18, 12, 12, 8}};
    };
  } // namespace internal

#endif // DOXYGEN
} // namespace PatchDistributors



DEAL_II_NAMESPACE_CLOSE

#endif // dealii__patch_distributors
