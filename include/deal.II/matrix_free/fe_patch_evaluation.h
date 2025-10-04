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


#ifndef dealii__fe_patch_evaluation_h
#define dealii__fe_patch_evaluation_h


#include <deal.II/base/utilities.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/patch_storage.h>

#include <algorithm>
#include <vector>


DEAL_II_NAMESPACE_OPEN


enum VectorizationType
{
  over_patches,
  within_patch
};


#ifndef DOXYGEN
namespace internal
{
  // Taken from:
  // https://stackoverflow.com/questions/17923683/why-does-stdarray-not-have-an-constructor-that-takes-a-value-for-the-array-to
  template <typename T, std::size_t... Is>
  std::array<T, sizeof...(Is)>
  make_array(const T &value, std::index_sequence<Is...>)
  {
    return {{(static_cast<void>(Is), value)...}};
  }

  template <std::size_t N, typename T>
  std::array<T, N>
  make_array(const T &value)
  {
    return make_array(value, std::make_index_sequence<N>());
  }


  template <template <int, int, int, int, typename, typename> typename T,
            int dim,
            int degree,
            int n_components,
            int n_q_points,
            typename NumberType,
            typename VectorizedNumberType>
  constexpr int
  extract_degree(const T<dim,
                         degree,
                         n_components,
                         n_q_points,
                         NumberType,
                         VectorizedNumberType> *)
  {
    return degree;
  }

  template <unsigned int      NLanes,
            unsigned int      NCells,
            VectorizationType vectorization>
  constexpr unsigned int
  n_evaluators_v()
  {
    static_assert(NLanes > 0 && NCells > 0, "n_lanes and n_cells must be > 0");
    if constexpr (vectorization == within_patch)
      {
        if constexpr (NLanes >= NCells)
          return 1;

        static_assert(
          NCells % NLanes == 0,
          "n_cells must be divisible by n_lanes for within_patch vectorization");
        return NCells / NLanes;
      }
    else
      {
        return NCells;
      }
  }


} // namespace internal

#endif // DOXYGEN

/**
 * @brief A class that provides evaluation capabilities for finite element
 * data on patches of cells.
 *
 * This class builds upon FEEvaluation to handle data associated with patches
 * of cells, which are typically used in algorithms like geometric multigrid or
 * patch smoothers. It manages the distribution and gathering of degrees of
 * freedom (DoFs) between the patch representation and the local cell-based
 * representation used by FEEvaluation.
 *
 * The Distributor must provide the following interface (prototypes only):
 *
 * @code
 * unsigned int n_patch_dofs() const;
 *
 * template <typename RegularOp, typename OverlapOp, typename SkippedOp>
 * void loop(RegularOp regular_op,
 *           OverlapOp overlap_op,
 *           SkippedOp skipped_op) const;
 *
 * // Reinitialize distributor mapping for a new patch.
 * void reinit(...);
 * @endcode
 *
 * Semantics expected by FEPatchEvaluation:
 *  - distribute_patch_to_local:
 *      * regular_op(patch_index, cell, cell_index) writes the patch value to
 *        the corresponding local slot.
 *      * overlap_op(patch_index, cell, cell_index) is used for DoFs shared by
 *        multiple cells; FEPatchEvaluation either writes the value or zeros
 *        the slot depending on copy_duplicates.
 *      * skipped_op(cell, cell_index) zeros out padded/unused local slots.
 *
 *  - gather_local_to_patch:
 *      * regular_op(patch_index, cell, cell_index) reads the local slot into
 *        the patch vector.
 *      * overlap_op(patch_index, cell, cell_index) accumulates contributions
 *        from multiple local slots when sum_overlapping == true.
 *      * skipped slots are ignored.
 *
 * Reinitialization note:
 *  - FEPatchEvaluation::reinit(...) will invoke the distributor's reinit
 *    to update its mapping for the new patch. The exact signature of that
 *    call and the form of any orientation data passed is not finalized yet.
 *    For the time being, all available distributors assume the standard
 *    (canonical) cell orientation and therefore do not require an explicit
 *    orientation argument, and the reinit call is a just a placeholder.
 *
 * Note: See PatchDistributors::Lookup for a concrete example implementation
 * based on precomputed lookup tables.
 *
 * @tparam FEEval      The underlying FEEvaluation type used for cell-local evaluations.
 * @tparam Distributor A class responsible for mapping between patch DoFs and local cell DoFs.
 * @tparam vectorization Specifies the vectorization strategy. Currently only
 *                       within_patch is supported.
 */
template <typename FEEval,
          typename Distributor,
          VectorizationType vectorization = within_patch>
class FEPatchEvaluation
{
public:
  using FEEvaluationType            = FEEval;
  using VectorizedNumber            = typename FEEvaluationType::NumberType;
  using Number                      = typename FEEvaluationType::number_type;
  const static unsigned int dim     = FEEvaluationType::dimension;
  const static unsigned int n_lanes = FEEvaluationType::n_lanes;

  const static unsigned int n_cells = 1 << dim;

  static constexpr int fe_degree =
    internal::extract_degree(static_cast<FEEvaluationType *>(nullptr));
  const static unsigned int n_dofs_per_cell =
    Utilities::fixed_power<dim>(fe_degree + 1);


  const static unsigned int n_evaluators =
    internal::n_evaluators_v<n_lanes, n_cells, vectorization>();

  using PatchStorageType = PatchStorage<MatrixFree<dim, Number>>;
  using CellDofsViewRaw  = StridedArrayView<Number, n_lanes>;

  /**
   * @brief Constructor for FEPatchEvaluation.
   *
   * Initializes the FEPatchEvaluation object with the provided patch storage
   * and a base FEEvaluation object. It sets up the internal data structures,
   * including the array of FEEvaluation instances (currently only one is
   * supported) and the views into the cell DoF data. It also verifies that
   * the MatrixFree object associated with the storage and the FEEvaluation
   * object are the same.
   *
   * @param storage A constant reference to the PatchStorage object that
   *   manages the patch data.
   * @param fe_eval A constant reference to the FEEvaluation object to be used
   *   as a template for the internal evaluators.
   */
  FEPatchEvaluation(const PatchStorageType &storage, const FEEval &fe_eval);

  /**
   * Construct a new FEPatchEvaluation by taking ownership of the internal
   * resources held by @p other. All internal handles, buffers and cached
   * metadata are transferred to the new object rather than deep-copied.
   * After the move, @p other is left in a valid but unspecified state and may
   * be safely destroyed or assigned to.
   *
   * Thread safety: Not safe to call concurrently with other non-const
   * operations
   * on @p other.
   */
  FEPatchEvaluation(FEPatchEvaluation &&other) noexcept;

  /**
   * @brief Reinitializes the FEPatchEvaluation for a specific patch.
   *
   * This function prepares the FEPatchEvaluation object to work with the
   * patch identified by `patch_index`. It reinitializes the underlying
   * FEEvaluation object(s) with the cells belonging to the specified patch
   * and updates the internal state, including the `current_patch_index`.
   *
   * @param patch_index The index of the patch to reinitialize for.
   */
  void
  reinit(const std::size_t patch_index);

  /**
   * @brief Returns the total number of unique degrees of freedom associated
   * with the current patch.
   *
   * This value is determined by the Distributor object.
   *
   * @return The number of patch degrees of freedom.
   */
  constexpr unsigned int
  n_patch_dofs() const;

  /**
   * @brief Returns the index of the patch currently associated with this
   * FEPatchEvaluation object.
   *
   * This index corresponds to the `patch_index` passed to the last call to
   * `reinit()`.
   *
   * @return A constant reference to the current patch index.
   */
  const unsigned int &
  get_current_patch_index() const;

  /**
   * @brief Distributes values from a patch vector to the local cell DoF
   * storage.
   *
   * This function takes a vector representing the DoF values for the entire
   * patch and distributes them into the internal cell-based storage managed
   * by the underlying FEEvaluation object(s). The Distributor handles the
   * mapping.
   *
   * @tparam NumberType The data type of the values in the patch vector.
   * @param patch_vector An ArrayView containing the patch DoF values. Its
   *   size must match `n_patch_dofs()`.
   * @param copy_duplicates A flag indicating how to handle DoFs that
   *   are shared by multiple cells within the patch. If true, the value is
   *   copied to all corresponding local DoF slots. False implies writing
   *   once.
   */
  template <typename NumberType>
  void
  distribute_patch_to_local(const ArrayView<const NumberType> &patch_vector,
                            const bool                         copy_duplicates);

  /**
   * @brief Gathers values from the local cell DoF storage into a patch
   * vector.
   *
   * This function collects the DoF values from the internal cell-based
   * storage and assembles them into a vector representing the DoFs for the
   * entire patch. The Distributor handles the mapping.
   *
   * @tparam NumberType The data type of the values in the patch vector.
   * @param patch_vector An ArrayView where the gathered patch DoF values will
   *   be stored. Its size must match `n_patch_dofs()`.
   * @param sum_overlapping A flag indicating how to handle DoFs
   *   shared by multiple cells. If true, contributions from different cells
   *   to the same patch DoF are summed.
   */
  template <typename NumberType>
  void
  gather_local_to_patch(const ArrayView<NumberType> &patch_vector,
                        const bool                   sum_overlapping) const;

  /**
   * @brief Reads DoF values from a global vector into the local cell storage
   * for the current patch.
   *
   * This function delegates to the `read_dof_values` method of the
   * underlying FEEvaluation object(s), effectively loading the relevant DoF
   * values from a global data structure (like a solution vector) into the
   * local storage used for computations on the current patch. Requires that
   * `reinit` has been called.
   *
   * @tparam VECTOR The type of the source global vector (e.g.,
   *   LinearAlgebra::distributed::Vector).
   * @param src The global vector from which to read DoF values.
   */
  template <typename VECTOR>
  void
  read_dof_values(const VECTOR &src);

  /**
   * @brief Distributes (adds) local DoF values from the cell storage to a
   * global vector.
   *
   * This function delegates to the `distribute_local_to_global` method of the
   * underlying FEEvaluation object(s). It takes the computed local DoF
   * values (e.g., residuals or updates) stored internally and adds them to
   * the corresponding entries in a global data structure. Requires that
   * `reinit` has been called.
   *
   * @tparam VECTOR The type of the destination global vector (e.g.,
   *   LinearAlgebra::distributed::Vector).
   * @param dst The global vector to which the local DoF values will be added.
   */
  template <typename VECTOR>
  void
  distribute_local_to_global(VECTOR &dst) const;

  /**
   * @brief Provides mutable access to a specific DoF value within the local
   * cell storage.
   *
   * Allows direct modification of the value associated with a specific DoF
   * (`dof`) on a specific cell (`cell`) within the patch's local data
   * representation. Assumes `vectorization == within_patch`.
   *
   * @param cell The local index of the cell within the patch (0 to n_cells-1).
   * @param dof The local index of the DoF within the cell.
   * @return A mutable reference to the DoF value.
   */
  auto &
  get_dof_value(const unsigned int &cell, const unsigned int &dof);

  /**
   * @brief Provides constant access to a specific DoF value within the local
   * cell storage.
   *
   * Allows reading the value associated with a specific DoF (`dof`) on a
   * specific cell (`cell`) within the patch's local data representation.
   * Assumes `vectorization == within_patch`.
   *
   * @param cell The local index of the cell within the patch (0 to n_cells-1).
   * @param dof The local index of the DoF within the cell.
   * @return A constant reference to the DoF value.
   */
  const auto &
  get_dof_value(const unsigned int &cell, const unsigned int &dof) const;

  /**
   * @brief Array containing the underlying FEEvaluation objects.
   *
   * Currently, only one evaluator (`n_evaluators == 1`) is supported. This
   * member provides access to the FEEvaluation instance used for cell-local
   * computations.
   */
  std::array<FEEvaluationType, n_evaluators> fe_evaluations;

private:
  std::array<CellDofsViewRaw, n_cells> cell_dofs_view_raw;

  const PatchStorageType &storage;
  unsigned int            current_patch_index;

  Distributor distributor;
};

#ifndef DOXYGEN

namespace internal
{
  // Helper using index_sequence to initialize std::array elements in-place.
  // This version accepts a pointer to the first FEEvaluation in an array of
  // evaluators (fe_evals). It supports the "within_patch" vectorization:
  // if there are more cells than lanes, multiple evaluators are used and
  // the views are filled by iterating lanes of evaluator 0, then lanes of
  // evaluator 1, etc.
  template <typename FEEval,
            typename Number,
            unsigned int NLanes,
            std::size_t... Is>
  std::array<StridedArrayView<Number, NLanes>, sizeof...(Is)>
  create_cell_dofs_view_array_impl(const ArrayView<FEEval> &fe_evals,
                                   unsigned int             NDofsPerCell,
                                   std::index_sequence<Is...>)
  {
    // Number of cells we are creating views for.
    constexpr std::size_t NCells = sizeof...(Is);

    // Number of evaluators required for within_patch vectorization.
    constexpr unsigned int n_evaluators_required =
      n_evaluators_v<NLanes, static_cast<unsigned int>(NCells), within_patch>();

    Assert(fe_evals.size() >= n_evaluators_required,
           ExcMessage("Wrong size of FEEvaluation objects provided!"));

    // Helper to create a single StridedArrayView for cell index `cell_index`.
    auto make_view =
      [&fe_evals, NDofsPerCell](
        std::size_t cell_index) -> StridedArrayView<Number, NLanes> {
      const unsigned int eval_idx =
        static_cast<unsigned int>(cell_index / NLanes);
      const unsigned int lane = static_cast<unsigned int>(cell_index % NLanes);

      // Create view pointing to the lane-th entry of the evaluator eval_idx.
      // We assume fe_evals references at least `n_evaluators_required`
      // elements.
      return StridedArrayView<Number, NLanes>(
        &(fe_evals[eval_idx].begin_dof_values()[0][lane]), NDofsPerCell);
    };

    // Build the std::array by expanding over the index sequence.
    return {{make_view(Is)...}};
  }

  // Top-level factory. Caller should pass an ArrayView referencing the
  // fe_evaluations container (e.g., ArrayView(fe_evaluations)).
  template <typename FEEval,
            unsigned int NCells,
            typename Number,
            unsigned int NLanes>
  std::array<StridedArrayView<Number, NLanes>, NCells>
  create_cell_dofs_view_array(const ArrayView<FEEval> &fe_evals,
                              unsigned int             NDofsPerCell)
  {
    static_assert(NCells > 0, "NCells must be > 0");
    static_assert(NLanes > 0, "NLanes must be > 0");

    return create_cell_dofs_view_array_impl<FEEval, Number, NLanes>(
      fe_evals, NDofsPerCell, std::make_index_sequence<NCells>());
  }

} // namespace internal


template <typename FEEval,
          typename Distributor,
          VectorizationType vectorization>
FEPatchEvaluation<FEEval, Distributor, vectorization>::FEPatchEvaluation(
  const PatchStorageType &storage,
  const FEEval           &fe_eval)
  : fe_evaluations(
      internal::make_array<n_evaluators, FEEvaluationType>(fe_eval))
  , cell_dofs_view_raw(
      internal::
        create_cell_dofs_view_array<FEEvaluationType, n_cells, Number, n_lanes>(
          ArrayView<FEEvaluationType>(fe_evaluations),
          n_dofs_per_cell))
  , storage(storage)
  , current_patch_index(numbers::invalid_unsigned_int)
{
  Assert(storage.get_matrix_free().get() == &fe_eval.get_matrix_free(),
         ExcMessage("MatrixFree objects do not match!"));
}

template <typename FEEval,
          typename Distributor,
          VectorizationType vectorization>
FEPatchEvaluation<FEEval, Distributor, vectorization>::FEPatchEvaluation(
  FEPatchEvaluation &&other) noexcept
  : fe_evaluations(std::move(other.fe_evaluations))
  , cell_dofs_view_raw(
      internal::
        create_cell_dofs_view_array<FEEvaluationType, n_cells, Number, n_lanes>(
          ArrayView<FEEvaluationType>(fe_evaluations),
          n_dofs_per_cell))
  , storage(other.storage)
  , current_patch_index(numbers::invalid_unsigned_int)
  , distributor(std::move(other.distributor))
{
  other.current_patch_index = numbers::invalid_unsigned_int;
}


template <typename FEEval,
          typename Distributor,
          VectorizationType vectorization>
void
FEPatchEvaluation<FEEval, Distributor, vectorization>::reinit(
  const std::size_t patch_index)
{
  // static_assert(n_evaluators == 1, "Only one evaluator is supported for
  // now");
  static_assert(vectorization == within_patch,
                "Only vectorization within_patch is supported");

  if constexpr (n_evaluators == 1 && vectorization == within_patch)
    {
      const auto &patch_cells =
        storage.get_regular_patch(patch_index).get_cells();

      if constexpr (n_lanes > n_cells)
        {
          std::array<typename PatchStorageType::CellIndex, n_lanes>
            padded_cells;
          for (unsigned int i = 0; i < n_cells; ++i)
            padded_cells[i] = patch_cells[i];
          for (unsigned int i = n_cells; i < n_lanes; ++i)
            padded_cells[i] = numbers::invalid_unsigned_int;

          fe_evaluations[0].reinit(padded_cells);
        }
      else
        {
          fe_evaluations[0].reinit(patch_cells);
        }

      current_patch_index = patch_index;

      // Check consistency of cell_dofs_view_raw with fe_evaluations.
      // This is a debug-only check.
      for (unsigned int i = 0; i < n_cells; ++i)
        {
          // Assert that the view points to the correct memory location by
          // checking the address of the first element accessible through the
          // view's operator[].
          Assert(
            &cell_dofs_view_raw[i][0] ==
              &(fe_evaluations[0].begin_dof_values()[0][i]),
            ExcInternalError(
              "cell_dofs_view_raw does not point to the expected memory location after reinit."));
        }
      // The lexicographic view uses the updated cell_dofs_view_raw, so it
      // should reflect the changes.
      return;
    }

  if constexpr (n_evaluators > 1 && vectorization == within_patch)
    {
      {
        const auto &patch_cells =
          storage.get_regular_patch(patch_index).get_cells();

        // Partition patch_cells into chunks of size n_lanes and reinit each
        // evaluator.
        for (unsigned int eval = 0; eval < n_evaluators; ++eval)
          {
            std::array<typename PatchStorageType::CellIndex, n_lanes> chunk;
            const unsigned int offset = eval * n_lanes;
            for (unsigned int lane = 0; lane < n_lanes; ++lane)
              chunk[lane] = patch_cells[offset + lane];

            fe_evaluations[eval].reinit(chunk);
          }

        current_patch_index = patch_index;

        // Debug-only consistency check of cell_dofs_view_raw against
        // fe_evaluations.
        for (unsigned int i = 0; i < n_cells; ++i)
          {
            const unsigned int eval_idx = i / n_lanes;
            const unsigned int lane     = i % n_lanes;
            Assert(
              &cell_dofs_view_raw[i][0] ==
                &(fe_evaluations[eval_idx].begin_dof_values()[0][lane]),
              ExcInternalError(
                "cell_dofs_view_raw does not point to the expected memory location after reinit."));
          }

        return;
      }
    }
}


template <typename FEEval,
          typename Distributor,
          VectorizationType vectorization>
inline constexpr unsigned int
FEPatchEvaluation<FEEval, Distributor, vectorization>::n_patch_dofs() const
{
  return distributor.n_patch_dofs();
}

template <typename FEEval,
          typename Distributor,
          VectorizationType vectorization>
const unsigned int &
FEPatchEvaluation<FEEval, Distributor, vectorization>::get_current_patch_index()
  const
{
  return current_patch_index;
}


template <typename FEEval,
          typename Distributor,
          VectorizationType vectorization>
auto &
FEPatchEvaluation<FEEval, Distributor, vectorization>::get_dof_value(
  const unsigned int &cell,
  const unsigned int &dof)
{
  static_assert(vectorization == within_patch,
                "Only one evaluator is supported for now");

  return cell_dofs_view_raw[cell][dof];
}


template <typename FEEval,
          typename Distributor,
          VectorizationType vectorization>
const auto &
FEPatchEvaluation<FEEval, Distributor, vectorization>::get_dof_value(
  const unsigned int &cell,
  const unsigned int &dof) const
{
  static_assert(vectorization == within_patch,
                "Only one evaluator is supported for now");

  return cell_dofs_view_raw[cell][dof];
}

template <typename FEEval,
          typename Distributor,
          VectorizationType vectorization>
template <typename NumberType>
inline void
FEPatchEvaluation<FEEval, Distributor, vectorization>::
  distribute_patch_to_local(const ArrayView<const NumberType> &patch_vector,
                            const bool                         copy_duplicates)
{
  AssertDimension(patch_vector.size(), distributor.n_patch_dofs());
  Assert(current_patch_index != numbers::invalid_unsigned_int,
         ExcNotInitialized());

  const auto regular_operation = [&](const unsigned int patch_index,
                                     const unsigned int cell,
                                     const unsigned int cell_index) {
    cell_dofs_view_raw[cell][cell_index] = patch_vector[patch_index];
  };

  const auto overlap_operation = [&](const unsigned int patch_index,
                                     const unsigned int cell,
                                     const unsigned int cell_index) {
    if (copy_duplicates)
      cell_dofs_view_raw[cell][cell_index] = patch_vector[patch_index];
    else
      cell_dofs_view_raw[cell][cell_index] = 0;
  };

  const auto skipped_operation = [&](const unsigned int cell,
                                     const unsigned int cell_index) {
    cell_dofs_view_raw[cell][cell_index] = 0;
  };
  distributor.loop(regular_operation, overlap_operation, skipped_operation);
}

template <typename FEEval,
          typename Distributor,
          VectorizationType vectorization>
template <typename NumberType>
inline void
FEPatchEvaluation<FEEval, Distributor, vectorization>::gather_local_to_patch(
  const ArrayView<NumberType> &patch_vector,
  const bool                   sum_overlapping) const
{
  AssertDimension(patch_vector.size(), distributor.n_patch_dofs());
  Assert(current_patch_index != numbers::invalid_unsigned_int,
         ExcNotInitialized());


  const auto regular_operation = [&](const unsigned int patch_index,
                                     const unsigned int cell,
                                     const unsigned int cell_index) {
    patch_vector[patch_index] = cell_dofs_view_raw[cell][cell_index];
  };

  const auto overlap_operation = [&](const unsigned int patch_index,
                                     const unsigned int cell,
                                     const unsigned int cell_index) {
    if (sum_overlapping)
      patch_vector[patch_index] += cell_dofs_view_raw[cell][cell_index];
  };

  distributor.loop(regular_operation, overlap_operation, [](const auto &...) {
  });
}

template <typename FEEval,
          typename Distributor,
          VectorizationType vectorization>
template <typename VECTOR>
void
FEPatchEvaluation<FEEval, Distributor, vectorization>::read_dof_values(
  const VECTOR &src)
{
  Assert(current_patch_index != numbers::invalid_unsigned_int,
         ExcNotInitialized());
  for (auto &fe_eval : fe_evaluations)
    fe_eval.read_dof_values(src);
}

template <typename FEEval,
          typename Distributor,
          VectorizationType vectorization>
template <typename VECTOR>
void
FEPatchEvaluation<FEEval, Distributor, vectorization>::
  distribute_local_to_global(VECTOR &dst) const
{
  Assert(current_patch_index != numbers::invalid_unsigned_int,
         ExcNotInitialized());
  for (auto &fe_eval : fe_evaluations)
    fe_eval.distribute_local_to_global(dst);
}

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif // dealii__fe_patch_evaluation_h
