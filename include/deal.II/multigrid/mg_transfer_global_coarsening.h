// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_mg_transfer_global_coarsening_h
#define dealii_mg_transfer_global_coarsening_h

#include <deal.II/base/memory_space.h>
#include <deal.II/base/mg_level_object.h>
#include <deal.II/base/mpi_remote_point_evaluation.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/matrix_free/constraint_info.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/shape_info.h>

#include <deal.II/multigrid/mg_base.h>
#include <deal.II/multigrid/mg_transfer_matrix_free.h>

#include <deal.II/non_matching/mapping_info.h>



DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
namespace internal
{
  class MGTwoLevelTransferImplementation;
}

namespace RepartitioningPolicyTools
{
  template <int dim, int spacedim>
  class Base;
}

template <int dim,
          typename Number,
          typename MemorySpace = ::dealii::MemorySpace::Host>
class MGTransferMF;
#endif


/**
 * @addtogroup mg
 * @{
 */

namespace mg
{
  /**
   * A structure with boost::signal objects for optional processing in a
   * non-nested multigrid solver.
   *
   * Similarly to mg::Signals, each signal is called twice: once before and
   * once after the action is performed. The two calls only differ in the
   * booleanargument @p before, which is true the first time and false the
   * second.
   *
   */
  struct SignalsNonNested
  {
    /**
     * This signal is triggered before and after the call to the actual
     * evaluation function inside RemotePointEvaluation::evaluate_and_process()
     * during prolongation.
     */
    boost::signals2::signal<void(const bool before)> prolongation_cell_loop;

    /**
     * This signal is triggered before and after the call to the actual
     * evaluation function inside RemotePointEvaluation::process_and_evaluate()
     * during restriction.
     */
    boost::signals2::signal<void(const bool before)> restriction_cell_loop;

    /**
     * This signal is triggered before and after the call to
     * RemotePointEvaluation::evaluate_and_process() used in
     * MGTwoLevelTransferNonNested::prolongate_and_add(). The difference
     * with the @p prolongation_cell_loop signal is that also the
     * communication phase is included.
     */
    boost::signals2::signal<void(const bool before)> prolongation;

    /**
     * This signal is triggered before and after the call to
     * RemotePointEvaluation::process_and_evaluate() used in
     * MGTwoLevelTransferNonNested::restrict_and_add(). Similarly to
     * the @p restriction_cell_loop signal, also the communication phase is
     * included.
     */
    boost::signals2::signal<void(const bool before)> restriction;
  };
} // namespace mg

/**
 * Global coarsening utility functions.
 */
namespace MGTransferGlobalCoarseningTools
{
  /**
   * Common polynomial coarsening sequences.
   *
   * @note These polynomial coarsening sequences up to a degree of 9 are
   *   precompiled in MGTwoLevelTransfer. See also:
   *   MGTwoLevelTransfer::fast_polynomial_transfer_supported()
   */
  enum class PolynomialCoarseningSequenceType
  {
    /**
     * Half polynomial degree by integer division. For example, for degree=7
     * the following sequence would be obtained:: 7 -> 3 -> 1
     */
    bisect,
    /**
     * Decrease the polynomial degree by one. E.g., for degree=7 following
     * sequence would result: 7 -> 6 -> 5 -> 4 -> 3 -> 2 -> 1
     */
    decrease_by_one,
    /**
     * Decrease the polynomial degree to one. E.g., for degree=7 following
     * sequence would result: 7 -> 1
     */
    go_to_one
  };

  /**
   * For a given @p degree and polynomial coarsening sequence @p p_sequence,
   * determine the next coarser degree.
   */
  unsigned int
  create_next_polynomial_coarsening_degree(
    const unsigned int                      degree,
    const PolynomialCoarseningSequenceType &p_sequence);

  /**
   * For a given @p max_degree and polynomial coarsening sequence @p p_sequence,
   * determine the full sequence of polynomial degrees, sorted in ascending
   * order.
   */
  std::vector<unsigned int>
  create_polynomial_coarsening_sequence(
    const unsigned int                      max_degree,
    const PolynomialCoarseningSequenceType &p_sequence);

  /**
   * For a given triangulation @p tria, determine the geometric coarsening
   * sequence by repeated global coarsening of the provided triangulation.
   *
   * @note For convenience, a reference to the input triangulation is stored in
   *   the last entry of the return vector.
   * @note Currently, not implemented for parallel::fullydistributed::Triangulation.
   * @note The type of the returned triangulations is the same as of the input
   *   triangulation.
   */
  template <int dim, int spacedim>
  std::vector<std::shared_ptr<const Triangulation<dim, spacedim>>>
  create_geometric_coarsening_sequence(
    const Triangulation<dim, spacedim> &tria);

  /**
   * Similar to the above function but also taking a @p policy for
   * repartitioning the triangulations on the coarser levels. If
   * @p preserve_fine_triangulation is set, the input triangulation is not
   * altered,
   * else the triangulation is coarsened. If @p repartition_fine_triangulation
   * is set, the triangulation on the finest level is repartitioned as well. If
   * the flags are set to true/false, the input triangulation is simply used as
   * the finest triangulation.
   *
   * @note For convenience, a reference to the input triangulation is stored in
   *   the last entry of the return vector.
   * @note The type of the returned triangulations is
   *   parallel::fullydistributed::Triangulation.
   * @note Currently, only implemented for parallel::distributed::Triangulation.
   */
  template <int dim, int spacedim>
  std::vector<std::shared_ptr<const Triangulation<dim, spacedim>>>
  create_geometric_coarsening_sequence(
    Triangulation<dim, spacedim>                         &tria,
    const RepartitioningPolicyTools::Base<dim, spacedim> &policy,
    const bool preserve_fine_triangulation,
    const bool repartition_fine_triangulation);

  /**
   * Similar to the above function but taking in a constant version of
   * @p tria. As a consequence, it can not be used for coarsening directly,
   * so a temporary copy will be created internally.
   */
  template <int dim, int spacedim>
  std::vector<std::shared_ptr<const Triangulation<dim, spacedim>>>
  create_geometric_coarsening_sequence(
    const Triangulation<dim, spacedim>                   &tria,
    const RepartitioningPolicyTools::Base<dim, spacedim> &policy,
    const bool repartition_fine_triangulation = false);

} // namespace MGTransferGlobalCoarseningTools



/**
 * An abstract base class for transfer operators between two multigrid levels.
 * The implementation of
 * restriction and prolongation between levels is delegated to derived classes,
 * which implement prolongate_and_add_internal() and restrict_and_add_internal()
 * accordingly.
 */
template <typename VectorType>
class MGTwoLevelTransferBase : public EnableObserverPointer
{
public:
  static_assert(
    std::is_same_v<
      VectorType,
      LinearAlgebra::distributed::Vector<typename VectorType::value_type,
                                         MemorySpace::Host>> ||
      std::is_same_v<
        VectorType,
        LinearAlgebra::distributed::Vector<typename VectorType::value_type,
                                           MemorySpace::Default>>,
    "This class is currently only implemented for vectors of "
    "type LinearAlgebra::distributed::Vector.");

  /**
   * The scalar type used by the vector-type template argument.
   */
  using Number = typename VectorType::value_type;

  /**
   * Default constructor.
   */
  MGTwoLevelTransferBase();

  /**
   * Perform prolongation on a solution vector.
   */
  void
  prolongate_and_add(VectorType &dst, const VectorType &src) const;

  /**
   * Perform restriction on a residual vector.
   */
  void
  restrict_and_add(VectorType &dst, const VectorType &src) const;

  /**
   * Perform interpolation of a solution vector from the fine level to the
   * coarse level. This function is different from restriction, where a
   * weighted residual is transferred to a coarser level (transposition of
   * prolongation matrix). In other words, restriction acts on right hand
   * side vectors, whereas interpolation acts on solution vectors.
   */
  virtual void
  interpolate(VectorType &dst, const VectorType &src) const = 0;

  /**
   * Enable inplace vector operations if external and internal vectors
   * are compatible. The returned pair indicates if the operation
   * was successful on the coarse and the fine level.
   */
  virtual std::pair<bool, bool>
  enable_inplace_operations_if_possible(
    const std::shared_ptr<const Utilities::MPI::Partitioner>
      &partitioner_coarse,
    const std::shared_ptr<const Utilities::MPI::Partitioner>
      &partitioner_fine) = 0;

  /**
   * Return the memory consumption of the allocated memory in this class.
   */
  virtual std::size_t
  memory_consumption() const = 0;

protected:
  /**
   * Perform prolongation on vectors with correct ghosting.
   */
  virtual void
  prolongate_and_add_internal(VectorType &dst, const VectorType &src) const = 0;

  /**
   * Perform restriction on vectors with correct ghosting.
   */
  virtual void
  restrict_and_add_internal(VectorType &dst, const VectorType &src) const = 0;

  /**
   * A wrapper around update_ghost_values() optimized in case the
   * present vector has the same parallel layout of one of the external
   * partitioners.
   */
  void
  update_ghost_values(const VectorType &vec) const;

  /**
   * A wrapper around compress() optimized in case the
   * present vector has the same parallel layout of one of the external
   * partitioners.
   */
  void
  compress(VectorType &vec, const VectorOperation::values op) const;

  /**
   * A wrapper around zero_out_ghost_values() optimized in case the
   * present vector has the same parallel layout of one of the external
   * partitioners.
   */
  void
  zero_out_ghost_values(const VectorType &vec) const;

  /**
   * Enable inplace vector operations if external and internal vectors
   * are compatible.
   */
  template <int dim, std::size_t width, typename IndexType>
  std::pair<bool, bool>
  internal_enable_inplace_operations_if_possible(
    const std::shared_ptr<const Utilities::MPI::Partitioner>
      &partitioner_coarse,
    const std::shared_ptr<const Utilities::MPI::Partitioner> &partitioner_fine,
    bool &vec_fine_needs_ghost_update,
    internal::MatrixFreeFunctions::
      ConstraintInfo<dim, VectorizedArray<Number, width>, IndexType>
                              &constraint_info_coarse,
    std::vector<unsigned int> &dof_indices_fine);

  /**
   * Flag if the finite elements on the fine cells are continuous. If yes,
   * the multiplicity of DoF sharing a vertex/line as well as constraints have
   * to be taken into account via weights.
   */
  bool fine_element_is_continuous;

public:
  /**
   * Partitioner needed by the intermediate vector.
   */
  std::shared_ptr<const Utilities::MPI::Partitioner> partitioner_coarse;

  /**
   * Partitioner needed by the intermediate vector.
   */
  std::shared_ptr<const Utilities::MPI::Partitioner> partitioner_fine;

protected:
  /**
   * Internal vector on which the actual prolongation/restriction is performed.
   */
  mutable VectorType vec_coarse;

  /**
   * Internal vector needed for collecting all degrees of freedom of the fine
   * cells. It is only initialized if the fine-level DoF indices touch DoFs
   * other than the locally active ones (which we always assume can be
   * accessed by the given vectors in the prolongate/restrict functions),
   * otherwise it is left at size zero.
   */
  mutable VectorType vec_fine;

  /**
   * Bool indicating whether fine vector has relevant ghost values.
   */
  bool vec_fine_needs_ghost_update;

  /**
   * Embedded partitioner for efficient communication if locally relevant DoFs
   * are a subset of an external Partitioner object.
   */
  std::shared_ptr<const Utilities::MPI::Partitioner>
    partitioner_coarse_embedded;

  /**
   * Embedded partitioner for efficient communication if locally relevant DoFs
   * are a subset of an external Partitioner object.
   */
  std::shared_ptr<const Utilities::MPI::Partitioner> partitioner_fine_embedded;

  /**
   * Buffer for efficient communication if locally relevant DoFs
   * are a subset of an external Partitioner object.
   */
  mutable AlignedVector<Number> buffer_coarse_embedded;

  /**
   * Buffer for efficient communication if locally relevant DoFs
   * are a subset of an external Partitioner object.
   */
  mutable AlignedVector<Number> buffer_fine_embedded;
};



/**
 * Class for transfer between two multigrid levels for p- or global
 * coarsening. It relies on a list of DoF indices associated with the cells on
 * the coarse and fine side of the transfer, and implements a cell-by-cell
 * (matrix-free) interpolation setup with the reference-cell embedding
 * matrices.
 *
 * The implementation of this class is explained in detail in @cite munch2022gc.
 *
 * There are two possible ways to use this class. In the first option, the
 * transfer is built from the underlying DoFHandler and AffineConstraints
 * objects on the coarse and fine side, collecting an explicit copy of all
 * indices on both sides. This works for a relatively wide set of
 * FiniteElement combinations, including p-adaptive schemes using
 * hp::FECollection. The second, more setup-efficient approach is to build the
 * transfer between two multigrid levels for polynomial coarsening
 * (p-coarsening) from two MatrixFree objects that might already exist from
 * other parts of the code. In this case, we require that both objects share
 * the same triangulation (but differ through their DoFHandler descriptions)
 * and are described by the respective DoFHandler/AffineConstraints pair. This
 * second variant is more efficient because no queries to the DoFHandler need
 * to be made, reducing both the setup time and the overall memory
 * consumption. Note that not all options are supported for the second entry
 * point, and we fall back to the first option in such a case.
 */
template <int dim, typename VectorType>
class MGTwoLevelTransfer : public MGTwoLevelTransferBase<VectorType>
{
public:
  static_assert(
    std::is_same_v<
      VectorType,
      LinearAlgebra::distributed::Vector<typename VectorType::value_type,
                                         MemorySpace::Host>> ||
      std::is_same_v<
        VectorType,
        LinearAlgebra::distributed::Vector<typename VectorType::value_type,
                                           MemorySpace::Default>>,
    "This class is currently only implemented for vectors of "
    "type LinearAlgebra::distributed::Vector.");

  /**
   * The scalar type used by the vector-type template argument.
   */
  using Number = typename VectorType::value_type;

  /**
   * A data type representing a vectorized array of the same kind of objects
   * stored in the `VectorType`.
   */
  using VectorizedArrayType = VectorizedArray<Number>;

  /**
   * Set up global coarsening between the given DoFHandler objects (
   * @p dof_handler_fine and @p dof_handler_coarse). The transfer
   * can be only performed on active levels.
   */
  void
  reinit_geometric_transfer(
    const DoFHandler<dim>           &dof_handler_fine,
    const DoFHandler<dim>           &dof_handler_coarse,
    const AffineConstraints<Number> &constraint_fine =
      AffineConstraints<Number>(),
    const AffineConstraints<Number> &constraint_coarse =
      AffineConstraints<Number>(),
    const unsigned int mg_level_fine   = numbers::invalid_unsigned_int,
    const unsigned int mg_level_coarse = numbers::invalid_unsigned_int);

  /**
   * Set up polynomial coarsening between the given DoFHandler objects (
   * @p dof_handler_fine and @p dof_handler_coarse). Polynomial transfers
   * can be only performed on active levels (`numbers::invalid_unsigned_int`)
   * or on coarse-grid levels, i.e., levels without hanging nodes.
   *
   * @note The function polynomial_transfer_supported() can be used to
   *   check if the given polynomial coarsening strategy is supported.
   */
  void
  reinit_polynomial_transfer(
    const DoFHandler<dim>           &dof_handler_fine,
    const DoFHandler<dim>           &dof_handler_coarse,
    const AffineConstraints<Number> &constraint_fine =
      AffineConstraints<Number>(),
    const AffineConstraints<Number> &constraint_coarse =
      AffineConstraints<Number>(),
    const unsigned int mg_level_fine   = numbers::invalid_unsigned_int,
    const unsigned int mg_level_coarse = numbers::invalid_unsigned_int);

  /**
   * Set up the transfer operator between the given DoFHandler objects (
   * @p dof_handler_fine and @p dof_handler_coarse). Depending on the
   * underlying Triangulation objects polynomial or geometrical global
   * coarsening is performed.
   *
   * @note While geometric transfer can be only performed on active levels
   *   (`numbers::invalid_unsigned_int`), polynomial transfers can also be
   *   performed on coarse-grid levels, i.e., levels without hanging nodes.
   *
   * @note The function polynomial_transfer_supported() can be used to
   *   check if the given polynomial coarsening strategy is supported.
   */
  void
  reinit(const DoFHandler<dim>           &dof_handler_fine,
         const DoFHandler<dim>           &dof_handler_coarse,
         const AffineConstraints<Number> &constraint_fine =
           AffineConstraints<Number>(),
         const AffineConstraints<Number> &constraint_coarse =
           AffineConstraints<Number>(),
         const unsigned int mg_level_fine   = numbers::invalid_unsigned_int,
         const unsigned int mg_level_coarse = numbers::invalid_unsigned_int);

  /**
   * Set up polynomial coarsening between the DoFHandler objects underlying
   * two MatrixFree objects and the respective numbers for the DoFHandler
   * objects within MatrixFree. This reinit() function allows for a more
   * efficient setup of the transfer operator and reduces the overall memory
   * consumption of a multigrid cycle in case the same MatrixFree objects are
   * also used for smoothers and residual evaluation on the two involved
   * levels.
   */
  void
  reinit(const MatrixFree<dim, Number> &matrix_free_fine,
         const unsigned int             dof_no_fine,
         const MatrixFree<dim, Number> &matrix_free_coarse,
         const unsigned int             dof_no_coarse);

  /**
   * Check if a fast templated version of the polynomial transfer between
   * @p fe_degree_fine and @p fe_degree_coarse is available.
   *
   * @note Currently, the polynomial coarsening strategies: 1) go-to-one,
   *   2) bisect, and 3) decrease-by-one are precompiled with templates for
   *   degrees up to 9.
   */
  static bool
  fast_polynomial_transfer_supported(const unsigned int fe_degree_fine,
                                     const unsigned int fe_degree_coarse);

  /**
   * @copydoc MGTwoLevelTransferBase::interpolate
   */
  void
  interpolate(VectorType &dst, const VectorType &src) const override;

  /**
   * Enable inplace vector operations if external and internal vectors
   * are compatible.
   */
  std::pair<bool, bool>
  enable_inplace_operations_if_possible(
    const std::shared_ptr<const Utilities::MPI::Partitioner>
      &partitioner_coarse,
    const std::shared_ptr<const Utilities::MPI::Partitioner> &partitioner_fine)
    override;

  /**
   * Return the memory consumption of the allocated memory in this class.
   */
  std::size_t
  memory_consumption() const override;

protected:
  void
  prolongate_and_add_internal(VectorType       &dst,
                              const VectorType &src) const override;

  void
  restrict_and_add_internal(VectorType       &dst,
                            const VectorType &src) const override;

private:
  /**
   * A multigrid transfer scheme. A multrigrid transfer class can have different
   * transfer schemes to enable p-adaptivity (one transfer scheme per
   * polynomial degree pair) and to enable global coarsening (one transfer
   * scheme for transfer between children and parent cells, as well as, one
   * transfer scheme for cells that are not refined).
   */
  struct MGTransferScheme
  {
    /**
     * Number of coarse cells.
     */
    unsigned int n_coarse_cells;

    /**
     * Number of degrees of freedom of a coarse cell.
     *
     * @note For tensor-product elements, the value equals
     *   `n_components * (degree_coarse + 1)^dim`.
     */
    unsigned int n_dofs_per_cell_coarse;

    /**
     * Number of degrees of freedom of fine cell.
     *
     * @note For tensor-product elements, the value equals
     *   `n_components * (n_dofs_per_cell_fine + 1)^dim`.
     */
    unsigned int n_dofs_per_cell_fine;

    /**
     * Polynomial degree of the finite element of a coarse cell.
     */
    unsigned int degree_coarse;

    /**
     * "Polynomial degree" of the finite element of the union of all children
     * of a coarse cell, i.e., actually `degree_fine * 2 + 1` if a cell is
     * refined.
     */
    unsigned int degree_fine;

    /**
     * Prolongation matrix used for the prolongate_and_add() and
     * restrict_and_add() functions.
     */
    AlignedVector<double> prolongation_matrix;

    /**
     * Restriction matrix used for the interpolate() function.
     */
    AlignedVector<double> restriction_matrix;

    /**
     * ShapeInfo description of the coarse cell. Needed during the
     * fast application of hanging-node constraints.
     */
    internal::MatrixFreeFunctions::ShapeInfo<double> shape_info_coarse;
  };

  /**
   * Transfer schemes.
   */
  std::vector<MGTransferScheme> schemes;

  /**
   * Helper class for reading from and writing to global coarse vectors and for
   * applying constraints.
   */
  internal::MatrixFreeFunctions::
    ConstraintInfo<dim, VectorizedArrayType, types::global_dof_index>
      constraint_info_coarse;

  /**
   * Helper class for reading from and writing to global fine vectors.
   */
  internal::MatrixFreeFunctions::
    ConstraintInfo<dim, VectorizedArrayType, types::global_dof_index>
      constraint_info_fine;

  struct MatrixFreeRelatedData
  {
    /**
     * Matrix-free object on the fine side.
     */
    ObserverPointer<const MatrixFree<dim, Number>> matrix_free_fine;

    /**
     * Index within the list of DoFHandler objects in the matrix_free_fine
     * object.
     */
    unsigned int dof_handler_index_fine;

    /**
     * Matrix-free object on the coarse side.
     */
    ObserverPointer<const MatrixFree<dim, Number>> matrix_free_coarse;

    /**
     * Index within the list of DoFHandler objects in the matrix_free_coarse
     * object.
     */
    unsigned int dof_handler_index_coarse;

    /**
     * The two matrix-free objects will in general not agree on the order the
     * cells are traversed. Thus, the loop will be run by the matrix-free object
     * on the fine side, and the coarse side will adapt to those cell indices.
     */
    std::vector<std::array<unsigned int, VectorizedArrayType::size()>>
      cell_list_fine_to_coarse;
  };

  /**
   * In case this class is built with MatrixFree objects (see the respective
   * reinit() function), we set up this data structure and skip the other
   * fields of the class.
   */
  std::unique_ptr<MatrixFreeRelatedData> matrix_free_data;

  /**
   * CRS-like pointer to the start into the weights array, as that array can
   * be compressed or in full format.
   */
  std::vector<unsigned int> weights_start;

  /**
   * Weights for continuous elements, either in full format or compressed into
   * 3^dim doubles per cell if possible.
   */
  AlignedVector<VectorizedArrayType> weights;

  /**
   * Store whether the weights are in compressed format or not, in the
   * ordering of the weights_start array.
   */
  std::vector<unsigned char> weights_are_compressed;

  /**
   * Number of components.
   */
  unsigned int n_components;

  /**
   * Pointer to the DoFHandler object used during initialization.
   */
  ObserverPointer<const DoFHandler<dim>> dof_handler_fine;

  /**
   * Multigrid level used during initialization.
   */
  unsigned int mg_level_fine;

  friend class internal::MGTwoLevelTransferImplementation;

  friend class MGTransferMF<dim, Number, typename VectorType::memory_space>;
};



/**
 * Class for transfer between two non-nested multigrid levels.
 */
template <int dim, typename VectorType>
class MGTwoLevelTransferNonNested : public MGTwoLevelTransferBase<VectorType>
{
private:
  static_assert(
    std::is_same_v<
      VectorType,
      LinearAlgebra::distributed::Vector<typename VectorType::value_type>>,
    "This class is currently only implemented for vectors of "
    "type LinearAlgebra::distributed::Vector.");

  using Number              = typename VectorType::value_type;
  using VectorizedArrayType = VectorizedArray<Number, 1>;

  mg::SignalsNonNested signals_non_nested;

public:
  /**
   * AdditionalData structure that can be used to tweak parameters
   * related to the search procedure (used internally by RemotePointEvaluation)
   * or, in the future, transfer operators needed by the non-nested multigrid
   * algorithm.
   */
  struct AdditionalData
  {
    /**
     * Constructor. By default, the @p tolerance and @p rtree_level parameters
     * are set to the default values used in the constructor of
     * RemotePointEvaluation, i.e. 1e-6 and 0, respectively. The last Boolean
     * parameter @p enforce_all_points_found is true by default and checks
     * that all points submitted internally to RemotePointEvaluation::reinit()
     * have been found.
     *
     */
    AdditionalData(const double       tolerance                = 1e-6,
                   const unsigned int rtree_level              = 0,
                   const bool         enforce_all_points_found = true)
      : tolerance(tolerance)
      , rtree_level(rtree_level)
      , enforce_all_points_found(enforce_all_points_found)
    {}

    /**
     * Tolerance parameter. See the constructor of RemotePointEvaluation for
     * more details.
     */
    double tolerance;

    /**
     * RTree level parameter. See the constructor of RemotePointEvaluation for
     * more details.
     *
     */
    unsigned int rtree_level;

    /**
     * If set to true, it checks if RemotePointEvaluation::all_points_found()
     * evaluates to true internally during the each call to reinit() from one
     * level to the next one, ensuring that all submitted points have been found
     * inside the domain.
     *
     */
    bool enforce_all_points_found;
  };

  /**
   * Constructor.
   */
  MGTwoLevelTransferNonNested(const AdditionalData &data = AdditionalData());

  /**
   * Set up transfer operator between the given DoFHandler objects (
   * @p dof_handler_fine and @p dof_handler_coarse).
   */
  void
  reinit(const DoFHandler<dim>           &dof_handler_fine,
         const DoFHandler<dim>           &dof_handler_coarse,
         const Mapping<dim>              &mapping_fine,
         const Mapping<dim>              &mapping_coarse,
         const AffineConstraints<Number> &constraint_fine =
           AffineConstraints<Number>(),
         const AffineConstraints<Number> &constraint_coarse =
           AffineConstraints<Number>());

  /**
   * Perform interpolation of a solution vector from the fine level to the
   * coarse level. This function is different from restriction, where a
   * weighted residual is transferred to a coarser level (transposition of
   * prolongation matrix).
   */
  void
  interpolate(VectorType &dst, const VectorType &src) const override;

  /**
   * Enable inplace vector operations if external and internal vectors
   * are compatible.
   */
  std::pair<bool, bool>
  enable_inplace_operations_if_possible(
    const std::shared_ptr<const Utilities::MPI::Partitioner>
      &partitioner_coarse,
    const std::shared_ptr<const Utilities::MPI::Partitioner> &partitioner_fine)
    override;

  /**
   * Return the memory consumption of the allocated memory in this class.
   */
  std::size_t
  memory_consumption() const override;

  /**
   * Connect a function to mg::SignalsNonNested::prolongation_cell_loop.
   */
  boost::signals2::connection
  connect_prolongation_cell_loop(const std::function<void(const bool)> &slot);

  /**
   * Connect a function to mg::SignalsNonNested::restriction_cell_loop.
   */
  boost::signals2::connection
  connect_restriction_cell_loop(const std::function<void(const bool)> &slot);

  /**
   * Connect a function to mg::SignalsNonNested::prolongation.
   */
  boost::signals2::connection
  connect_prolongation(const std::function<void(const bool)> &slot);

  /**
   * Connect a function to mg::SignalsNonNested::restriction.
   */
  boost::signals2::connection
  connect_restriction(const std::function<void(const bool)> &slot);

protected:
  AdditionalData additional_data;
  /**
   * Perform prolongation.
   */
  void
  prolongate_and_add_internal(VectorType       &dst,
                              const VectorType &src) const override;

  /**
   * Perform restriction.
   */
  void
  restrict_and_add_internal(VectorType       &dst,
                            const VectorType &src) const override;

private:
  /**
   * Perform prolongation for correct number of components.
   */
  template <int n_components>
  void
  prolongate_and_add_internal_comp(VectorType       &dst,
                                   const VectorType &src) const;

  /**
   * Perform restriction for correct number of components.
   */
  template <int n_components>
  void
  restrict_and_add_internal_comp(VectorType &dst, const VectorType &src) const;

  /**
   * Pointer to the DoFHandler object used during initialization.
   */
  ObserverPointer<const DoFHandler<dim>> dof_handler_fine;

  /**
   * Multigrid level used during initialization.
   */
  unsigned int mg_level_fine;

  /**
   * Object to evaluate shape functions on one mesh on visited support points of
   * the other mesh.
   */
  Utilities::MPI::RemotePointEvaluation<dim> rpe;

  /**
   * Vectors for input/output for rpe.
   */
  mutable std::vector<Number> rpe_input_output;

  /**
   * Buffers to be reused by rpe.
   */
  mutable std::vector<Number> rpe_buffer;

  /**
   * MappingInfo object needed as Mapping argument by FEPointEvaluation.
   */
  std::shared_ptr<NonMatching::MappingInfo<dim, dim, Number>> mapping_info;

  /**
   * Helper class for reading from and writing to global vectors and for
   * applying constraints.
   */
  internal::MatrixFreeFunctions::
    ConstraintInfo<dim, VectorizedArrayType, unsigned int>
      constraint_info;

  /**
   * Finite element of the coarse DoFHandler passed to reinit().
   */
  std::unique_ptr<FiniteElement<dim>> fe_coarse;

  /**
   * DoF indices of the fine cells, expressed in indices local to the MPI
   * rank.
   */
  std::vector<unsigned int> level_dof_indices_fine;

  /**
   * CRS like structure which points to DoFs associated with the same support
   * point. The vector stays empty if only one DoF corresponds to one support
   * point.
   */
  std::vector<unsigned int> level_dof_indices_fine_ptrs;

  friend class MGTransferMF<dim, Number, typename VectorType::memory_space>;
};



/**
 * Implementation of the MGTransferBase. In contrast to
 * other multigrid transfer operators, the user can provide separate
 * transfer operators of type MGTwoLevelTransfer between each level. The
 * sequence of functions calls for setup is:
 * @code
 * MGTransferGlobalCoarsening mg_transfer;
 * mg_transfer.initialize_two_level_transfers(two_level_transfers);
 * mg_transfer.build(partitioners);
 * @endcode
 *
 * Alternatively, this class can also be set up as in the case of
 * MGTransferMatrixFree:
 * @code
 * MGTransferGlobalCoarsening mg_transfer;
 * mg_transfer.initialize_constraints(mg_constrained_dofs);
 * mg_transfer.build(dof_handler, partitioners);
 * @endcode
 * However, this is way to set up is currently only working for globally
 * refined meshes.
 *
 * This class currently only works for the tensor-product finite elements
 * FE_Q, FE_DGQ, and FE_DGP and simplex elements FE_SimplexP and FE_SimplexDGP
 * as well as for systems involving multiple components of one of these
 * elements. Other elements are currently not implemented.
 *
 * The implementation of this class is explained in detail in @cite munch2022gc.
 */
template <int dim, typename Number, typename MemorySpace>
class MGTransferMF : public dealii::MGLevelGlobalTransfer<
                       LinearAlgebra::distributed::Vector<Number, MemorySpace>>
{
public:
  /**
   * Value type.
   */
  using VectorType = LinearAlgebra::distributed::Vector<Number, MemorySpace>;

  /**
   * Default constructor.
   *
   * @note See also MGTransferMatrixFree.
   */
  MGTransferMF();

  /**
   * @name Global coarsening.
   */
  /** @{ */

  /**
   * Constructor taking a collection of transfer operators (with the coarsest
   * level kept empty in @p transfer) and an optional function that initializes the
   * internal level vectors within the function call copy_to_mg() if used in the
   * context of PreconditionMG. The template parameter @p MGTwoLevelTransferObject should derive from
   * MGTwoLevelTransferBase and implement the transfer operation (see for
   * instance MGTwoLevelTransfer). It can also be a std::shared_ptr or
   * std::unique_ptr to the actual transfer operator.
   */
  template <typename MGTwoLevelTransferObject>
  MGTransferMF(const MGLevelObject<MGTwoLevelTransferObject> &transfer,
               const std::function<void(const unsigned int, VectorType &)>
                 &initialize_dof_vector = {});

  /**
   * Set two-level transfers.
   */
  template <typename MGTwoLevelTransferObject>
  void
  initialize_two_level_transfers(
    const MGLevelObject<MGTwoLevelTransferObject> &transfer);

  /**
   * Similar function to MGTransferMatrixFree::build() with the difference that
   * the information for the prolongation for each level has already been built.
   * So this function only tries to optimize the data structures of the
   * two-level transfer operators, e.g., by enabling inplace vector operations,
   * by checking if @p external_partitioners and the internal ones are
   * compatible.
   */
  void
  build(const std::vector<std::shared_ptr<const Utilities::MPI::Partitioner>>
          &external_partitioners = {});

  /**
   * Same as above but taking a lambda for initializing vector instead of
   * partitioners.
   */
  void
  build(const std::function<void(const unsigned int, VectorType &)>
          &initialize_dof_vector);

  /** @} */

  /**
   * @name Local smoothing.
   */
  /** @{ */

  /**
   * Constructor with constraints. Equivalent to the default constructor
   * followed by initialize_constraints().
   *
   * @note See also MGTransferMatrixFree.
   */
  MGTransferMF(const MGConstrainedDoFs &mg_constrained_dofs);

  /**
   * Initialize the constraints to be used in build().
   *
   * @note See also MGTransferMatrixFree.
   */
  void
  initialize_constraints(const MGConstrainedDoFs &mg_constrained_dofs);

  /** @} */

  /**
   * @name Global coarsening and local smoothing.
   */
  /** @{ */

  /**
   * Actually build the information for the prolongation for each level.
   *
   * @note In the case of global coarsening, you can pass a @p dof_handler
   * with different DoF numbering as the one used within the provided
   * two-level transfer objects  into this function. In this case, vector
   * entries are permuted during copy_to_mg(), copy_from_mg(), and
   * interpolate_to_mg().
   */
  void
  build(const DoFHandler<dim> &dof_handler,
        const std::vector<std::shared_ptr<const Utilities::MPI::Partitioner>>
          &external_partitioners = {});

  /**
   * Same as above but taking a lambda for initializing vector instead of
   * partitioners.
   *
   * @note In the case of global coarsening, you can pass into this function
   * a @p dof_handler with different DoF numbering as the one used within the
   * provided two-level transfer objects. In this case, vector entries are
   * permuted during copy_to_mg(), copy_from_mg(), and interpolate_to_mg().
   */
  void
  build(const DoFHandler<dim> &dof_handler,
        const std::function<void(const unsigned int, VectorType &)>
          &initialize_dof_vector);

  /** @} */

  /**
   * @name Transfer functions.
   */
  /** @{ */

  /**
   * Perform prolongation.
   */
  void
  prolongate(const unsigned int to_level,
             VectorType        &dst,
             const VectorType  &src) const override;

  /**
   * Perform prolongation.
   */
  void
  prolongate_and_add(const unsigned int to_level,
                     VectorType        &dst,
                     const VectorType  &src) const override;

  /**
   * Perform restriction.
   */
  virtual void
  restrict_and_add(const unsigned int from_level,
                   VectorType        &dst,
                   const VectorType  &src) const override;

  /**
   * Initialize internal vectors and copy @p src vector
   * (associated to @p dof_handler) to the finest multigrid level.
   *
   * @note The @p dof_handler object needs to be the same as the DoFHandler passed
   * directly to the function build() or indirectly to the function
   * initialize_two_level_transfers(). Alternatively, the numbering of the DoFs
   * need to be same.
   */
  template <class InVector>
  void
  copy_to_mg(const DoFHandler<dim>     &dof_handler,
             MGLevelObject<VectorType> &dst,
             const InVector            &src) const;

  /**
   * Copy the values on the finest multigrid level to @p dst
   * vector (associated to @p dof_handler).
   *
   * @note The @p dof_handler object needs to be the same as the DoFHandler passed
   * directly to the function build() or indirectly to the function
   * initialize_two_level_transfers(). Alternatively, the numbering of the DoFs
   * need to be same.
   */
  template <class OutVector>
  void
  copy_from_mg(const DoFHandler<dim>           &dof_handler,
               OutVector                       &dst,
               const MGLevelObject<VectorType> &src) const;

  /**
   * Interpolate fine-mesh field @p src (associated to @p dof_handler)
   * to each multigrid level and
   * store the result in @p dst. This function is different from
   * restriction, where a weighted residual is transferred to a coarser
   * level (transposition of prolongation matrix).
   *
   * The argument @p dst has to be initialized with the correct size according
   * to the number of levels of the triangulation.
   *
   * If an inner vector of @p dst is empty or has incorrect locally owned size,
   * it will be resized to locally relevant degrees of freedom on each level.
   *
   * @note The @p dof_handler object needs to be the same as the DoFHandler passed
   * directly to the function build() or indirectly to the function
   * initialize_two_level_transfers(). Alternatively, the numbering of the DoFs
   * need to be same.
   */
  template <class InVector>
  void
  interpolate_to_mg(const DoFHandler<dim>     &dof_handler,
                    MGLevelObject<VectorType> &dst,
                    const InVector            &src) const;

  /**
   * Interpolate fine-mesh field @p src to each multigrid level and
   * store the result in @p dst.
   *
   * @note In contrast to the last function, no @p dof_handler object needs passed.
   * This function used the DoFHandler passed directly to the function build()
   * or indirectly to the function initialize_two_level_transfers().
   */
  template <class InVector>
  void
  interpolate_to_mg(MGLevelObject<VectorType> &dst, const InVector &src) const;

  /** @} */

  /**
   * @name Utility functions.
   */
  /** @{ */

  /**
   * Return the memory consumption of the allocated memory in this class.
   *
   * @note Counts also the memory consumption of the underlying two-level
   *   transfer operators.
   */
  std::size_t
  memory_consumption() const;

  /**
   * Minimum level.
   */
  unsigned int
  min_level() const;

  /**
   * Maximum level.
   */
  unsigned int
  max_level() const;

  /**
   * Clear all data fields and brings the class into a condition similar
   * to after having called the default constructor.
   */
  void
  clear();

  /** @} */

private:
  /**
   * Initial internal transfer operator.
   *
   * @note See also MGTransferMatrixFree.
   */
  void
  initialize_internal_transfer(
    const DoFHandler<dim>                          &dof_handler,
    const ObserverPointer<const MGConstrainedDoFs> &mg_constrained_dofs);

  /**
   * Retrieve finest DoFHandler from two-level transfer objects.
   */
  std::pair<const DoFHandler<dim> *, unsigned int>
  get_dof_handler_fine() const;

  /**
   * Initialize copy indices for MGTransferMF::copy_to_mg(),
   * MGTransferMF::copy_to_mg(), and MGTransferMF::interpolate_to_mg()
   * in the case of global coarsening.
   */
  void
  fill_and_communicate_copy_indices_global_coarsening(
    const DoFHandler<dim> &dof_handler);

  /**
   * Set references to two-level transfer operators to be used.
   */
  template <typename MGTwoLevelTransferObject>
  void
  initialize_transfer_references(
    const MGLevelObject<MGTwoLevelTransferObject> &transfer);

  /**
   * Function to initialize internal level vectors.
   */
  template <class InVector>
  void
  initialize_dof_vector(const unsigned int level,
                        VectorType        &vector,
                        const InVector    &vector_reference,
                        const bool         omit_zeroing_entries) const;

  /**
   * Check that the internal DoFHandler is compatible with the external one
   * provided by copy_to_mg(), copy_from_mg() and interpolate_to_mg()
   * used, e.g., by PreconditionMG.
   */
  void
  assert_dof_handler(const DoFHandler<dim> &dof_handler_out) const;

  /**
   * Internal transfer operator.
   *
   * @note See also MGTransferMatrixFree.
   */
  MGLevelObject<MGTwoLevelTransfer<dim, VectorType>> internal_transfer;

  /**
   * Collection of the two-level transfer operators.
   */
  MGLevelObject<ObserverPointer<MGTwoLevelTransferBase<VectorType>>> transfer;

  /**
   * External partitioners used during initialize_dof_vector().
   */
  std::vector<std::shared_ptr<const Utilities::MPI::Partitioner>>
    external_partitioners;
};



/**
 * This class works with LinearAlgebra::distributed::BlockVector and
 * performs exactly the same transfer operations for each block as
 * MGTransferMF.
 */
template <int dim, typename Number>
class MGTransferBlockMF
  : public MGTransferBlockMatrixFreeBase<
      dim,
      Number,
      MGTransferMF<dim, Number, ::dealii::MemorySpace::Host>>
{
public:
  /**
   * Constructor.
   */
  MGTransferBlockMF(const MGTransferMF<dim, Number, ::dealii::MemorySpace::Host>
                      &transfer_operator);

  /**
   * Constructor.
   *
   * @note See also MGTransferBlockMatrixFree.
   */
  MGTransferBlockMF() = default;

  /**
   * Constructor.
   *
   * @note See also MGTransferBlockMatrixFree.
   */
  MGTransferBlockMF(const MGConstrainedDoFs &mg_constrained_dofs);

  /**
   * Constructor.
   *
   * @note See also MGTransferBlockMatrixFree.
   */
  MGTransferBlockMF(const std::vector<MGConstrainedDoFs> &mg_constrained_dofs);

  /**
   * Initialize the constraints to be used in build().
   *
   * @note See also MGTransferBlockMatrixFree.
   */
  void
  initialize_constraints(const MGConstrainedDoFs &mg_constrained_dofs);

  /**
   * Same as above for the case that each block has its own DoFHandler.
   *
   * @note See also MGTransferBlockMatrixFree.
   */
  void
  initialize_constraints(
    const std::vector<MGConstrainedDoFs> &mg_constrained_dofs);

  /**
   * Actually build the information for the prolongation for each level.
   *
   * @note See also MGTransferBlockMatrixFree.
   */
  void
  build(const DoFHandler<dim> &dof_handler);

  /**
   * Same as above for the case that each block has its own DoFHandler.
   *
   * @note See also MGTransferBlockMatrixFree.
   */
  void
  build(const std::vector<const DoFHandler<dim> *> &dof_handler);

protected:
  const MGTransferMF<dim, Number, ::dealii::MemorySpace::Host> &
  get_matrix_free_transfer(const unsigned int b) const override;

private:
  /**
   * Internal non-block version of transfer operation.
   */
  std::vector<MGTransferMF<dim, Number, ::dealii::MemorySpace::Host>>
    transfer_operators_internal;

  /**
   * Non-block version of transfer operation.
   */
  std::vector<ObserverPointer<
    const MGTransferMF<dim, Number, ::dealii::MemorySpace::Host>>>
    transfer_operators;
};



template <int dim, typename VectorType>
using MGTransferGlobalCoarsening = MGTransferMF<dim,
                                                typename VectorType::value_type,
                                                ::dealii::MemorySpace::Host>;

template <int dim, typename VectorType>
using MGTransferBlockGlobalCoarsening =
  MGTransferBlockMF<dim, typename VectorType::value_type>;


/** @} */

#ifndef DOXYGEN

/* ----------------------- Inline functions --------------------------------- */



template <int dim, typename Number, typename MemorySpace>
template <typename MGTwoLevelTransferObject>
MGTransferMF<dim, Number, MemorySpace>::MGTransferMF(
  const MGLevelObject<MGTwoLevelTransferObject> &transfer,
  const std::function<void(const unsigned int, VectorType &)>
    &initialize_dof_vector)
{
  this->transfer.clear();
  this->internal_transfer.clear();

  this->initialize_transfer_references(transfer);
  this->build(initialize_dof_vector);
}



template <int dim, typename Number, typename MemorySpace>
template <typename MGTwoLevelTransferObject>
void
MGTransferMF<dim, Number, MemorySpace>::initialize_two_level_transfers(
  const MGLevelObject<MGTwoLevelTransferObject> &transfer)
{
  this->initialize_transfer_references(transfer);
}



template <int dim, typename Number, typename MemorySpace>
template <typename MGTwoLevelTransferObject>
void
MGTransferMF<dim, Number, MemorySpace>::initialize_transfer_references(
  const MGLevelObject<MGTwoLevelTransferObject> &transfer)
{
  const unsigned int min_level = transfer.min_level();
  const unsigned int max_level = transfer.max_level();

  this->transfer.resize(min_level, max_level);

  // Note that transfer[min_level] is empty and never used:
  for (unsigned int l = min_level + 1; l <= max_level; ++l)
    this->transfer[l] = &const_cast<MGTwoLevelTransferBase<VectorType> &>(
      static_cast<const MGTwoLevelTransferBase<VectorType> &>(
        Utilities::get_underlying_value(transfer[l])));
}



template <int dim, typename Number, typename MemorySpace>
template <class InVector>
void
MGTransferMF<dim, Number, MemorySpace>::initialize_dof_vector(
  const unsigned int level,
  VectorType        &vec,
  const InVector    &vec_reference,
  const bool         omit_zeroing_entries) const
{
  std::shared_ptr<const Utilities::MPI::Partitioner> partitioner;

  if (external_partitioners.empty())
    {
      partitioner = vec_reference.get_partitioner();
    }
  else
    {
      Assert(transfer.min_level() <= level && level <= transfer.max_level(),
             ExcInternalError());

      partitioner = external_partitioners[level - transfer.min_level()];
    }

  // check if vectors are already correctly initialized

  // yes: same partitioners are used
  if (vec.get_partitioner().get() == partitioner.get())
    {
      if (omit_zeroing_entries == false)
        vec = 0;
      return; // nothing to do
    }

  // yes: vectors are compatible
  if (vec.size() == partitioner->size() &&
      vec.locally_owned_size() == partitioner->locally_owned_size())
    {
      if (omit_zeroing_entries == false)
        vec = 0;
      return; // nothing to do
    }

  // no
  vec.reinit(partitioner, omit_zeroing_entries);
}



template <int dim, typename Number, typename MemorySpace>
template <class InVector>
void
MGTransferMF<dim, Number, MemorySpace>::copy_to_mg(
  const DoFHandler<dim>     &dof_handler,
  MGLevelObject<VectorType> &dst,
  const InVector            &src) const
{
  assert_dof_handler(dof_handler);

  for (unsigned int level = dst.min_level(); level <= dst.max_level(); ++level)
    {
      const bool zero_out_values =
        (this->perform_plain_copy == false &&
         this->perform_renumbered_plain_copy == false) ||
        level != dst.max_level();

      this->initialize_dof_vector(level, dst[level], src, !zero_out_values);
    }

  if (this->perform_plain_copy)
    {
      dst[dst.max_level()].copy_locally_owned_data_from(src);
    }
  else if (this->perform_renumbered_plain_copy)
    {
      auto &dst_level = dst[dst.max_level()];

      for (unsigned int i = 0; i < this->copy_indices.back().n_cols(); ++i)
        dst_level.local_element(this->copy_indices.back()(1, i)) =
          src.local_element(i);
    }
  else
    {
      this->ghosted_global_vector = src;
      this->ghosted_global_vector.update_ghost_values();

      for (unsigned int l = dst.max_level() + 1; l != dst.min_level();)
        {
          --l;

          auto &dst_level = dst[l];

          const auto copy_unknowns = [&](const auto &indices) {
            for (unsigned int i = 0; i < indices.n_cols(); ++i)
              dst_level.local_element(indices(1, i)) =
                this->ghosted_global_vector.local_element(indices(0, i));
          };

          copy_unknowns(this->copy_indices[l]);
          copy_unknowns(this->copy_indices_level_mine[l]);

          dst_level.compress(VectorOperation::insert);
        }
    }
}



template <int dim, typename Number, typename MemorySpace>
template <class OutVector>
void
MGTransferMF<dim, Number, MemorySpace>::copy_from_mg(
  const DoFHandler<dim>           &dof_handler,
  OutVector                       &dst,
  const MGLevelObject<VectorType> &src) const
{
  assert_dof_handler(dof_handler);

  if (this->perform_plain_copy)
    {
      dst.zero_out_ghost_values();
      dst.copy_locally_owned_data_from(src[src.max_level()]);
    }
  else if (this->perform_renumbered_plain_copy)
    {
      const auto &src_level = src[src.max_level()];
      dst.zero_out_ghost_values();
      for (unsigned int i = 0; i < this->copy_indices.back().n_cols(); ++i)
        dst.local_element(i) =
          src_level.local_element(this->copy_indices.back()(1, i));
    }
  else
    {
      dst = 0;
      for (unsigned int l = src.min_level(); l <= src.max_level(); ++l)
        {
          auto &ghosted_vector = this->ghosted_level_vector[l];

          if (this->ghosted_level_vector[l].size() > 0)
            ghosted_vector = src[l];

          const auto *const ghosted_vector_ptr =
            (this->ghosted_level_vector[l].size() > 0) ? &ghosted_vector :
                                                         &src[l];

          ghosted_vector_ptr->update_ghost_values();

          const auto copy_unknowns = [&](const auto &indices) {
            for (unsigned int i = 0; i < indices.n_cols(); ++i)
              dst.local_element(indices(0, i)) =
                ghosted_vector_ptr->local_element(indices(1, i));
          };

          copy_unknowns(this->copy_indices[l]);
          copy_unknowns(this->copy_indices_global_mine[l]);
        }
      dst.compress(VectorOperation::insert);
    }
}



template <int dim, typename Number, typename MemorySpace>
template <class InVector>
void
MGTransferMF<dim, Number, MemorySpace>::interpolate_to_mg(
  MGLevelObject<VectorType> &dst,
  const InVector            &src) const
{
  DoFHandler<dim> dof_handler_dummy;

  this->interpolate_to_mg(dof_handler_dummy, dst, src);
}



template <int dim, typename Number, typename MemorySpace>
template <class InVector>
void
MGTransferMF<dim, Number, MemorySpace>::interpolate_to_mg(
  const DoFHandler<dim>     &dof_handler,
  MGLevelObject<VectorType> &dst,
  const InVector            &src) const
{
  assert_dof_handler(dof_handler);

  const unsigned int min_level = transfer.min_level();
  const unsigned int max_level = transfer.max_level();

  AssertDimension(min_level, dst.min_level());
  AssertDimension(max_level, dst.max_level());

  for (unsigned int level = min_level; level <= max_level; ++level)
    {
      const bool zero_out_values = false;
      this->initialize_dof_vector(level, dst[level], src, !zero_out_values);
    }

  if (this->perform_plain_copy)
    {
      dst[max_level].copy_locally_owned_data_from(src);

      for (unsigned int l = max_level; l > min_level; --l)
        this->transfer[l]->interpolate(dst[l - 1], dst[l]);
    }
  else if (this->perform_renumbered_plain_copy)
    {
      auto &dst_level = dst[max_level];

      for (unsigned int i = 0; i < this->solution_copy_indices.back().n_cols();
           ++i)
        dst_level.local_element(this->solution_copy_indices.back()(1, i)) =
          src.local_element(i);

      for (unsigned int l = max_level; l > min_level; --l)
        this->transfer[l]->interpolate(dst[l - 1], dst[l]);
    }
  else
    {
      this->solution_ghosted_global_vector = src;
      this->solution_ghosted_global_vector.update_ghost_values();

      for (unsigned int l = max_level + 1; l != min_level;)
        {
          --l;

          auto &dst_level = dst[l];

          const auto copy_unknowns = [&](const auto &indices) {
            for (unsigned int i = 0; i < indices.n_cols(); ++i)
              dst_level.local_element(indices(1, i)) =
                this->solution_ghosted_global_vector.local_element(
                  indices(0, i));
          };

          copy_unknowns(this->solution_copy_indices[l]);
          copy_unknowns(this->solution_copy_indices_level_mine[l]);

          dst_level.compress(VectorOperation::insert);

          if (l != min_level)
            this->transfer[l]->interpolate(dst[l - 1], dst[l]);
        }
    }
}

#endif

DEAL_II_NAMESPACE_CLOSE

#endif
