// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2023 by the deal.II authors
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

#ifndef dealii_mg_transfer_matrix_free_h
#define dealii_mg_transfer_matrix_free_h

#include <deal.II/base/config.h>

#include <deal.II/base/mg_level_object.h>
#include <deal.II/base/mpi_remote_point_evaluation.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/matrix_free/constraint_info.h>
#include <deal.II/matrix_free/shape_info.h>

#include <deal.II/multigrid/mg_base.h>
#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/multigrid/mg_transfer.h>
#include <deal.II/multigrid/mg_transfer_internal.h>
#include <deal.II/multigrid/mg_transfer_matrix_free.h>

#include <deal.II/non_matching/mapping_info.h>


DEAL_II_NAMESPACE_OPEN


// Forward declarations
#ifndef DOXYGEN
namespace internal
{
  class MGTwoLevelTransferImplementation;
}

template <int dim, typename Number>
class MGTransferMatrixFree;
#endif


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
 * Abstract base class for transfer operators between two multigrid levels.
 */
template <typename VectorType>
class MGTwoLevelTransferBase : public Subscriptor
{
public:
  /**
   * Perform prolongation.
   */
  virtual void
  prolongate_and_add(VectorType &dst, const VectorType &src) const = 0;

  /**
   * Perform restriction.
   */
  virtual void
  restrict_and_add(VectorType &dst, const VectorType &src) const = 0;

  /**
   * Perform interpolation of a solution vector from the fine level to the
   * coarse level. This function is different from restriction, where a
   * weighted residual is transferred to a coarser level (transposition of
   * prolongation matrix).
   */
  virtual void
  interpolate(VectorType &dst, const VectorType &src) const = 0;

  /**
   * Enable inplace vector operations if external and internal vectors
   * are compatible.
   */
  virtual void
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
};



/**
 * Base class for transfer operators between two multigrid levels.
 * Specialization for LinearAlgebra::distributed::Vector. The implementation of
 * restriction and prolongation between levels is delegated to derived classes,
 * which implement prolongate_and_add_internal() and restrict_and_add_internal()
 * accordingly.
 */
template <typename Number>
class MGTwoLevelTransferBase<LinearAlgebra::distributed::Vector<Number>>
  : public Subscriptor
{
public:
  using VectorType = LinearAlgebra::distributed::Vector<Number>;

  /**
   * Default constructor.
   */
  MGTwoLevelTransferBase();

  /**
   * Perform prolongation.
   */
  void
  prolongate_and_add(VectorType &dst, const VectorType &src) const;

  /**
   * Perform restriction.
   */
  void
  restrict_and_add(VectorType &dst, const VectorType &src) const;

  /**
   * Perform interpolation of a solution vector from the fine level to the
   * coarse level.
   */
  virtual void
  interpolate(VectorType &dst, const VectorType &src) const = 0;

  /**
   * Enable inplace vector operations if external and internal vectors
   * are compatible.
   */
  virtual void
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
  prolongate_and_add_internal(
    LinearAlgebra::distributed::Vector<Number>       &dst,
    const LinearAlgebra::distributed::Vector<Number> &src) const = 0;

  /**
   * Perform restriction on vectors with correct ghosting.
   */
  virtual void
  restrict_and_add_internal(
    LinearAlgebra::distributed::Vector<Number>       &dst,
    const LinearAlgebra::distributed::Vector<Number> &src) const = 0;

  /**
   * A wrapper around update_ghost_values() optimized in case the
   * present vector has the same parallel layout of one of the external
   * partitioners.
   */
  void
  update_ghost_values(
    const LinearAlgebra::distributed::Vector<Number> &vec) const;

  /**
   * A wrapper around compress() optimized in case the
   * present vector has the same parallel layout of one of the external
   * partitioners.
   */
  void
  compress(LinearAlgebra::distributed::Vector<Number> &vec,
           const VectorOperation::values               op) const;

  /**
   * A wrapper around zero_out_ghost_values() optimized in case the
   * present vector has the same parallel layout of one of the external
   * partitioners.
   */
  void
  zero_out_ghost_values(
    const LinearAlgebra::distributed::Vector<Number> &vec) const;

  /**
   * Enable inplace vector operations if external and internal vectors
   * are compatible.
   */
  template <int dim, std::size_t width>
  void
  internal_enable_inplace_operations_if_possible(
    const std::shared_ptr<const Utilities::MPI::Partitioner>
      &partitioner_coarse,
    const std::shared_ptr<const Utilities::MPI::Partitioner> &partitioner_fine,
    bool &vec_fine_needs_ghost_update,
    internal::MatrixFreeFunctions::ConstraintInfo<
      dim,
      VectorizedArray<Number, width>> &constraint_info_coarse,
    std::vector<unsigned int>         &dof_indices_fine);

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
   * Internal vector on that the actual prolongation/restriction is performed.
   */
  mutable LinearAlgebra::distributed::Vector<Number> vec_coarse;

  /**
   * Internal vector needed for collecting all degrees of freedom of the fine
   * cells. It is only initialized if the fine-level DoF indices touch DoFs
   * other than the locally active ones (which we always assume can be
   * accessed by the given vectors in the prolongate/restrict functions),
   * otherwise it is left at size zero.
   */
  mutable LinearAlgebra::distributed::Vector<Number> vec_fine;

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
 * Class for transfer between two multigrid levels for p- or global coarsening.
 *
 * The implementation of this class is explained in detail in @cite munch2022gc.
 */
template <int dim, typename VectorType>
class MGTwoLevelTransfer : public MGTwoLevelTransferBase<VectorType>
{
public:
  /**
   * Perform prolongation.
   */
  void
  prolongate_and_add(VectorType &dst, const VectorType &src) const override;

  /**
   * Perform restriction.
   */
  void
  restrict_and_add(VectorType &dst, const VectorType &src) const override;

  /**
   * Perform interpolation of a solution vector from the fine level to the
   * coarse level.
   */
  void
  interpolate(VectorType &dst, const VectorType &src) const override;

  /**
   * Enable inplace vector operations if external and internal vectors
   * are compatible.
   */
  void
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
};



/**
 * Class for transfer between two multigrid levels for p- or global coarsening.
 * Specialization for LinearAlgebra::distributed::Vector.
 *
 * The implementation of this class is explained in detail in @cite munch2022gc.
 */
template <int dim, typename Number>
class MGTwoLevelTransfer<dim, LinearAlgebra::distributed::Vector<Number>>
  : public MGTwoLevelTransferBase<LinearAlgebra::distributed::Vector<Number>>
{
  using VectorizedArrayType = VectorizedArray<Number>;

public:
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
   * Set up transfer operator between the given DoFHandler objects (
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
   * Perform interpolation of a solution vector from the fine level to the
   * coarse level.
   */
  void
  interpolate(
    LinearAlgebra::distributed::Vector<Number>       &dst,
    const LinearAlgebra::distributed::Vector<Number> &src) const override;

  /**
   * Enable inplace vector operations if external and internal vectors
   * are compatible.
   */
  void
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
  prolongate_and_add_internal(
    LinearAlgebra::distributed::Vector<Number>       &dst,
    const LinearAlgebra::distributed::Vector<Number> &src) const override;

  void
  restrict_and_add_internal(
    LinearAlgebra::distributed::Vector<Number>       &dst,
    const LinearAlgebra::distributed::Vector<Number> &src) const override;

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
     * Prolongation matrix for non-tensor-product elements.
     */
    AlignedVector<VectorizedArrayType> prolongation_matrix;

    /**
     * 1d prolongation matrix for tensor-product elements.
     */
    AlignedVector<VectorizedArrayType> prolongation_matrix_1d;

    /**
     * Restriction matrix for non-tensor-product elements.
     */
    AlignedVector<VectorizedArrayType> restriction_matrix;

    /**
     * 1d restriction matrix for tensor-product elements.
     */
    AlignedVector<VectorizedArrayType> restriction_matrix_1d;

    /**
     * ShapeInfo description of the coarse cell. Needed during the
     * fast application of hanging-node constraints.
     */
    internal::MatrixFreeFunctions::ShapeInfo<VectorizedArrayType>
      shape_info_coarse;
  };

  /**
   * Transfer schemes.
   */
  std::vector<MGTransferScheme> schemes;

  /**
   * Helper class for reading from and writing to global coarse vectors and for
   * applying constraints.
   */
  internal::MatrixFreeFunctions::ConstraintInfo<dim, VectorizedArrayType>
    constraint_info_coarse;

  /**
   * Helper class for reading from and writing to global fine vectors.
   */
  internal::MatrixFreeFunctions::ConstraintInfo<dim, VectorizedArrayType>
    constraint_info_fine;

  /**
   * Weights for continuous elements.
   */
  std::vector<Number> weights; // TODO: vectorize

  /**
   * Weights for continuous elements, compressed into 3^dim doubles per
   * cell if possible.
   */
  AlignedVector<VectorizedArrayType> weights_compressed;

  /**
   * Number of components.
   */
  unsigned int n_components;

  /**
   * Pointer to the DoFHandler object used during initialization.
   */
  SmartPointer<const DoFHandler<dim>> dof_handler_fine;

  /**
   * Muligird level used during initialization.
   */
  unsigned int mg_level_fine;

  friend class internal::MGTwoLevelTransferImplementation;

  friend class MGTransferMatrixFree<dim, Number>;
};



/**
 * Class for transfer between two non-nested multigrid levels.
 *
 */
template <int dim, typename VectorType>
class MGTwoLevelTransferNonNested : public MGTwoLevelTransferBase<VectorType>
{
public:
  /**
   * Perform prolongation.
   */
  void
  prolongate_and_add(VectorType &dst, const VectorType &src) const override;

  /**
   * Perform restriction.
   */
  void
  restrict_and_add(VectorType &dst, const VectorType &src) const override;

  /**
   * Perform interpolation of a solution vector from the fine level to the
   * coarse level. This function is different from restriction, where a
   * weighted residual is transferred to a coarser level (transposition of
   * prolongation matrix).
   */
  void
  interpolate(VectorType &dst, const VectorType &src) const override;

  /**
   * Return the memory consumption of the allocated memory in this class.
   */
  std::size_t
  memory_consumption() const override;
};



/**
 * Class for transfer between two non-nested multigrid levels.
 *
 * Specialization for LinearAlgebra::distributed::Vector.
 *
 */
template <int dim, typename Number>
class MGTwoLevelTransferNonNested<dim,
                                  LinearAlgebra::distributed::Vector<Number>>
  : public MGTwoLevelTransferBase<LinearAlgebra::distributed::Vector<Number>>
{
private:
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
  interpolate(
    LinearAlgebra::distributed::Vector<Number>       &dst,
    const LinearAlgebra::distributed::Vector<Number> &src) const override;

  /**
   * Enable inplace vector operations if external and internal vectors
   * are compatible.
   */
  void
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
  prolongate_and_add_internal(
    LinearAlgebra::distributed::Vector<Number>       &dst,
    const LinearAlgebra::distributed::Vector<Number> &src) const override;

  /**
   * Perform restriction.
   */
  void
  restrict_and_add_internal(
    LinearAlgebra::distributed::Vector<Number>       &dst,
    const LinearAlgebra::distributed::Vector<Number> &src) const override;

private:
  /**
   * Perform prolongation for correct number of components.
   */
  template <int n_components>
  void
  prolongate_and_add_internal_comp(
    LinearAlgebra::distributed::Vector<Number>       &dst,
    const LinearAlgebra::distributed::Vector<Number> &src) const;

  /**
   * Perform restriction for correct number of components.
   */
  template <int n_components>
  void
  restrict_and_add_internal_comp(
    LinearAlgebra::distributed::Vector<Number>       &dst,
    const LinearAlgebra::distributed::Vector<Number> &src) const;

  /**
   * Pointer to the DoFHandler object used during initialization.
   */
  SmartPointer<const DoFHandler<dim>> dof_handler_fine;

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
   * MappingInfo object needed as Mapping argument by FEPointEvaluation.
   */
  std::shared_ptr<NonMatching::MappingInfo<dim, dim, Number>> mapping_info;

  /**
   * Helper class for reading from and writing to global vectors and for
   * applying constraints.
   */
  internal::MatrixFreeFunctions::ConstraintInfo<dim, VectorizedArrayType>
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

  friend class MGTransferMatrixFree<dim, Number>;
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
 * FE_Q and FE_DGQ and simplex elements FE_SimplexP and FE_SimplexDGP as well as
 * for systems involving multiple components of one of these elements. Other
 * elements are currently not implemented.
 *
 * The implementation of this class is explained in detail in @cite munch2022gc.
 */
template <int dim, typename Number>
class MGTransferMatrixFree : public dealii::MGLevelGlobalTransfer<
                               LinearAlgebra::distributed::Vector<Number>>
{
public:
  /**
   * Value type.
   */
  using VectorType = LinearAlgebra::distributed::Vector<Number>;

  /**
   * Default constructor.
   *
   * @note See also MGTransferMatrixFree.
   */
  MGTransferMatrixFree();

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
  MGTransferMatrixFree(
    const MGLevelObject<MGTwoLevelTransferObject> &transfer,
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
  MGTransferMatrixFree(const MGConstrainedDoFs &mg_constrained_dofs);

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

  /**
   * Indicate if global coarsening or local smoothing is performed. Note
   * that in the case local smoothing copy_to_mg() and copy_from_mg()
   * are more complex and edge constraints need to be considered during
   * the multigrid algorithm.
   */
  bool
  performs_global_coarsening() const;

  /** @} */

private:
  /**
   * Initial internal transfer operator.
   *
   * @note See also MGTransferMatrixFree.
   */
  void
  initialize_internal_transfer(
    const DoFHandler<dim>                       &dof_handler,
    const SmartPointer<const MGConstrainedDoFs> &mg_constrained_dofs);

  /**
   * Retrieve finest DoFHandler from two-level transfer objects.
   */
  std::pair<const DoFHandler<dim> *, unsigned int>
  get_dof_handler_fine() const;

  /**
   * Initialize copy indices for MGTransferMatrixFree::copy_to_mg(),
   * MGTransferMatrixFree::copy_to_mg(), and
   * MGTransferMatrixFree::interpolate_to_mg() in the case of global coarsening.
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
  MGLevelObject<SmartPointer<MGTwoLevelTransferBase<VectorType>>> transfer;

  /**
   * External partitioners used during initialize_dof_vector().
   */
  std::vector<std::shared_ptr<const Utilities::MPI::Partitioner>>
    external_partitioners;
};



/**
 * Base class of MGTransferBlockMatrixFree. While MGTransferBlockMatrixFree
 * contains all the setup routines of the transfer operators for the blocks,
 * this class simply applies them, e.g., for restricting and prolongating.
 */
template <int dim, typename Number, typename TransferType>
class MGTransferBlockMatrixFreeBase
  : public MGTransferBase<LinearAlgebra::distributed::BlockVector<Number>>
{
public:
  MGTransferBlockMatrixFreeBase(const bool same_for_all)
    : same_for_all(same_for_all)
  {}

  /**
   * Prolongate a vector from level <tt>to_level-1</tt> to level
   * <tt>to_level</tt> using the embedding matrices of the underlying finite
   * element. The previous content of <tt>dst</tt> is overwritten.
   *
   * @param to_level The index of the level to prolongate to, which is the
   * level of @p dst.
   *
   * @param src is a vector with as many elements as there are degrees of
   * freedom on the coarser level involved.
   *
   * @param dst has as many elements as there are degrees of freedom on the
   * finer level.
   */
  virtual void
  prolongate(
    const unsigned int                                     to_level,
    LinearAlgebra::distributed::BlockVector<Number>       &dst,
    const LinearAlgebra::distributed::BlockVector<Number> &src) const override;

  virtual void
  prolongate_and_add(
    const unsigned int                                     to_level,
    LinearAlgebra::distributed::BlockVector<Number>       &dst,
    const LinearAlgebra::distributed::BlockVector<Number> &src) const override;

  /**
   * Restrict a vector from level <tt>from_level</tt> to level
   * <tt>from_level-1</tt> using the transpose operation of the prolongate()
   * method. If the region covered by cells on level <tt>from_level</tt> is
   * smaller than that of level <tt>from_level-1</tt> (local refinement), then
   * some degrees of freedom in <tt>dst</tt> are active and will not be
   * altered. For the other degrees of freedom, the result of the restriction
   * is added.
   *
   * @param from_level The index of the level to restrict from, which is the
   * level of @p src.
   *
   * @param src is a vector with as many elements as there are degrees of
   * freedom on the finer level involved.
   *
   * @param dst has as many elements as there are degrees of freedom on the
   * coarser level.
   */
  virtual void
  restrict_and_add(
    const unsigned int                                     from_level,
    LinearAlgebra::distributed::BlockVector<Number>       &dst,
    const LinearAlgebra::distributed::BlockVector<Number> &src) const override;

  /**
   * Transfer from a block-vector on the global grid to block-vectors defined
   * on each of the levels separately for active degrees of freedom.
   * In particular, for a globally refined mesh only the finest level in @p dst
   * is filled as a plain copy of @p src. All the other level objects are left
   * untouched.
   *
   * This function will initialize @p dst accordingly if needed as required by
   * the Multigrid class.
   */
  template <typename BlockVectorType2>
  void
  copy_to_mg(
    const DoFHandler<dim>                                          &dof_handler,
    MGLevelObject<LinearAlgebra::distributed::BlockVector<Number>> &dst,
    const BlockVectorType2                                         &src) const;

  /**
   * Same as above for the case that each block has its own DoFHandler.
   */
  template <typename BlockVectorType2>
  void
  copy_to_mg(
    const std::vector<const DoFHandler<dim> *>                     &dof_handler,
    MGLevelObject<LinearAlgebra::distributed::BlockVector<Number>> &dst,
    const BlockVectorType2                                         &src) const;

  /**
   * Transfer from multi-level block-vector to normal vector.
   */
  template <typename BlockVectorType2>
  void
  copy_from_mg(
    const DoFHandler<dim> &dof_handler,
    BlockVectorType2      &dst,
    const MGLevelObject<LinearAlgebra::distributed::BlockVector<Number>> &src)
    const;

  /**
   * Same as above for the case that each block has its own DoFHandler.
   */
  template <typename BlockVectorType2>
  void
  copy_from_mg(
    const std::vector<const DoFHandler<dim> *> &dof_handler,
    BlockVectorType2                           &dst,
    const MGLevelObject<LinearAlgebra::distributed::BlockVector<Number>> &src)
    const;

  /**
   * This class can both be used with a single DoFHandler
   * or a separate DoFHandler for each block.
   */
  static const bool supports_dof_handler_vector = true;

protected:
  /**
   * Return the right non-block transfer operator. Has to be implemented by
   * the derived class.
   */
  virtual const TransferType &
  get_matrix_free_transfer(const unsigned int b) const = 0;

  /**
   * A flag to indicate whether the same DoFHandler is used for all
   * the components or if each block has its own DoFHandler.
   */
  const bool same_for_all;
};



/**
 * This class works with LinearAlgebra::distributed::BlockVector and
 * performs exactly the same transfer operations for each block as
 * MGTransferMatrixFree.
 */
template <int dim, typename Number>
class MGTransferBlockMatrixFree
  : public MGTransferBlockMatrixFreeBase<dim,
                                         Number,
                                         MGTransferMatrixFree<dim, Number>>
{
public:
  /**
   * Constructor.
   */
  MGTransferBlockMatrixFree(
    const MGTransferMatrixFree<dim, Number> &transfer_operator);

  /**
   * Constructor.
   *
   * @note See also MGTransferBlockMatrixFree.
   */
  MGTransferBlockMatrixFree() = default;

  /**
   * Constructor.
   *
   * @note See also MGTransferBlockMatrixFree.
   */
  MGTransferBlockMatrixFree(const MGConstrainedDoFs &mg_constrained_dofs);

  /**
   * Constructor.
   *
   * @note See also MGTransferBlockMatrixFree.
   */
  MGTransferBlockMatrixFree(
    const std::vector<MGConstrainedDoFs> &mg_constrained_dofs);

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
  const MGTransferMatrixFree<dim, Number> &
  get_matrix_free_transfer(const unsigned int b) const override;

private:
  /**
   * Internal non-block version of transfer operation.
   */
  std::vector<MGTransferMatrixFree<dim, Number>> transfer_operators_internal;

  /**
   * Non-block version of transfer operation.
   */
  std::vector<SmartPointer<const MGTransferMatrixFree<dim, Number>>>
    transfer_operators;
};



#ifndef DOXYGEN

/* ----------------------- Inline functions --------------------------------- */



template <int dim, typename Number>
template <typename MGTwoLevelTransferObject>
MGTransferMatrixFree<dim, Number>::MGTransferMatrixFree(
  const MGLevelObject<MGTwoLevelTransferObject> &transfer,
  const std::function<void(const unsigned int, VectorType &)>
    &initialize_dof_vector)
{
  this->transfer.clear();
  this->internal_transfer.clear();

  this->initialize_transfer_references(transfer);
  this->build(initialize_dof_vector);
}



template <int dim, typename Number>
template <typename MGTwoLevelTransferObject>
void
MGTransferMatrixFree<dim, Number>::initialize_two_level_transfers(
  const MGLevelObject<MGTwoLevelTransferObject> &transfer)
{
  this->initialize_transfer_references(transfer);
}



template <int dim, typename Number>
template <typename MGTwoLevelTransferObject>
void
MGTransferMatrixFree<dim, Number>::initialize_transfer_references(
  const MGLevelObject<MGTwoLevelTransferObject> &transfer)
{
  const unsigned int min_level = transfer.min_level();
  const unsigned int max_level = transfer.max_level();

  this->transfer.resize(min_level, max_level);

  for (unsigned int l = min_level; l <= max_level; ++l)
    this->transfer[l] = &const_cast<MGTwoLevelTransferBase<VectorType> &>(
      static_cast<const MGTwoLevelTransferBase<VectorType> &>(
        Utilities::get_underlying_value(transfer[l])));
}



template <int dim, typename Number>
template <class InVector>
void
MGTransferMatrixFree<dim, Number>::initialize_dof_vector(
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



template <int dim, typename Number>
template <class InVector>
void
MGTransferMatrixFree<dim, Number>::copy_to_mg(
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



template <int dim, typename Number>
template <class OutVector>
void
MGTransferMatrixFree<dim, Number>::copy_from_mg(
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



template <int dim, typename Number>
template <class InVector>
void
MGTransferMatrixFree<dim, Number>::interpolate_to_mg(
  MGLevelObject<VectorType> &dst,
  const InVector            &src) const
{
  DoFHandler<dim> dof_handler_dummy;

  this->interpolate_to_mg(dof_handler_dummy, dst, src);
}



template <int dim, typename Number>
template <class InVector>
void
MGTransferMatrixFree<dim, Number>::interpolate_to_mg(
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

          const bool has_ghost_elements = dst_level.has_ghost_elements();

          if (has_ghost_elements)
            dst_level.zero_out_ghost_values();

          const auto copy_unknowns = [&](const auto &indices) {
            for (unsigned int i = 0; i < indices.n_cols(); ++i)
              dst_level.local_element(indices(1, i)) =
                this->solution_ghosted_global_vector.local_element(
                  indices(0, i));
          };

          copy_unknowns(this->solution_copy_indices[l]);
          copy_unknowns(this->solution_copy_indices_level_mine[l]);

          dst_level.compress(VectorOperation::insert);

          if (has_ghost_elements)
            dst_level.update_ghost_values();

          if (l != min_level)
            this->transfer[l]->interpolate(dst[l - 1], dst[l]);
        }
    }
}



template <int dim, typename Number, typename TransferType>
template <typename BlockVectorType2>
void
MGTransferBlockMatrixFreeBase<dim, Number, TransferType>::copy_to_mg(
  const DoFHandler<dim>                                          &dof_handler,
  MGLevelObject<LinearAlgebra::distributed::BlockVector<Number>> &dst,
  const BlockVectorType2                                         &src) const
{
  Assert(same_for_all,
         ExcMessage(
           "This object was initialized with support for usage with one "
           "DoFHandler for each block, but this method assumes that "
           "the same DoFHandler is used for all the blocks!"));
  const std::vector<const DoFHandler<dim> *> mg_dofs(src.n_blocks(),
                                                     &dof_handler);

  copy_to_mg(mg_dofs, dst, src);
}



template <int dim, typename Number, typename TransferType>
template <typename BlockVectorType2>
void
MGTransferBlockMatrixFreeBase<dim, Number, TransferType>::copy_to_mg(
  const std::vector<const DoFHandler<dim> *>                     &dof_handler,
  MGLevelObject<LinearAlgebra::distributed::BlockVector<Number>> &dst,
  const BlockVectorType2                                         &src) const
{
  const unsigned int n_blocks = src.n_blocks();
  AssertDimension(dof_handler.size(), n_blocks);

  if (n_blocks == 0)
    return;

  const unsigned int min_level = dst.min_level();
  const unsigned int max_level = dst.max_level();

  for (unsigned int level = min_level; level <= max_level; ++level)
    if (dst[level].n_blocks() != n_blocks)
      dst[level].reinit(n_blocks);

  // FIXME: this a quite ugly as we need a temporary object:
  MGLevelObject<LinearAlgebra::distributed::Vector<Number>> dst_non_block(
    min_level, max_level);

  for (unsigned int b = 0; b < n_blocks; ++b)
    {
      const unsigned int data_block = same_for_all ? 0 : b;
      get_matrix_free_transfer(data_block)
        .copy_to_mg(*dof_handler[b], dst_non_block, src.block(b));

      for (unsigned int l = min_level; l <= max_level; ++l)
        dst[l].block(b) = dst_non_block[l];
    }

  for (unsigned int level = min_level; level <= max_level; ++level)
    dst[level].collect_sizes();
}

template <int dim, typename Number, typename TransferType>
template <typename BlockVectorType2>
void
MGTransferBlockMatrixFreeBase<dim, Number, TransferType>::copy_from_mg(
  const DoFHandler<dim> &dof_handler,
  BlockVectorType2      &dst,
  const MGLevelObject<LinearAlgebra::distributed::BlockVector<Number>> &src)
  const
{
  Assert(same_for_all,
         ExcMessage(
           "This object was initialized with support for usage with one "
           "DoFHandler for each block, but this method assumes that "
           "the same DoFHandler is used for all the blocks!"));
  const std::vector<const DoFHandler<dim> *> mg_dofs(dst.n_blocks(),
                                                     &dof_handler);

  copy_from_mg(mg_dofs, dst, src);
}

template <int dim, typename Number, typename TransferType>
template <typename BlockVectorType2>
void
MGTransferBlockMatrixFreeBase<dim, Number, TransferType>::copy_from_mg(
  const std::vector<const DoFHandler<dim> *> &dof_handler,
  BlockVectorType2                           &dst,
  const MGLevelObject<LinearAlgebra::distributed::BlockVector<Number>> &src)
  const
{
  const unsigned int n_blocks = dst.n_blocks();
  AssertDimension(dof_handler.size(), n_blocks);

  if (n_blocks == 0)
    return;

  const unsigned int min_level = src.min_level();
  const unsigned int max_level = src.max_level();

  for (unsigned int l = min_level; l <= max_level; ++l)
    AssertDimension(src[l].n_blocks(), dst.n_blocks());

  // FIXME: this a quite ugly as we need a temporary object:
  MGLevelObject<LinearAlgebra::distributed::Vector<Number>> src_non_block(
    min_level, max_level);

  for (unsigned int b = 0; b < n_blocks; ++b)
    {
      for (unsigned int l = min_level; l <= max_level; ++l)
        {
          src_non_block[l].reinit(src[l].block(b));
          src_non_block[l] = src[l].block(b);
        }
      const unsigned int data_block = same_for_all ? 0 : b;
      get_matrix_free_transfer(data_block)
        .copy_from_mg(*dof_handler[b], dst.block(b), src_non_block);
    }
}

#endif


DEAL_II_NAMESPACE_CLOSE

#endif
