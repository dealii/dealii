// ---------------------------------------------------------------------
//
// Copyright (C) 2020 - 2023 by the deal.II authors
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

#ifndef dealii_mg_transfer_global_coarsening_h
#define dealii_mg_transfer_global_coarsening_h

#include <deal.II/base/mg_level_object.h>
#include <deal.II/base/mpi_remote_point_evaluation.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/matrix_free/constraint_info.h>
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
#endif



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
    Triangulation<dim, spacedim> &                        tria,
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
    const Triangulation<dim, spacedim> &                  tria,
    const RepartitioningPolicyTools::Base<dim, spacedim> &policy,
    const bool repartition_fine_triangulation = false);

} // namespace MGTransferGlobalCoarseningTools


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
   * Perform prolongation.
   */
  virtual void
  prolongate_and_add(VectorType &dst, const VectorType &src) const;

  /**
   * Perform restriction.
   */
  virtual void
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
    LinearAlgebra::distributed::Vector<Number> &      dst,
    const LinearAlgebra::distributed::Vector<Number> &src) const = 0;

  /**
   * Perform restriction on vectors with correct ghosting.
   */
  virtual void
  restrict_and_add_internal(
    LinearAlgebra::distributed::Vector<Number> &      dst,
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
    internal::MatrixFreeFunctions::ConstraintInfo<
      dim,
      VectorizedArray<Number, width>> &constraint_info_coarse,
    std::vector<unsigned int> &        dof_indices_fine);

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
   * Internal vector needed for collecting all degrees of freedom of the fine
   * cells. It is only initialized if the fine-level DoF indices touch DoFs
   * other than the locally active ones (which we always assume can be
   * accessed by the given vectors in the prolongate/restrict functions),
   * otherwise it is left at size zero.
   */
  mutable LinearAlgebra::distributed::Vector<Number> vec_fine;

  /**
   * Internal vector on that the actual prolongation/restriction is performed.
   */
  mutable LinearAlgebra::distributed::Vector<Number> vec_coarse;

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
    const DoFHandler<dim> &          dof_handler_fine,
    const DoFHandler<dim> &          dof_handler_coarse,
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
    const DoFHandler<dim> &          dof_handler_fine,
    const DoFHandler<dim> &          dof_handler_coarse,
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
  reinit(const DoFHandler<dim> &          dof_handler_fine,
         const DoFHandler<dim> &          dof_handler_coarse,
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
    LinearAlgebra::distributed::Vector<Number> &      dst,
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
    LinearAlgebra::distributed::Vector<Number> &      dst,
    const LinearAlgebra::distributed::Vector<Number> &src) const override;

  void
  restrict_and_add_internal(
    LinearAlgebra::distributed::Vector<Number> &      dst,
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

  friend class internal::MGTwoLevelTransferImplementation;
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

public:
  /**
   * Set up transfer operator between the given DoFHandler objects (
   * @p dof_handler_fine and @p dof_handler_coarse).
   */
  void
  reinit(const DoFHandler<dim> &          dof_handler_fine,
         const DoFHandler<dim> &          dof_handler_coarse,
         const Mapping<dim> &             mapping_fine,
         const Mapping<dim> &             mapping_coarse,
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
    LinearAlgebra::distributed::Vector<Number> &      dst,
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
  /**
   * Perform prolongation.
   */
  void
  prolongate_and_add_internal(
    LinearAlgebra::distributed::Vector<Number> &      dst,
    const LinearAlgebra::distributed::Vector<Number> &src) const override;

  /**
   * Perform restriction.
   */
  void
  restrict_and_add_internal(
    LinearAlgebra::distributed::Vector<Number> &      dst,
    const LinearAlgebra::distributed::Vector<Number> &src) const override;

private:
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
};



/**
 * Implementation of the MGTransferBase. In contrast to
 * other multigrid transfer operators, the user can provide separate
 * transfer operators of type MGTwoLevelTransfer between each level. The
 * sequence of functions calls for setup is:
 * @code
 * MGTransferGlobalCoarsening mg_transfer;
 * mg_transfer.intitialize_two_level_transfers(two_level_transfers);
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
class MGTransferMF : public dealii::MGLevelGlobalTransfer<
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
  MGTransferMF() = default;

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
  intitialize_two_level_transfers(
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

  /**
   * Actually build the information for the prolongation for each level.
   *
   * @note See also MGTransferMatrixFree.
   */
  void
  build(const DoFHandler<dim> &dof_handler,
        const std::vector<std::shared_ptr<const Utilities::MPI::Partitioner>>
          &external_partitioners = {});

  /**
   * Same as above but taking a lambda for initializing vector instead of
   * partitioners.
   *
   * @note See also MGTransferMatrixFree.
   */
  void
  build(const DoFHandler<dim> &dof_handler,
        const std::function<void(const unsigned int, VectorType &)>
          &initialize_dof_vector);

  /** @} */

  /**
   * @name Tranfer functions.
   */
  /** @{ */

  /**
   * Perform prolongation.
   */
  void
  prolongate(const unsigned int to_level,
             VectorType &       dst,
             const VectorType & src) const override;

  /**
   * Perform prolongation.
   */
  void
  prolongate_and_add(const unsigned int to_level,
                     VectorType &       dst,
                     const VectorType & src) const override;

  /**
   * Perform restriction.
   */
  virtual void
  restrict_and_add(const unsigned int from_level,
                   VectorType &       dst,
                   const VectorType & src) const override;

  /**
   * Initialize internal vectors and copy @p src vector to the finest
   * multigrid level.
   *
   * @note DoFHandler is not needed here, but is required by the interface.
   */
  template <class InVector>
  void
  copy_to_mg(const DoFHandler<dim> &    dof_handler,
             MGLevelObject<VectorType> &dst,
             const InVector &           src) const;

  /**
   * Initialize internal vectors and copy the values on the finest
   * multigrid level to @p dst vector.
   *
   * @note DoFHandler is not needed here, but is required by the interface.
   */
  template <class OutVector>
  void
  copy_from_mg(const DoFHandler<dim> &          dof_handler,
               OutVector &                      dst,
               const MGLevelObject<VectorType> &src) const;

  /**
   * Interpolate fine-mesh field @p src to each multigrid level in
   * @p dof_handler and store the result in @p dst. This function is different
   * from restriction, where a weighted residual is
   * transferred to a coarser level (transposition of prolongation matrix).
   *
   * The argument @p dst has to be initialized with the correct size according
   * to the number of levels of the triangulation.
   *
   * If an inner vector of @p dst is empty or has incorrect locally owned size,
   * it will be resized to locally relevant degrees of freedom on each level.
   */
  template <class InVector>
  void
  interpolate_to_mg(MGLevelObject<VectorType> &dst, const InVector &src) const;

  /**
   * Like the above function but with a user-provided DoFHandler as
   * additional argument. However, this DoFHandler is not used internally, but
   * is required to be able to use MGTransferMF and
   * MGTransferMatrixFree as template argument.
   */
  template <class InVector>
  void
  interpolate_to_mg(const DoFHandler<dim> &    dof_handler,
                    MGLevelObject<VectorType> &dst,
                    const InVector &           src) const;

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
  intitialize_internal_transfer(
    const DoFHandler<dim> &                      dof_handler,
    const SmartPointer<const MGConstrainedDoFs> &mg_constrained_dofs);

  /**
   * Set references to two-level transfer operators to be used.
   */
  template <typename MGTwoLevelTransferObject>
  void
  intitialize_transfer_references(
    const MGLevelObject<MGTwoLevelTransferObject> &transfer);

  /**
   * Function to initialize internal level vectors.
   */
  template <class InVector>
  void
  initialize_dof_vector(const unsigned int level,
                        VectorType &       vector,
                        const InVector &   vector_reference,
                        const bool         omit_zeroing_entries) const;

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
 * This class works with LinearAlgebra::distributed::BlockVector and
 * performs exactly the same transfer operations for each block as
 * MGTransferMF.
 */
template <int dim, typename Number>
class MGTransferBlockMF
  : public MGTransferBlockMatrixFreeBase<dim, Number, MGTransferMF<dim, Number>>
{
public:
  /**
   * Constructor.
   */
  MGTransferBlockMF(const MGTransferMF<dim, Number> &transfer_operator);

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
  const MGTransferMF<dim, Number> &
  get_matrix_free_transfer(const unsigned int b) const override;

private:
  /**
   * Internal non-block version of transfer operation.
   */
  std::vector<MGTransferMF<dim, Number>> transfer_operators_internal;

  /**
   * Non-block version of transfer operation.
   */
  std::vector<SmartPointer<const MGTransferMF<dim, Number>>> transfer_operators;
};



template <int dim, typename VectorType>
using MGTransferGlobalCoarsening =
  MGTransferMF<dim, typename VectorType::value_type>;

template <int dim, typename VectorType>
using MGTransferBlockGlobalCoarsening =
  MGTransferBlockMF<dim, typename VectorType::value_type>;



#ifndef DOXYGEN

/* ----------------------- Inline functions --------------------------------- */



template <int dim, typename Number>
template <typename MGTwoLevelTransferObject>
MGTransferMF<dim, Number>::MGTransferMF(
  const MGLevelObject<MGTwoLevelTransferObject> &transfer,
  const std::function<void(const unsigned int, VectorType &)>
    &initialize_dof_vector)
{
  this->intitialize_transfer_references(transfer);
  this->build(initialize_dof_vector);
}



template <int dim, typename Number>
template <typename MGTwoLevelTransferObject>
void
MGTransferMF<dim, Number>::intitialize_two_level_transfers(
  const MGLevelObject<MGTwoLevelTransferObject> &transfer)
{
  this->intitialize_transfer_references(transfer);
}



template <int dim, typename Number>
template <typename MGTwoLevelTransferObject>
void
MGTransferMF<dim, Number>::intitialize_transfer_references(
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
MGTransferMF<dim, Number>::initialize_dof_vector(
  const unsigned int level,
  VectorType &       vec,
  const InVector &   vec_reference,
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

  // check if vectors are already correctly initalized

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
MGTransferMF<dim, Number>::copy_to_mg(const DoFHandler<dim> &    dof_handler,
                                      MGLevelObject<VectorType> &dst,
                                      const InVector &           src) const
{
  (void)dof_handler;

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
MGTransferMF<dim, Number>::copy_from_mg(
  const DoFHandler<dim> &          dof_handler,
  OutVector &                      dst,
  const MGLevelObject<VectorType> &src) const
{
  (void)dof_handler;

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
MGTransferMF<dim, Number>::interpolate_to_mg(MGLevelObject<VectorType> &dst,
                                             const InVector &src) const
{
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



template <int dim, typename Number>
template <class InVector>
void
MGTransferMF<dim, Number>::interpolate_to_mg(const DoFHandler<dim> &dof_handler,
                                             MGLevelObject<VectorType> &dst,
                                             const InVector &src) const
{
  (void)dof_handler;

  this->interpolate_to_mg(dst, src);
}

#endif

DEAL_II_NAMESPACE_CLOSE

#endif
