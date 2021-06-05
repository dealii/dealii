// ---------------------------------------------------------------------
//
// Copyright (C) 2020 - 2021 by the deal.II authors
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

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/multigrid/mg_base.h>

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
   * @p tria and as a consequence not allowing to directly use it for
   * coarsening, requiring that internally a temporal copy is created.
   */
  template <int dim, int spacedim>
  std::vector<std::shared_ptr<const Triangulation<dim, spacedim>>>
  create_geometric_coarsening_sequence(
    const Triangulation<dim, spacedim> &                  tria,
    const RepartitioningPolicyTools::Base<dim, spacedim> &policy,
    const bool repartition_fine_triangulation = false);

} // namespace MGTransferGlobalCoarseningTools



/**
 * Class for transfer between two multigrid levels for p- or global coarsening.
 */
template <int dim, typename VectorType>
class MGTwoLevelTransfer
{
public:
  /**
   * Perform prolongation.
   */
  void
  prolongate(VectorType &dst, const VectorType &src) const;

  /**
   * Perform restriction.
   */
  void
  restrict_and_add(VectorType &dst, const VectorType &src) const;

  /**
   * Perform interpolation of a solution vector from the fine level to the
   * coarse level. This function is different from restriction, where a
   * weighted residual is transferred to a coarser level (transposition of
   * prolongation matrix).
   */
  void
  interpolate(VectorType &dst, const VectorType &src) const;
};



/**
 * Class for transfer between two multigrid levels for p- or global coarsening.
 * Specialization for LinearAlgebra::distributed::Vector.
 */
template <int dim, typename Number>
class MGTwoLevelTransfer<dim, LinearAlgebra::distributed::Vector<Number>>
{
public:
  /**
   * Set up global coarsening between the given DoFHandler objects (
   * @p dof_handler_fine and @p dof_handler_coarse). The transfer
   * can be only performed on active levels.
   */
  void
  reinit_geometric_transfer(const DoFHandler<dim> &          dof_handler_fine,
                            const DoFHandler<dim> &          dof_handler_coarse,
                            const AffineConstraints<Number> &constraint_fine,
                            const AffineConstraints<Number> &constraint_coarse);

  /**
   * Set up polynomial coarsening between the given DoFHandler objects (
   * @p dof_handler_fine and @p dof_handler_coarse). Polynomial transfers
   * can be only performed on active levels (`numbers::invalid_unsigned_int`)
   * or on coarse-grid levels.
   *
   * @note The function polynomial_transfer_supported() can be used to
   *   check if the given polynomial coarsening strategy is supported.
   */
  void
  reinit_polynomial_transfer(
    const DoFHandler<dim> &          dof_handler_fine,
    const DoFHandler<dim> &          dof_handler_coarse,
    const AffineConstraints<Number> &constraint_fine,
    const AffineConstraints<Number> &constraint_coarse,
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
   * Perform prolongation.
   */
  void
  prolongate(LinearAlgebra::distributed::Vector<Number> &      dst,
             const LinearAlgebra::distributed::Vector<Number> &src) const;

  /**
   * Perform restriction.
   */
  void
  restrict_and_add(LinearAlgebra::distributed::Vector<Number> &      dst,
                   const LinearAlgebra::distributed::Vector<Number> &src) const;

  /**
   * Perform interpolation of a solution vector from the fine level to the
   * coarse level. This function is different from restriction, where a
   * weighted residual is transferred to a coarser level (transposition of
   * prolongation matrix).
   */
  void
  interpolate(LinearAlgebra::distributed::Vector<Number> &      dst,
              const LinearAlgebra::distributed::Vector<Number> &src) const;

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
     */
    unsigned int dofs_per_cell_coarse;

    /**
     * Number of degrees of freedom of fine cell.
     */
    unsigned int dofs_per_cell_fine;

    /**
     * Polynomial degree of the finite element of the coarse cells.
     */
    unsigned int degree_coarse;

    /**
     * Polynomial degree of the finite element of the fine cells.
     */
    unsigned int degree_fine;

    /**
     * Weights for continuous elements.
     */
    std::vector<Number> weights;

    /**
     * Prolongation matrix for non-tensor-product elements.
     */
    AlignedVector<VectorizedArray<Number>> prolongation_matrix;

    /**
     * 1D prolongation matrix for tensor-product elements.
     */
    AlignedVector<VectorizedArray<Number>> prolongation_matrix_1d;

    /**
     * Restriction matrix for non-tensor-product elements.
     */
    AlignedVector<VectorizedArray<Number>> restriction_matrix;

    /**
     * 1D restriction matrix for tensor-product elements.
     */
    AlignedVector<VectorizedArray<Number>> restriction_matrix_1d;

    /**
     * DoF indices of the coarse cells, expressed in indices local to the MPI
     * rank.
     */
    std::vector<unsigned int> level_dof_indices_coarse;

    /**
     * DoF indices of the fine cells, expressed in indices local to the MPI
     * rank.
     */
    std::vector<unsigned int> level_dof_indices_fine;
  };

  /**
   * Transfer schemes.
   */
  std::vector<MGTransferScheme> schemes;

  /**
   * Flag if the finite elements on the fine cells are continuous. If yes,
   * the multiplicity of DoF sharing a vertex/line as well as constraints have
   * to be taken into account via weights.
   */
  bool fine_element_is_continuous;

  /**
   * Partitioner needed by the intermediate vector.
   */
  std::shared_ptr<const Utilities::MPI::Partitioner> partitioner_fine;

  /**
   * Partitioner needed by the intermediate vector.
   */
  std::shared_ptr<const Utilities::MPI::Partitioner> partitioner_coarse;

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
   * Internal vector for performing manual constraint_coarse.distribute(), which
   * is needed for acceptable performance.
   */
  mutable LinearAlgebra::distributed::Vector<Number> vec_coarse_constraints;

  /**
   * Constraint-entry indices for manually performing
   * constraint_coarse.distribute() in MPI-local indices (for performance
   * reasons).
   */
  std::vector<unsigned int> constraint_coarse_distribute_indices;

  /**
   * Constraint-entry values for manually performing
   * constraint_coarse.distribute() in MPI-local indices (for performance
   * reasons).
   */
  std::vector<Number> constraint_coarse_distribute_values;

  /**
   * Pointers to the constraint entries for performing manual
   * constraint_coarse.distribute().
   */
  std::vector<unsigned int> constraint_coarse_distribute_ptr;

  /**
   * Constraint-entry indices for performing manual
   * constraint_coarse.distribute_local_to_global().
   */
  std::vector<unsigned int> distribute_local_to_global_indices;

  /**
   * Constraint-entry values for performing manual
   * constraint_coarse.distribute_local_to_global().
   */
  std::vector<Number> distribute_local_to_global_values;

  /**
   * Pointers to the constraint entries for performing manual
   * constraint_coarse.distribute_local_to_global().
   */
  std::vector<unsigned int> distribute_local_to_global_ptr;

  /**
   * Number of components.
   */
  unsigned int n_components;

  friend class internal::MGTwoLevelTransferImplementation;
};



/**
 * Implementation of the MGTransferBase. In contrast to
 * other multigrid transfer operators, the user can provide separate
 * transfer operators of type MGTwoLevelTransfer between each level.
 *
 * This class currently only works for tensor-product finite elements based on
 * FE_Q and FE_DGQ elements. Systems involving multiple components of
 * one of these element, as well as, systems with different elements or other
 * elements are currently not implemented.
 */
template <int dim, typename VectorType>
class MGTransferGlobalCoarsening : public dealii::MGTransferBase<VectorType>
{
public:
  /**
   * Value type.
   */
  using Number = typename VectorType::value_type;

  /**
   * Constructor taking a collection of transfer operators (with the coarsest
   * level kept
   * empty in @p transfer) and an optional function that initializes the
   * internal level vectors within the function call copy_to_mg() if used in the
   * context of PreconditionMG.
   */
  MGTransferGlobalCoarsening(
    const MGLevelObject<MGTwoLevelTransfer<dim, VectorType>> &transfer,
    const std::function<void(const unsigned int, VectorType &)>
      &initialize_dof_vector = {});

  /**
   * Perform prolongation.
   */
  void
  prolongate(const unsigned int to_level,
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
  template <class InVector, int spacedim>
  void
  copy_to_mg(const DoFHandler<dim, spacedim> &dof_handler,
             MGLevelObject<VectorType> &      dst,
             const InVector &                 src) const;

  /**
   * Initialize internal vectors and copy the values on the finest
   * multigrid level to @p dst vector.
   *
   * @note DoFHandler is not needed here, but is required by the interface.
   */
  template <class OutVector, int spacedim>
  void
  copy_from_mg(const DoFHandler<dim, spacedim> &dof_handler,
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
   *
   * @note DoFHandler is not needed here, but is required by the interface.
   */
  template <class InVector, int spacedim>
  void
  interpolate_to_mg(const DoFHandler<dim, spacedim> &dof_handler,
                    MGLevelObject<VectorType> &      dst,
                    const InVector &                 src) const;

private:
  /**
   * Collection of the two-level transfer operators.
   */
  const MGLevelObject<MGTwoLevelTransfer<dim, VectorType>> &transfer;

  /**
   * %Function to initialize internal level vectors.
   */
  const std::function<void(const unsigned int, VectorType &)>
    initialize_dof_vector;
};



#ifndef DOXYGEN

/* ----------------------- Inline functions --------------------------------- */



template <int dim, typename VectorType>
MGTransferGlobalCoarsening<dim, VectorType>::MGTransferGlobalCoarsening(
  const MGLevelObject<MGTwoLevelTransfer<dim, VectorType>> &transfer,
  const std::function<void(const unsigned int, VectorType &)>
    &initialize_dof_vector)
  : transfer(transfer)
  , initialize_dof_vector(initialize_dof_vector)
{}



template <int dim, typename VectorType>
void
MGTransferGlobalCoarsening<dim, VectorType>::prolongate(
  const unsigned int to_level,
  VectorType &       dst,
  const VectorType & src) const
{
  this->transfer[to_level].prolongate(dst, src);
}



template <int dim, typename VectorType>
void
MGTransferGlobalCoarsening<dim, VectorType>::restrict_and_add(
  const unsigned int from_level,
  VectorType &       dst,
  const VectorType & src) const
{
  this->transfer[from_level].restrict_and_add(dst, src);
}



template <int dim, typename VectorType>
template <class InVector, int spacedim>
void
MGTransferGlobalCoarsening<dim, VectorType>::copy_to_mg(
  const DoFHandler<dim, spacedim> &dof_handler,
  MGLevelObject<VectorType> &      dst,
  const InVector &                 src) const
{
  (void)dof_handler;

  Assert(
    initialize_dof_vector,
    ExcMessage(
      "To be able to use this function, a function to initialize an internal "
      "DoF vector has to be provided in the constructor of "
      "MGTransferGlobalCoarsening."));

  for (unsigned int level = dst.min_level(); level <= dst.max_level(); ++level)
    initialize_dof_vector(level, dst[level]);

  dst[dst.max_level()].copy_locally_owned_data_from(src);
}



template <int dim, typename VectorType>
template <class OutVector, int spacedim>
void
MGTransferGlobalCoarsening<dim, VectorType>::copy_from_mg(
  const DoFHandler<dim, spacedim> &dof_handler,
  OutVector &                      dst,
  const MGLevelObject<VectorType> &src) const
{
  (void)dof_handler;

  dst.copy_locally_owned_data_from(src[src.max_level()]);
}



template <int dim, typename VectorType>
template <class InVector, int spacedim>
void
MGTransferGlobalCoarsening<dim, VectorType>::interpolate_to_mg(
  const DoFHandler<dim, spacedim> &dof_handler,
  MGLevelObject<VectorType> &      dst,
  const InVector &                 src) const
{
  (void)dof_handler;

  Assert(
    initialize_dof_vector,
    ExcMessage(
      "To be able to use this function, a function to initialize an internal "
      "DoF vector has to be provided in the constructor of "
      "MGTransferGlobalCoarsening."));

  const unsigned int min_level = transfer.min_level();
  const unsigned int max_level = transfer.max_level();

  AssertDimension(min_level, dst.min_level());
  AssertDimension(max_level, dst.max_level());

  for (unsigned int level = min_level; level <= max_level; ++level)
    initialize_dof_vector(level, dst[level]);

  dst[transfer.max_level()].copy_locally_owned_data_from(src);

  for (unsigned int l = max_level; l > min_level; --l)
    this->transfer[l].interpolate(dst[l - 1], dst[l]);
}

#endif

DEAL_II_NAMESPACE_CLOSE

#endif
