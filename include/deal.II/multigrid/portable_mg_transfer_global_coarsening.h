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

#ifndef dealii_portable_mg_transfer_global_coarsening_h
#define dealii_portable_mg_transfer_global_coarsening_h

#include <deal.II/base/memory_space.h>
#include <deal.II/base/mpi_remote_point_evaluation.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/matrix_free/constraint_info.h>
#include <deal.II/matrix_free/portable_matrix_free.h>
#include <deal.II/matrix_free/shape_info.h>

#include <deal.II/multigrid/mg_base.h>
#include <deal.II/multigrid/mg_transfer_matrix_free.h>
#include <deal.II/multigrid/portable_mg_transfer_global_coarsening.h>

DEAL_II_NAMESPACE_OPEN

namespace Portable
{

  // Forward declarations
#ifndef DOXYGEN
  namespace internal
  {
    class MGTwoLevelTransferImplementation;

    template <typename Number>
    struct MGTransferSharedData;

  } // namespace internal
#endif

  /**
   * @addtogroup mg
   * @{
   */

  /**
   * Class for transfer between two multigrid levels for p- or global
   * coarsening. It relies on a list of DoF indices associated with the cells on
   * the coarse and fine side of the transfer, and implements a cell-by-cell
   * (matrix-free) interpolation setup with the reference-cell embedding
   * matrices.
   *
   * The implementation of this class is similar to MGTwoLevelTransfer which is
   * explained in detail in @cite munch2022gc.
   *
   * The transfer is built from the underlying DoFHandler and AffineConstraints
   * objects on the coarse and fine side, collecting an explicit copy of all
   * indices on both sides. For the moment, this works with FE_Q and a single
   * DoFHandler. The intergration of other types of finite elements as well as
   * polynomial transfer is underway.
   */
  template <int dim, typename VectorType>
  class MGTwoLevelTransfer
    : public dealii::internal::MGTwoLevelTransferCore<dim, VectorType>
  {
  public:
    static_assert(
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
     * Set up polynomial coarsening between the DoFHandler objects underlying
     * two MatrixFree objects and the respective numbers for the DoFHandler
     * objects within MatrixFree. This reinit() function allows for a more
     * efficient setup of the transfer operator and reduces the overall memory
     * consumption of a multigrid cycle in case the same MatrixFree objects are
     * also used for smoothers and residual evaluation on the two involved
     * levels.
     */
    void
    reinit_polynomial_transfer(
      const MatrixFree<dim, Number> &matrix_free_fine,
      const unsigned int             dof_handler_index_fine,
      const MatrixFree<dim, Number> &matrix_free_coarse,
      const unsigned int             dof_handler_index_coarse);

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
      const std::shared_ptr<const Utilities::MPI::Partitioner>
        &partitioner_fine) override;

    /**
     * Return the memory consumption of the allocated memory in this class.
     */
    std::size_t
    memory_consumption() const override;

    /**
     * Return the DoFHandler associated with the fine level and the level index.
     */
    std::pair<const DoFHandler<dim> *, unsigned int>
    get_dof_handler_fine() const override;

    /**
     * A multigrid transfer scheme. A multrigrid transfer class can have
     * different transfer schemes to enable p-adaptivity (one transfer scheme
     * per polynomial degree pair) and to enable global coarsening (one transfer
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
      Kokkos::View<Number *, MemorySpace::Default::kokkos_space>
        prolongation_matrix;

      /**
       * DoF indices coarse.
       */
      Kokkos::View<unsigned int **, MemorySpace::Default::kokkos_space>
        dof_indices_coarse;

      /**
       * DoF indices fine.
       */
      Kokkos::View<unsigned int **, MemorySpace::Default::kokkos_space>
        dof_indices_fine;

      /**
       * ShapeInfo description of the coarse cell. Needed during the
       * fast application of hanging-node constraints.
       */
      dealii::internal::MatrixFreeFunctions::ShapeInfo<double>
        shape_info_coarse;

      /**
       * Weights for continuous elements.
       */
      Kokkos::View<Number **, MemorySpace::Default::kokkos_space> weights;
    };


    /**
     * Struct used for copying the data to the execution kernel handled by
     * internal::ApplyCellKernel. This is similar to MatrixFree<dim,
     * Number>::Data.
     */
    struct TransferCellData
    {
      using TeamHandle = Kokkos::TeamPolicy<
        MemorySpace::Default::kokkos_space::execution_space>::member_type;

      using SharedViewValues =
        Kokkos::View<Number *,
                     MemorySpace::Default::kokkos_space::execution_space::
                       scratch_memory_space,
                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

      TeamHandle team_member;

      const int cell_index;

      const MGTransferScheme &transfer_data;

      const SharedViewValues &prolongation_matrix_device;

      /**
       * Memory for coarse dof values.
       */
      SharedViewValues &values_coarse;

      /**
       * Memory for fine dof values.
       */
      SharedViewValues &values_fine;

      /**
       * Memory for temporary arrays required by kernel evaluation.
       */
      SharedViewValues &scratch_pad;
    };



  protected:
    void
    prolongate_and_add_internal(VectorType       &dst,
                                const VectorType &src) const override;

    void
    restrict_and_add_internal(VectorType       &dst,
                              const VectorType &src) const override;

  private:
    /**
     * Transfer schemes.
     */
    std::vector<MGTransferScheme> schemes;

    /**
     * Helper class for reading from and writing to global coarse vectors and
     * for applying constraints.
     */
    dealii::internal::MatrixFreeFunctions::
      ConstraintInfo<dim, VectorizedArray<Number, 1>, types::global_dof_index>
        constraint_info_coarse;

    /**
     * Helper class for reading from and writing to global fine vectors.
     */
    dealii::internal::MatrixFreeFunctions::
      ConstraintInfo<dim, VectorizedArray<Number, 1>, types::global_dof_index>
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
    };

    /**
     * In case this class is built with MatrixFree objects (see the respective
     * reinit() function), we set up this data structure and skip the other
     * fields of the class.
     */
    std::unique_ptr<MatrixFreeRelatedData> matrix_free_data;

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
  };

} // namespace Portable

DEAL_II_NAMESPACE_CLOSE

#endif
