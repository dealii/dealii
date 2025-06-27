// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_matrix_free_tools_h
#define dealii_matrix_free_tools_h

#include <deal.II/base/config.h>

#include <deal.II/grid/tria.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/portable_fe_evaluation.h>
#include <deal.II/matrix_free/portable_matrix_free.h>
#include <deal.II/matrix_free/vector_access_internal.h>

#include <Kokkos_Core.hpp>


DEAL_II_NAMESPACE_OPEN

/**
 * A namespace for utility functions in the context of matrix-free operator
 * evaluation.
 */
namespace MatrixFreeTools
{
  namespace internal
  {
    template <int dim, typename Number, bool is_face_>
    class ComputeMatrixScratchData
    {
    public:
      using FEEvalType = FEEvaluationData<dim, Number, is_face_>;

      std::vector<unsigned int> dof_numbers;
      std::vector<unsigned int> quad_numbers;
      std::vector<unsigned int> n_components;
      std::vector<unsigned int> first_selected_components;
      std::vector<unsigned int> batch_type;
      static const bool         is_face = is_face_;

      std::function<std::vector<std::unique_ptr<FEEvalType>>(
        const std::pair<unsigned int, unsigned int> &)>
        op_create;
      std::function<void(std::vector<std::unique_ptr<FEEvalType>> &,
                         const unsigned int)>
        op_reinit;
      std::function<void(std::vector<std::unique_ptr<FEEvalType>> &)>
        op_compute;
    };
  } // namespace internal

  /**
   * Modify @p additional_data so that cells are categorized
   * according to their boundary IDs, making face integrals in the case of
   * cell-centric loop simpler.
   */
  template <int dim, typename AdditionalData>
  void
  categorize_by_boundary_ids(const Triangulation<dim> &tria,
                             AdditionalData           &additional_data);



  /**
   * Compute the diagonal of a linear operator (@p diagonal_global), given
   * @p matrix_free and the local cell integral operation @p cell_operation. The
   * vector is initialized to the right size in the function.
   *
   * The parameters @p dof_no, @p quad_no, and @p first_selected_component are
   * passed to the constructor of the FEEvaluation that is internally set up.
   *
   * The parameter @p first_vector_component is used to select the right
   * starting block in a block vector.
   */
  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components,
            typename Number,
            typename VectorizedArrayType,
            typename VectorType>
  void
  compute_diagonal(
    const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free,
    VectorType                                         &diagonal_global,
    const std::function<void(FEEvaluation<dim,
                                          fe_degree,
                                          n_q_points_1d,
                                          n_components,
                                          Number,
                                          VectorizedArrayType> &)>
                      &cell_operation,
    const unsigned int dof_no                   = 0,
    const unsigned int quad_no                  = 0,
    const unsigned int first_selected_component = 0,
    const unsigned int first_vector_component   = 0);

  /**
   * Same as above but for Portable::MatrixFree.
   */
  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components,
            typename Number,
            typename MemorySpace,
            typename QuadOperation>
  void
  compute_diagonal(
    const Portable::MatrixFree<dim, Number>                 &matrix_free,
    LinearAlgebra::distributed::Vector<Number, MemorySpace> &diagonal_global,
    const QuadOperation                                     &quad_operation,
    EvaluationFlags::EvaluationFlags                         evaluation_flags,
    EvaluationFlags::EvaluationFlags                         integration_flags,
    const unsigned int                                       dof_no  = 0,
    const unsigned int                                       quad_no = 0,
    const unsigned int first_selected_component                      = 0,
    const unsigned int first_vector_component                        = 0);

  /**
   * Same as above but with a class and a function pointer.
   */
  template <typename CLASS,
            int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components,
            typename Number,
            typename VectorizedArrayType,
            typename VectorType>
  void
  compute_diagonal(
    const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free,
    VectorType                                         &diagonal_global,
    void (CLASS::*cell_operation)(FEEvaluation<dim,
                                               fe_degree,
                                               n_q_points_1d,
                                               n_components,
                                               Number,
                                               VectorizedArrayType> &) const,
    const CLASS       *owning_class,
    const unsigned int dof_no                   = 0,
    const unsigned int quad_no                  = 0,
    const unsigned int first_selected_component = 0,
    const unsigned int first_vector_component   = 0);



  /**
   * Compute the diagonal of a linear operator (@p diagonal_global), given
   * @p matrix_free and the local cell integral operation @p cell_operation,
   * interior face integral operation @p face_operation, and boundary face
   * integral operation @p boundary_operation. The
   * vector is initialized to the right size in the function.
   *
   * The parameters @p dof_no, @p quad_no, and @p first_selected_component are
   * passed to the constructor of the FEEvaluation that is internally set up.
   */
  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components,
            typename Number,
            typename VectorizedArrayType,
            typename VectorType>
  void
  compute_diagonal(
    const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free,
    VectorType                                         &diagonal_global,
    const std::function<void(FEEvaluation<dim,
                                          fe_degree,
                                          n_q_points_1d,
                                          n_components,
                                          Number,
                                          VectorizedArrayType> &)>
      &cell_operation,
    const std::function<void(FEFaceEvaluation<dim,
                                              fe_degree,
                                              n_q_points_1d,
                                              n_components,
                                              Number,
                                              VectorizedArrayType> &,
                             FEFaceEvaluation<dim,
                                              fe_degree,
                                              n_q_points_1d,
                                              n_components,
                                              Number,
                                              VectorizedArrayType> &)>
      &face_operation,
    const std::function<void(FEFaceEvaluation<dim,
                                              fe_degree,
                                              n_q_points_1d,
                                              n_components,
                                              Number,
                                              VectorizedArrayType> &)>
                      &boundary_operation,
    const unsigned int dof_no                   = 0,
    const unsigned int quad_no                  = 0,
    const unsigned int first_selected_component = 0,
    const unsigned int first_vector_component   = 0);



  /**
   * Same as above but with a class and a function pointer.
   */
  template <typename CLASS,
            int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components,
            typename Number,
            typename VectorizedArrayType,
            typename VectorType>
  void
  compute_diagonal(
    const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free,
    VectorType                                         &diagonal_global,
    void (CLASS::*cell_operation)(FEEvaluation<dim,
                                               fe_degree,
                                               n_q_points_1d,
                                               n_components,
                                               Number,
                                               VectorizedArrayType> &) const,
    void (CLASS::*face_operation)(FEFaceEvaluation<dim,
                                                   fe_degree,
                                                   n_q_points_1d,
                                                   n_components,
                                                   Number,
                                                   VectorizedArrayType> &,
                                  FEFaceEvaluation<dim,
                                                   fe_degree,
                                                   n_q_points_1d,
                                                   n_components,
                                                   Number,
                                                   VectorizedArrayType> &)
      const,
    void (CLASS::*boundary_operation)(FEFaceEvaluation<dim,
                                                       fe_degree,
                                                       n_q_points_1d,
                                                       n_components,
                                                       Number,
                                                       VectorizedArrayType> &)
      const,
    const CLASS       *owning_class,
    const unsigned int dof_no                   = 0,
    const unsigned int quad_no                  = 0,
    const unsigned int first_selected_component = 0,
    const unsigned int first_vector_component   = 0);



  /**
   * Compute the matrix representation of a linear operator (@p matrix), given
   * @p matrix_free and the local cell integral operation @p cell_operation.
   * Constrained entries on the diagonal are set to one.
   *
   * The parameters @p dof_no, @p quad_no, and @p first_selected_component are
   * passed to the constructor of the FEEvaluation that is internally set up.
   */
  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components,
            typename Number,
            typename VectorizedArrayType,
            typename MatrixType>
  void
  compute_matrix(
    const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free,
    const AffineConstraints<Number>                    &constraints,
    MatrixType                                         &matrix,
    const std::function<void(FEEvaluation<dim,
                                          fe_degree,
                                          n_q_points_1d,
                                          n_components,
                                          Number,
                                          VectorizedArrayType> &)>
                      &cell_operation,
    const unsigned int dof_no                   = 0,
    const unsigned int quad_no                  = 0,
    const unsigned int first_selected_component = 0);



  /**
   * Same as above but with a class and a function pointer.
   */
  template <typename CLASS,
            int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components,
            typename Number,
            typename VectorizedArrayType,
            typename MatrixType>
  void
  compute_matrix(
    const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free,
    const AffineConstraints<Number>                    &constraints,
    MatrixType                                         &matrix,
    void (CLASS::*cell_operation)(FEEvaluation<dim,
                                               fe_degree,
                                               n_q_points_1d,
                                               n_components,
                                               Number,
                                               VectorizedArrayType> &) const,
    const CLASS       *owning_class,
    const unsigned int dof_no                   = 0,
    const unsigned int quad_no                  = 0,
    const unsigned int first_selected_component = 0);


  namespace internal
  {
    /**
     * Compute the diagonal of a linear operator (@p diagonal_global), given
     * @p matrix_free and the local cell, face and boundary integral operation.
     */
    template <int dim,
              typename Number,
              typename VectorizedArrayType,
              typename VectorType,
              typename VectorType2>
    void
    compute_diagonal(
      const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free,
      const internal::ComputeMatrixScratchData<dim, VectorizedArrayType, false>
        &data_cell,
      const internal::ComputeMatrixScratchData<dim, VectorizedArrayType, true>
        &data_face,
      const internal::ComputeMatrixScratchData<dim, VectorizedArrayType, true>
                                 &data_boundary,
      VectorType                 &diagonal_global,
      std::vector<VectorType2 *> &diagonal_global_components);

    /**
     * Compute the matrix representation of a linear operator (@p matrix), given
     * @p matrix_free and the local cell, face and boundary integral operation.
     */
    template <int dim,
              typename Number,
              typename VectorizedArrayType,
              typename MatrixType>
    void
    compute_matrix(
      const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free,
      const AffineConstraints<Number>                    &constraints,
      const internal::ComputeMatrixScratchData<dim, VectorizedArrayType, false>
        &cell_operation,
      const internal::ComputeMatrixScratchData<dim, VectorizedArrayType, true>
        &face_operation,
      const internal::ComputeMatrixScratchData<dim, VectorizedArrayType, true>
                 &boundary_operation,
      MatrixType &matrix);
  } // namespace internal



  /**
   * Compute the matrix representation of a linear operator (@p matrix), given
   * @p matrix_free and the local cell integral operation @p cell_operation,
   * interior face integral operation @p face_operation, and boundary face
   * integal operation @p boundary_operation.
   *
   * The parameters @p dof_no, @p quad_no, and @p first_selected_component are
   * passed to the constructor of the FEEvaluation that is internally set up.
   */
  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components,
            typename Number,
            typename VectorizedArrayType,
            typename MatrixType>
  void
  compute_matrix(
    const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free,
    const AffineConstraints<Number>                    &constraints,
    MatrixType                                         &matrix,
    const std::function<void(FEEvaluation<dim,
                                          fe_degree,
                                          n_q_points_1d,
                                          n_components,
                                          Number,
                                          VectorizedArrayType> &)>
      &cell_operation,
    const std::function<void(FEFaceEvaluation<dim,
                                              fe_degree,
                                              n_q_points_1d,
                                              n_components,
                                              Number,
                                              VectorizedArrayType> &,
                             FEFaceEvaluation<dim,
                                              fe_degree,
                                              n_q_points_1d,
                                              n_components,
                                              Number,
                                              VectorizedArrayType> &)>
      &face_operation,
    const std::function<void(FEFaceEvaluation<dim,
                                              fe_degree,
                                              n_q_points_1d,
                                              n_components,
                                              Number,
                                              VectorizedArrayType> &)>
                      &boundary_operation,
    const unsigned int dof_no                   = 0,
    const unsigned int quad_no                  = 0,
    const unsigned int first_selected_component = 0);



  /**
   * Same as above but with a class and a function pointer.
   */
  template <typename CLASS,
            int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components,
            typename Number,
            typename VectorizedArrayType,
            typename MatrixType>
  void
  compute_matrix(
    const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free,
    const AffineConstraints<Number>                    &constraints,
    MatrixType                                         &matrix,
    void (CLASS::*cell_operation)(FEEvaluation<dim,
                                               fe_degree,
                                               n_q_points_1d,
                                               n_components,
                                               Number,
                                               VectorizedArrayType> &) const,
    void (CLASS::*face_operation)(FEFaceEvaluation<dim,
                                                   fe_degree,
                                                   n_q_points_1d,
                                                   n_components,
                                                   Number,
                                                   VectorizedArrayType> &,
                                  FEFaceEvaluation<dim,
                                                   fe_degree,
                                                   n_q_points_1d,
                                                   n_components,
                                                   Number,
                                                   VectorizedArrayType> &)
      const,
    void (CLASS::*boundary_operation)(FEFaceEvaluation<dim,
                                                       fe_degree,
                                                       n_q_points_1d,
                                                       n_components,
                                                       Number,
                                                       VectorizedArrayType> &)
      const,
    const CLASS       *owning_class,
    const unsigned int dof_no                   = 0,
    const unsigned int quad_no                  = 0,
    const unsigned int first_selected_component = 0);



  /**
   * A wrapper around MatrixFree to help users to deal with DoFHandler
   * objects involving cells without degrees of freedom, i.e.,
   * cells using FE_Nothing as element type.  In the following we call such
   * cells deactivated. All other cells are activated. In contrast to
   * MatrixFree, this class skips deactivated cells and faces between activated
   * and deactivated cells are treated as boundary faces.
   */
  template <int dim,
            typename Number,
            typename VectorizedArrayType = VectorizedArray<Number>>
  class ElementActivationAndDeactivationMatrixFree
  {
  public:
    /**
     * Struct that helps to configure
     * ElementActivationAndDeactivationMatrixFree.
     */
    struct AdditionalData
    {
      /**
       * Constructor.
       */
      AdditionalData(const unsigned int dof_index = 0)
        : dof_index(dof_index)
      {}

      /**
       * Index of the DoFHandler within MatrixFree to be used.
       */
      unsigned int dof_index;
    };

    /**
     * Reinitialize class based on a given MatrixFree instance.
     * Particularly, the index of the valid FE is determined.
     *
     * @note At the moment, only DoFHandler objects are accepted
     * that are initialized with FECollection objects with at most
     * two finite elements (i.e., `FE_Nothing` and another finite
     * element).
     */
    void
    reinit(const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free,
           const AdditionalData &additional_data = AdditionalData())
    {
      this->matrix_free = &matrix_free;

      std::vector<unsigned int> valid_fe_indices;

      const auto &fe_collection =
        matrix_free.get_dof_handler(additional_data.dof_index)
          .get_fe_collection();

      for (unsigned int i = 0; i < fe_collection.size(); ++i)
        if (fe_collection[i].n_dofs_per_cell() > 0)
          valid_fe_indices.push_back(i);

      // TODO: relax this so that arbitrary number of valid
      // FEs can be accepted
      AssertDimension(valid_fe_indices.size(), 1);

      fe_index_valid = *valid_fe_indices.begin();
    }

    /**
     * Loop over all activated cells.
     *
     * For the meaning of the parameters see MatrixFree::cell_loop().
     */
    template <typename VectorTypeOut, typename VectorTypeIn>
    void
    cell_loop(const std::function<void(
                const MatrixFree<dim, Number, VectorizedArrayType> &,
                VectorTypeOut &,
                const VectorTypeIn &,
                const std::pair<unsigned int, unsigned int> &)> &cell_operation,
              VectorTypeOut                                     &dst,
              const VectorTypeIn                                &src,
              const bool zero_dst_vector = false) const
    {
      const auto ebd_cell_operation = [&](const auto &matrix_free,
                                          auto       &dst,
                                          const auto &src,
                                          const auto &range) {
        const auto category = matrix_free.get_cell_range_category(range);

        if (category != fe_index_valid)
          return;

        cell_operation(matrix_free, dst, src, range);
      };

      matrix_free->template cell_loop<VectorTypeOut, VectorTypeIn>(
        ebd_cell_operation, dst, src, zero_dst_vector);
    }

    /**
     * Loop over all activated cells and faces. Faces between activated and
     * deactivated cells are treated as boundary faces with the boundary ID
     * numbers::internal_face_boundary_id.
     *
     * For the meaning of the parameters see MatrixFree::cell_loop().
     */
    template <typename VectorTypeOut, typename VectorTypeIn>
    void
    loop(const std::function<
           void(const MatrixFree<dim, Number, VectorizedArrayType> &,
                VectorTypeOut &,
                const VectorTypeIn &,
                const std::pair<unsigned int, unsigned int> &)> &cell_operation,
         const std::function<
           void(const MatrixFree<dim, Number, VectorizedArrayType> &,
                VectorTypeOut &,
                const VectorTypeIn &,
                const std::pair<unsigned int, unsigned int> &)> &face_operation,
         const std::function<
           void(const MatrixFree<dim, Number, VectorizedArrayType> &,
                VectorTypeOut &,
                const VectorTypeIn &,
                const std::pair<unsigned int, unsigned int> &,
                const bool)> &boundary_operation,
         VectorTypeOut       &dst,
         const VectorTypeIn  &src,
         const bool           zero_dst_vector = false) const
    {
      const auto ebd_cell_operation = [&](const auto &matrix_free,
                                          auto       &dst,
                                          const auto &src,
                                          const auto &range) {
        const auto category = matrix_free.get_cell_range_category(range);

        if (category != fe_index_valid)
          return;

        cell_operation(matrix_free, dst, src, range);
      };

      const auto ebd_internal_or_boundary_face_operation =
        [&](const auto &matrix_free,
            auto       &dst,
            const auto &src,
            const auto &range) {
          const auto category = matrix_free.get_face_range_category(range);

          const unsigned int type =
            static_cast<unsigned int>(category.first == fe_index_valid) +
            static_cast<unsigned int>(category.second == fe_index_valid);

          if (type == 0)
            return; // deactivated face -> nothing to do

          if (type == 1) // boundary face
            boundary_operation(
              matrix_free, dst, src, range, category.first == fe_index_valid);
          else if (type == 2) // internal face
            face_operation(matrix_free, dst, src, range);
        };

      matrix_free->template loop<VectorTypeOut, VectorTypeIn>(
        ebd_cell_operation,
        ebd_internal_or_boundary_face_operation,
        ebd_internal_or_boundary_face_operation,
        dst,
        src,
        zero_dst_vector);
    }

  private:
    /**
     * Reference to the underlying MatrixFree object.
     */
    ObserverPointer<const MatrixFree<dim, Number, VectorizedArrayType>>
      matrix_free;

    /**
     * Index of the valid FE. Currently, only a single one is supported.
     */
    unsigned int fe_index_valid;
  };

  // implementations

#ifndef DOXYGEN

  template <int dim, typename AdditionalData>
  void
  categorize_by_boundary_ids(const Triangulation<dim> &tria,
                             AdditionalData           &additional_data)
  {
    // ... determine if we are on an active or a multigrid level
    const unsigned int level = additional_data.mg_level;
    const bool         is_mg = (level != numbers::invalid_unsigned_int);

    // ... create empty list for the category of each cell
    if (is_mg)
      additional_data.cell_vectorization_category.assign(
        std::distance(tria.begin(level), tria.end(level)), 0);
    else
      additional_data.cell_vectorization_category.assign(tria.n_active_cells(),
                                                         0);

    // ... set up scaling factor
    std::vector<unsigned int> factors(GeometryInfo<dim>::faces_per_cell);

    auto bids = tria.get_boundary_ids();
    std::sort(bids.begin(), bids.end());

    {
      unsigned int n_bids = bids.size() + 1;
      int          offset = 1;
      for (unsigned int i = 0; i < GeometryInfo<dim>::faces_per_cell;
           i++, offset = offset * n_bids)
        factors[i] = offset;
    }

    const auto to_category = [&](const auto &cell) {
      unsigned int c_num = 0;
      for (unsigned int i = 0; i < GeometryInfo<dim>::faces_per_cell; ++i)
        {
          const auto face = cell->face(i);
          if (face->at_boundary() && !cell->has_periodic_neighbor(i))
            c_num +=
              factors[i] * (1 + std::distance(bids.begin(),
                                              std::find(bids.begin(),
                                                        bids.end(),
                                                        face->boundary_id())));
        }
      return c_num;
    };

    if (!is_mg)
      {
        for (auto cell = tria.begin_active(); cell != tria.end(); ++cell)
          {
            if (cell->is_locally_owned())
              additional_data
                .cell_vectorization_category[cell->active_cell_index()] =
                to_category(cell);
          }
      }
    else
      {
        for (auto cell = tria.begin(level); cell != tria.end(level); ++cell)
          {
            if (cell->is_locally_owned_on_level())
              additional_data.cell_vectorization_category[cell->index()] =
                to_category(cell);
          }
      }

    // ... finalize set up of matrix_free
    additional_data.hold_all_faces_to_owned_cells        = true;
    additional_data.cell_vectorization_categories_strict = true;
    additional_data.mapping_update_flags_faces_by_cells =
      additional_data.mapping_update_flags_inner_faces |
      additional_data.mapping_update_flags_boundary_faces;
  }

  namespace internal
  {
    template <typename Number>
    struct LocalCSR
    {
      std::vector<unsigned int> row_lid_to_gid;
      std::vector<unsigned int> row;
      std::vector<unsigned int> col;
      std::vector<Number>       val;

      std::vector<unsigned int>                          inverse_lookup_rows;
      std::vector<std::pair<unsigned int, unsigned int>> inverse_lookup_origins;
    };

    template <int dim, typename VectorizedArrayType, bool is_face>
    class ComputeDiagonalHelper
    {
    public:
      using FEEvaluationType =
        FEEvaluationData<dim, VectorizedArrayType, is_face>;

      using Number = typename VectorizedArrayType::value_type;
      static const unsigned int n_lanes = VectorizedArrayType::size();

      ComputeDiagonalHelper()
        : phi(nullptr)
        , matrix_free(nullptr)
        , dofs_per_component(0)
        , n_components(0)
      {}

      ComputeDiagonalHelper(const ComputeDiagonalHelper &)
        : phi(nullptr)
        , matrix_free(nullptr)
        , dofs_per_component(0)
        , n_components(0)
      {}

      void
      initialize(
        FEEvaluationType                                   &phi,
        const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free,
        const unsigned int                                  n_components)
      {
        // if we are in hp mode and the number of unknowns changed, we must
        // clear the map of entries
        if (dofs_per_component !=
            phi.get_shape_info().dofs_per_component_on_cell)
          {
            locally_relevant_constraints_hn_map.clear();
            dofs_per_component =
              phi.get_shape_info().dofs_per_component_on_cell;
          }
        this->n_components  = n_components;
        this->dofs_per_cell = n_components * dofs_per_component;
        this->phi           = &phi;
        this->matrix_free   = &matrix_free;
      }

      void
      reinit(const unsigned int cell)
      {
        // STEP 1: get relevant information from FEEvaluation
        const auto        &matrix_free     = *this->matrix_free;
        const auto        &dof_info        = phi->get_dof_info();
        const unsigned int n_fe_components = dof_info.start_components.back();

        // if we have a block vector with components with the same DoFHandler,
        // each component is described with same set of constraints and
        // we consider the shift in components only during access of the global
        // vector
        const unsigned int first_selected_component =
          n_fe_components == 1 ? 0 : phi->get_first_selected_component();

        this->n_lanes_filled =
          is_face ? matrix_free.n_active_entries_per_face_batch(cell) :
                    matrix_free.n_active_entries_per_cell_batch(cell);

        // STEP 2: setup CSR storage of transposed locally-relevant
        //   constraint matrix

        // (constrained local index, global index of dof
        // constraints, weight)
        std::vector<std::tuple<unsigned int, unsigned int, Number>>
          locally_relevant_constraints, locally_relevant_constraints_tmp;
        locally_relevant_constraints.reserve(dofs_per_cell);
        std::vector<unsigned int>  constraint_position;
        std::vector<unsigned char> is_constrained_hn;

        const std::array<unsigned int, n_lanes> &cells =
          this->phi->get_cell_ids();

        for (unsigned int v = 0; v < n_lanes_filled; ++v)
          {
            Assert(cells[v] != numbers::invalid_unsigned_int,
                   ExcInternalError());

            const unsigned int *dof_indices;
            unsigned int        index_indicators, next_index_indicators;

            const unsigned int start =
              cells[v] * n_fe_components + first_selected_component;
            dof_indices =
              dof_info.dof_indices.data() + dof_info.row_starts[start].first;
            index_indicators      = dof_info.row_starts[start].second;
            next_index_indicators = dof_info.row_starts[start + 1].second;

            // STEP 2a: setup locally-relevant constraint matrix in a
            //   coordinate list (COO)
            locally_relevant_constraints.clear();

            if (n_components == 1 || n_fe_components == 1)
              {
                unsigned int ind_local = 0;
                for (; index_indicators != next_index_indicators;
                     ++index_indicators, ++ind_local)
                  {
                    const std::pair<unsigned short, unsigned short> indicator =
                      dof_info.constraint_indicator[index_indicators];

                    for (unsigned int j = 0; j < indicator.first;
                         ++j, ++ind_local)
                      locally_relevant_constraints.emplace_back(ind_local,
                                                                dof_indices[j],
                                                                1.0);

                    dof_indices += indicator.first;

                    const Number *data_val =
                      matrix_free.constraint_pool_begin(indicator.second);
                    const Number *end_pool =
                      matrix_free.constraint_pool_end(indicator.second);

                    for (; data_val != end_pool; ++data_val, ++dof_indices)
                      locally_relevant_constraints.emplace_back(ind_local,
                                                                *dof_indices,
                                                                *data_val);
                  }

                AssertIndexRange(ind_local, dofs_per_component + 1);

                for (; ind_local < dofs_per_component;
                     ++dof_indices, ++ind_local)
                  locally_relevant_constraints.emplace_back(ind_local,
                                                            *dof_indices,
                                                            1.0);
              }
            else
              {
                // case with vector-valued finite elements where all
                // components are included in one single vector. Assumption:
                // first come all entries to the first component, then all
                // entries to the second one, and so on. This is ensured by
                // the way MatrixFree reads out the indices.
                for (unsigned int comp = 0; comp < n_components; ++comp)
                  {
                    unsigned int ind_local = 0;

                    // check whether there is any constraint on the current
                    // cell
                    for (; index_indicators != next_index_indicators;
                         ++index_indicators, ++ind_local)
                      {
                        const std::pair<unsigned short, unsigned short>
                          indicator =
                            dof_info.constraint_indicator[index_indicators];

                        // run through values up to next constraint
                        for (unsigned int j = 0; j < indicator.first;
                             ++j, ++ind_local)
                          locally_relevant_constraints.emplace_back(
                            comp * dofs_per_component + ind_local,
                            dof_indices[j],
                            1.0);
                        dof_indices += indicator.first;

                        const Number *data_val =
                          matrix_free.constraint_pool_begin(indicator.second);
                        const Number *end_pool =
                          matrix_free.constraint_pool_end(indicator.second);

                        for (; data_val != end_pool; ++data_val, ++dof_indices)
                          locally_relevant_constraints.emplace_back(
                            comp * dofs_per_component + ind_local,
                            *dof_indices,
                            *data_val);
                      }

                    AssertIndexRange(ind_local, dofs_per_component + 1);

                    // get the dof values past the last constraint
                    for (; ind_local < dofs_per_component;
                         ++dof_indices, ++ind_local)
                      locally_relevant_constraints.emplace_back(
                        comp * dofs_per_component + ind_local,
                        *dof_indices,
                        1.0);

                    if (comp + 1 < n_components)
                      next_index_indicators =
                        dof_info.row_starts[start + comp + 2].second;
                  }
              }

            // we only need partial sortedness for the algorithm below in that
            // all entries for a particular row must be adjacent. this is
            // ensured by the way we fill the field, but check it again
            for (unsigned int i = 1; i < locally_relevant_constraints.size();
                 ++i)
              Assert(std::get<0>(locally_relevant_constraints[i]) >=
                       std::get<0>(locally_relevant_constraints[i - 1]),
                     ExcInternalError());

            // STEP 2c: apply hanging-node constraints
            if (dof_info.hanging_node_constraint_masks.size() > 0 &&
                dof_info.hanging_node_constraint_masks_comp.size() > 0 &&
                dof_info.hanging_node_constraint_masks_comp
                  [phi->get_active_fe_index()][first_selected_component])
              {
                const auto mask =
                  dof_info.hanging_node_constraint_masks[cells[v]];

                // cell has hanging nodes
                if (mask != dealii::internal::MatrixFreeFunctions::
                              unconstrained_compressed_constraint_kind)
                  {
                    // check if hanging node internpolation matrix has been set
                    // up
                    if (locally_relevant_constraints_hn_map.find(mask) ==
                        locally_relevant_constraints_hn_map.end())
                      fill_constraint_type_into_map(mask);

                    const auto &locally_relevant_constraints_hn =
                      locally_relevant_constraints_hn_map[mask];

                    locally_relevant_constraints_tmp.clear();
                    if (locally_relevant_constraints_tmp.capacity() <
                        locally_relevant_constraints.size())
                      locally_relevant_constraints_tmp.reserve(
                        locally_relevant_constraints.size() +
                        locally_relevant_constraints_hn.size());

                    // combine with other constraints: to avoid binary
                    // searches, we first build a list of where constraints
                    // are pointing to, and then merge the two lists
                    constraint_position.assign(dofs_per_cell,
                                               numbers::invalid_unsigned_int);
                    for (auto &a : locally_relevant_constraints)
                      if (constraint_position[std::get<0>(a)] ==
                          numbers::invalid_unsigned_int)
                        constraint_position[std::get<0>(a)] =
                          std::distance(locally_relevant_constraints.data(),
                                        &a);
                    is_constrained_hn.assign(dofs_per_cell, false);
                    for (auto &hn : locally_relevant_constraints_hn)
                      is_constrained_hn[std::get<0>(hn)] = 1;

                    // not constrained from hanging nodes
                    for (const auto &a : locally_relevant_constraints)
                      if (is_constrained_hn[std::get<0>(a)] == 0)
                        locally_relevant_constraints_tmp.push_back(a);

                    // dof is constrained by hanging nodes: build transitive
                    // closure
                    for (const auto &hn : locally_relevant_constraints_hn)
                      if (constraint_position[std::get<1>(hn)] !=
                          numbers::invalid_unsigned_int)
                        {
                          AssertIndexRange(constraint_position[std::get<1>(hn)],
                                           locally_relevant_constraints.size());
                          auto other = locally_relevant_constraints.begin() +
                                       constraint_position[std::get<1>(hn)];
                          AssertDimension(std::get<0>(*other), std::get<1>(hn));

                          for (; other != locally_relevant_constraints.end() &&
                                 std::get<0>(*other) == std::get<1>(hn);
                               ++other)
                            locally_relevant_constraints_tmp.emplace_back(
                              std::get<0>(hn),
                              std::get<1>(*other),
                              std::get<2>(hn) * std::get<2>(*other));
                        }

                    std::swap(locally_relevant_constraints,
                              locally_relevant_constraints_tmp);
                  }
              }

            // STEP 2d: transpose COO
            std::sort(locally_relevant_constraints.begin(),
                      locally_relevant_constraints.end(),
                      [](const auto &a, const auto &b) {
                        if (std::get<1>(a) < std::get<1>(b))
                          return true;
                        return (std::get<1>(a) == std::get<1>(b)) &&
                               (std::get<0>(a) < std::get<0>(b));
                      });

            // STEP 2e: translate COO to CRS
            auto &c_pool = c_pools[v];
            {
              c_pool.row_lid_to_gid.clear();
              c_pool.row.clear();
              c_pool.row.push_back(0);
              c_pool.col.clear();
              c_pool.val.clear();

              if (locally_relevant_constraints.size() > 0)
                c_pool.row_lid_to_gid.emplace_back(
                  std::get<1>(locally_relevant_constraints.front()));
              for (const auto &j : locally_relevant_constraints)
                {
                  if (c_pool.row_lid_to_gid.back() != std::get<1>(j))
                    {
                      c_pool.row_lid_to_gid.push_back(std::get<1>(j));
                      c_pool.row.push_back(c_pool.val.size());
                    }

                  c_pool.col.emplace_back(std::get<0>(j));
                  c_pool.val.emplace_back(std::get<2>(j));
                }

              if (c_pool.val.size() > 0)
                c_pool.row.push_back(c_pool.val.size());

              c_pool.inverse_lookup_rows.clear();
              c_pool.inverse_lookup_rows.resize(1 + dofs_per_cell);
              for (const unsigned int i : c_pool.col)
                c_pool.inverse_lookup_rows[1 + i]++;
              // transform to offsets
              std::partial_sum(c_pool.inverse_lookup_rows.begin(),
                               c_pool.inverse_lookup_rows.end(),
                               c_pool.inverse_lookup_rows.begin());
              AssertDimension(c_pool.inverse_lookup_rows.back(),
                              c_pool.col.size());

              c_pool.inverse_lookup_origins.resize(c_pool.col.size());
              std::vector<unsigned int> inverse_lookup_count(dofs_per_cell);
              for (unsigned int row = 0; row < c_pool.row.size() - 1; ++row)
                for (unsigned int col = c_pool.row[row];
                     col < c_pool.row[row + 1];
                     ++col)
                  {
                    const unsigned int index = c_pool.col[col];
                    c_pool.inverse_lookup_origins
                      [c_pool.inverse_lookup_rows[index] +
                       inverse_lookup_count[index]] = std::make_pair(row, col);
                    ++inverse_lookup_count[index];
                  }
            }
          }

        // STEP 3: compute element matrix A_e, apply
        //   locally-relevant constraints C_e^T * A_e * C_e, and get the
        //   the diagonal entry
        //     (C_e^T * A_e * C_e)(i,i)
        //   or
        //     C_e^T(i,:) * A_e * C_e(:,i).
        //
        //   Since, we compute the element matrix column-by-column and as a
        //   result never actually have the full element matrix, we actually
        //   perform following steps:
        //    1) loop over all columns of the element matrix
        //     a) compute column i
        //     b) compute for each j (rows of C_e^T):
        //          (C_e^T(j,:) * A_e(:,i)) * C_e(i,j)
        //       or
        //          (C_e^T(j,:) * A_e(:,i)) * C_e^T(j,i)
        //       This gives a contribution the j-th entry of the
        //       locally-relevant diagonal and comprises the multiplication
        //       by the locally-relevant constraint matrix from the left and
        //       the right. There is no contribution to the j-th vector
        //       entry if the j-th row of C_e^T is empty or C_e^T(j,i) is
        //       zero.

        // set size locally-relevant diagonal
        for (unsigned int v = 0; v < n_lanes_filled; ++v)
          diagonals_local_constrained[v].assign(
            c_pools[v].row_lid_to_gid.size() *
              (n_fe_components == 1 ? n_components : 1),
            Number(0.0));

        // check if fast path can be taken via FEEvaluation
        bool use_fast_path = true;

        for (unsigned int v = 0; v < n_lanes_filled; ++v)
          {
            auto &c_pool = c_pools[v];

            for (unsigned int i = 0; i < c_pool.row.size() - 1; ++i)
              {
                if ((c_pool.row[i + 1] - c_pool.row[i]) > 1)
                  {
                    use_fast_path = false;
                    break;
                  }
                else if (((c_pool.row[i + 1] - c_pool.row[i]) == 1) &&
                         (c_pool.val[c_pool.row[i]] != 1.0))
                  {
                    use_fast_path = false;
                    break;
                  }
              }

            if (use_fast_path == false)
              break;
          }

        this->has_simple_constraints_ = use_fast_path;
      }

      void
      fill_constraint_type_into_map(
        const dealii::internal::MatrixFreeFunctions::compressed_constraint_kind
          mask)
      {
        auto &constraints_hn = locally_relevant_constraints_hn_map[mask];

        // assume that we constrain one face, i.e., (fe_degree + 1)^(dim-1)
        // unknowns - we might have more or less entries, but this is a good
        // first guess
        const unsigned int degree =
          phi->get_shape_info().data.front().fe_degree;
        constraints_hn.reserve(Utilities::pow(degree + 1, dim - 1));

        // 1) collect hanging-node constraints for cell assuming
        // scalar finite element
        values_dofs.resize(dofs_per_component);
        std::array<
          dealii::internal::MatrixFreeFunctions::compressed_constraint_kind,
          VectorizedArrayType::size()>
          constraint_mask;
        constraint_mask.fill(dealii::internal::MatrixFreeFunctions::
                               unconstrained_compressed_constraint_kind);
        constraint_mask[0] = mask;

        for (unsigned int i = 0; i < dofs_per_component; ++i)
          {
            for (unsigned int j = 0; j < dofs_per_component; ++j)
              values_dofs[j] = VectorizedArrayType();
            values_dofs[i] = Number(1);

            dealii::internal::FEEvaluationHangingNodesFactory<
              dim,
              Number,
              VectorizedArrayType>::apply(1,
                                          degree,
                                          phi->get_shape_info(),
                                          false,
                                          constraint_mask,
                                          values_dofs.data());

            const Number tolerance =
              std::max(Number(1e-12),
                       std::numeric_limits<Number>::epsilon() * 16);
            for (unsigned int j = 0; j < dofs_per_component; ++j)
              if (std::abs(values_dofs[j][0]) > tolerance &&
                  (j != i ||
                   std::abs(values_dofs[j][0] - Number(1)) > tolerance))
                constraints_hn.emplace_back(j, i, values_dofs[j][0]);
          }

        // 2) extend for multiple components
        const unsigned int n_hn_constraints = constraints_hn.size();
        constraints_hn.resize(n_hn_constraints * n_components);

        for (unsigned int c = 1; c < n_components; ++c)
          for (unsigned int i = 0; i < n_hn_constraints; ++i)
            constraints_hn[c * n_hn_constraints + i] = std::make_tuple(
              std::get<0>(constraints_hn[i]) + c * dofs_per_component,
              std::get<1>(constraints_hn[i]) + c * dofs_per_component,
              std::get<2>(constraints_hn[i]));
      }

      void
      prepare_basis_vector(const unsigned int i)
      {
        this->i = i;

        // compute i-th column of element stiffness matrix:
        // this could be simply performed as done at the moment with
        // matrix-free operator evaluation applied to a ith-basis vector
        VectorizedArrayType *dof_values = phi->begin_dof_values();
        for (unsigned int j = 0; j < dofs_per_cell; ++j)
          dof_values[j] = VectorizedArrayType();
        dof_values[i] = Number(1);
      }

      void
      zero_basis_vector()
      {
        VectorizedArrayType *dof_values = phi->begin_dof_values();
        for (unsigned int j = 0; j < dofs_per_cell; ++j)
          dof_values[j] = VectorizedArrayType();
      }

      void
      submit()
      {
        // if we have a block vector with components with the same DoFHandler,
        // we need to figure out which component and which DoF within the
        // component are we currently considering
        const unsigned int n_fe_components =
          phi->get_dof_info().start_components.back();
        const unsigned int comp =
          n_fe_components == 1 ? i / dofs_per_component : 0;
        const unsigned int i_comp =
          n_fe_components == 1 ? (i % dofs_per_component) : i;

        // apply local constraint matrix from left and from right:
        // loop over all rows of transposed constrained matrix
        for (unsigned int v = 0; v < n_lanes_filled; ++v)
          {
            const auto &c_pool = c_pools[v];

            for (unsigned int jj = c_pool.inverse_lookup_rows[i_comp];
                 jj < c_pool.inverse_lookup_rows[i_comp + 1];
                 ++jj)
              {
                const unsigned int j = c_pool.inverse_lookup_origins[jj].first;
                // apply constraint matrix from the left
                Number temp = 0.0;
                for (unsigned int k = c_pool.row[j]; k < c_pool.row[j + 1]; ++k)
                  temp += c_pool.val[k] *
                          phi->begin_dof_values()[comp * dofs_per_component +
                                                  c_pool.col[k]][v];

                // apply constraint matrix from the right
                diagonals_local_constrained
                  [v][j + comp * c_pools[v].row_lid_to_gid.size()] +=
                  temp * c_pool.val[c_pool.inverse_lookup_origins[jj].second];
              }
          }
      }

      template <typename VectorType>
      inline void
      distribute_local_to_global(std::vector<VectorType *> &diagonal_global)
      {
        // STEP 4: assembly results: add into global vector
        const unsigned int n_fe_components =
          phi->get_dof_info().start_components.back();

        if (n_fe_components == 1)
          AssertDimension(diagonal_global.size(), n_components);

        for (unsigned int v = 0; v < n_lanes_filled; ++v)
          // if we have a block vector with components with the same
          // DoFHandler, we need to loop over all components manually and
          // need to apply the correct shift
          for (unsigned int j = 0; j < c_pools[v].row.size() - 1; ++j)
            for (unsigned int comp = 0;
                 comp < (n_fe_components == 1 ?
                           static_cast<unsigned int>(n_components) :
                           1);
                 ++comp)
              ::dealii::internal::vector_access_add(
                *diagonal_global[n_fe_components == 1 ? comp : 0],
                c_pools[v].row_lid_to_gid[j],
                diagonals_local_constrained
                  [v][j + comp * c_pools[v].row_lid_to_gid.size()]);
      }

      bool
      has_simple_constraints() const
      {
        return has_simple_constraints_;
      }

    private:
      FEEvaluationType                                   *phi;
      const MatrixFree<dim, Number, VectorizedArrayType> *matrix_free;

      unsigned int dofs_per_component;
      unsigned int dofs_per_cell;
      unsigned int n_components;

      unsigned int i;

      unsigned int n_lanes_filled;

      std::array<internal::LocalCSR<Number>, n_lanes> c_pools;

      // local storage: buffer so that we access the global vector once
      // note: may be larger then dofs_per_cell in the presence of
      // constraints!
      std::array<std::vector<Number>, n_lanes> diagonals_local_constrained;

      std::map<
        dealii::internal::MatrixFreeFunctions::compressed_constraint_kind,
        std::vector<std::tuple<unsigned int, unsigned int, Number>>>
        locally_relevant_constraints_hn_map;

      // scratch array
      AlignedVector<VectorizedArrayType> values_dofs;

      bool has_simple_constraints_;
    };

    template <bool is_face,
              int  dim,
              typename Number,
              typename VectorizedArrayType>
    bool
    is_fe_nothing(
      const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free,
      const std::pair<unsigned int, unsigned int>        &range,
      const unsigned int                                  dof_no,
      const unsigned int                                  quad_no,
      const unsigned int first_selected_component,
      const unsigned int fe_degree,
      const unsigned int n_q_points_1d,
      const bool         is_interior_face = true)
    {
      const unsigned int static_n_q_points =
        is_face ? Utilities::pow(n_q_points_1d, dim - 1) :
                  Utilities::pow(n_q_points_1d, dim);

      unsigned int active_fe_index = 0;
      if (!is_face)
        active_fe_index = matrix_free.get_cell_active_fe_index(range, dof_no);
      else if (is_interior_face)
        active_fe_index =
          matrix_free.get_face_range_category(range, dof_no).first;
      else
        active_fe_index =
          matrix_free.get_face_range_category(range, dof_no).second;

      const auto init_data = dealii::internal::
        extract_initialization_data<is_face, dim, Number, VectorizedArrayType>(
          matrix_free,
          dof_no,
          first_selected_component,
          quad_no,
          fe_degree,
          static_n_q_points,
          active_fe_index,
          numbers::invalid_unsigned_int /*active_quad_index*/,
          numbers::invalid_unsigned_int /*face_type*/);

      return init_data.shape_info->dofs_per_component_on_cell == 0;
    }



    /**
     * Portable compute kernel for MatrixFreeTools::compute_diagonal().
     */
    template <int dim,
              int fe_degree,
              int n_q_points_1d,
              int n_components,
              typename Number,
              typename QuadOperation>
    class ComputeDiagonalCellAction
    {
    public:
      ComputeDiagonalCellAction(
        const QuadOperation                   &quad_operation,
        const EvaluationFlags::EvaluationFlags evaluation_flags,
        const EvaluationFlags::EvaluationFlags integration_flags)
        : m_quad_operation(quad_operation)
        , m_evaluation_flags(evaluation_flags)
        , m_integration_flags(integration_flags)
      {}

      KOKKOS_FUNCTION void
      operator()(const typename Portable::MatrixFree<dim, Number>::Data *data,
                 const Portable::DeviceVector<Number> &,
                 Portable::DeviceVector<Number> &dst) const
      {
        Portable::
          FEEvaluation<dim, fe_degree, n_q_points_1d, n_components, Number>
            fe_eval(data);
        const typename Portable::MatrixFree<dim, Number>::PrecomputedData
                 *gpu_data = data->precomputed_data;
        const int cell     = data->cell_index;

        constexpr int dofs_per_cell = decltype(fe_eval)::tensor_dofs_per_cell;
        typename decltype(fe_eval)::value_type
          diagonal[dofs_per_cell / n_components] = {};
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          {
            const auto c = i % n_components;

            Kokkos::parallel_for(
              Kokkos::TeamThreadRange(data->team_member,
                                      dofs_per_cell / n_components),
              [&](unsigned int j) {
                typename decltype(fe_eval)::value_type val = {};

                if constexpr (n_components == 1)
                  {
                    val = (i == j) ? 1 : 0;
                  }
                else
                  {
                    val[c] = (i / n_components == j) ? 1 : 0;
                  }

                fe_eval.submit_dof_value(val, j);
              });

            data->team_member.team_barrier();

            Portable::internal::
              resolve_hanging_nodes<dim, fe_degree, false, Number>(
                data->team_member,
                gpu_data->constraint_weights,
                gpu_data->constraint_mask(cell * n_components + c),
                Kokkos::subview(data->shared_data->values, Kokkos::ALL, c));

            fe_eval.evaluate(m_evaluation_flags);
            fe_eval.apply_for_each_quad_point(m_quad_operation);
            fe_eval.integrate(m_integration_flags);

            Portable::internal::
              resolve_hanging_nodes<dim, fe_degree, true, Number>(
                data->team_member,
                gpu_data->constraint_weights,
                gpu_data->constraint_mask(cell * n_components + c),
                Kokkos::subview(data->shared_data->values, Kokkos::ALL, c));

            Kokkos::single(Kokkos::PerTeam(data->team_member), [&] {
              if constexpr (n_components == 1)
                diagonal[i] = fe_eval.get_dof_value(i);
              else
                diagonal[i / n_components][i % n_components] =
                  fe_eval.get_dof_value(i / n_components)[i % n_components];
            });

            data->team_member.team_barrier();
          }

        Kokkos::single(Kokkos::PerTeam(data->team_member), [&] {
          for (unsigned int i = 0; i < dofs_per_cell / n_components; ++i)
            fe_eval.submit_dof_value(diagonal[i], i);
        });

        data->team_member.team_barrier();

        // We need to do the same as distribute_local_to_global but without
        // constraints since we have already taken care of them earlier
        if (gpu_data->use_coloring)
          {
            Kokkos::parallel_for(
              Kokkos::TeamThreadRange(data->team_member, dofs_per_cell),
              [&](const int &i) {
                dst[gpu_data->local_to_global(i, cell)] +=
                  data->shared_data->values(i % (dofs_per_cell / n_components),
                                            i / (dofs_per_cell / n_components));
              });
          }
        else
          {
            Kokkos::parallel_for(
              Kokkos::TeamThreadRange(data->team_member, dofs_per_cell),
              [&](const int &i) {
                Kokkos::atomic_add(&dst[gpu_data->local_to_global(i, cell)],
                                   data->shared_data->values(
                                     i % (dofs_per_cell / n_components),
                                     i / (dofs_per_cell / n_components)));
              });
          }
      };

      static constexpr unsigned int n_local_dofs =
        dealii::Utilities::pow(fe_degree + 1, dim) * n_components;
      static constexpr unsigned int n_q_points =
        dealii::Utilities::pow(n_q_points_1d, dim);

    private:
      const QuadOperation                    m_quad_operation;
      const EvaluationFlags::EvaluationFlags m_evaluation_flags;
      const EvaluationFlags::EvaluationFlags m_integration_flags;
    };

  } // namespace internal



  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components,
            typename Number,
            typename MemorySpace,
            typename QuadOperation>
  void
  compute_diagonal(
    const Portable::MatrixFree<dim, Number>                 &matrix_free,
    LinearAlgebra::distributed::Vector<Number, MemorySpace> &diagonal_global,
    const QuadOperation                                     &quad_operation,
    EvaluationFlags::EvaluationFlags                         evaluation_flags,
    EvaluationFlags::EvaluationFlags                         integration_flags,
    const unsigned int                                       dof_no,
    const unsigned int                                       quad_no,
    const unsigned int first_selected_component,
    const unsigned int first_vector_component)
  {
    Assert(dof_no == 0, ExcNotImplemented());
    Assert(quad_no == 0, ExcNotImplemented());
    Assert(first_selected_component == 0, ExcNotImplemented());
    Assert(first_vector_component == 0, ExcNotImplemented());

    matrix_free.initialize_dof_vector(diagonal_global);


    internal::ComputeDiagonalCellAction<dim,
                                        fe_degree,
                                        n_q_points_1d,
                                        n_components,
                                        Number,
                                        QuadOperation>
      cell_action(quad_operation, evaluation_flags, integration_flags);
    LinearAlgebra::distributed::Vector<Number, MemorySpace> dummy;
    matrix_free.cell_loop(cell_action, dummy, diagonal_global);

    matrix_free.set_constrained_values(Number(1.), diagonal_global);
  }

  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components,
            typename Number,
            typename VectorizedArrayType,
            typename VectorType>
  void
  compute_diagonal(
    const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free,
    VectorType                                         &diagonal_global,
    const std::function<void(FEEvaluation<dim,
                                          fe_degree,
                                          n_q_points_1d,
                                          n_components,
                                          Number,
                                          VectorizedArrayType> &)>
                      &cell_operation,
    const unsigned int dof_no,
    const unsigned int quad_no,
    const unsigned int first_selected_component,
    const unsigned int first_vector_component)
  {
    compute_diagonal<dim,
                     fe_degree,
                     n_q_points_1d,
                     n_components,
                     Number,
                     VectorizedArrayType,
                     VectorType>(matrix_free,
                                 diagonal_global,
                                 cell_operation,
                                 {},
                                 {},
                                 dof_no,
                                 quad_no,
                                 first_selected_component,
                                 first_vector_component);
  }

  template <typename CLASS,
            int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components,
            typename Number,
            typename VectorizedArrayType,
            typename VectorType>
  void
  compute_diagonal(
    const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free,
    VectorType                                         &diagonal_global,
    void (CLASS::*cell_operation)(FEEvaluation<dim,
                                               fe_degree,
                                               n_q_points_1d,
                                               n_components,
                                               Number,
                                               VectorizedArrayType> &) const,
    const CLASS       *owning_class,
    const unsigned int dof_no,
    const unsigned int quad_no,
    const unsigned int first_selected_component,
    const unsigned int first_vector_component)
  {
    compute_diagonal<dim,
                     fe_degree,
                     n_q_points_1d,
                     n_components,
                     Number,
                     VectorizedArrayType,
                     VectorType>(
      matrix_free,
      diagonal_global,
      [&](auto &phi) { (owning_class->*cell_operation)(phi); },
      dof_no,
      quad_no,
      first_selected_component,
      first_vector_component);
  }

  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components,
            typename Number,
            typename VectorizedArrayType,
            typename VectorType>
  void
  compute_diagonal(
    const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free,
    VectorType                                         &diagonal_global,
    const std::function<void(FEEvaluation<dim,
                                          fe_degree,
                                          n_q_points_1d,
                                          n_components,
                                          Number,
                                          VectorizedArrayType> &)>
      &cell_operation,
    const std::function<void(FEFaceEvaluation<dim,
                                              fe_degree,
                                              n_q_points_1d,
                                              n_components,
                                              Number,
                                              VectorizedArrayType> &,
                             FEFaceEvaluation<dim,
                                              fe_degree,
                                              n_q_points_1d,
                                              n_components,
                                              Number,
                                              VectorizedArrayType> &)>
      &face_operation,
    const std::function<void(FEFaceEvaluation<dim,
                                              fe_degree,
                                              n_q_points_1d,
                                              n_components,
                                              Number,
                                              VectorizedArrayType> &)>
                      &boundary_operation,
    const unsigned int dof_no,
    const unsigned int quad_no,
    const unsigned int first_selected_component,
    const unsigned int first_vector_component)
  {
    std::vector<typename dealii::internal::BlockVectorSelector<
      VectorType,
      IsBlockVector<VectorType>::value>::BaseVectorType *>
      diagonal_global_components(n_components);

    for (unsigned int d = 0; d < n_components; ++d)
      diagonal_global_components[d] = dealii::internal::
        BlockVectorSelector<VectorType, IsBlockVector<VectorType>::value>::
          get_vector_component(diagonal_global, d + first_vector_component);

    const auto &dof_info = matrix_free.get_dof_info(dof_no);

    if (dof_info.start_components.back() == 1)
      for (unsigned int comp = 0; comp < n_components; ++comp)
        {
          Assert(diagonal_global_components[comp] != nullptr,
                 ExcMessage("The finite element underlying this FEEvaluation "
                            "object is scalar, but you requested " +
                            std::to_string(n_components) +
                            " components via the template argument in "
                            "FEEvaluation. In that case, you must pass an "
                            "std::vector<VectorType> or a BlockVector to " +
                            "read_dof_values and distribute_local_to_global."));
          dealii::internal::check_vector_compatibility(
            *diagonal_global_components[comp], matrix_free, dof_info);
        }
    else
      {
        dealii::internal::check_vector_compatibility(
          *diagonal_global_components[0], matrix_free, dof_info);
      }

    using FEEvalType = FEEvaluation<dim,
                                    fe_degree,
                                    n_q_points_1d,
                                    n_components,
                                    Number,
                                    VectorizedArrayType>;

    using FEFaceEvalType = FEFaceEvaluation<dim,
                                            fe_degree,
                                            n_q_points_1d,
                                            n_components,
                                            Number,
                                            VectorizedArrayType>;

    internal::ComputeMatrixScratchData<dim, VectorizedArrayType, false>
      data_cell;

    data_cell.dof_numbers               = {dof_no};
    data_cell.quad_numbers              = {quad_no};
    data_cell.n_components              = {n_components};
    data_cell.first_selected_components = {first_selected_component};
    data_cell.batch_type                = {0};

    data_cell.op_create =
      [&](const std::pair<unsigned int, unsigned int> &range) {
        std::vector<
          std::unique_ptr<FEEvaluationData<dim, VectorizedArrayType, false>>>
          phi;

        if (!internal::is_fe_nothing<false>(matrix_free,
                                            range,
                                            dof_no,
                                            quad_no,
                                            first_selected_component,
                                            fe_degree,
                                            n_q_points_1d))
          phi.emplace_back(std::make_unique<FEEvalType>(
            matrix_free, range, dof_no, quad_no, first_selected_component));

        return phi;
      };

    data_cell.op_reinit = [](auto &phi, const unsigned batch) {
      if (phi.size() == 1)
        static_cast<FEEvalType &>(*phi[0]).reinit(batch);
    };

    if (cell_operation)
      data_cell.op_compute = [&](auto &phi) {
        cell_operation(static_cast<FEEvalType &>(*phi[0]));
      };

    internal::ComputeMatrixScratchData<dim, VectorizedArrayType, true>
      data_face;

    data_face.dof_numbers               = {dof_no, dof_no};
    data_face.quad_numbers              = {quad_no, quad_no};
    data_face.n_components              = {n_components, n_components};
    data_face.first_selected_components = {first_selected_component,
                                           first_selected_component};
    data_face.batch_type                = {1, 2};

    data_face.op_create =
      [&](const std::pair<unsigned int, unsigned int> &range) {
        std::vector<
          std::unique_ptr<FEEvaluationData<dim, VectorizedArrayType, true>>>
          phi;

        if (!internal::is_fe_nothing<true>(matrix_free,
                                           range,
                                           dof_no,
                                           quad_no,
                                           first_selected_component,
                                           fe_degree,
                                           n_q_points_1d,
                                           true) &&
            !internal::is_fe_nothing<true>(matrix_free,
                                           range,
                                           dof_no,
                                           quad_no,
                                           first_selected_component,
                                           fe_degree,
                                           n_q_points_1d,
                                           false))
          {
            phi.emplace_back(
              std::make_unique<FEFaceEvalType>(matrix_free,
                                               range,
                                               true,
                                               dof_no,
                                               quad_no,
                                               first_selected_component));
            phi.emplace_back(
              std::make_unique<FEFaceEvalType>(matrix_free,
                                               range,
                                               false,
                                               dof_no,
                                               quad_no,
                                               first_selected_component));
          }

        return phi;
      };

    data_face.op_reinit = [](auto &phi, const unsigned batch) {
      if (phi.size() == 2)
        {
          static_cast<FEFaceEvalType &>(*phi[0]).reinit(batch);
          static_cast<FEFaceEvalType &>(*phi[1]).reinit(batch);
        }
    };

    if (face_operation)
      data_face.op_compute = [&](auto &phi) {
        face_operation(static_cast<FEFaceEvalType &>(*phi[0]),
                       static_cast<FEFaceEvalType &>(*phi[1]));
      };

    internal::ComputeMatrixScratchData<dim, VectorizedArrayType, true>
      data_boundary;

    data_boundary.dof_numbers               = {dof_no};
    data_boundary.quad_numbers              = {quad_no};
    data_boundary.n_components              = {n_components};
    data_boundary.first_selected_components = {first_selected_component};
    data_boundary.batch_type                = {1};

    data_boundary
      .op_create = [&](const std::pair<unsigned int, unsigned int> &range) {
      std::vector<
        std::unique_ptr<FEEvaluationData<dim, VectorizedArrayType, true>>>
        phi;

      if (!internal::is_fe_nothing<true>(matrix_free,
                                         range,
                                         dof_no,
                                         quad_no,
                                         first_selected_component,
                                         fe_degree,
                                         n_q_points_1d,
                                         true))
        phi.emplace_back(std::make_unique<FEFaceEvalType>(
          matrix_free, range, true, dof_no, quad_no, first_selected_component));

      return phi;
    };

    data_boundary.op_reinit = [](auto &phi, const unsigned batch) {
      if (phi.size() == 1)
        static_cast<FEFaceEvalType &>(*phi[0]).reinit(batch);
    };

    if (boundary_operation)
      data_boundary.op_compute = [&](auto &phi) {
        boundary_operation(static_cast<FEFaceEvalType &>(*phi[0]));
      };

    internal::compute_diagonal(matrix_free,
                               data_cell,
                               data_face,
                               data_boundary,
                               diagonal_global,
                               diagonal_global_components);
  }

  namespace internal
  {
    template <int dim,
              typename Number,
              typename VectorizedArrayType,
              typename VectorType,
              typename VectorType2>
    void
    compute_diagonal(
      const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free,
      const internal::ComputeMatrixScratchData<dim, VectorizedArrayType, false>
        &data_cell,
      const internal::ComputeMatrixScratchData<dim, VectorizedArrayType, true>
        &data_face,
      const internal::ComputeMatrixScratchData<dim, VectorizedArrayType, true>
                                 &data_boundary,
      VectorType                 &diagonal_global,
      std::vector<VectorType2 *> &diagonal_global_components)
    {
      // TODO: can we remove diagonal_global_components as argument?

      int dummy = 0;

      using Helper =
        internal::ComputeDiagonalHelper<dim, VectorizedArrayType, false>;

      using HelperFace =
        internal::ComputeDiagonalHelper<dim, VectorizedArrayType, true>;

      Threads::ThreadLocalStorage<std::vector<Helper>> scratch_data;
      Threads::ThreadLocalStorage<std::vector<HelperFace>>
        scratch_data_internal;
      Threads::ThreadLocalStorage<std::vector<HelperFace>> scratch_data_bc;

      const auto batch_operation =
        [&](auto                                        &data,
            auto                                        &scratch_data,
            const std::pair<unsigned int, unsigned int> &range) {
          if (!data.op_compute)
            return; // nothing to do

          auto phi = data.op_create(range);

          const unsigned int n_blocks = phi.size();

          auto &helpers = scratch_data.get();
          helpers.resize(n_blocks);

          for (unsigned int b = 0; b < n_blocks; ++b)
            helpers[b].initialize(*phi[b], matrix_free, data.n_components[b]);

          for (unsigned int batch = range.first; batch < range.second; ++batch)
            {
              data.op_reinit(phi, batch);

              for (unsigned int b = 0; b < n_blocks; ++b)
                helpers[b].reinit(batch);

              if (n_blocks > 1)
                {
                  Assert(std::all_of(helpers.begin(),
                                     helpers.end(),
                                     [](const auto &helper) {
                                       return helper.has_simple_constraints();
                                     }),
                         ExcNotImplemented());
                }

              for (unsigned int b = 0; b < n_blocks; ++b)
                {
                  for (unsigned int i = 0;
                       i < phi[b]->get_shape_info().dofs_per_component_on_cell *
                             data.n_components[b];
                       ++i)
                    {
                      for (unsigned int bb = 0; bb < n_blocks; ++bb)
                        if (b == bb)
                          helpers[bb].prepare_basis_vector(i);
                        else
                          helpers[bb].zero_basis_vector();

                      data.op_compute(phi);
                      helpers[b].submit();
                    }

                  helpers[b].distribute_local_to_global(
                    diagonal_global_components);
                }
            }
        };

      const auto cell_operation_wrapped =
        [&](const auto &, auto &, const auto &, const auto range) {
          batch_operation(data_cell, scratch_data, range);
        };

      const auto face_operation_wrapped =
        [&](const auto &, auto &, const auto &, const auto range) {
          batch_operation(data_face, scratch_data_internal, range);
        };

      const auto boundary_operation_wrapped =
        [&](const auto &, auto &, const auto &, const auto range) {
          batch_operation(data_boundary, scratch_data_bc, range);
        };

      if (data_face.op_compute || data_boundary.op_compute)
        matrix_free.template loop<VectorType, int>(cell_operation_wrapped,
                                                   face_operation_wrapped,
                                                   boundary_operation_wrapped,
                                                   diagonal_global,
                                                   dummy,
                                                   false);
      else
        matrix_free.template cell_loop<VectorType, int>(cell_operation_wrapped,
                                                        diagonal_global,
                                                        dummy,
                                                        false);
    }
  } // namespace internal

  template <typename CLASS,
            int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components,
            typename Number,
            typename VectorizedArrayType,
            typename VectorType>
  void
  compute_diagonal(
    const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free,
    VectorType                                         &diagonal_global,
    void (CLASS::*cell_operation)(FEEvaluation<dim,
                                               fe_degree,
                                               n_q_points_1d,
                                               n_components,
                                               Number,
                                               VectorizedArrayType> &) const,
    void (CLASS::*face_operation)(FEFaceEvaluation<dim,
                                                   fe_degree,
                                                   n_q_points_1d,
                                                   n_components,
                                                   Number,
                                                   VectorizedArrayType> &,
                                  FEFaceEvaluation<dim,
                                                   fe_degree,
                                                   n_q_points_1d,
                                                   n_components,
                                                   Number,
                                                   VectorizedArrayType> &)
      const,
    void (CLASS::*boundary_operation)(FEFaceEvaluation<dim,
                                                       fe_degree,
                                                       n_q_points_1d,
                                                       n_components,
                                                       Number,
                                                       VectorizedArrayType> &)
      const,
    const CLASS       *owning_class,
    const unsigned int dof_no,
    const unsigned int quad_no,
    const unsigned int first_selected_component,
    const unsigned int first_vector_component)
  {
    compute_diagonal<dim,
                     fe_degree,
                     n_q_points_1d,
                     n_components,
                     Number,
                     VectorizedArrayType,
                     VectorType>(
      matrix_free,
      diagonal_global,
      [&](auto &phi) { (owning_class->*cell_operation)(phi); },
      [&](auto &phi_m, auto &phi_p) {
        (owning_class->*face_operation)(phi_m, phi_p);
      },
      [&](auto &phi) { (owning_class->*boundary_operation)(phi); },
      dof_no,
      quad_no,
      first_selected_component,
      first_vector_component);
  }

  namespace internal
  {
    /**
     * If value type of matrix and constrains equals, return a reference
     * to the given AffineConstraint instance.
     */
    template <
      typename MatrixType,
      typename Number,
      std::enable_if_t<std::is_same_v<
        std::remove_const_t<
          std::remove_reference_t<typename MatrixType::value_type>>,
        std::remove_const_t<std::remove_reference_t<Number>>>> * = nullptr>
    const AffineConstraints<typename MatrixType::value_type> &
    create_new_affine_constraints_if_needed(
      const MatrixType &,
      const AffineConstraints<Number> &constraints,
      std::unique_ptr<AffineConstraints<typename MatrixType::value_type>> &)
    {
      return constraints;
    }

    /**
     * If value type of matrix and constrains do not equal, a new
     * AffineConstraint instance with the value type of the matrix is
     * created and a reference to it is returned.
     */
    template <
      typename MatrixType,
      typename Number,
      std::enable_if_t<!std::is_same_v<
        std::remove_const_t<
          std::remove_reference_t<typename MatrixType::value_type>>,
        std::remove_const_t<std::remove_reference_t<Number>>>> * = nullptr>
    const AffineConstraints<typename MatrixType::value_type> &
    create_new_affine_constraints_if_needed(
      const MatrixType &,
      const AffineConstraints<Number> &constraints,
      std::unique_ptr<AffineConstraints<typename MatrixType::value_type>>
        &new_constraints)
    {
      new_constraints =
        std::make_unique<AffineConstraints<typename MatrixType::value_type>>();
      new_constraints->copy_from(constraints);

      return *new_constraints;
    }
  } // namespace internal

  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components,
            typename Number,
            typename VectorizedArrayType,
            typename MatrixType>
  void
  compute_matrix(
    const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free,
    const AffineConstraints<Number>                    &constraints_in,
    MatrixType                                         &matrix,
    const std::function<void(FEEvaluation<dim,
                                          fe_degree,
                                          n_q_points_1d,
                                          n_components,
                                          Number,
                                          VectorizedArrayType> &)>
                      &cell_operation,
    const unsigned int dof_no,
    const unsigned int quad_no,
    const unsigned int first_selected_component)
  {
    compute_matrix<dim,
                   fe_degree,
                   n_q_points_1d,
                   n_components,
                   Number,
                   VectorizedArrayType,
                   MatrixType>(matrix_free,
                               constraints_in,
                               matrix,
                               cell_operation,
                               {},
                               {},
                               dof_no,
                               quad_no,
                               first_selected_component);
  }

  template <typename CLASS,
            int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components,
            typename Number,
            typename VectorizedArrayType,
            typename MatrixType>
  void
  compute_matrix(
    const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free,
    const AffineConstraints<Number>                    &constraints,
    MatrixType                                         &matrix,
    void (CLASS::*cell_operation)(FEEvaluation<dim,
                                               fe_degree,
                                               n_q_points_1d,
                                               n_components,
                                               Number,
                                               VectorizedArrayType> &) const,
    const CLASS       *owning_class,
    const unsigned int dof_no,
    const unsigned int quad_no,
    const unsigned int first_selected_component)
  {
    compute_matrix<dim,
                   fe_degree,
                   n_q_points_1d,
                   n_components,
                   Number,
                   VectorizedArrayType,
                   MatrixType>(
      matrix_free,
      constraints,
      matrix,
      [&](auto &phi) { (owning_class->*cell_operation)(phi); },
      dof_no,
      quad_no,
      first_selected_component);
  }

  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components,
            typename Number,
            typename VectorizedArrayType,
            typename MatrixType>
  void
  compute_matrix(
    const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free,
    const AffineConstraints<Number>                    &constraints_in,
    MatrixType                                         &matrix,
    const std::function<void(FEEvaluation<dim,
                                          fe_degree,
                                          n_q_points_1d,
                                          n_components,
                                          Number,
                                          VectorizedArrayType> &)>
      &cell_operation,
    const std::function<void(FEFaceEvaluation<dim,
                                              fe_degree,
                                              n_q_points_1d,
                                              n_components,
                                              Number,
                                              VectorizedArrayType> &,
                             FEFaceEvaluation<dim,
                                              fe_degree,
                                              n_q_points_1d,
                                              n_components,
                                              Number,
                                              VectorizedArrayType> &)>
      &face_operation,
    const std::function<void(FEFaceEvaluation<dim,
                                              fe_degree,
                                              n_q_points_1d,
                                              n_components,
                                              Number,
                                              VectorizedArrayType> &)>
                      &boundary_operation,
    const unsigned int dof_no,
    const unsigned int quad_no,
    const unsigned int first_selected_component)
  {
    using FEEvalType = FEEvaluation<dim,
                                    fe_degree,
                                    n_q_points_1d,
                                    n_components,
                                    Number,
                                    VectorizedArrayType>;

    using FEFaceEvalType = FEFaceEvaluation<dim,
                                            fe_degree,
                                            n_q_points_1d,
                                            n_components,
                                            Number,
                                            VectorizedArrayType>;

    internal::ComputeMatrixScratchData<dim, VectorizedArrayType, false>
      data_cell;

    data_cell.dof_numbers               = {dof_no};
    data_cell.quad_numbers              = {quad_no};
    data_cell.n_components              = {n_components};
    data_cell.first_selected_components = {first_selected_component};
    data_cell.batch_type                = {0};

    data_cell.op_create =
      [&](const std::pair<unsigned int, unsigned int> &range) {
        std::vector<
          std::unique_ptr<FEEvaluationData<dim, VectorizedArrayType, false>>>
          phi;

        if (!internal::is_fe_nothing<false>(matrix_free,
                                            range,
                                            dof_no,
                                            quad_no,
                                            first_selected_component,
                                            fe_degree,
                                            n_q_points_1d))
          phi.emplace_back(std::make_unique<FEEvalType>(
            matrix_free, range, dof_no, quad_no, first_selected_component));

        return phi;
      };

    data_cell.op_reinit = [](auto &phi, const unsigned batch) {
      if (phi.size() == 1)
        static_cast<FEEvalType &>(*phi[0]).reinit(batch);
    };

    if (cell_operation)
      data_cell.op_compute = [&](auto &phi) {
        cell_operation(static_cast<FEEvalType &>(*phi[0]));
      };

    internal::ComputeMatrixScratchData<dim, VectorizedArrayType, true>
      data_face;

    data_face.dof_numbers               = {dof_no, dof_no};
    data_face.quad_numbers              = {quad_no, quad_no};
    data_face.n_components              = {n_components, n_components};
    data_face.first_selected_components = {first_selected_component,
                                           first_selected_component};
    data_face.batch_type                = {1, 2};

    data_face.op_create =
      [&](const std::pair<unsigned int, unsigned int> &range) {
        std::vector<
          std::unique_ptr<FEEvaluationData<dim, VectorizedArrayType, true>>>
          phi;

        if (!internal::is_fe_nothing<true>(matrix_free,
                                           range,
                                           dof_no,
                                           quad_no,
                                           first_selected_component,
                                           fe_degree,
                                           n_q_points_1d,
                                           true) &&
            !internal::is_fe_nothing<true>(matrix_free,
                                           range,
                                           dof_no,
                                           quad_no,
                                           first_selected_component,
                                           fe_degree,
                                           n_q_points_1d,
                                           false))
          {
            phi.emplace_back(
              std::make_unique<FEFaceEvalType>(matrix_free,
                                               range,
                                               true,
                                               dof_no,
                                               quad_no,
                                               first_selected_component));
            phi.emplace_back(
              std::make_unique<FEFaceEvalType>(matrix_free,
                                               range,
                                               false,
                                               dof_no,
                                               quad_no,
                                               first_selected_component));
          }

        return phi;
      };

    data_face.op_reinit = [](auto &phi, const unsigned batch) {
      if (phi.size() == 2)
        {
          static_cast<FEFaceEvalType &>(*phi[0]).reinit(batch);
          static_cast<FEFaceEvalType &>(*phi[1]).reinit(batch);
        }
    };

    if (face_operation)
      data_face.op_compute = [&](auto &phi) {
        face_operation(static_cast<FEFaceEvalType &>(*phi[0]),
                       static_cast<FEFaceEvalType &>(*phi[1]));
      };

    internal::ComputeMatrixScratchData<dim, VectorizedArrayType, true>
      data_boundary;

    data_boundary.dof_numbers               = {dof_no};
    data_boundary.quad_numbers              = {quad_no};
    data_boundary.n_components              = {n_components};
    data_boundary.first_selected_components = {first_selected_component};
    data_boundary.batch_type                = {1};

    data_boundary
      .op_create = [&](const std::pair<unsigned int, unsigned int> &range) {
      std::vector<
        std::unique_ptr<FEEvaluationData<dim, VectorizedArrayType, true>>>
        phi;

      if (!internal::is_fe_nothing<true>(matrix_free,
                                         range,
                                         dof_no,
                                         quad_no,
                                         first_selected_component,
                                         fe_degree,
                                         n_q_points_1d,
                                         true))
        phi.emplace_back(std::make_unique<FEFaceEvalType>(
          matrix_free, range, true, dof_no, quad_no, first_selected_component));

      return phi;
    };

    data_boundary.op_reinit = [](auto &phi, const unsigned batch) {
      if (phi.size() == 1)
        static_cast<FEFaceEvalType &>(*phi[0]).reinit(batch);
    };

    if (boundary_operation)
      data_boundary.op_compute = [&](auto &phi) {
        boundary_operation(static_cast<FEFaceEvalType &>(*phi[0]));
      };

    internal::compute_matrix(
      matrix_free, constraints_in, data_cell, data_face, data_boundary, matrix);
  }

  namespace internal
  {
    template <int dim,
              typename Number,
              typename VectorizedArrayType,
              typename MatrixType>
    void
    compute_matrix(
      const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free,
      const AffineConstraints<Number>                    &constraints_in,
      const internal::ComputeMatrixScratchData<dim, VectorizedArrayType, false>
        &data_cell,
      const internal::ComputeMatrixScratchData<dim, VectorizedArrayType, true>
        &data_face,
      const internal::ComputeMatrixScratchData<dim, VectorizedArrayType, true>
                 &data_boundary,
      MatrixType &matrix)
    {
      std::unique_ptr<AffineConstraints<typename MatrixType::value_type>>
        constraints_for_matrix;
      const AffineConstraints<typename MatrixType::value_type> &constraints =
        internal::create_new_affine_constraints_if_needed(
          matrix, constraints_in, constraints_for_matrix);

      const auto batch_operation =
        [&matrix_free, &constraints, &matrix](
          auto &data, const std::pair<unsigned int, unsigned int> &range) {
          if (!data.op_compute)
            return; // nothing to do

          auto phi = data.op_create(range);

          const unsigned int n_blocks = phi.size();

          if (n_blocks == 0)
            return;

          Table<1, unsigned int> dofs_per_cell(n_blocks);

          Table<1, std::vector<types::global_dof_index>> dof_indices(n_blocks);
          Table<2, std::vector<types::global_dof_index>> dof_indices_mf(
            n_blocks, VectorizedArrayType::size());
          Table<1, std::vector<unsigned int>> lexicographic_numbering(n_blocks);
          Table<2,
                std::array<FullMatrix<typename MatrixType::value_type>,
                           VectorizedArrayType::size()>>
            matrices(n_blocks, n_blocks);

          for (unsigned int b = 0; b < n_blocks; ++b)
            {
              const auto &fe = matrix_free.get_dof_handler(data.dof_numbers[b])
                                 .get_fe(phi[b]->get_active_fe_index());

              const auto component_base =
                matrix_free.get_dof_info(data.dof_numbers[b])
                  .component_to_base_index[data.first_selected_components[b]];
              const auto component_in_base =
                data.first_selected_components[b] -
                matrix_free.get_dof_info(data.dof_numbers[b])
                  .start_components[component_base];

              const auto &shape_info = matrix_free.get_shape_info(
                data.dof_numbers[b],
                data.quad_numbers[b],
                component_base,
                phi[b]->get_active_fe_index(),
                phi[b]->get_active_quadrature_index());

              dofs_per_cell[b] =
                shape_info.dofs_per_component_on_cell * data.n_components[b];

              dof_indices[b].resize(fe.n_dofs_per_cell());

              for (unsigned int v = 0; v < VectorizedArrayType::size(); ++v)
                dof_indices_mf[b][v].resize(dofs_per_cell[b]);

              lexicographic_numbering[b].insert(
                lexicographic_numbering[b].begin(),
                shape_info.lexicographic_numbering.begin() +
                  component_in_base * shape_info.dofs_per_component_on_cell,
                shape_info.lexicographic_numbering.begin() +
                  (component_in_base + data.n_components[b]) *
                    shape_info.dofs_per_component_on_cell);
            }

          for (unsigned int bj = 0; bj < n_blocks; ++bj)
            for (unsigned int bi = 0; bi < n_blocks; ++bi)
              std::fill_n(matrices[bi][bj].begin(),
                          VectorizedArrayType::size(),
                          FullMatrix<typename MatrixType::value_type>(
                            dofs_per_cell[bi], dofs_per_cell[bj]));

          for (auto batch = range.first; batch < range.second; ++batch)
            {
              data.op_reinit(phi, batch);

              const unsigned int n_filled_lanes =
                data.is_face ?
                  matrix_free.n_active_entries_per_face_batch(batch) :
                  matrix_free.n_active_entries_per_cell_batch(batch);

              for (unsigned int v = 0; v < n_filled_lanes; ++v)
                for (unsigned int b = 0; b < n_blocks; ++b)
                  {
                    unsigned int const cell_index =
                      (data.batch_type[b] == 0) ?
                        (batch * VectorizedArrayType::size() + v) :
                        ((data.batch_type[b] == 1) ?
                           matrix_free.get_face_info(batch).cells_interior[v] :
                           matrix_free.get_face_info(batch).cells_exterior[v]);

                    const auto cell_iterator = matrix_free.get_cell_iterator(
                      cell_index / VectorizedArrayType::size(),
                      cell_index % VectorizedArrayType::size(),
                      data.dof_numbers[b]);

                    if (matrix_free.get_mg_level() !=
                        numbers::invalid_unsigned_int)
                      cell_iterator->get_mg_dof_indices(dof_indices[b]);
                    else
                      cell_iterator->get_dof_indices(dof_indices[b]);

                    for (unsigned int j = 0; j < dofs_per_cell[b]; ++j)
                      dof_indices_mf[b][v][j] =
                        dof_indices[b][lexicographic_numbering[b][j]];
                  }

              for (unsigned int bj = 0; bj < n_blocks; ++bj)
                {
                  for (unsigned int j = 0; j < dofs_per_cell[bj]; ++j)
                    {
                      for (unsigned int bi = 0; bi < n_blocks; ++bi)
                        for (unsigned int i = 0; i < dofs_per_cell[bi]; ++i)
                          phi[bi]->begin_dof_values()[i] =
                            (bj == bi) ? static_cast<Number>(i == j) : 0.0;

                      data.op_compute(phi);

                      for (unsigned int bi = 0; bi < n_blocks; ++bi)
                        for (unsigned int i = 0; i < dofs_per_cell[bi]; ++i)
                          for (unsigned int v = 0; v < n_filled_lanes; ++v)
                            matrices[bi][bj][v](i, j) =
                              phi[bi]->begin_dof_values()[i][v];
                    }

                  for (unsigned int v = 0; v < n_filled_lanes; ++v)
                    for (unsigned int bi = 0; bi < n_blocks; ++bi)
                      if (bi == bj)
                        // specialization for blocks on the diagonal
                        // to writing into diagonal elements of the
                        // matrix if the corresponding degree of freedom
                        // is constrained, see also the documentation
                        // of AffineConstraints::distribute_local_to_global()
                        constraints.distribute_local_to_global(
                          matrices[bi][bi][v], dof_indices_mf[bi][v], matrix);
                      else
                        constraints.distribute_local_to_global(
                          matrices[bi][bj][v],
                          dof_indices_mf[bi][v],
                          dof_indices_mf[bj][v],
                          matrix);
                }
            }
        };

      const auto cell_operation_wrapped =
        [&](const auto &, auto &, const auto &, const auto range) {
          batch_operation(data_cell, range);
        };

      const auto face_operation_wrapped =
        [&](const auto &, auto &, const auto &, const auto range) {
          batch_operation(data_face, range);
        };

      const auto boundary_operation_wrapped =
        [&](const auto &, auto &, const auto &, const auto range) {
          batch_operation(data_boundary, range);
        };

      if (data_face.op_compute || data_boundary.op_compute)
        {
          matrix_free.template loop<MatrixType, MatrixType>(
            cell_operation_wrapped,
            face_operation_wrapped,
            boundary_operation_wrapped,
            matrix,
            matrix);
        }
      else
        matrix_free.template cell_loop<MatrixType, MatrixType>(
          cell_operation_wrapped, matrix, matrix);

      matrix.compress(VectorOperation::add);
    }
  } // namespace internal

  template <typename CLASS,
            int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components,
            typename Number,
            typename VectorizedArrayType,
            typename MatrixType>
  void
  compute_matrix(
    const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free,
    const AffineConstraints<Number>                    &constraints,
    MatrixType                                         &matrix,
    void (CLASS::*cell_operation)(FEEvaluation<dim,
                                               fe_degree,
                                               n_q_points_1d,
                                               n_components,
                                               Number,
                                               VectorizedArrayType> &) const,
    void (CLASS::*face_operation)(FEFaceEvaluation<dim,
                                                   fe_degree,
                                                   n_q_points_1d,
                                                   n_components,
                                                   Number,
                                                   VectorizedArrayType> &,
                                  FEFaceEvaluation<dim,
                                                   fe_degree,
                                                   n_q_points_1d,
                                                   n_components,
                                                   Number,
                                                   VectorizedArrayType> &)
      const,
    void (CLASS::*boundary_operation)(FEFaceEvaluation<dim,
                                                       fe_degree,
                                                       n_q_points_1d,
                                                       n_components,
                                                       Number,
                                                       VectorizedArrayType> &)
      const,
    const CLASS       *owning_class,
    const unsigned int dof_no,
    const unsigned int quad_no,
    const unsigned int first_selected_component)
  {
    compute_matrix<dim,
                   fe_degree,
                   n_q_points_1d,
                   n_components,
                   Number,
                   VectorizedArrayType,
                   MatrixType>(
      matrix_free,
      constraints,
      matrix,
      [&](auto &phi) { (owning_class->*cell_operation)(phi); },
      [&](auto &phi_m, auto &phi_p) {
        (owning_class->*face_operation)(phi_m, phi_p);
      },
      [&](auto &phi) { (owning_class->*boundary_operation)(phi); },
      dof_no,
      quad_no,
      first_selected_component);
  }

#endif // DOXYGEN

} // namespace MatrixFreeTools

DEAL_II_NAMESPACE_CLOSE


#endif
