/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2022 - 2024 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Part of the source code is dual licensed under Apache-2.0 WITH
 * LLVM-exception OR LGPL-2.1-or-later. Detailed license information
 * governing the source code and code contributions can be found in
 * LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
 *
 * ------------------------------------------------------------------------
 *
 * Author: Maximilian Bergbauer, Technical University of Munich, 2024
 */

// @sect3{Include files}

// The first include files have all been treated in previous examples.

#include <deal.II/base/function.h>

#include <deal.II/base/convergence_table.h>
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/timer.h>

#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/hp/fe_collection.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/vector.h>

#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/fe_evaluation.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <vector>

// This header contains some common level set functions.
// For example, the spherical geometry that we use here (the same as in
// step-85).
#include <deal.II/base/function_signed_distance.h>

// We also need 2 already used headers in step-85 from the NonMatching
// namespace.
#include <deal.II/non_matching/mesh_classifier.h>
#include <deal.II/non_matching/quadrature_generator.h>

// Compared to the matrix-based tutorial step-85 we need this new header for
// precomputation of mapping information.
#include <deal.II/non_matching/mapping_info.h>

// And most important the header for flexible matrix-free evaluation.
#include <deal.II/matrix_free/fe_point_evaluation.h>

// @sect3{Start of the program}
namespace Step95
{
  using namespace dealii;

  // Helper functions to determine if a cell (batch) is inside or intersected
  // depending on the active_fe_index.
  inline bool is_inside(unsigned int active_fe_index)
  {
    return active_fe_index ==
           (unsigned int)NonMatching::LocationToLevelSet::inside;
  }

  inline bool is_intersected(unsigned int active_fe_index)
  {
    return active_fe_index ==
           (unsigned int)NonMatching::LocationToLevelSet::intersected;
  }

  // @sect3{The PoissonOperator class Template}
  // This operator applies the matrix-free operator evaluation with the vmult()
  // function.
  template <int dim>
  class PoissonOperator
  {
    // We define 3 types frequently used throughout the tutorial keep the
    // definition of the Number type at a central place.
    using Number              = double;
    using VectorizedArrayType = VectorizedArray<Number>;
    using VectorType          = LinearAlgebra::distributed::Vector<Number>;

    // We also define types for the matrix-free evaluation classes for cells and
    // faces and for structured and unstructured quadrature to keep the
    // definition of the template arguments for those at a central place.
    using CellIntegrator =
      FEEvaluation<dim, -1, 0, 1, Number, VectorizedArrayType>;
    using FaceIntegrator =
      FEFaceEvaluation<dim, -1, 0, 1, Number, VectorizedArrayType>;
    using GenericCellIntegrator =
      FEPointEvaluation<1, dim, dim, VectorizedArrayType>;
    using GenericFaceIntegrator =
      FEFacePointEvaluation<1, dim, dim, VectorizedArrayType>;

    // And finally, we introduce a shortcut for the width of the detected SIMD
    // vectorization.
    static constexpr unsigned int n_lanes = VectorizedArrayType::size();

  public:
    // In the reinit() function the PoissonOperator receives the necessary
    // information for matrix-free evaluation, i.e. the MatrixFree object for
    // inside cells and the NonMatching::MappingInfo objects for intersected
    // cells.
    void
    reinit(const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free_in,
           const NonMatching::MappingInfo<dim, dim, VectorizedArrayType>
             *mapping_info_cell_in,
           const NonMatching::MappingInfo<dim, dim, VectorizedArrayType>
             *mapping_info_surface_in,
           const NonMatching::MappingInfo<dim, dim, VectorizedArrayType>
                             *mapping_info_faces_in,
           const bool         is_dg_in,
           const unsigned int degree)
    {
      matrix_free          = &matrix_free_in;
      mapping_info_cell    = mapping_info_cell_in;
      mapping_info_surface = mapping_info_surface_in;
      mapping_info_faces   = mapping_info_faces_in;
      is_dg                = is_dg_in;

      fe_dgq = std::make_unique<FE_DGQ<dim>>(degree);

      evaluator_cell =
        std::make_unique<GenericCellIntegrator>(*mapping_info_cell, *fe_dgq);
      evaluator_surface =
        std::make_unique<GenericCellIntegrator>(*mapping_info_surface, *fe_dgq);
      if (is_dg)
        {
          evaluator_face_m =
            std::make_unique<GenericFaceIntegrator>(*mapping_info_faces,
                                                    *fe_dgq,
                                                    true);
          evaluator_face_p =
            std::make_unique<GenericFaceIntegrator>(*mapping_info_faces,
                                                    *fe_dgq,
                                                    false);
        }

      matrix_free->initialize_cell_data_vector(cell_diameter);
      for (unsigned int cell_batch_index = 0;
           cell_batch_index <
           matrix_free->n_cell_batches() + matrix_free->n_ghost_cell_batches();
           ++cell_batch_index)
        {
          auto &diameter = cell_diameter[cell_batch_index];
          for (unsigned int v = 0;
               v <
               matrix_free->n_active_entries_per_cell_batch(cell_batch_index);
               ++v)
            {
              const auto cell_accessor_inside =
                matrix_free->get_cell_iterator(cell_batch_index, v);

              diameter[v] = cell_accessor_inside->minimum_vertex_distance();
            }
        }
    }

    // This function is the interface function for linear solvers that applies
    // the operator evaluation which is an abstraction of a matrix-vector
    // product by executing loops over cells and faces.
    void vmult(VectorType &dst, const VectorType &src) const
    {
      matrix_free->loop(&PoissonOperator::local_apply_cell,
                        &PoissonOperator::local_apply_face,
                        &PoissonOperator::local_apply_boundary_face,
                        this,
                        dst,
                        src,
                        true);
    }

    // The right-hand-side is assembled with a loop over the cells.
    void rhs(VectorType &rhs, const Function<dim> &rhs_function)
    {
      const unsigned int dummy = 0;
      matrix_free->template cell_loop<VectorType, unsigned int>(
        [&](const MatrixFree<dim, Number, VectorizedArrayType> &,
            VectorType &dst,
            const VectorType &,
            const std::pair<unsigned int, unsigned int> &cell_range) {
          CellIntegrator evaluator(*matrix_free, dof_index, quad_index);

          const auto cell_range_category =
            matrix_free->get_cell_range_category(cell_range);

          if (is_inside(cell_range_category))
            {
              for (unsigned int cell_batch_index = cell_range.first;
                   cell_batch_index < cell_range.second;
                   ++cell_batch_index)
                {
                  evaluator.reinit(cell_batch_index);
                  for (unsigned int q : evaluator.quadrature_point_indices())
                    do_rhs_cell_term(evaluator, rhs_function, q);
                  evaluator.integrate(EvaluationFlags::values);
                  evaluator.distribute_local_to_global(dst);
                }
            }
          else if (is_intersected(cell_range_category))
            {
              const auto dofs_per_cell = evaluator.dofs_per_cell;

              for (unsigned int cell_batch_index = cell_range.first;
                   cell_batch_index < cell_range.second;
                   ++cell_batch_index)
                {
                  evaluator.reinit(cell_batch_index);

                  for (unsigned int v = 0; v < n_lanes; ++v)
                    {
                      const unsigned int cell_index =
                        cell_batch_index * n_lanes + v;

                      evaluator_cell->reinit(cell_index);

                      for (const auto q :
                           evaluator_cell->quadrature_point_indices())
                        do_rhs_cell_term(*evaluator_cell, rhs_function, q);

                      evaluator_cell->integrate(
                        StridedArrayView<Number, n_lanes>(
                          &evaluator.begin_dof_values()[0][v], dofs_per_cell),
                        EvaluationFlags::values);
                    }

                  evaluator.distribute_local_to_global(dst);
                }
            }
        },
        rhs,
        dummy);
    }

    // This functions partitions a vector with matrix-free functionality.
    void initialize_dof_vector(VectorType &vec) const
    {
      matrix_free->initialize_dof_vector(vec);
    }

    // With this function we can assemble the diagonal of the operator instead
    // of applying the matrix-vector product. It uses the same function than the
    // operator evalaution for the local operation on cells and faces but sets
    // the template argument assemble to true.
    void compute_diagonal(VectorType &diagonal) const
    {
      matrix_free->loop(&PoissonOperator::local_apply_cell<true>,
                        &PoissonOperator::local_apply_face<true>,
                        &PoissonOperator::local_apply_boundary_face,
                        this,
                        diagonal,
                        diagonal);
    }

  private:
    // We define two helper functions needed for matrix-free assembly
    // operations.
    template <typename Integrator>
    void create_zero_basis(Integrator &evaluator) const
    {
      for (unsigned int i = 0; i < evaluator.dofs_per_cell; ++i)
        evaluator.begin_dof_values()[i] = VectorizedArrayType(0.);
    }

    template <typename Integrator>
    void create_standard_basis(const unsigned int j,
                               Integrator        &evaluator) const
    {
      create_zero_basis(evaluator);
      evaluator.begin_dof_values()[j] = VectorizedArrayType(1.);
    }

    // This function implements the local cell operation.
    template <bool assemble = false>
    void local_apply_cell(
      const MatrixFree<dim, Number, VectorizedArrayType> &,
      VectorType                                  &dst,
      const VectorType                            &src,
      const std::pair<unsigned int, unsigned int> &cell_range) const
    {
      CellIntegrator evaluator(*matrix_free, dof_index, quad_index);

      const auto dofs_per_cell = evaluator.dofs_per_cell;

      AlignedVector<VectorizedArrayType> dof_buffer(assemble ? dofs_per_cell :
                                                               0);

      const auto cell_range_category =
        matrix_free->get_cell_range_category(cell_range);

      // We define the cell_operation for inside cells as the standard Poisson
      // term.
      auto inside_cell_operation = [&](CellIntegrator &evaluator) {
        evaluator.evaluate(EvaluationFlags::gradients);
        for (unsigned int q : evaluator.quadrature_point_indices())
          do_poisson_cell_term(evaluator, q);
        evaluator.integrate(EvaluationFlags::gradients);
      };

      // We define the cell_operation for intersected cells as the standard
      // Poisson term in the cut volume and the Nitsche term for weak Dirichlet
      // boundary enforcement on the cut surface.
      auto intersected_cell_operation = [&](CellIntegrator &evaluator) {
        const unsigned int cell_batch_index =
          evaluator.get_cell_or_face_batch_id();
        const VectorizedArrayType diameter = cell_diameter[cell_batch_index];

        const auto tau = compute_interior_penalty_parameter(diameter);

        for (unsigned int v = 0; v < n_lanes; ++v)
          {
            const unsigned int cell_index = cell_batch_index * n_lanes + v;

            evaluator_cell->reinit(cell_index);
            evaluator_surface->reinit(cell_index);

            evaluator_cell->evaluate(StridedArrayView<const Number, n_lanes>(
                                       &evaluator.begin_dof_values()[0][v],
                                       dofs_per_cell),
                                     EvaluationFlags::gradients);

            evaluator_surface->evaluate(StridedArrayView<const Number, n_lanes>(
                                          &evaluator.begin_dof_values()[0][v],
                                          dofs_per_cell),
                                        EvaluationFlags::values |
                                          EvaluationFlags::gradients);

            for (const auto q : evaluator_cell->quadrature_point_indices())
              do_poisson_cell_term(*evaluator_cell, q);

            for (const auto q : evaluator_surface->quadrature_point_indices())
              do_boundary_flux_term_homogeneous(*evaluator_surface, tau[v], q);

            evaluator_cell->integrate(StridedArrayView<Number, n_lanes>(
                                        &evaluator.begin_dof_values()[0][v],
                                        dofs_per_cell),
                                      EvaluationFlags::gradients);

            evaluator_surface->integrate(StridedArrayView<Number, n_lanes>(
                                           &evaluator.begin_dof_values()[0][v],
                                           dofs_per_cell),
                                         EvaluationFlags::values |
                                           EvaluationFlags::gradients,
                                         true);
          }
      };

      // Depending on the active_fe_index of the cell_range we select the
      // cell_operation to execute.
      std::function<void(CellIntegrator &)> cell_operation;
      if (is_inside(cell_range_category))
        cell_operation = inside_cell_operation;
      else if (is_intersected(cell_range_category))
        cell_operation = intersected_cell_operation;
      else
        return;

      // We loop over the cell batches of the current cell_range and apply the
      // cell_operation. For a vmult() we only apply once, for diagonal assembly
      // we apply the cell_operation for every local cell DoF on a cell DoF unit
      // vector.
      for (unsigned int cell_batch_index = cell_range.first;
           cell_batch_index < cell_range.second;
           ++cell_batch_index)
        {
          evaluator.reinit(cell_batch_index);

          if (assemble)
            {
              for (unsigned int j = 0; j < evaluator.dofs_per_cell; ++j)
                {
                  create_standard_basis(j, evaluator);

                  cell_operation(evaluator);

                  dof_buffer[j] = evaluator.begin_dof_values()[j];
                }

              for (unsigned int j = 0; j < evaluator.dofs_per_cell; ++j)
                evaluator.begin_dof_values()[j] = dof_buffer[j];
            }
          else
            {
              evaluator.read_dof_values(src);

              cell_operation(evaluator);
            }

          evaluator.distribute_local_to_global(dst);
        }
    }

    // This function implements the local face operation.
    template <bool assemble = false>
    void local_apply_face(
      const MatrixFree<dim, Number, VectorizedArrayType> &,
      VectorType                                  &dst,
      const VectorType                            &src,
      const std::pair<unsigned int, unsigned int> &face_range) const
    {
      FaceIntegrator evaluator_m(*matrix_free, true, dof_index, quad_index);
      FaceIntegrator evaluator_p(*matrix_free, false, dof_index, quad_index);

      AlignedVector<VectorizedArrayType> dof_buffer_m(
        evaluator_m.dofs_per_cell);
      AlignedVector<VectorizedArrayType> dof_buffer_p(
        evaluator_p.dofs_per_cell);

      AlignedVector<VectorizedArrayType> local_diagonal_m(
        assemble ? evaluator_m.dofs_per_cell : 0);
      AlignedVector<VectorizedArrayType> local_diagonal_p(
        assemble ? evaluator_p.dofs_per_cell : 0);

      const auto face_range_category =
        matrix_free->get_face_range_category(face_range);

      // We start with the face operations for the DG case.
      // We define the face_operation for inside face as the SIPG term.
      auto inside_face_operation_dg = [&](FaceIntegrator &evaluator_m,
                                          FaceIntegrator &evaluator_p) {
        const unsigned int face_batch_index =
          evaluator_m.get_cell_or_face_batch_id();

        const auto diameter =
          compute_diameter_of_inner_face_batch(face_batch_index);

        do_local_apply_sipg_term(evaluator_m,
                                 evaluator_p,
                                 evaluator_m.begin_dof_values(),
                                 evaluator_p.begin_dof_values(),
                                 diameter);
      };

      // We define the face_operation for mixed face (between an inside and an
      // intersected cell) as the SIPG term plus the ghost penalty term.
      auto mixed_face_operation_dg = [&](FaceIntegrator &evaluator_m,
                                         FaceIntegrator &evaluator_p) {
        const unsigned int face_batch_index =
          evaluator_m.get_cell_or_face_batch_id();

        const auto diameter =
          compute_diameter_of_inner_face_batch(face_batch_index);

        do_local_apply_sipg_term(evaluator_m,
                                 evaluator_p,
                                 dof_buffer_m.data(),
                                 dof_buffer_p.data(),
                                 diameter);

        do_local_apply_gp_face_term<true>(evaluator_m,
                                          evaluator_p,
                                          dof_buffer_m.data(),
                                          dof_buffer_p.data(),
                                          diameter,
                                          true);
      };

      // We define the face_operation for intersected face (between two
      // intersected cells) as the SIPG term plus the ghost penalty term.
      auto intersected_face_operation_dg = [&](FaceIntegrator &evaluator_m,
                                               FaceIntegrator &evaluator_p) {
        const unsigned int face_batch_index =
          evaluator_m.get_cell_or_face_batch_id();

        const auto diameter =
          compute_diameter_of_inner_face_batch(face_batch_index);
        const auto tau = compute_interior_penalty_parameter(diameter);

        evaluator_m.project_to_face(dof_buffer_m.data(),
                                    EvaluationFlags::values |
                                      EvaluationFlags::gradients);
        evaluator_p.project_to_face(dof_buffer_p.data(),
                                    EvaluationFlags::values |
                                      EvaluationFlags::gradients);

        for (unsigned int v = 0; v < n_lanes; ++v)
          {
            const unsigned int face_index = face_batch_index * n_lanes + v;

            evaluator_face_m->reinit(face_index);
            evaluator_face_p->reinit(face_index);

            evaluator_face_m->evaluate_in_face(
              &evaluator_m.get_scratch_data().begin()[0][v],
              EvaluationFlags::values | EvaluationFlags::gradients);
            evaluator_face_p->evaluate_in_face(
              &evaluator_p.get_scratch_data().begin()[0][v],
              EvaluationFlags::values | EvaluationFlags::gradients);

            for (const auto q : evaluator_face_m->quadrature_point_indices())
              do_flux_term(*evaluator_face_m, *evaluator_face_p, tau[v], q);

            evaluator_face_m->integrate_in_face(
              &evaluator_m.get_scratch_data().begin()[0][v],
              EvaluationFlags::values | EvaluationFlags::gradients);
            evaluator_face_p->integrate_in_face(
              &evaluator_p.get_scratch_data().begin()[0][v],
              EvaluationFlags::values | EvaluationFlags::gradients);
          }

        evaluator_m.collect_from_face(EvaluationFlags::values |
                                      EvaluationFlags::gradients);
        evaluator_p.collect_from_face(EvaluationFlags::values |
                                      EvaluationFlags::gradients);

        do_local_apply_gp_face_term<true>(evaluator_m,
                                          evaluator_p,
                                          dof_buffer_m.data(),
                                          dof_buffer_p.data(),
                                          diameter,
                                          true);
      };

      // For the CG case we only need to define the ghost penalty term for mixed
      // and intersected faces.
      auto intersected_mixed_face_operation_cg =
        [&](FaceIntegrator &evaluator_m, FaceIntegrator &evaluator_p) {
          const unsigned int face_batch_index =
            evaluator_m.get_cell_or_face_batch_id();

          const VectorizedArrayType diameter =
            compute_diameter_of_inner_face_batch(face_batch_index);

          do_local_apply_gp_face_term<false>(evaluator_m,
                                             evaluator_p,
                                             evaluator_m.begin_dof_values(),
                                             evaluator_p.begin_dof_values(),
                                             diameter,
                                             false);
        };

      // Depending on the active_fe_index of the two cells sharing the current
      // face we select the face_operation.
      std::function<void(FaceIntegrator &, FaceIntegrator &)> face_operation;
      bool buffer_dof_values = false;
      if (is_dg)
        {
          if (is_inside_face(face_range_category))
            face_operation = inside_face_operation_dg;
          else if (is_mixed_face(face_range_category))
            {
              face_operation    = mixed_face_operation_dg;
              buffer_dof_values = true;
            }
          else if (is_intersected_face(face_range_category))
            {
              face_operation    = intersected_face_operation_dg;
              buffer_dof_values = true;
            }
          else
            return;
        }
      else
        {
          if (is_mixed_face(face_range_category) ||
              is_intersected_face(face_range_category))
            face_operation = intersected_mixed_face_operation_cg;
          else
            return;
        }

      // We loop over the face batches of the current face_range and apply the
      // face_operation. Again, for the vmult() the face_operation is executed
      // once, for diagonal assembly for every unit DoF vector of the two
      // face-sharing cells.
      for (unsigned int face_batch_index = face_range.first;
           face_batch_index < face_range.second;
           ++face_batch_index)
        {
          evaluator_m.reinit(face_batch_index);
          evaluator_p.reinit(face_batch_index);

          if (assemble)
            {
              for (unsigned int j = 0; j < evaluator_m.dofs_per_cell; ++j)
                for (unsigned int p = 0; p < 2; ++p)
                  {
                    if (p == 0)
                      {
                        create_standard_basis(j, evaluator_m);
                        create_zero_basis(evaluator_p);
                      }
                    else
                      {
                        create_zero_basis(evaluator_m);
                        create_standard_basis(j, evaluator_p);
                      }

                    if (buffer_dof_values)
                      for (unsigned int i = 0; i < evaluator_m.dofs_per_cell;
                           ++i)
                        {
                          dof_buffer_m[i] = evaluator_m.begin_dof_values()[i];
                          dof_buffer_p[i] = evaluator_p.begin_dof_values()[i];
                        }

                    face_operation(evaluator_m, evaluator_p);

                    if (p == 0)
                      local_diagonal_m[j] = evaluator_m.begin_dof_values()[j];
                    else
                      local_diagonal_p[j] = evaluator_p.begin_dof_values()[j];
                  }

              for (unsigned int j = 0; j < evaluator_m.dofs_per_cell; ++j)
                {
                  evaluator_m.begin_dof_values()[j] = local_diagonal_m[j];
                  evaluator_p.begin_dof_values()[j] = local_diagonal_p[j];
                }
            }
          else
            {
              evaluator_m.read_dof_values(src);
              evaluator_p.read_dof_values(src);

              if (buffer_dof_values)
                for (unsigned int i = 0; i < evaluator_m.dofs_per_cell; ++i)
                  {
                    dof_buffer_m[i] = evaluator_m.begin_dof_values()[i];
                    dof_buffer_p[i] = evaluator_p.begin_dof_values()[i];
                  }

              face_operation(evaluator_m, evaluator_p);
            }

          evaluator_m.distribute_local_to_global(dst);
          evaluator_p.distribute_local_to_global(dst);
        }
    }

    // We do not need a local operation of the fitted boundary for this tutorial
    // as the whole Dirichlet boundary is immersed in the volume of the domain.
    // However, this tutorial is easily extendable to the case with fitted and
    // unfitted boundaries.
    void local_apply_boundary_face(
      const MatrixFree<dim, Number, VectorizedArrayType> &,
      VectorType &,
      const VectorType &,
      const std::pair<unsigned int, unsigned int> &) const
    {}

    // This is the actual implementation of the quadrature point operation of
    // the Poisson term in the weak form. It is templated over the Integrator
    // type to be usable by FEEvaluation as well as FEPointEvaluation.
    template <typename Integrator>
    inline void do_poisson_cell_term(Integrator        &evaluator,
                                     const unsigned int q) const
    {
      evaluator.submit_gradient(evaluator.get_gradient(q), q);
    }

    // The implementation for the SIPG term (needed for DG). Again, templated
    // over the Integrator type to be usable by FEFaceEvaluation and
    // FEFacePointEvaluation.
    template <typename Integrator, typename Number2>
    inline void do_flux_term(Integrator        &evaluator_m,
                             Integrator        &evaluator_p,
                             const Number2     &tau,
                             const unsigned int q) const
    {
      const auto normal_gradient_m = evaluator_m.get_normal_derivative(q);
      const auto normal_gradient_p = evaluator_p.get_normal_derivative(q);

      const auto value_m = evaluator_m.get_value(q);
      const auto value_p = evaluator_p.get_value(q);

      const auto jump_value = value_m - value_p;

      const auto central_flux_gradient =
        0.5 * (normal_gradient_m + normal_gradient_p);

      const auto value_terms = central_flux_gradient - tau * jump_value;

      evaluator_m.submit_value(-value_terms, q);
      evaluator_p.submit_value(value_terms, q);

      const auto gradient_terms = -0.5 * jump_value;

      evaluator_m.submit_normal_derivative(gradient_terms, q);
      evaluator_p.submit_normal_derivative(gradient_terms, q);
    }

    // The implementation of the Nitsche term.
    template <typename Integrator, typename Number2>
    inline void do_boundary_flux_term_homogeneous(Integrator    &evaluator_m,
                                                  const Number2 &tau,
                                                  const unsigned int q) const
    {
      const auto value           = evaluator_m.get_value(q);
      const auto normal_gradient = evaluator_m.get_normal_derivative(q);

      const auto value_term           = 2. * tau * value - normal_gradient;
      const auto normal_gradient_term = -value;

      evaluator_m.submit_value(value_term, q);
      evaluator_m.submit_normal_derivative(normal_gradient_term, q);
    }

    // The implementation of the face-based ghost penalty term (up to degree 2 /
    // normal hessians).
    template <bool do_normal_hessians, bool do_values, typename Integrator>
    inline void do_gp_face_term(
      Integrator                            &evaluator_m,
      Integrator                            &evaluator_p,
      const typename Integrator::NumberType &masked_factor_value,
      const typename Integrator::NumberType &masked_factor_gradient,
      const typename Integrator::NumberType &masked_factor_hessian,
      const unsigned int                     q) const
    {
      if (do_values)
        {
          const auto value_m = evaluator_m.get_value(q);
          const auto value_p = evaluator_p.get_value(q);

          const auto jump_value = value_m - value_p;

          const auto value_term = masked_factor_value * jump_value;

          evaluator_m.submit_value(value_term, q);
          evaluator_p.submit_value(-value_term, q);
        }

      {
        const auto normal_gradient_m = evaluator_m.get_normal_derivative(q);
        const auto normal_gradient_p = evaluator_p.get_normal_derivative(q);

        const auto jump_normal_gradient = normal_gradient_m - normal_gradient_p;

        const auto gradient_term =
          masked_factor_gradient * jump_normal_gradient;

        evaluator_m.submit_normal_derivative(gradient_term, q);
        evaluator_p.submit_normal_derivative(-gradient_term, q);
      }

      if (do_normal_hessians)
        {
          const auto normal_hessian_m = evaluator_m.get_normal_hessian(q);
          const auto normal_hessian_p = evaluator_p.get_normal_hessian(q);

          const auto jump_normal_hessian = normal_hessian_m - normal_hessian_p;

          const auto hessian_term = masked_factor_hessian * jump_normal_hessian;

          evaluator_m.submit_normal_hessian(hessian_term, q);
          evaluator_p.submit_normal_hessian(-hessian_term, q);
        }
    }

    // The implementation of the right-hand-side term evaluating the rhs
    // function (unfortunately we cannot evaluate a Function vectorized, so we
    // have to reshuffle the quadrature point data).
    template <typename Integrator>
    inline void do_rhs_cell_term(Integrator          &evaluator,
                                 const Function<dim> &rhs_function,
                                 const unsigned int   q) const
    {
      const auto q_points = evaluator.quadrature_point(q);

      VectorizedArrayType value = 0.;

      for (unsigned int v = 0; v < n_lanes; ++v)
        {
          Point<dim> q_point;
          for (unsigned int d = 0; d < dim; ++d)
            q_point[d] = q_points[d][v];

          value[v] = rhs_function.value(q_point);
        }

      evaluator.submit_value(value, q);
    }

    // This implements the face_operation of the SIPG term (setting values in
    // integrate()).
    void do_local_apply_sipg_term(FaceIntegrator            &evaluator_m,
                                  FaceIntegrator            &evaluator_p,
                                  const VectorizedArrayType *dof_ptr_m,
                                  const VectorizedArrayType *dof_ptr_p,
                                  const VectorizedArrayType &diameter) const
    {
      const auto tau = compute_interior_penalty_parameter(diameter);

      evaluator_m.evaluate(dof_ptr_m,
                           EvaluationFlags::values |
                             EvaluationFlags::gradients);
      evaluator_p.evaluate(dof_ptr_p,
                           EvaluationFlags::values |
                             EvaluationFlags::gradients);

      for (const auto q : evaluator_m.quadrature_point_indices())
        do_flux_term(evaluator_m, evaluator_p, tau, q);

      evaluator_m.integrate(EvaluationFlags::values |
                            EvaluationFlags::gradients);
      evaluator_p.integrate(EvaluationFlags::values |
                            EvaluationFlags::gradients);
    }

    // This implements the face_operation of the ghost penalty term (potentially
    // adding into the values in integrate() depending on sum_into_values).
    template <bool is_dg_>
    void do_local_apply_gp_face_term(FaceIntegrator            &evaluator_m,
                                     FaceIntegrator            &evaluator_p,
                                     const VectorizedArrayType *dof_ptr_m,
                                     const VectorizedArrayType *dof_ptr_p,
                                     const VectorizedArrayType &diameter,
                                     const bool sum_into_values) const
    {
      EvaluationFlags::EvaluationFlags evaluation_flags =
        EvaluationFlags::gradients;
      if (is_dg_)
        evaluation_flags |= EvaluationFlags::values;

      const unsigned int degree =
        matrix_free->get_dof_handler(dof_index).get_fe().degree;
      const bool do_hessians = degree > 1;
      if (do_hessians)
        evaluation_flags |= EvaluationFlags::hessians;

      Assert(degree <= 2,
             ExcMessage(
               "Face-based stabilization only implemented up to degree 2!"));

      evaluator_m.evaluate(dof_ptr_m, evaluation_flags);
      evaluator_p.evaluate(dof_ptr_p, evaluation_flags);

      const VectorizedArrayType factor_values   = 0.5 / diameter;
      const VectorizedArrayType factor_gradient = 0.5 * diameter;
      const VectorizedArrayType factor_hessians =
        0.5 * diameter * diameter * diameter;

      if (do_hessians)
        for (const auto q : evaluator_m.quadrature_point_indices())
          do_gp_face_term<true, is_dg_>(evaluator_m,
                                        evaluator_p,
                                        factor_values,
                                        factor_gradient,
                                        factor_hessians,
                                        q);
      else
        for (const auto q : evaluator_m.quadrature_point_indices())
          do_gp_face_term<false, is_dg_>(evaluator_m,
                                         evaluator_p,
                                         factor_values,
                                         factor_gradient,
                                         factor_hessians,
                                         q);

      evaluator_m.integrate(evaluation_flags, sum_into_values);
      evaluator_p.integrate(evaluation_flags, sum_into_values);
    }

    // Three helper functions to determine the category of face based on the
    // active_fe_index of the face-sharing cells.
    inline bool
    is_inside_face(std::pair<unsigned int, unsigned int> face_category) const
    {
      return is_inside(face_category.first) && is_inside(face_category.second);
    }

    inline bool
    is_mixed_face(std::pair<unsigned int, unsigned int> face_category) const
    {
      return (is_inside(face_category.first) &&
              is_intersected(face_category.second)) ||
             (is_intersected(face_category.first) &&
              is_inside(face_category.second));
    }

    inline bool is_intersected_face(
      std::pair<unsigned int, unsigned int> face_category) const
    {
      return is_intersected(face_category.first) &&
             is_intersected(face_category.second);
    }

    // Helper function to determine the relevant cell lengths of a face batch;
    VectorizedArrayType
    compute_diameter_of_inner_face_batch(unsigned int face_batch_index) const
    {
      const auto &face_info = matrix_free->get_face_info(face_batch_index);

      VectorizedArrayType diameter = 0.;
      for (unsigned int v = 0;
           v < matrix_free->n_active_entries_per_face_batch(face_batch_index);
           ++v)
        {
          const auto cell_batch_index_interior =
            face_info.cells_interior[v] / n_lanes;
          const auto cell_lane_index_interior =
            face_info.cells_interior[v] % n_lanes;
          const auto cell_batch_index_exterior =
            face_info.cells_exterior[v] / n_lanes;
          const auto cell_lane_index_exterior =
            face_info.cells_exterior[v] % n_lanes;

          diameter[v] = std::max(
            cell_diameter[cell_batch_index_interior][cell_lane_index_interior],
            cell_diameter[cell_batch_index_exterior][cell_lane_index_exterior]);
        }

      return diameter;
    }

    // Helper function which computes the interior penalty parameter for the
    // SIPG flux.
    VectorizedArrayType compute_interior_penalty_parameter(
      const VectorizedArrayType &diameter) const
    {
      return 5. / diameter;
    }

    // ObserverPointer of the MatrixFree and the NonMatching::MappingInfo
    // objects.
    ObserverPointer<const MatrixFree<dim, Number, VectorizedArrayType>>
      matrix_free;
    ObserverPointer<
      const NonMatching::MappingInfo<dim, dim, VectorizedArrayType>>
      mapping_info_cell;
    ObserverPointer<
      const NonMatching::MappingInfo<dim, dim, VectorizedArrayType>>
      mapping_info_surface;
    ObserverPointer<
      const NonMatching::MappingInfo<dim, dim, VectorizedArrayType>>
      mapping_info_faces;

    std::unique_ptr<FE_DGQ<dim>> fe_dgq;

    std::unique_ptr<GenericCellIntegrator> evaluator_cell;
    std::unique_ptr<GenericCellIntegrator> evaluator_surface;
    std::unique_ptr<GenericFaceIntegrator> evaluator_face_m;
    std::unique_ptr<GenericFaceIntegrator> evaluator_face_p;

    AlignedVector<VectorizedArrayType> cell_diameter;

    const unsigned int dof_index  = 0;
    const unsigned int quad_index = 0;

    bool is_dg = false;
  };

  // @sect3{The Jacobi preconditioner}
  // Assembly is done in the constructor by calling the compute_diagonal()
  // function of the operator.
  template <int dim>
  class JacobiPreconditioner
  {
  public:
    using VectorType = LinearAlgebra::distributed::Vector<double>;

    JacobiPreconditioner(const PoissonOperator<dim> &poisson_operator)
    {
      poisson_operator.initialize_dof_vector(inverse_diagonal);
      poisson_operator.compute_diagonal(inverse_diagonal);
      for (unsigned int i = 0; i < inverse_diagonal.locally_owned_size(); ++i)
        {
          if (std::abs(inverse_diagonal.local_element(i)) > 1.0e-10)
            inverse_diagonal.local_element(i) =
              1.0 / inverse_diagonal.local_element(i);
          else
            inverse_diagonal.local_element(i) = 1.0;
        }
    }

    void vmult(VectorType &dst, const VectorType &src) const
    {
      if (PointerComparison::equal(&dst, &src))
        {
          dst.scale(inverse_diagonal);
        }
      else
        {
          for (unsigned int i = 0; i < dst.locally_owned_size(); ++i)
            dst.local_element(i) =
              inverse_diagonal.local_element(i) * src.local_element(i);
        }
    }

  private:
    VectorType inverse_diagonal;
  };

  // @sect3{The PoissonSolver class template}
  template <int dim>
  class PoissonSolver
  {
    using Number              = double;
    using VectorizedArrayType = VectorizedArray<Number>;
    using VectorType          = LinearAlgebra::distributed::Vector<Number>;

    using CellIntegrator =
      FEEvaluation<dim, -1, 0, 1, Number, VectorizedArrayType>;
    using GenericCellIntegrator =
      FEPointEvaluation<1, dim, dim, VectorizedArrayType>;

    static constexpr unsigned int n_lanes = VectorizedArrayType::size();

  public:
    PoissonSolver();

    void run(bool is_dg, unsigned int fe_degree);

  private:
    void make_grid();

    void setup_discrete_level_set();

    void distribute_dofs();

    void setup_matrix_free();

    void setup_mapping_data();

    void solve();

    void output_results() const;

    double compute_L2_error() const;

    unsigned int fe_degree;

    const Functions::ConstantFunction<dim> rhs_function;

    parallel::distributed::Triangulation<dim> triangulation;

    ConditionalOStream pcout;

    // MatrixFree object
    MatrixFree<dim, double> matrix_free;

    PoissonOperator<dim> poisson_operator;

    // We need two separate DoFHandlers. The first manages the DoFs for the
    // discrete level set function that describes the geometry of the domain.
    std::unique_ptr<FE_Q<dim>> fe_level_set;
    DoFHandler<dim>            level_set_dof_handler;
    VectorType                 level_set;

    // The second DoFHandler manages the DoFs for the solution of the Poisson
    // equation.
    DoFHandler<dim> dof_handler;
    VectorType      solution;

    NonMatching::MeshClassifier<dim> mesh_classifier;

    VectorType rhs;

    const MappingQ<dim> mapping;

    // NonMatching::MappingInfo objects
    std::unique_ptr<NonMatching::MappingInfo<dim, dim, VectorizedArrayType>>
      mapping_info_cell;
    std::unique_ptr<NonMatching::MappingInfo<dim, dim, VectorizedArrayType>>
      mapping_info_surface;
    std::unique_ptr<NonMatching::MappingInfo<dim, dim, VectorizedArrayType>>
      mapping_info_faces;

    const unsigned int dof_index  = 0;
    const unsigned int quad_index = 0;

    bool is_dg;
  };



  template <int dim>
  PoissonSolver<dim>::PoissonSolver()
    : fe_degree(1)
    , rhs_function(4.0)
    , triangulation(MPI_COMM_WORLD)
    , pcout(std::cout,
            Utilities::MPI::this_mpi_process(
              triangulation.get_communicator()) == 0)
    , level_set_dof_handler(triangulation)
    , dof_handler(triangulation)
    , mesh_classifier(level_set_dof_handler, level_set)
    , mapping(1)
    , is_dg(false)
  {}



  // @sect3{Setting up the Background Mesh}
  // We generate a background mesh with perfectly Cartesian cells. Our domain is
  // a unit disc centered at the origin, so we need to make the background mesh
  // a bit larger than $[-1, 1]^{\text{dim}}$ to completely cover $\Omega$.
  template <int dim>
  void PoissonSolver<dim>::make_grid()
  {
    pcout << "Creating background mesh" << std::endl;

    triangulation.clear();
    GridGenerator::hyper_cube(triangulation, -1.21, 1.21);
    triangulation.refine_global(2);
  }



  // @sect3{Setting up the Discrete Level Set Function}
  // The discrete level set function is defined on the whole background mesh.
  // Thus, to set up the DoFHandler for the level set function, we distribute
  // DoFs over all elements in $\mathcal{T}_h$. We then set up the discrete
  // level set function by interpolating onto this finite element space.
  template <int dim>
  void PoissonSolver<dim>::setup_discrete_level_set()
  {
    pcout << "Setting up discrete level set function" << std::endl;

    fe_level_set = std::make_unique<FE_Q<dim>>(fe_degree);
    level_set_dof_handler.distribute_dofs(*fe_level_set);

    // We set up the level set vector with all locally relevant DoFs. This is
    // currently required by the NonMatching::MeshClassifier.
    level_set.reinit(level_set_dof_handler.locally_owned_dofs(),
                     DoFTools::extract_locally_relevant_dofs(
                       level_set_dof_handler),
                     triangulation.get_communicator());

    const Functions::SignedDistance::Sphere<dim> signed_distance_sphere;
    VectorTools::interpolate(level_set_dof_handler,
                             signed_distance_sphere,
                             level_set);

    level_set.update_ghost_values();
  }



  // @sect3{Distributing the DoFs}
  // We then use the NonMatching::MeshClassifier to check
  // NonMatching::LocationToLevelSet for each cell in the mesh and tell the
  // DoFHandler to use FE_Q or FE_DGQ on elements that are inside or
  // intersected, and FE_Nothing on the elements that are outside.
  template <int dim>
  void PoissonSolver<dim>::distribute_dofs()
  {
    pcout << "Distributing degrees of freedom" << std::endl;

    std::unique_ptr<FE_Poly<dim>> fe;
    if (is_dg)
      fe = std::make_unique<FE_DGQ<dim>>(fe_degree);
    else
      fe = std::make_unique<FE_Q<dim>>(fe_degree);

    hp::FECollection<dim> fe_collection;
    fe_collection.push_back(*fe);               // inside
    fe_collection.push_back(FE_Nothing<dim>()); // outside
    fe_collection.push_back(*fe);               // intersected

    for (const auto &cell : dof_handler.active_cell_iterators() |
                              IteratorFilters::LocallyOwnedCell())
      {
        const NonMatching::LocationToLevelSet cell_location =
          mesh_classifier.location_to_level_set(cell);

        if (cell_location == NonMatching::LocationToLevelSet::inside)
          cell->set_active_fe_index(
            (unsigned int)NonMatching::LocationToLevelSet::inside);
        else if (cell_location == NonMatching::LocationToLevelSet::outside)
          cell->set_active_fe_index(
            (unsigned int)NonMatching::LocationToLevelSet::outside);
        else if (cell_location == NonMatching::LocationToLevelSet::intersected)
          cell->set_active_fe_index(
            (unsigned int)NonMatching::LocationToLevelSet::intersected);
        else
          cell->set_active_fe_index(
            (unsigned int)NonMatching::LocationToLevelSet::outside);
      }

    dof_handler.distribute_dofs(fe_collection);
  }



  // @sect3{Setting up the MatrixFree object}
  template <int dim>
  void PoissonSolver<dim>::setup_matrix_free()
  {
    pcout << "Setting up matrix-free" << std::endl;
    QGauss<1>         quadrature(fe_degree + 1);
    AffineConstraints affine_constraints;
    affine_constraints.close();

    // setup MatrixFree::AdditionalData
    typename MatrixFree<dim, Number, VectorizedArrayType>::AdditionalData
      additional_data;
    additional_data.mapping_update_flags =
      update_gradients | update_values | update_quadrature_points;
    additional_data.mapping_update_flags_inner_faces =
      update_values | update_gradients;
    additional_data.mapping_update_flags_boundary_faces =
      update_values | update_gradients | update_quadrature_points;
    if (dof_handler.get_fe().degree > 1)
      additional_data.mapping_update_flags_inner_faces |= update_hessians;

    // setup MatrixFree object
    matrix_free.reinit(
      mapping, dof_handler, affine_constraints, quadrature, additional_data);
  }


  // @sect3{Setting up the NonMatching::MappingInfo objects}
  template <int dim>
  void PoissonSolver<dim>::setup_mapping_data()
  {
    pcout << "Setting up non matching mapping info" << std::endl;
    auto is_intersected_cell =
      [&](const TriaIterator<CellAccessor<dim, dim>> &cell) {
        return mesh_classifier.location_to_level_set(cell) ==
               NonMatching::LocationToLevelSet::intersected;
      };

    // We start by filling the containers in matrix-free ordering,
    std::vector<Quadrature<dim>> quad_vec_cells;
    quad_vec_cells.reserve(
      (matrix_free.n_cell_batches() + matrix_free.n_ghost_cell_batches()) *
      n_lanes);

    std::vector<NonMatching::ImmersedSurfaceQuadrature<dim>> quad_vec_surface;
    quad_vec_surface.reserve(
      (matrix_free.n_cell_batches() + matrix_free.n_ghost_cell_batches()) *
      n_lanes * n_lanes);

    hp::QCollection<1> q_collection1D(QGauss<1>(fe_degree + 1));

    NonMatching::DiscreteQuadratureGenerator<dim> quadrature_generator(
      q_collection1D, level_set_dof_handler, level_set);

    std::vector<typename DoFHandler<dim>::cell_iterator> vector_accessors;
    vector_accessors.reserve(
      (matrix_free.n_cell_batches() + matrix_free.n_ghost_cell_batches()) *
      n_lanes);
    for (unsigned int cell_batch = 0;
         cell_batch <
         matrix_free.n_cell_batches() + matrix_free.n_ghost_cell_batches();
         ++cell_batch)
      for (unsigned int v = 0; v < n_lanes; ++v)
        {
          if (v < matrix_free.n_active_entries_per_cell_batch(cell_batch))
            vector_accessors.push_back(
              matrix_free.get_cell_iterator(cell_batch, v));
          else
            vector_accessors.push_back(
              matrix_free.get_cell_iterator(cell_batch, 0));

          const auto &cell = vector_accessors.back();

          if (is_intersected_cell(cell))
            {
              quadrature_generator.generate(cell);

              quad_vec_cells.push_back(
                quadrature_generator.get_inside_quadrature());
              quad_vec_surface.push_back(
                quadrature_generator.get_surface_quadrature());
            }
          else
            {
              quad_vec_cells.emplace_back();
              quad_vec_surface.emplace_back();
            }
        }

    // then we initialize the NonMatching::MappingInfo objects to precompute
    // mapping information for cells and surface quadrature points.
    mapping_info_cell = std::make_unique<
      NonMatching::MappingInfo<dim, dim, VectorizedArray<Number>>>(
      mapping, update_values | update_gradients | update_JxW_values);
    mapping_info_cell->reinit_cells(vector_accessors, quad_vec_cells);

    mapping_info_surface = std::make_unique<
      NonMatching::MappingInfo<dim, dim, VectorizedArray<Number>>>(
      mapping,
      update_values | update_gradients | update_JxW_values |
        update_normal_vectors);
    mapping_info_surface->reinit_surface(vector_accessors, quad_vec_surface);

    // In case of DG, we also have to compute mapping data for cut faces, so we
    // again fill the containers in matrix-free ordering.
    if (is_dg)
      {
        NonMatching::DiscreteFaceQuadratureGenerator<dim>
          face_quadrature_generator(q_collection1D,
                                    level_set_dof_handler,
                                    level_set);

        std::vector<Quadrature<dim - 1>> quad_vec_faces;
        quad_vec_faces.reserve((matrix_free.n_inner_face_batches() +
                                matrix_free.n_boundary_face_batches() +
                                matrix_free.n_ghost_inner_face_batches()) *
                               n_lanes);
        std::vector<
          std::pair<typename DoFHandler<dim>::cell_iterator, unsigned int>>
          vector_face_accessors_m;
        vector_face_accessors_m.reserve(
          (matrix_free.n_inner_face_batches() +
           matrix_free.n_boundary_face_batches() +
           matrix_free.n_ghost_inner_face_batches()) *
          n_lanes);
        // Fill container for inner face batches,
        unsigned int face_batch = 0;
        for (; face_batch < matrix_free.n_inner_face_batches(); ++face_batch)
          {
            for (unsigned int v = 0; v < n_lanes; ++v)
              {
                if (v < matrix_free.n_active_entries_per_face_batch(face_batch))
                  vector_face_accessors_m.push_back(
                    matrix_free.get_face_iterator(face_batch, v, true));
                else
                  vector_face_accessors_m.push_back(
                    matrix_free.get_face_iterator(face_batch, 0, true));

                const auto &cell_m = vector_face_accessors_m.back().first;

                const unsigned int f = vector_face_accessors_m.back().second;

                if (is_intersected_cell(cell_m))
                  {
                    face_quadrature_generator.generate(cell_m, f);
                    quad_vec_faces.push_back(
                      face_quadrature_generator.get_inside_quadrature());
                  }
                else
                  quad_vec_faces.emplace_back();
              }
          }
        // then for boundary face batches,
        for (; face_batch < (matrix_free.n_inner_face_batches() +
                             matrix_free.n_boundary_face_batches());
             ++face_batch)
          {
            for (unsigned int v = 0; v < n_lanes; ++v)
              {
                if (v < matrix_free.n_active_entries_per_face_batch(face_batch))
                  vector_face_accessors_m.push_back(
                    matrix_free.get_face_iterator(face_batch, v, true));
                else
                  vector_face_accessors_m.push_back(
                    matrix_free.get_face_iterator(face_batch, 0, true));

                const auto &cell_m = vector_face_accessors_m.back().first;

                const unsigned int f = vector_face_accessors_m.back().second;

                if (is_intersected_cell(cell_m))
                  {
                    face_quadrature_generator.generate(cell_m, f);
                    quad_vec_faces.push_back(
                      face_quadrature_generator.get_inside_quadrature());
                  }
                else
                  quad_vec_faces.emplace_back();
              }
          }
        // and finally for ghost inner face batches.
        for (; face_batch < (matrix_free.n_inner_face_batches() +
                             matrix_free.n_boundary_face_batches() +
                             matrix_free.n_ghost_inner_face_batches());
             ++face_batch)
          {
            for (unsigned int v = 0; v < n_lanes; ++v)
              {
                if (v < matrix_free.n_active_entries_per_face_batch(face_batch))
                  vector_face_accessors_m.push_back(
                    matrix_free.get_face_iterator(face_batch, v, true));
                else
                  vector_face_accessors_m.push_back(
                    matrix_free.get_face_iterator(face_batch, 0, true));

                const auto &cell_m = vector_face_accessors_m.back().first;

                const unsigned int f = vector_face_accessors_m.back().second;

                if (is_intersected_cell(cell_m))
                  {
                    face_quadrature_generator.generate(cell_m, f);
                    quad_vec_faces.push_back(
                      face_quadrature_generator.get_inside_quadrature());
                  }
                else
                  quad_vec_faces.emplace_back();
              }
          }

        // And finally, initialize the NonMatching::MappingInfo object for the
        // cut faces.
        mapping_info_faces = std::make_unique<
          NonMatching::MappingInfo<dim, dim, VectorizedArray<Number>>>(
          mapping,
          update_values | update_gradients | update_JxW_values |
            update_normal_vectors);
        mapping_info_faces->reinit_faces(vector_face_accessors_m,
                                         quad_vec_faces);
      }
  }


  // @sect3{Solving the System}
  // Here we create the Jacobi preconditioner, which assembles the diagonal on
  // construction. Then, we call the preconditioned conjugate gradient solver.
  template <int dim>
  void PoissonSolver<dim>::solve()
  {
    pcout << "Solving system" << std::endl;

    JacobiPreconditioner<dim> jacobi_preconditioner(poisson_operator);

    const unsigned int   max_iterations = solution.size();
    SolverControl        solver_control(max_iterations);
    SolverCG<VectorType> solver(solver_control);
    solver.solve(poisson_operator, solution, rhs, jacobi_preconditioner);
  }



  // @sect3{Data Output}
  // Since both DoFHandler instances use the same triangulation, we can add both
  // the level set function and the solution to the same vtu-file. Further, we
  // do not want to output the cells that have NonMatching::LocationToLevelSet
  // value outside. To disregard them, we write a small lambda function and use
  // the set_cell_selection function of the DataOut class.
  template <int dim>
  void PoissonSolver<dim>::output_results() const
  {
    pcout << "Writing vtu file" << std::endl;

    DataOut<dim>          data_out;
    DataOutBase::VtkFlags flags;
    flags.write_higher_order_cells = true;
    data_out.set_flags(flags);
    data_out.add_data_vector(dof_handler, solution, "solution");
    data_out.add_data_vector(level_set_dof_handler, level_set, "level_set");

    data_out.set_cell_selection(
      [this](const typename Triangulation<dim>::cell_iterator &cell) {
        return cell->is_active() && cell->is_locally_owned() &&
               mesh_classifier.location_to_level_set(cell) !=
                 NonMatching::LocationToLevelSet::outside;
      });

    data_out.build_patches();
    data_out.write_vtu_with_pvtu_record("./",
                                        "step-95",
                                        0,
                                        triangulation.get_communicator());
  }



  // @sect3{L2-Error}
  // To test that the implementation works as expected, we want to compute the
  // error in the solution in the $L^2$-norm. The analytical solution to the
  // Poisson problem stated in the introduction reads
  // @f{align*}{
  //  u(x) = - \frac{2}{\text{dim}}(\| x \|^2 - 1) , \qquad x \in
  //  \overline{\Omega}.
  // @f}
  // We first create a function corresponding to the analytical solution:
  template <int dim>
  class AnalyticalSolution : public Function<dim>
  {
  public:
    double value(const Point<dim>  &point,
                 const unsigned int component = 0) const override;
  };



  template <int dim>
  double AnalyticalSolution<dim>::value(const Point<dim>  &point,
                                        const unsigned int component) const
  {
    AssertIndexRange(component, this->n_components);
    (void)component;

    return -2. / dim * (point.norm_square() - 1.);
  }



  // Of course, the analytical solution, and thus also the error, is only
  // defined in $\overline{\Omega}$. Thus, to compute the $L^2$-error we must
  // proceed in the same way as for the operator.
  template <int dim>
  double PoissonSolver<dim>::compute_L2_error() const
  {
    pcout << "Computing L2 error" << std::endl;

    const QGauss<1> quadrature_1D(fe_degree + 1);

    AnalyticalSolution<dim> analytical_solution;
    double                  error_L2_squared = 0;

    // The quadrature point operation to compute the $L^2$-error used for
    // FEEvaluation and FEPointEvaluation.
    auto l2_kernel = [](auto                &evaluator,
                        const Function<dim> &analytical_solution_function,
                        const unsigned int   q) {
      const auto q_points = evaluator.quadrature_point(q);

      VectorizedArrayType value = 0.;

      for (unsigned int v = 0; v < n_lanes; ++v)
        {
          Point<dim> q_point;
          for (unsigned int d = 0; d < dim; ++d)
            q_point[d] = q_points[d][v];

          value[v] = analytical_solution_function.value(q_point);
        }

      const auto difference = evaluator.get_value(q) - value;

      evaluator.submit_value(difference * difference, q);
    };

    // We then iterate over the cells that have NonMatching::LocationToLevelSet
    // value inside or intersected. For each quadrature point, we compute
    // the pointwise error and use this to compute the integral.
    unsigned int dummy = 0;
    matrix_free.template cell_loop<unsigned int, VectorType>(
      [&](const MatrixFree<dim, Number, VectorizedArrayType> &,
          unsigned int &,
          const VectorType                            &src,
          const std::pair<unsigned int, unsigned int> &cell_range) {
        CellIntegrator evaluator(matrix_free, dof_index, quad_index);

        FE_DGQ<dim> fe_dgq(fe_degree);

        GenericCellIntegrator evaluator_cell(*mapping_info_cell, fe_dgq);

        const auto cell_range_category =
          matrix_free.get_cell_range_category(cell_range);

        if (is_inside(cell_range_category))
          {
            for (unsigned int cell_batch_index = cell_range.first;
                 cell_batch_index < cell_range.second;
                 ++cell_batch_index)
              {
                evaluator.reinit(cell_batch_index);

                evaluator.read_dof_values(src);

                evaluator.evaluate(EvaluationFlags::values);

                for (unsigned int q : evaluator.quadrature_point_indices())
                  l2_kernel(evaluator, analytical_solution, q);

                for (unsigned int v = 0;
                     v < matrix_free.n_active_entries_per_cell_batch(
                           cell_batch_index);
                     ++v)
                  error_L2_squared += evaluator.integrate_value()[v];
              }
          }
        else if (is_intersected(cell_range_category))
          {
            const auto dofs_per_cell = evaluator.dofs_per_cell;

            for (unsigned int cell_batch_index = cell_range.first;
                 cell_batch_index < cell_range.second;
                 ++cell_batch_index)
              {
                evaluator.reinit(cell_batch_index);
                evaluator.read_dof_values(src);

                for (unsigned int v = 0;
                     v < matrix_free.n_active_entries_per_cell_batch(
                           cell_batch_index);
                     ++v)
                  {
                    const unsigned int cell_index =
                      cell_batch_index * n_lanes + v;

                    evaluator_cell.reinit(cell_index);

                    evaluator_cell.evaluate(
                      StridedArrayView<const Number, n_lanes>(
                        &evaluator.begin_dof_values()[0][v], dofs_per_cell),
                      EvaluationFlags::values);

                    for (const auto q :
                         evaluator_cell.quadrature_point_indices())
                      l2_kernel(evaluator_cell, analytical_solution, q);

                    error_L2_squared += evaluator_cell.integrate_value();
                  }
              }
          }
      },
      dummy,
      solution);

    return std::sqrt(
      Utilities::MPI::sum(error_L2_squared, triangulation.get_communicator()));
  }



  // @sect3{A Convergence Study}
  // Finally, we do a convergence study to check that the $L^2$-error decreases
  // with the expected rate. We refine the background mesh a few times. In each
  // refinement cycle, we solve the problem, compute the error, and add the
  // $L^2$-error and the mesh size to a ConvergenceTable.
  template <int dim>
  void PoissonSolver<dim>::run(const bool         is_dg_in,
                               const unsigned int fe_degree_in)
  {
    is_dg     = is_dg_in;
    fe_degree = fe_degree_in;

    if (is_dg)
      pcout << "Run DG convergence study with degree " << fe_degree
            << std::endl;
    else
      pcout << "Run CG convergence study with degree " << fe_degree
            << std::endl;

    dealii::Timer      timer;
    ConvergenceTable   convergence_table;
    const unsigned int n_refinements = 3;

    make_grid();
    for (unsigned int cycle = 0; cycle <= n_refinements; cycle++)
      {
        pcout << "Refinement cycle " << cycle << std::endl;
        triangulation.refine_global(1);
        setup_discrete_level_set();
        pcout << "Classifying cells" << std::endl;
        mesh_classifier.reclassify();
        distribute_dofs();
        setup_matrix_free();
        setup_mapping_data();

        poisson_operator.reinit(matrix_free,
                                mapping_info_cell.get(),
                                mapping_info_surface.get(),
                                mapping_info_faces.get(),
                                is_dg,
                                fe_degree);

        matrix_free.initialize_dof_vector(solution);
        matrix_free.initialize_dof_vector(rhs);

        poisson_operator.rhs(rhs, rhs_function);

        solve();
        if (cycle == 3)
          output_results();
        const double error_L2 = compute_L2_error();
        const double cell_side_length =
          triangulation.begin_active(triangulation.n_levels() - 1)
            ->minimum_vertex_distance();

        convergence_table.add_value("Cycle", cycle);
        convergence_table.add_value("Mesh size", cell_side_length);
        convergence_table.add_value("L2-Error", error_L2);

        convergence_table.evaluate_convergence_rates(
          "L2-Error", ConvergenceTable::reduction_rate_log2);
        convergence_table.set_scientific("L2-Error", true);

        pcout << std::endl;
        if (Utilities::MPI::this_mpi_process(
              triangulation.get_communicator()) == 0)
          convergence_table.write_text(pcout.get_stream());
        pcout << std::endl;
      }
    pcout << "wall time: " << timer.wall_time() << "\n" << std::endl;
  }

} // namespace Step95



// @sect3{The main() function}
int main(int argc, char **argv)
{
  dealii::Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);

  constexpr int      dim       = 2;
  const unsigned int fe_degree = 2;

  Step95::PoissonSolver<dim> poisson_solver;
  // run CG
  poisson_solver.run(false /* is_dg */, fe_degree);
  // run DG
  poisson_solver.run(true /* is_dg */, fe_degree);
}
