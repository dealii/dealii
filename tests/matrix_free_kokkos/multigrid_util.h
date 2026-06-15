// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#ifndef dealii_tests_multigrid_util_h
#define dealii_tests_multigrid_util_h

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/lac/diagonal_matrix.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>

#include <deal.II/matrix_free/portable_fe_evaluation.h>
#include <deal.II/matrix_free/portable_matrix_free.h>
#include <deal.II/matrix_free/tools.h>

#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/multigrid/mg_matrix.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/multigrid.h>
#include <deal.II/multigrid/portable_mg_transfer_global_coarsening.h>
#include <deal.II/multigrid/portable_mg_transfer_global_coarsening.templates.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


// using namespace dealii;

DEAL_II_NAMESPACE_OPEN

namespace Portable
{

  // --- Matrix-free operator utilities ---

  template <int dim, int fe_degree, typename Number = double>
  class LaplaceOperatorQuad
  {
  public:
    DEAL_II_HOST_DEVICE
    LaplaceOperatorQuad() = default;

    DEAL_II_HOST_DEVICE void
    operator()(FEEvaluation<dim, fe_degree, fe_degree + 1, 1, Number> *fe_eval,
               const int q_point) const
    {
      fe_eval->submit_gradient(fe_eval->get_gradient(q_point), q_point);
    }

    static const unsigned int n_q_points =
      dealii::Utilities::pow(fe_degree + 1, dim);
  };

  template <int dim, int fe_degree, typename Number = double>
  class LocalLaplaceOperator
  {
  public:
    LocalLaplaceOperator() = default;

    DEAL_II_HOST_DEVICE void
    operator()(const typename MatrixFree<dim, Number>::Data *data,
               const DeviceVector<Number>                   &src,
               DeviceVector<Number>                         &dst) const
    {
      FEEvaluation<dim, fe_degree, fe_degree + 1, 1, Number> fe_eval(data);
      fe_eval.read_dof_values(src);
      fe_eval.evaluate(EvaluationFlags::gradients);

      LaplaceOperatorQuad<dim, fe_degree> quad;
      data->for_each_quad_point(
        [&](const int &q_point) { quad(&fe_eval, q_point); });

      fe_eval.integrate(EvaluationFlags::gradients);
      fe_eval.distribute_local_to_global(dst);
    }

    static constexpr unsigned int n_q_points =
      Utilities::pow(fe_degree + 1, dim);
  };


  template <int dim, int fe_degree, typename Number = double>
  class LaplaceOperator : public EnableObserverPointer
  {
    using number = Number;

  public:
    LaplaceOperator() = default;

    void
    reinit(const DoFHandler<dim>           &dof_handler,
           const AffineConstraints<number> &constraints,
           const unsigned int mg_level = numbers::invalid_unsigned_int)
    {
      this->constraints = &constraints;

      const MappingQ<dim> mapping(fe_degree);

      typename MatrixFree<dim, number>::AdditionalData additional_data;

      additional_data.mapping_update_flags =
        update_gradients | update_JxW_values | update_quadrature_points;
      additional_data.overlap_communication_computation = false;

      additional_data.mg_level = mg_level;

      const QGauss<1> quadrature_1d(fe_degree + 1);
      this->matrix_free.reinit(
        mapping, dof_handler, constraints, quadrature_1d, additional_data);
    }

    void
    vmult(LinearAlgebra::distributed::Vector<number, MemorySpace::Default> &dst,
          const LinearAlgebra::distributed::Vector<number, MemorySpace::Default>
            &src) const
    {
      dst = 0.;
      LocalLaplaceOperator<dim, fe_degree, number> cell_operator;
      matrix_free.cell_loop(cell_operator, src, dst);
      matrix_free.copy_constrained_values(src, dst);
    }

    void
    Tvmult(
      LinearAlgebra::distributed::Vector<number, MemorySpace::Default> &dst,
      const LinearAlgebra::distributed::Vector<number, MemorySpace::Default>
        &src) const
    {
      this->vmult(dst, src);
    }

    void
    initialize_dof_vector(
      LinearAlgebra::distributed::Vector<number, MemorySpace::Default> &vec)
      const
    {
      this->matrix_free.initialize_dof_vector(vec);
    }

    void
    compute_diagonal()
    {
      this->inverse_diagonal_entries.reset(
        new DiagonalMatrix<
          LinearAlgebra::distributed::Vector<number, MemorySpace::Default>>());
      LinearAlgebra::distributed::Vector<number, MemorySpace::Default>
        &inverse_diagonal = inverse_diagonal_entries->get_vector();
      this->initialize_dof_vector(inverse_diagonal);

      LaplaceOperatorQuad<dim, fe_degree, number> operator_quad;

      MatrixFreeTools::
        compute_diagonal<dim, fe_degree, fe_degree + 1, 1, number>(
          matrix_free,
          inverse_diagonal,
          operator_quad,
          EvaluationFlags::gradients,
          EvaluationFlags::gradients);

      number *raw_diagonal = inverse_diagonal.get_values();

      Kokkos::parallel_for(
        inverse_diagonal.locally_owned_size(), KOKKOS_LAMBDA(int i) {
          Assert(
            raw_diagonal[i] > 0.,
            ExcMessage(
              "No diagonal entry in a positive definite operator should be zero"));
          raw_diagonal[i] = 1. / raw_diagonal[i];
        });
    }

    void
    rhs(LinearAlgebra::distributed::Vector<number, MemorySpace::Default> &dst,
        const AffineConstraints<number> &constraints) const
    {
      const auto &dof_handler = matrix_free.get_dof_handler();
      const auto &fe          = dof_handler.get_fe();

      LinearAlgebra::distributed::Vector<double, MemorySpace::Host> dst_host(
        dof_handler.locally_owned_dofs(),
        DoFTools::extract_locally_relevant_dofs(dof_handler),
        dof_handler.get_mpi_communicator());

      const QGauss<dim> quadrature_formula(fe_degree + 1);

      FEValues<dim> fe_values(fe,
                              quadrature_formula,
                              update_values | update_JxW_values);

      Vector<double> cell_rhs(n_local_dofs);

      std::vector<types::global_dof_index> local_dof_indices(n_local_dofs);

      for (const auto &cell : dof_handler.active_cell_iterators())
        {
          if (cell->is_locally_owned())
            {
              cell_rhs = 0;

              fe_values.reinit(cell);

              for (unsigned int q_index = 0; q_index < n_q_points; ++q_index)
                for (unsigned int i = 0; i < n_local_dofs; ++i)
                  cell_rhs(i) += (fe_values.shape_value(i, q_index) * 1.0 *
                                  fe_values.JxW(q_index));

              cell->get_dof_indices(local_dof_indices);
              constraints.distribute_local_to_global(cell_rhs,
                                                     local_dof_indices,
                                                     dst_host);
            }
        }

      dst_host.compress(VectorOperation::add);
      LinearAlgebra::ReadWriteVector<double> rw_vector(
        dof_handler.locally_owned_dofs());

      rw_vector.import_elements(dst_host, VectorOperation::insert);
      dst.import_elements(rw_vector, VectorOperation::insert);
    }

    types::global_dof_index
    m() const
    {
      if (matrix_free.get_vector_partitioner())
        return matrix_free.get_vector_partitioner()->size();
      else
        return matrix_free.get_dof_handler().n_dofs();
    }

    types::global_dof_index
    n() const
    {
      if (matrix_free.get_vector_partitioner())
        return matrix_free.get_vector_partitioner()->size();
      else
        return matrix_free.get_dof_handler().n_dofs();
    }

    number
    el(const types::global_dof_index row,
       const types::global_dof_index col) const
    {
      (void)col;
      Assert(row == col, ExcNotImplemented());
      Assert(inverse_diagonal_entries.get() != nullptr &&
               inverse_diagonal_entries->m() > 0,
             ExcNotInitialized());

      return 1.0 / (*inverse_diagonal_entries)(row, row);
    }

    std::shared_ptr<DiagonalMatrix<
      LinearAlgebra::distributed::Vector<number, MemorySpace::Default>>>
    get_matrix_diagonal_inverse() const
    {
      return inverse_diagonal_entries;
    }

    const MatrixFree<dim, number> &
    get_matrix_free() const
    {
      return matrix_free;
    }

    const std::shared_ptr<const Utilities::MPI::Partitioner> &
    get_vector_partitioner() const
    {
      return matrix_free.get_vector_partitioner();
    }

    static constexpr unsigned int n_local_dofs =
      Utilities::pow(fe_degree + 1, dim);
    static constexpr unsigned int n_q_points =
      Utilities::pow(fe_degree + 1, dim);

  private:
    MatrixFree<dim, number> matrix_free;

    ObserverPointer<const AffineConstraints<number>> constraints;

    std::shared_ptr<DiagonalMatrix<
      LinearAlgebra::distributed::Vector<number, MemorySpace::Default>>>
      inverse_diagonal_entries;

    std::vector<
      Kokkos::View<unsigned int **, MemorySpace::Default::kokkos_space>>
      dof_indices_per_color;
  };


  // --- Multigrid transfer utilities ---

  template <typename VectorType, typename SmootherType>
  class MGCoarseFromSmoother : public MGCoarseGridBase<VectorType>
  {
  public:
    MGCoarseFromSmoother(const SmootherType &mg_smoother, const bool is_empty)
      : smoother(mg_smoother)
      , is_empty(is_empty)
    {}

    virtual void
    operator()(const unsigned int level,
               VectorType        &dst,
               const VectorType  &src) const override
    {
      if (is_empty)
        return;
      smoother[level].vmult(dst, src);
    }

    const SmootherType &smoother;
    const bool          is_empty;
  };


  template <int dim,
            typename LevelMatrixType,
            typename VectorType,
            typename TransferType,
            typename SmootherType>
  class VCycleMultigrid : public EnableObserverPointer
  {
  public:
    using number = typename VectorType::value_type;

    VCycleMultigrid(const MGLevelObject<LevelMatrixType> &mg_matrices,
                    const MGLevelObject<TransferType>    &mg_transfers,
                    const MGLevelObject<SmootherType>    &mg_smoothers)
      : minlevel(mg_matrices.min_level())
      , maxlevel(mg_matrices.max_level())
      , mg_matrices(mg_matrices)
      , mg_transfers(mg_transfers)
      , mg_smoothers(mg_smoothers)
      , coarse(mg_smoothers, false)
      , solution(minlevel, maxlevel)
      , defect(minlevel, maxlevel)
      , t(minlevel, maxlevel)
    {
      for (unsigned int level = minlevel; level <= maxlevel; ++level)
        {
          mg_matrices[level].initialize_dof_vector(solution[level]);
          defect[level] = solution[level];
          t[level]      = solution[level];
        }
    }

    // Solve with the conjugate gradient method preconditioned by the V-cycle
    // (invoking this->vmult).
    void
    solve_cg(SolverControl    &solver_control,
             VectorType       &dst,
             const VectorType &src)
    {
      SolverCG<VectorType> solver_cg(solver_control);

      dst = 0;
      solver_cg.solve(mg_matrices.back(), dst, src, *this);

      solution[maxlevel] = dst;
    }

    /**
     * Implements vmult() for the CG solver.
     */
    void
    vmult(VectorType &dst, const VectorType &src) const
    {
      defect[maxlevel] = src;

      v_cycle(maxlevel);

      dst = solution[maxlevel];
    }

  private:
    /**
     * Implements the v-cycle
     */
    void
    v_cycle(const unsigned int level) const
    {
      if (level == minlevel)
        {
          // Accuracy on coarsest level should be comparable to overall level
          // accuracy (~1e-3)
          (coarse)(level, solution[level], defect[level]);

          return;
        }

      // Pre-smoothing
      mg_smoothers[level].vmult(solution[level], defect[level]);

      // Compute residual
      mg_matrices[level].vmult(t[level], solution[level]);
      t[level].sadd(-1.0, 1.0, defect[level]);

      // Restrict residual to the next coarser level
      defect[level - 1] = 0;
      mg_transfers[level].restrict_and_add(defect[level - 1], t[level]);

      // Recursive call to v_cycle on the coarser level
      v_cycle(level - 1);

      // Prolongate coarse correction and add to current solution
      mg_transfers[level].prolongate_and_add(solution[level],
                                             solution[level - 1]);

      // Post-smoothing
      mg_smoothers[level].step(solution[level], defect[level]);
    }

    /**
     * Lowest level of cells.
     */
    unsigned int minlevel;

    /**
     * Highest level of cells.
     */
    unsigned int maxlevel;


    /**
     * The matrix for each level.
     */
    const MGLevelObject<LevelMatrixType> &mg_matrices;


    /**
     * The transfer between each level.
     */
    const MGLevelObject<TransferType> &mg_transfers;


    /**
     * The smmother for each level.
     */

    const MGLevelObject<SmootherType> &mg_smoothers;

    /**
     * The coarse solver
     */
    MGCoarseFromSmoother<VectorType, MGLevelObject<SmootherType>> coarse;

    /**
     * The solution update after the multigrid step.
     */
    mutable MGLevelObject<VectorType> solution;

    /**
     * Input vector for the cycle. Contains the defect of the outer method
     * projected to the multilevel vectors.
     */
    mutable MGLevelObject<VectorType> defect;

    /**
     * Auxiliary vector.
     */
    mutable MGLevelObject<VectorType> t;
  };



} // namespace Portable

DEAL_II_NAMESPACE_CLOSE

#endif
