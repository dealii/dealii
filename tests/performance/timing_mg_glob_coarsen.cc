// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

//
// Description:
//
// A performance benchmark assessing a Poisson problem with the
// performance-oriented matrix-free framework. As opposed to the related
// timing_step_37 benchmark, this case uses the global-coarsening multigrid
// framework with p-multigrid and using a locally refined mesh with hanging
// nodes. It also uses a setup with multiple DoFHandler objects, imitating the
// projection from a related (higher-order) DG function space.
//
// Status: stable
//
// Note: this test is marked "stable" and used for performance
// instrumentation in our testsuite, https://dealii.org/performance_tests
//

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/timer.h>

#include <deal.II/distributed/repartitioning_policy_tools.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>
#include <deal.II/matrix_free/tools.h>

#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_matrix.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_transfer_global_coarsening.h>
#include <deal.II/multigrid/multigrid.h>

#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>
#include <memory>

#define ENABLE_MPI

#include "performance_test_driver.h"

using namespace dealii;

dealii::ConditionalOStream debug_output(std::cout, false);



template <int dim, typename number = double>
class LaplaceOperator : public EnableObserverPointer
{
public:
  using value_type = number;
  using VectorType = LinearAlgebra::distributed::Vector<number>;

  LaplaceOperator(){};

  void
  initialize(const Mapping<dim>              &mapping,
             const DoFHandler<dim>           &dof_handler,
             const AffineConstraints<number> &constraints,
             const DoFHandler<dim>           &dg_dof_handler)
  {
    const QGauss<1> quad(dof_handler.get_fe().degree + 1);
    const QGauss<1> dg_quad(dg_dof_handler.get_fe().degree + 1);
    typename MatrixFree<dim, number>::AdditionalData mf_data;
    mf_data.tasks_parallel_scheme =
      MatrixFree<dim, number>::AdditionalData::none;
    mf_data.mapping_update_flags |= update_quadrature_points;
    mf_data.mapping_update_flags_inner_faces =
      (update_gradients | update_JxW_values);
    mf_data.mapping_update_flags_boundary_faces =
      (update_gradients | update_JxW_values);
    AffineConstraints<number> dg_constraints;

    data.reinit(mapping,
                std::vector<const DoFHandler<dim> *>{
                  {&dof_handler, &dg_dof_handler}},
                std::vector<const AffineConstraints<number> *>{
                  {&constraints, &dg_constraints}},
                std::vector<Quadrature<1>>{{quad, dg_quad}},
                mf_data);
  }

  void
  initialize(const Mapping<dim>              &mapping,
             const DoFHandler<dim>           &dof_handler,
             const AffineConstraints<number> &constraints)
  {
    const QGauss<1> quad(dof_handler.get_fe().degree + 1);
    typename MatrixFree<dim, number>::AdditionalData mf_data;
    mf_data.tasks_parallel_scheme =
      MatrixFree<dim, number>::AdditionalData::none;
    Assert(dof_handler.get_fe().dofs_per_vertex > 0,
           ExcNotImplemented("Only continuous elements implemented"));

    data.reinit(mapping, dof_handler, constraints, quad, mf_data);
  }

  void
  vmult(VectorType &dst, const VectorType &src) const
  {
    data.cell_loop(&LaplaceOperator::local_apply, this, dst, src, true);
    for (const auto i : data.get_constrained_dofs())
      dst.local_element(i) = src.local_element(i);
  }

  void
  vmult(VectorType       &dst,
        const VectorType &src,
        const std::function<void(const unsigned int, const unsigned int)>
          &operation_before_loop,
        const std::function<void(const unsigned int, const unsigned int)>
          &operation_after_loop) const
  {
    data.cell_loop(&LaplaceOperator::local_apply,
                   this,
                   dst,
                   src,
                   operation_before_loop,
                   operation_after_loop);
  }

  void
  Tvmult(VectorType &dst, const VectorType &src) const
  {
    vmult(dst, src);
  }

  number
  el(const types::global_dof_index, const types::global_dof_index) const
  {
    AssertThrow(false, ExcNotImplemented());
    return number(0.);
  }

  types::global_dof_index
  m() const
  {
    return data.get_vector_partitioner()->size();
  }

  types::global_dof_index
  n() const
  {
    return data.get_vector_partitioner()->size();
  }

  void
  initialize_dof_vector(VectorType        &vector,
                        const unsigned int component = 0) const
  {
    data.initialize_dof_vector(vector, component);
  }

  void
  compute_inverse_diagonal()
  {
    inverse_diagonal_entries = std::make_shared<DiagonalMatrix<VectorType>>();
    data.initialize_dof_vector(inverse_diagonal_entries->get_vector());
    MatrixFreeTools::
      compute_diagonal<dim, -1, 0, 1, number, VectorizedArray<number>>(
        data, inverse_diagonal_entries->get_vector(), [](auto &eval) {
          eval.evaluate(EvaluationFlags::gradients);
          for (const unsigned int q : eval.quadrature_point_indices())
            eval.submit_gradient(eval.get_gradient(q), q);
          eval.integrate(EvaluationFlags::gradients);
        });

    for (number &entry : inverse_diagonal_entries->get_vector())
      if (std::abs(entry) > 1e-10)
        entry = 1. / entry;
      else
        entry = 1.;
  }

  const std::shared_ptr<DiagonalMatrix<VectorType>> &
  get_matrix_diagonal_inverse() const
  {
    return inverse_diagonal_entries;
  }

  const MatrixFree<dim, number> &
  get_matrix_free() const
  {
    return data;
  }

private:
  void
  local_apply(const MatrixFree<dim, number>               &data,
              VectorType                                  &dst,
              const VectorType                            &src,
              const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    FEEvaluation<dim, -1, 0, 1, number> eval(data);

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        eval.reinit(cell);
        eval.gather_evaluate(src, EvaluationFlags::gradients);
        for (const unsigned int q : eval.quadrature_point_indices())
          eval.submit_gradient(eval.get_gradient(q), q);
        eval.integrate_scatter(EvaluationFlags::gradients, dst);
      }
  }

  MatrixFree<dim, number>                     data;
  std::shared_ptr<DiagonalMatrix<VectorType>> inverse_diagonal_entries;
};



template <typename Number>
void
make_zero_mean(const std::vector<unsigned int>            &constrained_dofs,
               LinearAlgebra::distributed::Vector<Number> &vec)
{
  // set constrained entries to zero
  for (const unsigned int index : constrained_dofs)
    vec.local_element(index) = 0.;

  // rescale mean value computed among all vector entries to the vector size
  // without constraints
  const unsigned int n_unconstrained_dofs =
    vec.locally_owned_size() - constrained_dofs.size();
  vec.add(
    -vec.mean_value() * vec.size() /
    Utilities::MPI::sum(n_unconstrained_dofs, vec.get_mpi_communicator()));

  // set constrained entries to zero again, this should now have zero mean
  for (const unsigned int index : constrained_dofs)
    vec.local_element(index) = 0.;

  Assert(std::abs(vec.mean_value()) <
           std::numeric_limits<Number>::epsilon() * vec.size(),
         ExcInternalError());
}



// class to impose zero-mean constraint on coarse level
template <class VectorType = LinearAlgebra::distributed::Vector<double>>
class MGCoarseSolverSingular : public MGCoarseGridBase<VectorType>
{
public:
  void
  clear()
  {
    coarse_smooth = nullptr;
  }

  void
  initialize(const MGSmootherBase<VectorType> &coarse_smooth,
             const std::vector<unsigned int>  &constrained_dofs)
  {
    this->coarse_smooth    = &coarse_smooth;
    this->constrained_dofs = &constrained_dofs;
  }

  void
  operator()(const unsigned int level,
             VectorType        &dst,
             const VectorType  &src) const override
  {
    src_copy.reinit(src, true);
    src_copy.copy_locally_owned_data_from(src);
    make_zero_mean(*constrained_dofs, src_copy);
    coarse_smooth->apply(level, dst, src_copy);
    make_zero_mean(*constrained_dofs, dst);
  }

private:
  ObserverPointer<const MGSmootherBase<VectorType>> coarse_smooth;
  const std::vector<unsigned int>                  *constrained_dofs;

  mutable VectorType src_copy;
};


const Tensor<2, 3> deformation{
  {{1.05, 1e-3, 1e-2}, {1e-3, 1., -1e-3}, {1e-2, -1e-3, 0.95}}};

template <int dim>
class Solution : public Function<dim>
{
public:
  Solution()
    : Function<dim>(1)
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int = 0) const override
  {
    const Point<dim> x     = Point<dim>(invert(deformation) * p);
    double           value = 1.0;
    for (unsigned int d = 0; d < dim; ++d)
      value *= std::cos(8. * numbers::PI * x[d]);
    return value;
  }
};


template <int dim>
class RightHandSide : public Function<dim>
{
public:
  RightHandSide()
    : Function<dim>(1)
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int = 0) const override
  {
    Solution<dim> sol;
    return dim * 64. * numbers::PI * numbers::PI * sol.value(p);
  }
};



template <int dim>
class LaplaceProblem
{
public:
  LaplaceProblem(const unsigned int degree);

  Measurement
  run();

private:
  void
  setup_grid();
  void
  create_coarse_triangulations();
  void
  setup_dofs();
  void
  setup_matrix_free();
  void
  setup_smoother();
  void
  setup_transfer();
  void
  compute_rhs();
  void
  solve();
  void
  embed_solution_to_dg();

  parallel::distributed::Triangulation<dim>              triangulation;
  std::vector<std::shared_ptr<const Triangulation<dim>>> coarse_triangulations;
  MappingQ<dim>                                          mapping;
  FE_DGQ<dim>                                            dg_fe;
  DoFHandler<dim>                                        dg_dof_handler;
  MGLevelObject<std::unique_ptr<FE_Q<dim>>>              fes;
  MGLevelObject<DoFHandler<dim>>                         dof_handlers;

  LinearAlgebra::distributed::Vector<double> dg_rhs;
  LinearAlgebra::distributed::Vector<double> rhs;
  LinearAlgebra::distributed::Vector<double> solution;
  LinearAlgebra::distributed::Vector<double> dg_solution;

  LaplaceOperator<dim, double>               system_matrix;
  MGLevelObject<AffineConstraints<float>>    level_constraints;
  MGLevelObject<LaplaceOperator<dim, float>> level_matrices;
  using VectorTypeMG = LinearAlgebra::distributed::Vector<float>;

  using SmootherType =
    PreconditionChebyshev<LaplaceOperator<dim, float>, VectorTypeMG>;
  mg::SmootherRelaxation<SmootherType, VectorTypeMG> mg_smoother;

  MGLevelObject<std::unique_ptr<MGTwoLevelTransferBase<VectorTypeMG>>>
                                                                 mg_transfers;
  std::unique_ptr<MGTransferGlobalCoarsening<dim, VectorTypeMG>> mg_transfer;
};



template <int dim>
LaplaceProblem<dim>::LaplaceProblem(const unsigned int degree)
#ifdef DEAL_II_WITH_P4EST
  : triangulation(MPI_COMM_WORLD)
#else
  : triangulation()
#endif
  , mapping(1)
  , dg_fe(degree)
  , dg_dof_handler(triangulation)
{}



template <int dim>
void
LaplaceProblem<dim>::setup_grid()
{
  GridGenerator::hyper_cube(triangulation, 0., 1.);
  GridTools::transform([](const Point<dim> &p) { return deformation * p; },
                       triangulation);

  switch (get_testing_environment())
    {
      case TestingEnvironment::light:
        triangulation.refine_global(3);
        break;
      case TestingEnvironment::medium:
        triangulation.refine_global(4);
        break;
      case TestingEnvironment::heavy:
        triangulation.refine_global(5);
        break;
    }

  for (const auto &cell : triangulation.active_cell_iterators())
    if (cell->is_locally_owned() && cell->center().norm() < 1.1)
      cell->set_refine_flag();
  triangulation.execute_coarsening_and_refinement();
  for (const auto &cell : triangulation.active_cell_iterators())
    if (cell->is_locally_owned() &&
        cell->center().distance(Point<dim>(0.3, 0.3, 0.3)) < 0.5)
      cell->set_refine_flag();
  triangulation.execute_coarsening_and_refinement();
}


template <int dim>
void
LaplaceProblem<dim>::create_coarse_triangulations()
{
  coarse_triangulations =
    MGTransferGlobalCoarseningTools::create_geometric_coarsening_sequence(
      triangulation,
      RepartitioningPolicyTools::MinimalGranularityPolicy<dim>(16));
}


template <int dim>
void
LaplaceProblem<dim>::setup_dofs()
{
  dg_dof_handler.reinit(triangulation);
  dg_dof_handler.distribute_dofs(dg_fe);

  // the solver uses ph-multigrid according to
  // https://doi.org/10.1016/j.jcp.2020.109538 and
  // https://doi.org/10.1145/3580314

  // start by creating levels of continuous elements
  std::vector<unsigned int> p_levels({dg_fe.degree - 1});
  while (p_levels.back() > 2)
    p_levels.push_back(std::max(p_levels.back() - 2, 2u));
  fes.resize(0, p_levels.size() - 1);
  for (unsigned int level = 0; level < p_levels.size(); ++level)
    fes[level] =
      std::make_unique<FE_Q<dim>>(p_levels[p_levels.size() - 1 - level]);

  dof_handlers.resize(0, coarse_triangulations.size() - 1 + fes.max_level());
  level_constraints.resize(0, dof_handlers.max_level());
  for (unsigned int level = dof_handlers.min_level();
       level <= dof_handlers.max_level();
       ++level)
    {
      DoFHandler<dim> &dof_h = dof_handlers[level];
      dof_h.reinit(
        *coarse_triangulations[std::min(level,
                                        triangulation.n_global_levels() - 1)]);
      if (level < coarse_triangulations.size())
        dof_h.distribute_dofs(*fes[0]);
      else
        dof_h.distribute_dofs(*fes[level + 1 - coarse_triangulations.size()]);

      IndexSet relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_h);
      AffineConstraints<float> &constraints = level_constraints[level];
      constraints.reinit(dof_h.locally_owned_dofs(), relevant_dofs);
      DoFTools::make_hanging_node_constraints(dof_h, constraints);
      constraints.close();
      typename MatrixFree<dim, float>::AdditionalData additional_data;
      additional_data.tasks_parallel_scheme =
        MatrixFree<dim, float>::AdditionalData::none;

      DoFRenumbering::matrix_free_data_locality(dof_h,
                                                constraints,
                                                additional_data);

      // now create the final constraints object
      relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_h);
      constraints.clear();
      constraints.reinit(dof_h.locally_owned_dofs(), relevant_dofs);
      DoFTools::make_hanging_node_constraints(dof_h, constraints);
      constraints.close();
    }
}



template <int dim>
void
LaplaceProblem<dim>::setup_matrix_free()
{
  AffineConstraints<double> constraints_fine;
  constraints_fine.reinit(level_constraints.back().get_local_lines());
  constraints_fine.copy_from(level_constraints.back());
  system_matrix.initialize(mapping,
                           dof_handlers.back(),
                           constraints_fine,
                           dg_dof_handler);
  system_matrix.initialize_dof_vector(dg_rhs, 1);
  system_matrix.initialize_dof_vector(dg_solution, 1);
  system_matrix.initialize_dof_vector(rhs, 0);
  system_matrix.initialize_dof_vector(solution, 0);

  level_matrices.resize(0, dof_handlers.max_level());
  for (unsigned int level = dof_handlers.min_level();
       level <= dof_handlers.max_level();
       ++level)
    {
      level_matrices[level].initialize(mapping,
                                       dof_handlers[level],
                                       level_constraints[level]);
    }
}


template <int dim>
void
LaplaceProblem<dim>::setup_smoother()
{
  MGLevelObject<typename SmootherType::AdditionalData> smoother_data(
    0, dof_handlers.max_level());
  for (unsigned int level = dof_handlers.min_level();
       level <= dof_handlers.max_level();
       ++level)
    {
      level_matrices[level].compute_inverse_diagonal();

      // manually compute the eigenvalue estimate for Chebyshev because we
      // need to be careful with the constrained indices
      IterationNumberControl control(12, 1e-6, false, false);

      using VectorType = LinearAlgebra::distributed::Vector<float>;
      SolverCG<VectorType>        solver(control);
      internal::EigenvalueTracker eigenvalue_tracker;
      solver.connect_eigenvalues_slot(
        [&eigenvalue_tracker](const std::vector<double> &eigenvalues) {
          eigenvalue_tracker.slot(eigenvalues);
        });

      VectorType sol, rhs;
      level_matrices[level].initialize_dof_vector(sol);
      level_matrices[level].initialize_dof_vector(rhs);

      for (float &a : rhs)
        a = (double)rand() / RAND_MAX;
      make_zero_mean(
        level_matrices[level].get_matrix_free().get_constrained_dofs(), rhs);
      solver.solve(level_matrices[level],
                   sol,
                   rhs,
                   *level_matrices[level].get_matrix_diagonal_inverse());

      if (level > 0)
        {
          smoother_data[level].smoothing_range = 15.;
          smoother_data[level].degree          = 4;
        }
      else
        {
          // Coarse level: Use MG smoother as solver (should use p-multigrid
          // or AMG for complicated meshes)
          smoother_data[level].smoothing_range =
            eigenvalue_tracker.values.back() /
            eigenvalue_tracker.values.front();
          smoother_data[0].degree = numbers::invalid_unsigned_int;
        }
      smoother_data[level].max_eigenvalue = eigenvalue_tracker.values.back();
      smoother_data[level].eig_cg_n_iterations = 0;
      smoother_data[level].preconditioner =
        level_matrices[level].get_matrix_diagonal_inverse();
    }

  mg_smoother.initialize(level_matrices, smoother_data);
}


template <int dim>
void
LaplaceProblem<dim>::setup_transfer()
{
  mg_transfers.resize(0, dof_handlers.max_level());
  for (unsigned int level = 1; level <= dof_handlers.max_level(); ++level)
    {
      auto transfer = std::make_unique<MGTwoLevelTransfer<dim, VectorTypeMG>>();
      if (level < triangulation.n_global_levels())
        transfer->reinit(dof_handlers[level],
                         dof_handlers[level - 1],
                         level_constraints[level],
                         level_constraints[level - 1]);
      else
        transfer->reinit(level_matrices[level].get_matrix_free(),
                         0,
                         level_matrices[level - 1].get_matrix_free(),
                         0);

      mg_transfers[level] = std::move(transfer);
    }

  mg_transfer = std::make_unique<MGTransferGlobalCoarsening<dim, VectorTypeMG>>(
    mg_transfers, [&](const unsigned level, VectorTypeMG &vec) {
      level_matrices[level].initialize_dof_vector(vec);
    });
}


template <int dim>
void
LaplaceProblem<dim>::compute_rhs()
{
  // interpolate to nodes
  VectorTools::interpolate(mapping,
                           dg_dof_handler,
                           RightHandSide<dim>(),
                           dg_rhs);

  // do the interpolation 10 times to get better significance in the numbers
  for (unsigned int i = 0; i < 10; ++i)
    {
      rhs = 0.;
      FEEvaluation<dim, -1, 0, 1, double> dg_eval(
        system_matrix.get_matrix_free(), 1);
      FEEvaluation<dim, -1, 0, 1, double> eval(system_matrix.get_matrix_free(),
                                               0);
      for (unsigned int cell = 0;
           cell < system_matrix.get_matrix_free().n_cell_batches();
           ++cell)
        {
          eval.reinit(cell);
          dg_eval.reinit(cell);
          dg_eval.gather_evaluate(dg_rhs, EvaluationFlags::values);
          for (const unsigned int q : eval.quadrature_point_indices())
            eval.submit_value(dg_eval.get_value(q), q);
          eval.integrate_scatter(EvaluationFlags::values, rhs);
        }
      rhs.compress(VectorOperation::add);

      // since we use Neumann boundary conditions on the whole boundary, the
      // right hand side must have zero mean value to ensure a solvable system
      make_zero_mean(system_matrix.get_matrix_free().get_constrained_dofs(),
                     rhs);
    }
}



template <int dim>
void
LaplaceProblem<dim>::solve()
{
  MGCoarseSolverSingular<VectorTypeMG> mg_coarse;
  mg_coarse.initialize(
    mg_smoother, level_matrices[0].get_matrix_free().get_constrained_dofs());
  mg::Matrix<VectorTypeMG> mg_matrix(level_matrices);

  Multigrid<VectorTypeMG> mg(
    mg_matrix, mg_coarse, *mg_transfer, mg_smoother, mg_smoother);
  PreconditionMG<dim,
                 VectorTypeMG,
                 MGTransferGlobalCoarsening<dim, VectorTypeMG>>
    preconditioner(dof_handlers.back(), mg, *mg_transfer);

  SolverControl control(20, 1e-10 * rhs.l2_norm());
  SolverCG<LinearAlgebra::distributed::Vector<double>> solver(control);

  solver.solve(system_matrix, solution, rhs, preconditioner);
  AssertThrow(control.last_step() < 10,
              ExcMessage("Solve should converge in at most 10 iterations"));
}



template <int dim>
void
LaplaceProblem<dim>::embed_solution_to_dg()
{
  make_zero_mean(system_matrix.get_matrix_free().get_constrained_dofs(),
                 solution);

  FEEvaluation<dim, -1, 0, 1, double> dg_eval(system_matrix.get_matrix_free(),
                                              1);
  MatrixFreeOperators::CellwiseInverseMassMatrix<dim, -1> inverse_mass(dg_eval);
  FEEvaluation<dim, -1, 0, 1, double> eval(system_matrix.get_matrix_free(), 0);
  // to get better timings, run the evaluation 10 times
  for (unsigned int i = 0; i < 10; ++i)
    {
      solution.update_ghost_values();
      for (unsigned int cell = 0;
           cell < system_matrix.get_matrix_free().n_cell_batches();
           ++cell)
        {
          eval.reinit(cell);
          dg_eval.reinit(cell);
          eval.gather_evaluate(solution, EvaluationFlags::values);
          inverse_mass.transform_from_q_points_to_basis(
            1, eval.begin_values(), dg_eval.begin_dof_values());
          dg_eval.set_dof_values(dg_solution);
        }
      solution.zero_out_ghost_values();
    }

  // compute error
  double error = 0;
  for (unsigned int cell = 0;
       cell < system_matrix.get_matrix_free().n_cell_batches();
       ++cell)
    {
      dg_eval.reinit(cell);
      dg_eval.gather_evaluate(dg_solution, EvaluationFlags::values);
      Solution<dim> solution;
      double        local_error = 0;
      for (const unsigned int q : dg_eval.quadrature_point_indices())
        for (unsigned int v = 0;
             v <
             system_matrix.get_matrix_free().n_active_entries_per_cell_batch(
               cell);
             ++v)
          {
            Point<dim> quadrature_point;
            for (unsigned int d = 0; d < dim; ++d)
              quadrature_point[d] = dg_eval.quadrature_point(q)[d][v];
            local_error +=
              Utilities::fixed_power<2>(solution.value(quadrature_point) -
                                        dg_eval.get_value(q)[v]) *
              dg_eval.JxW(q)[v];
          }
      error += local_error;
    }
  error =
    std::sqrt(Utilities::MPI::sum(error, dg_solution.get_mpi_communicator()));
  // do to the deformed mesh, the chosen right hand side and solution match
  // only approximately - we request a tolerance of 1e-2
  AssertThrow(error < 1e-2, ExcMessage("Error should be less than 1e-2"));
}


template <int dim>
Measurement
LaplaceProblem<dim>::run()
{
  std::map<std::string, dealii::Timer> timer;

  timer["setup_grid"].start();
  setup_grid();
  timer["setup_grid"].stop();

  timer["setup_coarse_grids"].start();
  create_coarse_triangulations();
  timer["setup_coarse_grids"].stop();

  timer["setup_dofs"].start();
  setup_dofs();
  timer["setup_dofs"].stop();

  timer["setup_matrix_free"].start();
  setup_matrix_free();
  timer["setup_matrix_free"].stop();

  timer["setup_smoother"].start();
  setup_smoother();
  timer["setup_smoother"].stop();

  timer["setup_transfer"].start();
  setup_transfer();
  timer["setup_transfer"].stop();

  timer["compute_rhs"].start();
  compute_rhs();
  timer["compute_rhs"].stop();

  timer["solve"].start();
  solve();
  timer["solve"].stop();

  const unsigned int n_repeat = 50;
  timer["matvec_double"].start();
  for (unsigned int t = 0; t < n_repeat; ++t)
    system_matrix.vmult(rhs, solution);
  timer["matvec_double"].stop();

  LinearAlgebra::distributed::Vector<float> vec1, vec2;
  level_matrices[level_matrices.max_level()].initialize_dof_vector(vec1);
  vec2.reinit(vec1);
  timer["matvec_float"].start();
  for (unsigned int t = 0; t < n_repeat; ++t)
    level_matrices[level_matrices.max_level()].vmult(vec2, vec1);
  timer["matvec_float"].stop();

  timer["embed_dg_and_error"].start();
  embed_solution_to_dg();
  timer["embed_dg_and_error"].stop();

  debug_output << std::endl;
  return {timer["setup_grid"].wall_time(),
          timer["setup_coarse_grids"].wall_time(),
          timer["setup_dofs"].wall_time(),
          timer["setup_matrix_free"].wall_time(),
          timer["setup_smoother"].wall_time(),
          timer["setup_transfer"].wall_time(),
          timer["compute_rhs"].wall_time(),
          timer["solve"].wall_time(),
          timer["matvec_double"].wall_time(),
          timer["matvec_float"].wall_time(),
          timer["embed_dg_and_error"].wall_time()};
}


std::tuple<Metric, unsigned int, std::vector<std::string>>
describe_measurements()
{
  return {Metric::timing,
          5,
          {"setup_grid",
           "setup_coarse_grids",
           "setup_dofs",
           "setup_matrix_free",
           "setup_smoother",
           "setup_transfer",
           "compute_rhs",
           "solve",
           "matvec_double",
           "matvec_float",
           "embed_dg_and_error"}};
}


Measurement
perform_single_measurement()
{
  // run in 3d with degree 5, i.e., degree 4 for the FEM part making the
  // actual solve
  return LaplaceProblem<3>(5).run();
}
