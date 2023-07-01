// ---------------------------------------------------------------------
//
// Copyright (C) 2022 - 2023 by the deal.II authors
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

//
// Description:
//
// A performance benchmark based on step 37 (without variable coefficient)
// that measures timings for grid creation, setup of unknowns, multigrid
// levels and solve for a Poisson problem with the performance-oriented
// matrix-free framework.
//
// Status: stable
//
// Note: this test is marked "stable" and used for performance
// instrumentation in our testsuite, https://dealii.org/performance_tests
//

#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>

#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_matrix.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_transfer_matrix_free.h>
#include <deal.II/multigrid/multigrid.h>

#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>
#include <memory>

#define ENABLE_MPI

#include "performance_test_driver.h"

using namespace dealii;

dealii::ConditionalOStream debug_output(std::cout, false);

const unsigned int degree_finite_element = 3;



template <int dim, int fe_degree, typename number>
class LaplaceOperator
  : public MatrixFreeOperators::Base<dim,
                                     LinearAlgebra::distributed::Vector<number>>
{
public:
  using value_type = number;

  LaplaceOperator();

  void
  vmult(LinearAlgebra::distributed::Vector<number> &      dst,
        const LinearAlgebra::distributed::Vector<number> &src) const;

  void
  vmult(LinearAlgebra::distributed::Vector<number> &      dst,
        const LinearAlgebra::distributed::Vector<number> &src,
        const std::function<void(const unsigned int, const unsigned int)>
          &operation_before_loop,
        const std::function<void(const unsigned int, const unsigned int)>
          &operation_after_loop) const;

  virtual void
  compute_diagonal() override;

private:
  virtual void
  apply_add(
    LinearAlgebra::distributed::Vector<number> &      dst,
    const LinearAlgebra::distributed::Vector<number> &src) const override;

  void
  local_apply(const MatrixFree<dim, number> &                   data,
              LinearAlgebra::distributed::Vector<number> &      dst,
              const LinearAlgebra::distributed::Vector<number> &src,
              const std::pair<unsigned int, unsigned int> &cell_range) const;

  void
  local_compute_diagonal(
    const MatrixFree<dim, number> &              data,
    LinearAlgebra::distributed::Vector<number> & dst,
    const unsigned int &                         dummy,
    const std::pair<unsigned int, unsigned int> &cell_range) const;
};



template <int dim, int fe_degree, typename number>
LaplaceOperator<dim, fe_degree, number>::LaplaceOperator()
  : MatrixFreeOperators::Base<dim, LinearAlgebra::distributed::Vector<number>>()
{}



template <int dim, int fe_degree, typename number>
void
LaplaceOperator<dim, fe_degree, number>::local_apply(
  const MatrixFree<dim, number> &                   data,
  LinearAlgebra::distributed::Vector<number> &      dst,
  const LinearAlgebra::distributed::Vector<number> &src,
  const std::pair<unsigned int, unsigned int> &     cell_range) const
{
  FEEvaluation<dim, fe_degree, fe_degree + 1, 1, number> phi(data);

  for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
    {
      phi.reinit(cell);
      phi.gather_evaluate(src, EvaluationFlags::gradients);
      for (unsigned int q = 0; q < phi.n_q_points; ++q)
        phi.submit_gradient(phi.get_gradient(q), q);
      phi.integrate_scatter(EvaluationFlags::gradients, dst);
    }
}



template <int dim, int fe_degree, typename number>
void
LaplaceOperator<dim, fe_degree, number>::apply_add(
  LinearAlgebra::distributed::Vector<number> &      dst,
  const LinearAlgebra::distributed::Vector<number> &src) const
{
  this->data->cell_loop(&LaplaceOperator::local_apply, this, dst, src);
}



template <int dim, int fe_degree, typename number>
void
LaplaceOperator<dim, fe_degree, number>::vmult(
  LinearAlgebra::distributed::Vector<number> &      dst,
  const LinearAlgebra::distributed::Vector<number> &src) const
{
  this->data->cell_loop(&LaplaceOperator::local_apply, this, dst, src, true);
  for (const unsigned int i : this->data->get_constrained_dofs())
    dst.local_element(i) = src.local_element(i);
}



template <int dim, int fe_degree, typename number>
void
LaplaceOperator<dim, fe_degree, number>::vmult(
  LinearAlgebra::distributed::Vector<number> &      dst,
  const LinearAlgebra::distributed::Vector<number> &src,
  const std::function<void(const unsigned int, const unsigned int)>
    &operation_before_loop,
  const std::function<void(const unsigned int, const unsigned int)>
    &operation_after_loop) const
{
  this->data->cell_loop(&LaplaceOperator::local_apply,
                        this,
                        dst,
                        src,
                        operation_before_loop,
                        operation_after_loop);
  for (const unsigned int i : this->data->get_constrained_dofs())
    dst.local_element(i) = src.local_element(i);
}



template <int dim, int fe_degree, typename number>
void
LaplaceOperator<dim, fe_degree, number>::compute_diagonal()
{
  this->inverse_diagonal_entries.reset(
    new DiagonalMatrix<LinearAlgebra::distributed::Vector<number>>());
  LinearAlgebra::distributed::Vector<number> &inverse_diagonal =
    this->inverse_diagonal_entries->get_vector();
  this->data->initialize_dof_vector(inverse_diagonal);
  unsigned int dummy = 0;
  this->data->cell_loop(&LaplaceOperator::local_compute_diagonal,
                        this,
                        inverse_diagonal,
                        dummy);

  this->set_constrained_entries_to_one(inverse_diagonal);

  for (unsigned int i = 0; i < inverse_diagonal.locally_owned_size(); ++i)
    {
      Assert(inverse_diagonal.local_element(i) > 0.,
             ExcMessage("No diagonal entry in a positive definite operator "
                        "should be zero"));
      inverse_diagonal.local_element(i) =
        1. / inverse_diagonal.local_element(i);
    }
}



template <int dim, int fe_degree, typename number>
void
LaplaceOperator<dim, fe_degree, number>::local_compute_diagonal(
  const MatrixFree<dim, number> &             data,
  LinearAlgebra::distributed::Vector<number> &dst,
  const unsigned int &,
  const std::pair<unsigned int, unsigned int> &cell_range) const
{
  FEEvaluation<dim, fe_degree, fe_degree + 1, 1, number> phi(data);

  AlignedVector<VectorizedArray<number>> diagonal(phi.dofs_per_cell);

  for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
    {
      phi.reinit(cell);
      for (unsigned int i = 0; i < phi.dofs_per_cell; ++i)
        {
          for (unsigned int j = 0; j < phi.dofs_per_cell; ++j)
            phi.submit_dof_value(VectorizedArray<number>(), j);
          phi.submit_dof_value(make_vectorized_array<number>(1.), i);

          phi.evaluate(EvaluationFlags::gradients);
          for (unsigned int q = 0; q < phi.n_q_points; ++q)
            phi.submit_gradient(phi.get_gradient(q), q);
          phi.integrate(EvaluationFlags::gradients);
          diagonal[i] = phi.get_dof_value(i);
        }
      for (unsigned int i = 0; i < phi.dofs_per_cell; ++i)
        phi.submit_dof_value(diagonal[i], i);
      phi.distribute_local_to_global(dst);
    }
}



template <int dim, typename MatrixType>
class MGTransferMF
  : public MGTransferMatrixFree<dim, typename MatrixType::value_type>
{
public:
  void
  setup(const MGLevelObject<MatrixType> &laplace,
        const MGConstrainedDoFs &        mg_constrained_dofs)
  {
    this->MGTransferMatrixFree<dim, typename MatrixType::value_type>::
      initialize_constraints(mg_constrained_dofs);
    laplace_operator = &laplace;
  }

  template <class InVector, int spacedim>
  void
  copy_to_mg(
    const DoFHandler<dim, spacedim> &mg_dof,
    MGLevelObject<
      LinearAlgebra::distributed::Vector<typename MatrixType::value_type>> &dst,
    const InVector &src) const
  {
    for (unsigned int level = dst.min_level(); level <= dst.max_level();
         ++level)
      (*laplace_operator)[level].initialize_dof_vector(dst[level]);
    MGLevelGlobalTransfer<LinearAlgebra::distributed::Vector<
      typename MatrixType::value_type>>::copy_to_mg(mg_dof, dst, src);
  }

private:
  const MGLevelObject<MatrixType> *laplace_operator;
};



template <int dim>
class LaplaceProblem
{
public:
  LaplaceProblem();

  Measurement
  run();

private:
  void
  setup_grid();
  void
  setup_dofs();
  void
  setup_matrix_free();
  void
  assemble_rhs();
  void
  setup_transfer();
  void
  setup_smoother();
  void
  solve();

  parallel::distributed::Triangulation<dim> triangulation;

  FE_Q<dim>       fe;
  DoFHandler<dim> dof_handler;

  MappingQ1<dim> mapping;

  AffineConstraints<double> constraints;
  using SystemMatrixType = LaplaceOperator<dim, degree_finite_element, double>;
  SystemMatrixType system_matrix;

  MGConstrainedDoFs mg_constrained_dofs;
  using LevelMatrixType = LaplaceOperator<dim, degree_finite_element, float>;
  MGLevelObject<LevelMatrixType> mg_matrices;

  LinearAlgebra::distributed::Vector<double> solution;
  LinearAlgebra::distributed::Vector<double> system_rhs;

  MGTransferMF<dim, LevelMatrixType> mg_transfer;

  using SmootherType =
    PreconditionChebyshev<LevelMatrixType,
                          LinearAlgebra::distributed::Vector<float>>;
  mg::SmootherRelaxation<SmootherType,
                         LinearAlgebra::distributed::Vector<float>>
    mg_smoother;
};



template <int dim>
LaplaceProblem<dim>::LaplaceProblem()
#ifdef DEAL_II_WITH_P4EST
  : triangulation(
      MPI_COMM_WORLD,
      Triangulation<dim>::limit_level_difference_at_vertices,
      parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy)
#else
  : triangulation(Triangulation<dim>::limit_level_difference_at_vertices)
#endif
  , fe(degree_finite_element)
  , dof_handler(triangulation)
{}



template <int dim>
void
LaplaceProblem<dim>::setup_grid()
{
  GridGenerator::hyper_cube(triangulation, 0., 1.);
  switch (get_testing_environment())
    {
      case TestingEnvironment::light:
        triangulation.refine_global(5);
        break;
      case TestingEnvironment::medium:
        triangulation.refine_global(6);
        break;
      case TestingEnvironment::heavy:
        triangulation.refine_global(7);
        break;
    }
}



template <int dim>
void
LaplaceProblem<dim>::setup_dofs()
{
  system_matrix.clear();
  mg_matrices.clear_elements();

  dof_handler.distribute_dofs(fe);
  dof_handler.distribute_mg_dofs();

  debug_output << "Number of DoFs: " << dof_handler.n_dofs() << std::endl;

  IndexSet locally_relevant_dofs;
  DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);

  constraints.clear();
  constraints.reinit(locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  VectorTools::interpolate_boundary_values(
    mapping, dof_handler, 0, Functions::ZeroFunction<dim>(), constraints);
  constraints.close();

  // Renumber DoFs
  typename MatrixFree<dim, float>::AdditionalData additional_data;
  additional_data.tasks_parallel_scheme =
    MatrixFree<dim, float>::AdditionalData::none;

  const std::set<types::boundary_id> dirichlet_boundary = {0};
  mg_constrained_dofs.initialize(dof_handler);
  mg_constrained_dofs.make_zero_boundary_constraints(dof_handler,
                                                     dirichlet_boundary);

  for (unsigned int level = 0; level < triangulation.n_global_levels(); ++level)
    {
      IndexSet relevant_dofs;
      DoFTools::extract_locally_relevant_level_dofs(dof_handler,
                                                    level,
                                                    relevant_dofs);
      AffineConstraints<double> level_constraints;
      level_constraints.reinit(relevant_dofs);
      level_constraints.add_lines(
        mg_constrained_dofs.get_boundary_indices(level));
      level_constraints.close();
      additional_data.mg_level = level;

      DoFRenumbering::matrix_free_data_locality(dof_handler,
                                                level_constraints,
                                                additional_data);
    }
}



template <int dim>
void
LaplaceProblem<dim>::setup_matrix_free()
{
  typename MatrixFree<dim, double>::AdditionalData additional_data;
  additional_data.tasks_parallel_scheme =
    MatrixFree<dim, double>::AdditionalData::none;
  additional_data.mapping_update_flags =
    (update_gradients | update_JxW_values | update_quadrature_points);
  std::shared_ptr<MatrixFree<dim, double>> system_mf_storage(
    new MatrixFree<dim, double>());
  system_mf_storage->reinit(mapping,
                            dof_handler,
                            constraints,
                            QGauss<1>(fe.degree + 1),
                            additional_data);
  system_matrix.initialize(system_mf_storage);


  system_matrix.initialize_dof_vector(solution);
  system_matrix.initialize_dof_vector(system_rhs);

  const unsigned int nlevels = triangulation.n_global_levels();
  mg_matrices.resize(0, nlevels - 1);

  const std::set<types::boundary_id> dirichlet_boundary = {0};
  mg_constrained_dofs.initialize(dof_handler);
  mg_constrained_dofs.make_zero_boundary_constraints(dof_handler,
                                                     dirichlet_boundary);

  for (unsigned int level = 0; level < nlevels; ++level)
    {
      IndexSet relevant_dofs;
      DoFTools::extract_locally_relevant_level_dofs(dof_handler,
                                                    level,
                                                    relevant_dofs);
      AffineConstraints<double> level_constraints;
      level_constraints.reinit(relevant_dofs);
      level_constraints.add_lines(
        mg_constrained_dofs.get_boundary_indices(level));
      level_constraints.close();

      typename MatrixFree<dim, float>::AdditionalData additional_data;
      additional_data.tasks_parallel_scheme =
        MatrixFree<dim, float>::AdditionalData::none;
      additional_data.mapping_update_flags =
        (update_gradients | update_JxW_values | update_quadrature_points);
      additional_data.mg_level = level;
      std::shared_ptr<MatrixFree<dim, float>> mg_mf_storage_level(
        new MatrixFree<dim, float>());
      mg_mf_storage_level->reinit(mapping,
                                  dof_handler,
                                  level_constraints,
                                  QGauss<1>(fe.degree + 1),
                                  additional_data);

      mg_matrices[level].initialize(mg_mf_storage_level,
                                    mg_constrained_dofs,
                                    level);
    }
}



template <int dim>
void
LaplaceProblem<dim>::assemble_rhs()
{
  Timer time;

  system_rhs = 0;
  FEEvaluation<dim, degree_finite_element> phi(
    *system_matrix.get_matrix_free());
  for (unsigned int cell = 0;
       cell < system_matrix.get_matrix_free()->n_cell_batches();
       ++cell)
    {
      phi.reinit(cell);
      for (unsigned int q = 0; q < phi.n_q_points; ++q)
        phi.submit_value(make_vectorized_array<double>(1.0), q);
      phi.integrate(EvaluationFlags::values);
      phi.distribute_local_to_global(system_rhs);
    }
  system_rhs.compress(VectorOperation::add);
}



template <int dim>
void
LaplaceProblem<dim>::setup_transfer()
{
  mg_transfer.setup(mg_matrices, mg_constrained_dofs);
  std::vector<std::shared_ptr<const Utilities::MPI::Partitioner>> partitioners(
    dof_handler.get_triangulation().n_global_levels());
  for (unsigned int level = 0; level < partitioners.size(); ++level)
    {
      LinearAlgebra::distributed::Vector<float> vec;
      mg_matrices[level].initialize_dof_vector(vec);
      partitioners[level] = vec.get_partitioner();
    }
  mg_transfer.build(dof_handler, partitioners);
}



template <int dim>
void
LaplaceProblem<dim>::setup_smoother()
{
  MGLevelObject<typename SmootherType::AdditionalData> smoother_data;
  smoother_data.resize(0, triangulation.n_global_levels() - 1);
  for (unsigned int level = 0; level < triangulation.n_global_levels(); ++level)
    {
      if (level > 0)
        {
          smoother_data[level].smoothing_range     = 15.;
          smoother_data[level].degree              = 5;
          smoother_data[level].eig_cg_n_iterations = 10;
        }
      else
        {
          smoother_data[0].smoothing_range     = 1e-3;
          smoother_data[0].degree              = numbers::invalid_unsigned_int;
          smoother_data[0].eig_cg_n_iterations = mg_matrices[0].m();
        }
      mg_matrices[level].compute_diagonal();
      smoother_data[level].preconditioner =
        mg_matrices[level].get_matrix_diagonal_inverse();
    }
  mg_smoother.initialize(mg_matrices, smoother_data);
}


template <int dim>
void
LaplaceProblem<dim>::solve()
{
  MGCoarseGridApplySmoother<LinearAlgebra::distributed::Vector<float>>
    mg_coarse;
  mg_coarse.initialize(mg_smoother);

  mg::Matrix<LinearAlgebra::distributed::Vector<float>> mg_matrix(mg_matrices);

  MGLevelObject<MatrixFreeOperators::MGInterfaceOperator<LevelMatrixType>>
    mg_interface_matrices;
  mg_interface_matrices.resize(0, triangulation.n_global_levels() - 1);
  for (unsigned int level = 0; level < triangulation.n_global_levels(); ++level)
    mg_interface_matrices[level].initialize(mg_matrices[level]);
  mg::Matrix<LinearAlgebra::distributed::Vector<float>> mg_interface(
    mg_interface_matrices);

  Multigrid<LinearAlgebra::distributed::Vector<float>> mg(
    mg_matrix, mg_coarse, mg_transfer, mg_smoother, mg_smoother);
  mg.set_edge_matrices(mg_interface, mg_interface);

  PreconditionMG<dim,
                 LinearAlgebra::distributed::Vector<float>,
                 MGTransferMF<dim, LevelMatrixType>>
    preconditioner(dof_handler, mg, mg_transfer);


  SolverControl solver_control(100, 1e-12 * system_rhs.l2_norm());
  SolverCG<LinearAlgebra::distributed::Vector<double>> cg(solver_control);

  constraints.set_zero(solution);
  cg.solve(system_matrix, solution, system_rhs, preconditioner);
}



template <int dim>
Measurement
LaplaceProblem<dim>::run()
{
  std::map<std::string, dealii::Timer> timer;

  timer["setup_grid"].start();
  setup_grid();
  timer["setup_grid"].stop();

  timer["setup_dofs"].start();
  setup_dofs();
  timer["setup_dofs"].stop();

  timer["setup_matrix_free"].start();
  setup_matrix_free();
  timer["setup_matrix_free"].stop();

  timer["assemble_rhs"].start();
  assemble_rhs();
  timer["assemble_rhs"].stop();

  timer["setup_transfer"].start();
  setup_transfer();
  timer["setup_transfer"].stop();

  timer["setup_smoother"].start();
  setup_smoother();
  timer["setup_smoother"].stop();

  timer["solve"].start();
  solve();
  timer["solve"].stop();

  const unsigned int n_repeat = 50;
  timer["matvec_double"].start();
  for (unsigned int t = 0; t < n_repeat; ++t)
    system_matrix.vmult(system_rhs, solution);
  timer["matvec_double"].stop();

  LinearAlgebra::distributed::Vector<float> vec1, vec2;
  mg_matrices[mg_matrices.max_level()].initialize_dof_vector(vec1);
  vec2.reinit(vec1);
  timer["matvec_float"].start();
  for (unsigned int t = 0; t < n_repeat; ++t)
    mg_matrices[mg_matrices.max_level()].vmult(vec2, vec1);
  timer["matvec_float"].stop();

  mg_matrices[mg_matrices.max_level() - 1].initialize_dof_vector(vec1);
  vec2.reinit(vec1);
  timer["matvec_float_coarser"].start();
  for (unsigned int t = 0; t < n_repeat; ++t)
    mg_matrices[mg_matrices.max_level() - 1].vmult(vec2, vec1);
  timer["matvec_float_coarser"].stop();

  debug_output << std::endl;
  return {timer["setup_grid"].wall_time(),
          timer["setup_dofs"].wall_time(),
          timer["setup_matrix_free"].wall_time(),
          timer["assemble_rhs"].wall_time(),
          timer["setup_transfer"].wall_time(),
          timer["setup_smoother"].wall_time(),
          timer["solve"].wall_time(),
          timer["matvec_double"].wall_time(),
          timer["matvec_float"].wall_time(),
          timer["matvec_float_coarser"].wall_time()};
}


std::tuple<Metric, unsigned int, std::vector<std::string>>
describe_measurements()
{
  return {Metric::timing,
          4,
          {"setup_grid",
           "setup_dofs",
           "setup_matrix_free",
           "assemble_rhs",
           "setup_transfer",
           "setup_smoother",
           "solve",
           "matvec_double",
           "matvec_float",
           "matvec_float_coarser"}};
}


Measurement
perform_single_measurement()
{
  return LaplaceProblem<3>().run();
}
