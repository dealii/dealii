// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// same test as parallel_multigrid_adaptive_03 but using partition_color

#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_matrix.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_transfer_global_coarsening.h>
#include <deal.II/multigrid/multigrid.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


template <int dim,
          int fe_degree,
          int n_q_points_1d = fe_degree + 1,
          typename number   = double>
class LaplaceOperator : public EnableObserverPointer
{
public:
  using value_type = number;

  LaplaceOperator(){};


  void
  initialize(const Mapping<dim>      &mapping,
             const DoFHandler<dim>   &dof_handler,
             const MGConstrainedDoFs &mg_constrained_dofs,
             const std::map<types::boundary_id, const Function<dim> *>
                               &dirichlet_boundary,
             const unsigned int level,
             const bool         threaded)
  {
    const QGauss<1>                                  quad(n_q_points_1d);
    typename MatrixFree<dim, number>::AdditionalData addit_data;
    if (threaded)
      addit_data.tasks_parallel_scheme =
        MatrixFree<dim, number>::AdditionalData::partition_color;
    else
      addit_data.tasks_parallel_scheme =
        MatrixFree<dim, number>::AdditionalData::none;
    addit_data.tasks_block_size = 3;
    addit_data.mg_level         = level;
    AffineConstraints<double> constraints;
    if (level == numbers::invalid_unsigned_int)
      {
        const IndexSet relevant_dofs =
          DoFTools::extract_locally_relevant_dofs(dof_handler);
        constraints.reinit(dof_handler.locally_owned_dofs(), relevant_dofs);
        DoFTools::make_hanging_node_constraints(dof_handler, constraints);
        VectorTools::interpolate_boundary_values(dof_handler,
                                                 dirichlet_boundary,
                                                 constraints);
      }
    else
      {
        const IndexSet &locally_owned =
          dof_handler.locally_owned_mg_dofs(level);
        const IndexSet relevant_dofs =
          DoFTools::extract_locally_relevant_level_dofs(dof_handler, level);
        constraints.reinit(locally_owned, relevant_dofs);
        constraints.add_lines(mg_constrained_dofs.get_boundary_indices(level));

        const std::vector<types::global_dof_index> interface_indices =
          mg_constrained_dofs.get_refinement_edge_indices(level)
            .get_index_vector();
        edge_constrained_indices.clear();
        edge_constrained_indices.reserve(interface_indices.size());
        edge_constrained_values.resize(interface_indices.size());
        for (unsigned int i = 0; i < interface_indices.size(); ++i)
          if (locally_owned.is_element(interface_indices[i]))
            edge_constrained_indices.push_back(
              locally_owned.index_within_set(interface_indices[i]));
        have_interface_matrices =
          Utilities::MPI::max((unsigned int)edge_constrained_indices.size(),
                              MPI_COMM_WORLD) > 0;
      }
    constraints.close();

    data.reinit(mapping, dof_handler, constraints, quad, addit_data);

    if (level != numbers::invalid_unsigned_int)
      compute_inverse_diagonal();
  }

  void
  vmult(LinearAlgebra::distributed::Vector<number>       &dst,
        const LinearAlgebra::distributed::Vector<number> &src) const
  {
    dst = 0;
    vmult_add(dst, src);
  }

  void
  Tvmult(LinearAlgebra::distributed::Vector<number>       &dst,
         const LinearAlgebra::distributed::Vector<number> &src) const
  {
    dst = 0;
    vmult_add(dst, src);
  }

  void
  Tvmult_add(LinearAlgebra::distributed::Vector<number>       &dst,
             const LinearAlgebra::distributed::Vector<number> &src) const
  {
    vmult_add(dst, src);
  }

  void
  vmult_add(LinearAlgebra::distributed::Vector<number>       &dst,
            const LinearAlgebra::distributed::Vector<number> &src) const
  {
    Assert(src.partitioners_are_globally_compatible(
             *data.get_dof_info(0).vector_partitioner),
           ExcInternalError());
    Assert(dst.partitioners_are_globally_compatible(
             *data.get_dof_info(0).vector_partitioner),
           ExcInternalError());

    // set zero Dirichlet values on the input vector (and remember the src and
    // dst values because we need to reset them at the end)
    for (unsigned int i = 0; i < edge_constrained_indices.size(); ++i)
      {
        edge_constrained_values[i] = std::pair<number, number>(
          src.local_element(edge_constrained_indices[i]),
          dst.local_element(edge_constrained_indices[i]));
        const_cast<LinearAlgebra::distributed::Vector<number> &>(src)
          .local_element(edge_constrained_indices[i]) = 0.;
      }

    data.cell_loop(&LaplaceOperator::local_apply, this, dst, src);

    const std::vector<unsigned int> &constrained_dofs =
      data.get_constrained_dofs();
    for (unsigned int i = 0; i < constrained_dofs.size(); ++i)
      dst.local_element(constrained_dofs[i]) +=
        src.local_element(constrained_dofs[i]);

    // reset edge constrained values, multiply by unit matrix and add into
    // destination
    for (unsigned int i = 0; i < edge_constrained_indices.size(); ++i)
      {
        const_cast<LinearAlgebra::distributed::Vector<number> &>(src)
          .local_element(edge_constrained_indices[i]) =
          edge_constrained_values[i].first;
        dst.local_element(edge_constrained_indices[i]) =
          edge_constrained_values[i].second + edge_constrained_values[i].first;
      }
  }

  void
  vmult_interface_down(
    LinearAlgebra::distributed::Vector<number>       &dst,
    const LinearAlgebra::distributed::Vector<number> &src) const
  {
    Assert(src.partitioners_are_globally_compatible(
             *data.get_dof_info(0).vector_partitioner),
           ExcInternalError());
    Assert(dst.partitioners_are_globally_compatible(
             *data.get_dof_info(0).vector_partitioner),
           ExcInternalError());

    dst = 0;

    if (!have_interface_matrices)
      return;

    // set zero Dirichlet values on the input vector (and remember the src and
    // dst values because we need to reset them at the end)
    for (unsigned int i = 0; i < edge_constrained_indices.size(); ++i)
      {
        const double src_val = src.local_element(edge_constrained_indices[i]);
        const_cast<LinearAlgebra::distributed::Vector<number> &>(src)
          .local_element(edge_constrained_indices[i]) = 0.;
        edge_constrained_values[i] = std::pair<number, number>(
          src_val, dst.local_element(edge_constrained_indices[i]));
      }

    data.cell_loop(&LaplaceOperator::local_apply, this, dst, src);

    unsigned int c = 0;
    for (unsigned int i = 0; i < edge_constrained_indices.size(); ++i)
      {
        for (; c < edge_constrained_indices[i]; ++c)
          dst.local_element(c) = 0.;
        ++c;

        // reset the src values
        const_cast<LinearAlgebra::distributed::Vector<number> &>(src)
          .local_element(edge_constrained_indices[i]) =
          edge_constrained_values[i].first;
      }
    for (; c < dst.locally_owned_size(); ++c)
      dst.local_element(c) = 0.;
  }

  void
  vmult_interface_up(
    LinearAlgebra::distributed::Vector<number>       &dst,
    const LinearAlgebra::distributed::Vector<number> &src) const
  {
    Assert(src.partitioners_are_globally_compatible(
             *data.get_dof_info(0).vector_partitioner),
           ExcInternalError());
    Assert(dst.partitioners_are_globally_compatible(
             *data.get_dof_info(0).vector_partitioner),
           ExcInternalError());

    dst = 0;

    if (!have_interface_matrices)
      return;

    LinearAlgebra::distributed::Vector<number> src_cpy(src);
    unsigned int                               c = 0;
    for (unsigned int i = 0; i < edge_constrained_indices.size(); ++i)
      {
        for (; c < edge_constrained_indices[i]; ++c)
          src_cpy.local_element(c) = 0.;
        ++c;
      }
    for (; c < src_cpy.locally_owned_size(); ++c)
      src_cpy.local_element(c) = 0.;

    data.cell_loop(&LaplaceOperator::local_apply, this, dst, src_cpy);
    for (unsigned int i = 0; i < edge_constrained_indices.size(); ++i)
      {
        dst.local_element(edge_constrained_indices[i]) = 0.;
      }
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

  number
  el(const unsigned int row, const unsigned int col) const
  {
    AssertThrow(false,
                ExcMessage("Matrix-free does not allow for entry access"));
    return number();
  }

  void
  initialize_dof_vector(
    LinearAlgebra::distributed::Vector<number> &vector) const
  {
    if (!vector.partitioners_are_compatible(
          *data.get_dof_info(0).vector_partitioner))
      data.initialize_dof_vector(vector);
    Assert(vector.partitioners_are_globally_compatible(
             *data.get_dof_info(0).vector_partitioner),
           ExcInternalError());
  }

  const LinearAlgebra::distributed::Vector<number> &
  get_matrix_diagonal_inverse() const
  {
    Assert(inverse_diagonal_entries.size() > 0, ExcNotInitialized());
    return inverse_diagonal_entries;
  }

  const std::shared_ptr<const Utilities::MPI::Partitioner> &
  get_vector_partitioner() const
  {
    return data.get_vector_partitioner();
  }


private:
  void
  local_apply(const MatrixFree<dim, number>                    &data,
              LinearAlgebra::distributed::Vector<number>       &dst,
              const LinearAlgebra::distributed::Vector<number> &src,
              const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    FEEvaluation<dim, fe_degree, n_q_points_1d, 1, number> phi(data);

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        phi.reinit(cell);
        phi.read_dof_values(src);
        phi.evaluate(EvaluationFlags::gradients);
        for (unsigned int q = 0; q < phi.n_q_points; ++q)
          phi.submit_gradient(phi.get_gradient(q), q);
        phi.integrate(EvaluationFlags::gradients);
        phi.distribute_local_to_global(dst);
      }
  }

  void
  compute_inverse_diagonal()
  {
    data.initialize_dof_vector(inverse_diagonal_entries);
    unsigned int dummy = 0;
    data.cell_loop(&LaplaceOperator::local_diagonal_cell,
                   this,
                   inverse_diagonal_entries,
                   dummy);

    const std::vector<unsigned int> &constrained_dofs =
      data.get_constrained_dofs();
    for (unsigned int i = 0; i < constrained_dofs.size(); ++i)
      inverse_diagonal_entries.local_element(constrained_dofs[i]) = 1.;
    for (unsigned int i = 0; i < edge_constrained_indices.size(); ++i)
      {
        inverse_diagonal_entries.local_element(edge_constrained_indices[i]) =
          1.;
      }


    for (unsigned int i = 0; i < inverse_diagonal_entries.locally_owned_size();
         ++i)
      if (std::abs(inverse_diagonal_entries.local_element(i)) > 1e-10)
        inverse_diagonal_entries.local_element(i) =
          1. / inverse_diagonal_entries.local_element(i);
      else
        inverse_diagonal_entries.local_element(i) = 1.;
  }

  void
  local_diagonal_cell(
    const MatrixFree<dim, number>              &data,
    LinearAlgebra::distributed::Vector<number> &dst,
    const unsigned int &,
    const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    FEEvaluation<dim, fe_degree, n_q_points_1d, 1, number> phi(data);

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        phi.reinit(cell);

        VectorizedArray<number> local_diagonal_vector[phi.tensor_dofs_per_cell];
        for (unsigned int i = 0; i < phi.dofs_per_cell; ++i)
          {
            for (unsigned int j = 0; j < phi.dofs_per_cell; ++j)
              phi.begin_dof_values()[j] = VectorizedArray<number>();
            phi.begin_dof_values()[i] = 1.;
            phi.evaluate(EvaluationFlags::gradients);
            for (unsigned int q = 0; q < phi.n_q_points; ++q)
              phi.submit_gradient(phi.get_gradient(q), q);
            phi.integrate(EvaluationFlags::gradients);
            local_diagonal_vector[i] = phi.begin_dof_values()[i];
          }
        for (unsigned int i = 0; i < phi.tensor_dofs_per_cell; ++i)
          phi.begin_dof_values()[i] = local_diagonal_vector[i];
        phi.distribute_local_to_global(dst);
      }
  }

  MatrixFree<dim, number>                        data;
  LinearAlgebra::distributed::Vector<number>     inverse_diagonal_entries;
  std::vector<unsigned int>                      edge_constrained_indices;
  mutable std::vector<std::pair<number, number>> edge_constrained_values;
  bool                                           have_interface_matrices;
};



template <typename LAPLACEOPERATOR>
class MGInterfaceMatrix : public EnableObserverPointer
{
public:
  void
  initialize(const LAPLACEOPERATOR &laplace)
  {
    this->laplace = &laplace;
  }

  void
  vmult(LinearAlgebra::distributed::Vector<typename LAPLACEOPERATOR::value_type>
          &dst,
        const LinearAlgebra::distributed::Vector<
          typename LAPLACEOPERATOR::value_type> &src) const
  {
    laplace->vmult_interface_down(dst, src);
  }

  void
  Tvmult(
    LinearAlgebra::distributed::Vector<typename LAPLACEOPERATOR::value_type>
      &dst,
    const LinearAlgebra::distributed::Vector<
      typename LAPLACEOPERATOR::value_type> &src) const
  {
    laplace->vmult_interface_up(dst, src);
  }

private:
  ObserverPointer<const LAPLACEOPERATOR> laplace;
};



template <typename MatrixType, typename Number>
class MGCoarseIterative
  : public MGCoarseGridBase<LinearAlgebra::distributed::Vector<Number>>
{
public:
  MGCoarseIterative()
  {}

  void
  initialize(const MatrixType &matrix)
  {
    coarse_matrix = &matrix;
  }

  virtual void
  operator()(const unsigned int                                level,
             LinearAlgebra::distributed::Vector<Number>       &dst,
             const LinearAlgebra::distributed::Vector<Number> &src) const
  {
    ReductionControl solver_control(1e4, 1e-50, 1e-10);
    SolverCG<LinearAlgebra::distributed::Vector<Number>> solver_coarse(
      solver_control);
    solver_coarse.solve(*coarse_matrix, dst, src, PreconditionIdentity());
  }

  const MatrixType *coarse_matrix;
};



template <int dim, int fe_degree, int n_q_points_1d, typename number>
void
do_test(const DoFHandler<dim> &dof, const bool threaded)
{
  deallog << "Testing " << dof.get_fe().get_name();
  deallog << std::endl;
  deallog << "Number of degrees of freedom: " << dof.n_dofs() << std::endl;

  AffineConstraints<double> hanging_node_constraints;
  const IndexSet            locally_relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(dof);
  hanging_node_constraints.reinit(dof.locally_owned_dofs(),
                                  locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints(dof, hanging_node_constraints);
  hanging_node_constraints.close();

  MGConstrainedDoFs                                   mg_constrained_dofs;
  Functions::ZeroFunction<dim>                        zero_function;
  std::map<types::boundary_id, const Function<dim> *> dirichlet_boundary;
  dirichlet_boundary[0] = &zero_function;
  mg_constrained_dofs.initialize(dof);
  mg_constrained_dofs.make_zero_boundary_constraints(dof, {0});

  MappingQ<dim>                                          mapping(fe_degree + 1);
  LaplaceOperator<dim, fe_degree, n_q_points_1d, double> fine_matrix;
  fine_matrix.initialize(mapping,
                         dof,
                         mg_constrained_dofs,
                         dirichlet_boundary,
                         numbers::invalid_unsigned_int,
                         threaded);

  LinearAlgebra::distributed::Vector<double> in, sol;
  fine_matrix.initialize_dof_vector(in);
  fine_matrix.initialize_dof_vector(sol);

  // set constant rhs vector
  for (unsigned int i = 0; i < in.locally_owned_size(); ++i)
    if (!hanging_node_constraints.is_constrained(
          in.get_partitioner()->local_to_global(i)))
      in.local_element(i) = 1.;

  // set up multigrid in analogy to step-37
  using LevelMatrixType =
    LaplaceOperator<dim, fe_degree, n_q_points_1d, number>;

  MGLevelObject<LevelMatrixType> mg_matrices;
  mg_matrices.resize(0, dof.get_triangulation().n_global_levels() - 1);
  for (unsigned int level = 0;
       level < dof.get_triangulation().n_global_levels();
       ++level)
    {
      mg_matrices[level].initialize(
        mapping, dof, mg_constrained_dofs, dirichlet_boundary, level, threaded);
    }
  MGLevelObject<MGInterfaceMatrix<LevelMatrixType>> mg_interface_matrices;
  mg_interface_matrices.resize(0,
                               dof.get_triangulation().n_global_levels() - 1);
  for (unsigned int level = 0;
       level < dof.get_triangulation().n_global_levels();
       ++level)
    mg_interface_matrices[level].initialize(mg_matrices[level]);

  std::vector<std::shared_ptr<const Utilities::MPI::Partitioner>> partitioners;
  for (unsigned int level = mg_matrices.min_level();
       level <= mg_matrices.max_level();
       ++level)
    partitioners.push_back(mg_matrices[level].get_vector_partitioner());

  MGTransferMF<dim, double> mg_transfer(mg_constrained_dofs);
  mg_transfer.build(dof, partitioners);

  MGCoarseIterative<LevelMatrixType, number> mg_coarse;
  mg_coarse.initialize(mg_matrices[0]);

  using SMOOTHER =
    PreconditionChebyshev<LevelMatrixType,
                          LinearAlgebra::distributed::Vector<number>>;
  MGSmootherPrecondition<LevelMatrixType,
                         SMOOTHER,
                         LinearAlgebra::distributed::Vector<number>>
    mg_smoother;

  MGLevelObject<typename SMOOTHER::AdditionalData> smoother_data;
  smoother_data.resize(0, dof.get_triangulation().n_global_levels() - 1);
  for (unsigned int level = 0;
       level < dof.get_triangulation().n_global_levels();
       ++level)
    {
      smoother_data[level].smoothing_range     = 15.;
      smoother_data[level].degree              = 5;
      smoother_data[level].eig_cg_n_iterations = 15;
      smoother_data[level].preconditioner.reset(
        new DiagonalMatrix<LinearAlgebra::distributed::Vector<number>>());
      smoother_data[level].preconditioner->get_vector() =
        mg_matrices[level].get_matrix_diagonal_inverse();
    }

  mg_smoother.initialize(mg_matrices, smoother_data);

  mg::Matrix<LinearAlgebra::distributed::Vector<number>> mg_matrix(mg_matrices);
  mg::Matrix<LinearAlgebra::distributed::Vector<number>> mg_interface(
    mg_interface_matrices);

  Multigrid<LinearAlgebra::distributed::Vector<number>> mg(
    mg_matrix, mg_coarse, mg_transfer, mg_smoother, mg_smoother);
  mg.set_edge_matrices(mg_interface, mg_interface);
  PreconditionMG<dim,
                 LinearAlgebra::distributed::Vector<number>,
                 MGTransferMF<dim, double>>
    preconditioner(dof, mg, mg_transfer);

  {
    // avoid output from inner (coarse-level) solver
    deallog.depth_file(3);

    ReductionControl control(30, 1e-20, 1e-7);
    SolverCG<LinearAlgebra::distributed::Vector<double>> solver(control);
    solver.solve(fine_matrix, sol, in, preconditioner);
  }
}



template <int dim, int fe_degree, typename Number>
void
test()
{
  parallel::distributed::Triangulation<dim> tria(
    MPI_COMM_WORLD,
    Triangulation<dim>::limit_level_difference_at_vertices,
    parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(8 - 2 * dim);
  const unsigned int n_runs = fe_degree == 1 ? 6 - dim : 5 - dim;
  for (unsigned int i = 0; i < n_runs; ++i)
    {
      for (typename Triangulation<dim>::active_cell_iterator cell =
             tria.begin_active();
           cell != tria.end();
           ++cell)
        if (cell->is_locally_owned() &&
            ((cell->center().norm() < 0.5 &&
              (cell->level() < 5 || cell->center().norm() > 0.45)) ||
             (dim == 2 && cell->center().norm() > 1.2)))
          cell->set_refine_flag();
      tria.execute_coarsening_and_refinement();
      FE_Q<dim>       fe(fe_degree);
      DoFHandler<dim> dof(tria);
      dof.distribute_dofs(fe);
      dof.distribute_mg_dofs();

      deallog.push("threaded");
      do_test<dim, fe_degree, fe_degree + 1, Number>(dof, true);
      deallog.pop();
    }
}



int
main(int argc, char **argv)
{
  // The original issue with partition_color
  // is hit with 2 threads and 4 cores.
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 2);

  mpi_initlog();

  test<2, 1, double>();
}
