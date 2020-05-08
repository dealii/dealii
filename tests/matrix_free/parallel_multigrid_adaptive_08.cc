// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2020 by the deal.II authors
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



// similar to parallel_multigrid_adaptive_02 but without using a separate
// transfer class that builds the appropriate vectors. Furthermore, we want to
// avoid setting some vectors to zero in the cell loop. Rather,
// MatrixFreeOperators::LaplaceOperator::adjust_ghost_range_if_necessary()
// will do that - or rather a variant of that, given that we want to modify
// the cell loop of the LaplaceOperator class and provide our own. This forces
// us to reimplement a few things, but all ideas are the same as in the
// matrix-free operators.

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
#include <deal.II/matrix_free/operators.h>

#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_matrix.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_transfer_matrix_free.h>
#include <deal.II/multigrid/multigrid.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


using namespace dealii::MatrixFreeOperators;


template <int dim, int fe_degree, typename Number>
class MyLaplaceOperator : public MatrixFreeOperators::LaplaceOperator<
                            dim,
                            fe_degree,
                            fe_degree + 1,
                            1,
                            LinearAlgebra::distributed::Vector<Number>>
{
public:
  void
  initialize(std::shared_ptr<const MatrixFree<dim, Number>> data,
             const MGConstrainedDoFs &                      mg_constrained_dofs,
             const unsigned int                             level)
  {
    MatrixFreeOperators::Base<dim, LinearAlgebra::distributed::Vector<Number>>::
      initialize(data,
                 mg_constrained_dofs,
                 level,
                 std::vector<unsigned int>({0}));

    std::vector<types::global_dof_index> interface_indices;
    mg_constrained_dofs.get_refinement_edge_indices(level).fill_index_vector(
      interface_indices);
    vmult_edge_constrained_indices.clear();
    vmult_edge_constrained_indices.reserve(interface_indices.size());
    vmult_edge_constrained_values.resize(interface_indices.size());
    const IndexSet &locally_owned =
      this->data->get_dof_handler(0).locally_owned_mg_dofs(level);
    for (unsigned int i = 0; i < interface_indices.size(); ++i)
      if (locally_owned.is_element(interface_indices[i]))
        vmult_edge_constrained_indices.push_back(
          locally_owned.index_within_set(interface_indices[i]));
  }

  void
  initialize(std::shared_ptr<const MatrixFree<dim, Number>> data,
             const std::vector<unsigned int> &              mask = {})
  {
    MatrixFreeOperators::Base<dim, LinearAlgebra::distributed::Vector<Number>>::
      initialize(data, mask);
  }


  void
  vmult(LinearAlgebra::distributed::Vector<Number> &      dst,
        const LinearAlgebra::distributed::Vector<Number> &src) const
  {
    adjust_ghost_range_if_necessary(src);
    adjust_ghost_range_if_necessary(dst);

    for (unsigned int i = 0; i < vmult_edge_constrained_indices.size(); ++i)
      {
        vmult_edge_constrained_values[i] =
          src.local_element(vmult_edge_constrained_indices[i]);
        const_cast<LinearAlgebra::distributed::Vector<Number> &>(src)
          .local_element(vmult_edge_constrained_indices[i]) = 0.;
      }

    // zero dst within the loop
    this->data->cell_loop(
      &MyLaplaceOperator::local_apply, this, dst, src, true);

    for (auto i : this->data->get_constrained_dofs(0))
      dst.local_element(i) = src.local_element(i);
    for (unsigned int i = 0; i < vmult_edge_constrained_indices.size(); ++i)
      {
        dst.local_element(vmult_edge_constrained_indices[i]) =
          vmult_edge_constrained_values[i];
        const_cast<LinearAlgebra::distributed::Vector<Number> &>(src)
          .local_element(vmult_edge_constrained_indices[i]) =
          vmult_edge_constrained_values[i];
      }
  }

private:
  void
  local_apply(const MatrixFree<dim, Number> &                   data,
              LinearAlgebra::distributed::Vector<Number> &      dst,
              const LinearAlgebra::distributed::Vector<Number> &src,
              const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    FEEvaluation<dim, fe_degree, fe_degree + 1, 1, Number> phi(data);

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        phi.reinit(cell);
        phi.gather_evaluate(src, false, true);
        for (unsigned int q = 0; q < phi.n_q_points; ++q)
          phi.submit_gradient(phi.get_gradient(q), q);
        phi.integrate_scatter(false, true, dst);
      }
  }

  void
  adjust_ghost_range_if_necessary(
    const LinearAlgebra::distributed::Vector<Number> &vec) const
  {
    if (vec.get_partitioner().get() ==
        this->data->get_dof_info(0).vector_partitioner.get())
      return;

    Assert(vec.get_partitioner()->local_size() ==
             this->data->get_dof_info(0).vector_partitioner->local_size(),
           ExcMessage("The vector passed to the vmult() function does not have "
                      "the correct size for compatibility with MatrixFree."));
    LinearAlgebra::distributed::Vector<Number> copy_vec(vec);
    const_cast<LinearAlgebra::distributed::Vector<Number> &>(vec).reinit(
      this->data->get_dof_info(0).vector_partitioner);
    const_cast<LinearAlgebra::distributed::Vector<Number> &>(vec)
      .copy_locally_owned_data_from(copy_vec);
  }

  std::vector<unsigned int> vmult_edge_constrained_indices;

  mutable std::vector<double> vmult_edge_constrained_values;
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
             LinearAlgebra::distributed::Vector<Number> &      dst,
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
do_test(const DoFHandler<dim> &dof)
{
  deallog << "Testing " << dof.get_fe().get_name();
  deallog << std::endl;
  deallog << "Number of degrees of freedom: " << dof.n_dofs() << std::endl;

  IndexSet locally_relevant_dofs;
  DoFTools::extract_locally_relevant_dofs(dof, locally_relevant_dofs);

  // Dirichlet BC
  Functions::ZeroFunction<dim>                        zero_function;
  std::map<types::boundary_id, const Function<dim> *> dirichlet_boundary;
  dirichlet_boundary[0] = &zero_function;

  // fine-level constraints
  AffineConstraints<double> constraints;
  constraints.reinit(locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints(dof, constraints);
  VectorTools::interpolate_boundary_values(dof,
                                           dirichlet_boundary,
                                           constraints);
  constraints.close();

  // level constraints:
  MGConstrainedDoFs mg_constrained_dofs;
  mg_constrained_dofs.initialize(dof, dirichlet_boundary);

  MappingQ<dim> mapping(fe_degree + 1);

  MyLaplaceOperator<dim, fe_degree, double> fine_matrix;
  std::shared_ptr<MatrixFree<dim, double>>  fine_level_data(
     new MatrixFree<dim, double>());

  typename MatrixFree<dim, double>::AdditionalData fine_level_additional_data;
  fine_level_additional_data.tasks_parallel_scheme =
    MatrixFree<dim, double>::AdditionalData::none;
  fine_level_additional_data.tasks_block_size = 3;
  fine_level_data->reinit(mapping,
                          dof,
                          constraints,
                          QGauss<1>(n_q_points_1d),
                          fine_level_additional_data);

  fine_matrix.initialize(fine_level_data);
  fine_matrix.compute_diagonal();


  LinearAlgebra::distributed::Vector<double> in, sol;
  fine_matrix.initialize_dof_vector(in);
  fine_matrix.initialize_dof_vector(sol);

  // set constant rhs vector
  {
    // this is to make it consistent with parallel_multigrid_adaptive.cc
    AffineConstraints<double> hanging_node_constraints;
    hanging_node_constraints.reinit(locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints(dof, hanging_node_constraints);
    hanging_node_constraints.close();

    for (unsigned int i = 0; i < in.local_size(); ++i)
      if (!hanging_node_constraints.is_constrained(
            in.get_partitioner()->local_to_global(i)))
        in.local_element(i) = 1.;
  }

  // set up multigrid in analogy to step-37
  using LevelMatrixType = MyLaplaceOperator<dim, fe_degree, number>;

  MGLevelObject<LevelMatrixType>         mg_matrices;
  MGLevelObject<MatrixFree<dim, number>> mg_level_data;
  mg_matrices.resize(0, dof.get_triangulation().n_global_levels() - 1);
  mg_level_data.resize(0, dof.get_triangulation().n_global_levels() - 1);
  for (unsigned int level = 0;
       level < dof.get_triangulation().n_global_levels();
       ++level)
    {
      typename MatrixFree<dim, number>::AdditionalData mg_additional_data;
      mg_additional_data.tasks_parallel_scheme =
        MatrixFree<dim, number>::AdditionalData::none;
      mg_additional_data.tasks_block_size = 3;
      mg_additional_data.mg_level         = level;

      AffineConstraints<double> level_constraints;
      IndexSet                  relevant_dofs;
      DoFTools::extract_locally_relevant_level_dofs(dof, level, relevant_dofs);
      level_constraints.reinit(relevant_dofs);
      level_constraints.add_lines(
        mg_constrained_dofs.get_boundary_indices(level));
      level_constraints.close();

      mg_level_data[level].reinit(mapping,
                                  dof,
                                  level_constraints,
                                  QGauss<1>(n_q_points_1d),
                                  mg_additional_data);
      mg_matrices[level].initialize(std::make_shared<MatrixFree<dim, number>>(
                                      mg_level_data[level]),
                                    mg_constrained_dofs,
                                    level);
      mg_matrices[level].compute_diagonal();
    }
  MGLevelObject<MGInterfaceOperator<LevelMatrixType>> mg_interface_matrices;
  mg_interface_matrices.resize(0,
                               dof.get_triangulation().n_global_levels() - 1);
  for (unsigned int level = 0;
       level < dof.get_triangulation().n_global_levels();
       ++level)
    mg_interface_matrices[level].initialize(mg_matrices[level]);

  MGTransferMatrixFree<dim, number> mg_transfer(mg_constrained_dofs);
  mg_transfer.build(dof);

  MGCoarseIterative<LevelMatrixType, number> mg_coarse;
  mg_coarse.initialize(mg_matrices[0]);

  typedef PreconditionChebyshev<LevelMatrixType,
                                LinearAlgebra::distributed::Vector<number>>
    SMOOTHER;
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
      smoother_data[level].degree              = 6;
      smoother_data[level].eig_cg_n_iterations = 15;
      smoother_data[level].preconditioner =
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
                 MGTransferMatrixFree<dim, number>>
    preconditioner(dof, mg, mg_transfer);

  {
    // avoid output from inner (coarse-level) solver
    deallog.depth_file(3);

    ReductionControl control(30, 1e-20, 1e-7);
    SolverCG<LinearAlgebra::distributed::Vector<double>> solver(control);
    solver.solve(fine_matrix, sol, in, preconditioner);
  }

  fine_matrix.clear();
  for (unsigned int level = 0;
       level < dof.get_triangulation().n_global_levels();
       ++level)
    mg_matrices[level].clear();
}



template <int dim, int fe_degree>
void
test()
{
  parallel::distributed::Triangulation<dim> tria(
    MPI_COMM_WORLD,
    Triangulation<dim>::limit_level_difference_at_vertices,
    parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(6 - dim);
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

      do_test<dim, fe_degree, fe_degree + 1, double>(dof);
    }
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);
  mpi_initlog(true);

  test<2, 2>();
}
