// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// same as parallel_multigrid_adaptive_02, but with Jacobi smoother.
// The only reason for this test is to provide a reference for
// parallel_multigrid_adaptive_06

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
#include <deal.II/multigrid/mg_transfer_global_coarsening.h>
#include <deal.II/multigrid/multigrid.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


using namespace dealii::MatrixFreeOperators;

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
             LinearAlgebra::distributed::Vector<double>       &dst,
             const LinearAlgebra::distributed::Vector<double> &src) const
  {
    ReductionControl solver_control(1e4, 1e-50, 1e-10);
    SolverCG<LinearAlgebra::distributed::Vector<double>> solver_coarse(
      solver_control);
    solver_coarse.solve(*coarse_matrix, dst, src, PreconditionIdentity());
  }

  const MatrixType *coarse_matrix;
};



template <int dim, int fe_degree, int n_q_points_1d, typename number>
void
do_test(const DoFHandler<dim> &dof)
{
  if (std::is_same_v<number, float> == true)
    {
      deallog.push("float");
    }
  else
    {
    }

  deallog << "Testing " << dof.get_fe().get_name();
  deallog << std::endl;
  deallog << "Number of degrees of freedom: " << dof.n_dofs() << std::endl;

  const IndexSet locally_relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(dof);

  // Dirichlet BC
  Functions::ZeroFunction<dim>                        zero_function;
  std::map<types::boundary_id, const Function<dim> *> dirichlet_boundary;
  dirichlet_boundary[0] = &zero_function;

  // fine-level constraints
  AffineConstraints<double> constraints;
  constraints.reinit(dof.locally_owned_dofs(), locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints(dof, constraints);
  VectorTools::interpolate_boundary_values(dof,
                                           dirichlet_boundary,
                                           constraints);
  constraints.close();

  // level constraints:
  MGConstrainedDoFs mg_constrained_dofs;
  mg_constrained_dofs.initialize(dof);
  mg_constrained_dofs.make_zero_boundary_constraints(dof, {0});

  MappingQ<dim> mapping(fe_degree + 1);

  LaplaceOperator<dim,
                  fe_degree,
                  n_q_points_1d,
                  1,
                  LinearAlgebra::distributed::Vector<number>>
                                           fine_matrix;
  std::shared_ptr<MatrixFree<dim, number>> fine_level_data(
    new MatrixFree<dim, number>());

  typename MatrixFree<dim, number>::AdditionalData fine_level_additional_data;
  fine_level_additional_data.tasks_parallel_scheme =
    MatrixFree<dim, number>::AdditionalData::none;
  fine_level_additional_data.tasks_block_size = 3;
  fine_level_data->reinit(mapping,
                          dof,
                          constraints,
                          QGauss<1>(n_q_points_1d),
                          fine_level_additional_data);

  fine_matrix.initialize(fine_level_data);
  fine_matrix.compute_diagonal();


  LinearAlgebra::distributed::Vector<number> in, sol;
  fine_matrix.initialize_dof_vector(in);
  fine_matrix.initialize_dof_vector(sol);

  // set constant rhs vector
  {
    // this is to make it consistent with parallel_multigrid_adaptive.cc
    AffineConstraints<double> hanging_node_constraints;
    hanging_node_constraints.reinit(dof.locally_owned_dofs(),
                                    locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints(dof, hanging_node_constraints);
    hanging_node_constraints.close();

    for (unsigned int i = 0; i < in.locally_owned_size(); ++i)
      if (!hanging_node_constraints.is_constrained(
            in.get_partitioner()->local_to_global(i)))
        in.local_element(i) = 1.;
  }

  // set up multigrid in analogy to step-37
  using LevelMatrixType =
    LaplaceOperator<dim,
                    fe_degree,
                    n_q_points_1d,
                    1,
                    LinearAlgebra::distributed::Vector<number>>;

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
      const IndexSet            relevant_dofs =
        DoFTools::extract_locally_relevant_level_dofs(dof, level);
      level_constraints.reinit(dof.locally_owned_mg_dofs(level), relevant_dofs);
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

  std::vector<std::shared_ptr<const Utilities::MPI::Partitioner>> partitioners;
  for (unsigned int level = mg_matrices.min_level();
       level <= mg_matrices.max_level();
       ++level)
    partitioners.push_back(mg_level_data[level].get_vector_partitioner());

  MGTransferGlobalCoarsening<dim, LinearAlgebra::distributed::Vector<double>>
    mg_transfer(mg_constrained_dofs);
  mg_transfer.build(dof, partitioners);

  MGCoarseIterative<LevelMatrixType, number> mg_coarse;
  mg_coarse.initialize(mg_matrices[0]);

  using SMOOTHER = PreconditionJacobi<LevelMatrixType>;
  MGSmootherPrecondition<LevelMatrixType,
                         SMOOTHER,
                         LinearAlgebra::distributed::Vector<double>>
    mg_smoother;

  mg_smoother.initialize(mg_matrices, typename SMOOTHER::AdditionalData(0.8));

  mg::Matrix<LinearAlgebra::distributed::Vector<double>> mg_matrix(mg_matrices);
  mg::Matrix<LinearAlgebra::distributed::Vector<double>> mg_interface(
    mg_interface_matrices);

  Multigrid<LinearAlgebra::distributed::Vector<double>> mg(
    mg_matrix, mg_coarse, mg_transfer, mg_smoother, mg_smoother);
  mg.set_edge_matrices(mg_interface, mg_interface);
  PreconditionMG<
    dim,
    LinearAlgebra::distributed::Vector<double>,
    MGTransferGlobalCoarsening<dim, LinearAlgebra::distributed::Vector<double>>>
    preconditioner(dof, mg, mg_transfer);

  {
    // avoid output from inner (coarse-level) solver
    deallog.depth_file(3);

    ReductionControl control(30, 1e-20, 1e-7);
    SolverCG<LinearAlgebra::distributed::Vector<double>> solver(control);
    solver.solve(fine_matrix, sol, in, preconditioner);
  }

  if (std::is_same_v<number, float> == true)
    deallog.pop();

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
  mpi_initlog();
  deallog << std::setprecision(4);

  {
    deallog.push("2d");
    test<2, 1>();
    test<2, 3>();
    deallog.pop();
    deallog.push("3d");
    test<3, 1>();
    test<3, 2>();
    deallog.pop();
  }
}
