// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2021 by the deal.II authors
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


// same test as parallel_multigrid_mf but using the Laplace operator from
// MatrixFreeOperators::Laplace and handing the vector partitioning of those
// operators to the MG transfer

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

  MappingQ<dim> mapping(fe_degree + 1);
  MatrixFreeOperators::LaplaceOperator<
    dim,
    fe_degree,
    n_q_points_1d,
    1,
    LinearAlgebra::distributed::Vector<double>>
    fine_matrix;

  // Dirichlet BC
  Functions::ZeroFunction<dim>                        zero_function;
  std::map<types::boundary_id, const Function<dim> *> dirichlet_boundary;
  dirichlet_boundary[0] = &zero_function;

  // fine-level constraints
  IndexSet locally_relevant_dofs;
  DoFTools::extract_locally_relevant_dofs(dof, locally_relevant_dofs);
  AffineConstraints<double> constraints;
  constraints.reinit(locally_relevant_dofs);
  VectorTools::interpolate_boundary_values(dof,
                                           dirichlet_boundary,
                                           constraints);
  constraints.close();

  std::shared_ptr<MatrixFree<dim, double>> fine_level_data(
    new MatrixFree<dim, double>());

  typename MatrixFree<dim, double>::AdditionalData fine_level_additional_data;
  fine_level_data->reinit(mapping,
                          dof,
                          constraints,
                          QGauss<1>(n_q_points_1d),
                          fine_level_additional_data);

  fine_matrix.initialize(fine_level_data);

  LinearAlgebra::distributed::Vector<double> in, sol;
  fine_matrix.initialize_dof_vector(in);
  fine_matrix.initialize_dof_vector(sol);

  // set constant rhs vector
  in = 1.;

  // set up multigrid in analogy to step-37
  MGConstrainedDoFs mg_constrained_dofs;
  mg_constrained_dofs.initialize(dof);
  mg_constrained_dofs.make_zero_boundary_constraints(dof, {0});

  using LevelMatrixType = MatrixFreeOperators::LaplaceOperator<
    dim,
    fe_degree,
    n_q_points_1d,
    1,
    LinearAlgebra::distributed::Vector<number>>;

  MGLevelObject<LevelMatrixType>                          mg_matrices;
  MGLevelObject<std::shared_ptr<MatrixFree<dim, number>>> mg_level_data;
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

      mg_level_data[level].reset(new MatrixFree<dim, number>());
      mg_level_data[level]->reinit(mapping,
                                   dof,
                                   level_constraints,
                                   QGauss<1>(n_q_points_1d),
                                   mg_additional_data);
      mg_matrices[level].initialize(mg_level_data[level],
                                    mg_constrained_dofs,
                                    level);
      mg_matrices[level].compute_diagonal();
    }

  MGTransferMF<dim, number>                                       mg_transfer;
  std::vector<std::shared_ptr<const Utilities::MPI::Partitioner>> partitioners(
    dof.get_triangulation().n_global_levels());
  for (unsigned int level = 0; level < partitioners.size(); ++level)
    partitioners[level] =
      mg_level_data[level]->get_dof_info().vector_partitioner;
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
      smoother_data[level].preconditioner =
        mg_matrices[level].get_matrix_diagonal_inverse();
    }

  // temporarily disable deallog for the setup of the preconditioner that
  // involves a CG solver for eigenvalue estimation
  mg_smoother.initialize(mg_matrices, smoother_data);

  mg::Matrix<LinearAlgebra::distributed::Vector<number>> mg_matrix(mg_matrices);

  Multigrid<LinearAlgebra::distributed::Vector<number>> mg(
    mg_matrix, mg_coarse, mg_transfer, mg_smoother, mg_smoother);
  PreconditionMG<dim,
                 LinearAlgebra::distributed::Vector<number>,
                 MGTransferMF<dim, number>>
    preconditioner(dof, mg, mg_transfer);

  {
    // avoid output from inner (coarse-level) solver
    deallog.depth_file(2);
    ReductionControl control(30, 1e-20, 1e-7);
    SolverCG<LinearAlgebra::distributed::Vector<double>> solver(control);
    solver.solve(fine_matrix, sol, in, preconditioner);
  }
}



template <int dim, int fe_degree, typename number>
void
test()
{
  for (unsigned int i = 5; i < 7; ++i)
    {
      parallel::distributed::Triangulation<dim> tria(
        MPI_COMM_WORLD,
        Triangulation<dim>::limit_level_difference_at_vertices,
        parallel::distributed::Triangulation<
          dim>::construct_multigrid_hierarchy);
      GridGenerator::hyper_cube(tria);
      tria.refine_global(i - dim);

      FE_Q<dim>       fe(fe_degree);
      DoFHandler<dim> dof(tria);
      dof.distribute_dofs(fe);
      dof.distribute_mg_dofs();

      do_test<dim, fe_degree, fe_degree + 1, number>(dof);
    }
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);

  mpi_initlog();

  {
    test<2, 1, double>();
    test<2, 2, float>();

    test<3, 1, double>();
    test<3, 2, float>();
  }
}
