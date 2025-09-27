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



// This test is similar to parallel_multigrid_adaptive_06 but we also test
// for different polynomial degree in different blocks.
// We expect to have the same iteration numbers as in
// parallel_multigrid_adaptive_06 with respect to the highest polynomial
// degree used.


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
#include <deal.II/multigrid/mg_transfer.h>
#include <deal.II/multigrid/mg_transfer_matrix_free.h>
#include <deal.II/multigrid/multigrid.h>

#include <deal.II/numerics/vector_tools.h>

#include <algorithm>

#include "../tests.h"


using namespace dealii::MatrixFreeOperators;

/**
 * Block-Laplace operator with Block vector.
 */
template <int dim,
          int fe_degree_1,
          int fe_degree_2,
          int n_q_points_1d,
          typename BlockVectorType =
            LinearAlgebra::distributed::BlockVector<double>>
class BlockLaplace : public EnableObserverPointer
{
public:
  using value_type = typename BlockVectorType::value_type;
  using size_type  = typename BlockVectorType::size_type;

  BlockLaplace()
    : EnableObserverPointer()
  {}

  void
  initialize(std::shared_ptr<const MatrixFree<dim, value_type>> data)
  {
    laplace1.initialize(data, std::vector<unsigned int>(1, 0));
    laplace2.initialize(data, std::vector<unsigned int>(1, 1));
  }

  void
  initialize(std::shared_ptr<const MatrixFree<dim, value_type>> data,
             const std::vector<MGConstrainedDoFs> &mg_constrained_dofs,
             const unsigned int                    level)
  {
    laplace1.initialize(data,
                        mg_constrained_dofs[0],
                        level,
                        std::vector<unsigned int>(1, 0));
    laplace2.initialize(data,
                        mg_constrained_dofs[1],
                        level,
                        std::vector<unsigned int>(1, 1));
  }

  void
  vmult_interface_down(BlockVectorType &dst, const BlockVectorType &src) const
  {
    laplace1.vmult_interface_down(dst.block(0), src.block(0));
    laplace2.vmult_interface_down(dst.block(1), src.block(1));
  }

  void
  vmult_interface_up(BlockVectorType &dst, const BlockVectorType &src) const
  {
    laplace1.vmult_interface_up(dst.block(0), src.block(0));
    laplace2.vmult_interface_up(dst.block(1), src.block(1));
  }

  void
  vmult(BlockVectorType &dst, const BlockVectorType &src) const
  {
    laplace1.vmult(dst.block(0), src.block(0));
    laplace2.vmult(dst.block(1), src.block(1));
  }

  void
  Tvmult(BlockVectorType &dst, const BlockVectorType &src) const
  {
    laplace1.Tvmult(dst.block(0), src.block(0));
    laplace2.Tvmult(dst.block(1), src.block(1));
  }

  void
  vmult_add(BlockVectorType &dst, const BlockVectorType &src) const
  {
    laplace1.vmult_add(dst.block(0), src.block(0));
    laplace2.vmult_add(dst.block(1), src.block(1));
  }

  void
  Tvmult_add(BlockVectorType &dst, const BlockVectorType &src) const
  {
    laplace1.Tvmult_add(dst.block(0), src.block(0));
    laplace2.Tvmult_add(dst.block(1), src.block(1));
  }

  void
  precondition_Jacobi(BlockVectorType       &dst,
                      const BlockVectorType &src,
                      const value_type       omega) const
  {
    laplace1.precondition_Jacobi(dst.block(0), src.block(0), omega);
    laplace2.precondition_Jacobi(dst.block(1), src.block(1), omega);
  }

  void
  compute_diagonal()
  {
    laplace1.compute_diagonal();
    laplace2.compute_diagonal();
  }


  virtual void
  clear()
  {
    laplace1.clear();
    laplace2.clear();
  }

private:
  MatrixFreeOperators::LaplaceOperator<dim,
                                       fe_degree_1,
                                       n_q_points_1d,
                                       1,
                                       typename BlockVectorType::BlockType>
    laplace1;
  MatrixFreeOperators::LaplaceOperator<dim,
                                       fe_degree_2,
                                       n_q_points_1d,
                                       1,
                                       typename BlockVectorType::BlockType>
    laplace2;
};



template <typename MatrixType, typename Number>
class MGCoarseIterative
  : public MGCoarseGridBase<LinearAlgebra::distributed::BlockVector<Number>>
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
  operator()(const unsigned int                                     level,
             LinearAlgebra::distributed::BlockVector<Number>       &dst,
             const LinearAlgebra::distributed::BlockVector<Number> &src) const
  {
    ReductionControl solver_control(1e4, 1e-50, 1e-10);
    SolverCG<LinearAlgebra::distributed::BlockVector<Number>> solver_coarse(
      solver_control);
    solver_coarse.solve(*coarse_matrix, dst, src, PreconditionIdentity());
  }

  const MatrixType *coarse_matrix;
};



template <int dim,
          int fe_degree_1,
          int fe_degree_2,
          int n_q_points_1d,
          typename number>
void
do_test(const std::vector<const DoFHandler<dim> *> &dof)
{
  const unsigned int nb = 2;
  if (std::is_same_v<number, float> == true)
    {
      deallog.push("float");
    }
  else
    {
    }

  for (unsigned int i = 0; i < dof.size(); ++i)
    {
      deallog << "Testing " << dof[i]->get_fe().get_name();
      deallog << std::endl;
      deallog << "Number of degrees of freedom: " << dof[i]->n_dofs()
              << std::endl;
    }

  std::vector<IndexSet> locally_relevant_dofs(dof.size());
  for (unsigned int i = 0; i < dof.size(); ++i)
    locally_relevant_dofs[i] = DoFTools::extract_locally_relevant_dofs(*dof[i]);

  // Dirichlet BC
  Functions::ZeroFunction<dim>                        zero_function;
  std::map<types::boundary_id, const Function<dim> *> dirichlet_boundary;
  dirichlet_boundary[0] = &zero_function;

  // fine-level constraints
  std::vector<AffineConstraints<double>>         constraints(dof.size());
  std::vector<const AffineConstraints<double> *> constraints_ptrs(dof.size());
  for (unsigned int i = 0; i < dof.size(); ++i)
    {
      constraints[i].reinit(dof[i]->locally_owned_dofs(),
                            locally_relevant_dofs[i]);
      DoFTools::make_hanging_node_constraints(*dof[i], constraints[i]);
      VectorTools::interpolate_boundary_values(*dof[i],
                                               dirichlet_boundary,
                                               constraints[i]);
      constraints[i].close();
      constraints_ptrs[i] = &constraints[i];
    }
  QGauss<1>     quad(n_q_points_1d);
  constexpr int max_degree = std::max(fe_degree_1, fe_degree_2);
  MappingQ<dim> mapping(max_degree);

  typename MatrixFree<dim, number>::AdditionalData fine_level_additional_data;
  fine_level_additional_data.tasks_parallel_scheme =
    MatrixFree<dim, number>::AdditionalData::none;
  fine_level_additional_data.tasks_block_size = 3;

  std::shared_ptr<MatrixFree<dim, double>> fine_level_data(
    new MatrixFree<dim, double>());
  fine_level_data->reinit(
    mapping, dof, constraints_ptrs, quad, fine_level_additional_data);

  BlockLaplace<dim,
               fe_degree_1,
               fe_degree_2,
               n_q_points_1d,
               LinearAlgebra::distributed::BlockVector<number>>
    fine_matrix;

  fine_matrix.initialize(fine_level_data);
  fine_matrix.compute_diagonal();

  LinearAlgebra::distributed::BlockVector<number> in(nb), sol(nb);
  for (unsigned int b = 0; b < nb; ++b)
    {
      fine_level_data->initialize_dof_vector(in.block(b), b);
      fine_level_data->initialize_dof_vector(sol.block(b), b);
    }

  in.collect_sizes();
  sol.collect_sizes();

  // set constant rhs vector
  {
    for (unsigned int b = 0; b < nb; ++b)
      {
        // this is to make it consistent with parallel_multigrid_adaptive.cc
        AffineConstraints<double> hanging_node_constraints;

        hanging_node_constraints.reinit(dof[b]->locally_owned_dofs(),
                                        locally_relevant_dofs[b]);
        DoFTools::make_hanging_node_constraints(*dof[b],
                                                hanging_node_constraints);
        hanging_node_constraints.close();

        for (unsigned int i = 0; i < in.block(b).locally_owned_size(); ++i)
          if (!hanging_node_constraints.is_constrained(
                in.block(b).get_partitioner()->local_to_global(i)))
            in.block(b).local_element(i) = 1.;
      }
  }

  // level constraints:
  std::vector<MGConstrainedDoFs> mg_constrained_dofs(dof.size());
  for (unsigned int i = 0; i < dof.size(); ++i)
    {
      mg_constrained_dofs[i].initialize(*dof[i]);
      mg_constrained_dofs[i].make_zero_boundary_constraints(*dof[i], {0});
    }

  // set up multigrid in analogy to step-37
  using LevelMatrixType =
    BlockLaplace<dim,
                 fe_degree_1,
                 fe_degree_2,
                 n_q_points_1d,
                 LinearAlgebra::distributed::BlockVector<number>>;

  MGLevelObject<LevelMatrixType>         mg_matrices;
  MGLevelObject<MatrixFree<dim, number>> mg_level_data;
  mg_matrices.resize(0, dof[0]->get_triangulation().n_global_levels() - 1);
  mg_level_data.resize(0, dof[0]->get_triangulation().n_global_levels() - 1);
  for (unsigned int level = 0;
       level < dof[0]->get_triangulation().n_global_levels();
       ++level)
    {
      typename MatrixFree<dim, number>::AdditionalData mg_additional_data;
      mg_additional_data.tasks_parallel_scheme =
        MatrixFree<dim, number>::AdditionalData::none;
      mg_additional_data.tasks_block_size = 3;
      mg_additional_data.mg_level         = level;

      std::vector<AffineConstraints<double>> level_constraints(dof.size());
      std::vector<const AffineConstraints<double> *> level_constraints_ptrs(
        dof.size());
      for (unsigned int i = 0; i < dof.size(); ++i)
        {
          const IndexSet relevant_dofs =
            DoFTools::extract_locally_relevant_level_dofs(*dof[i], level);
          level_constraints[i].reinit(dof[i]->locally_owned_mg_dofs(level),
                                      relevant_dofs);
          level_constraints[i].add_lines(
            mg_constrained_dofs[i].get_boundary_indices(level));
          level_constraints[i].close();
          level_constraints_ptrs[i] = &level_constraints[i];
        }

      mg_level_data[level].reinit(
        mapping, dof, level_constraints_ptrs, quad, mg_additional_data);
      mg_matrices[level].initialize(std::make_shared<MatrixFree<dim, number>>(
                                      mg_level_data[level]),
                                    mg_constrained_dofs,
                                    level);
      mg_matrices[level].compute_diagonal();
    }

  MGLevelObject<MGInterfaceOperator<LevelMatrixType>> mg_interface_matrices;
  mg_interface_matrices.resize(0,
                               dof[0]->get_triangulation().n_global_levels() -
                                 1);
  for (unsigned int level = 0;
       level < dof[0]->get_triangulation().n_global_levels();
       ++level)
    mg_interface_matrices[level].initialize(mg_matrices[level]);

  MGTransferBlockMatrixFree<dim, number> mg_transfer(mg_constrained_dofs);
  mg_transfer.build(dof);

  MGCoarseIterative<LevelMatrixType, number> mg_coarse;
  mg_coarse.initialize(mg_matrices[0]);

  using SMOOTHER = PreconditionJacobi<LevelMatrixType>;
  MGSmootherPrecondition<LevelMatrixType,
                         SMOOTHER,
                         LinearAlgebra::distributed::BlockVector<number>>
    mg_smoother;

  mg_smoother.initialize(mg_matrices, typename SMOOTHER::AdditionalData(0.8));

  mg::Matrix<LinearAlgebra::distributed::BlockVector<number>> mg_matrix(
    mg_matrices);
  mg::Matrix<LinearAlgebra::distributed::BlockVector<number>> mg_interface(
    mg_interface_matrices);

  Multigrid<LinearAlgebra::distributed::BlockVector<number>> mg(
    mg_matrix, mg_coarse, mg_transfer, mg_smoother, mg_smoother);
  mg.set_edge_matrices(mg_interface, mg_interface);
  PreconditionMG<dim,
                 LinearAlgebra::distributed::BlockVector<number>,
                 MGTransferBlockMatrixFree<dim, number>>
    preconditioner(dof, mg, mg_transfer);

  {
    // avoid output from inner (coarse-level) solver
    deallog.depth_file(3);

    ReductionControl control(30, 1e-20, 1e-7);
    SolverCG<LinearAlgebra::distributed::BlockVector<number>> solver(control);
    solver.solve(fine_matrix, sol, in, preconditioner);
  }

  if (std::is_same_v<number, float> == true)
    deallog.pop();

  fine_matrix.clear();
  for (unsigned int level = 0;
       level < dof[0]->get_triangulation().n_global_levels();
       ++level)
    mg_matrices[level].clear();
}



template <int dim, int fe_degree_1, int fe_degree_2>
void
test()
{
  parallel::distributed::Triangulation<dim> tria(
    MPI_COMM_WORLD,
    Triangulation<dim>::limit_level_difference_at_vertices,
    parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(6 - dim);
  constexpr int      max_degree = std::max(fe_degree_1, fe_degree_2);
  const unsigned int n_runs     = max_degree == 1 ? 6 - dim : 5 - dim;
  for (unsigned int i = 0; i < n_runs; ++i)
    {
      for (typename Triangulation<dim>::active_cell_iterator cell =
             tria.begin_active();
           cell != tria.end();
           ++cell)
        if (cell->is_locally_owned() &&
            (((cell->center().norm() < 0.5 &&
               (cell->level() < 5 || cell->center().norm() > 0.45)) ||
              (dim == 2 && cell->center().norm() > 1.2))))
          cell->set_refine_flag();
      tria.execute_coarsening_and_refinement();
      FE_Q<dim>       fe_1(fe_degree_1);
      DoFHandler<dim> dof_1(tria);
      dof_1.distribute_dofs(fe_1);
      dof_1.distribute_mg_dofs();

      FE_Q<dim>       fe_2(fe_degree_2);
      DoFHandler<dim> dof_2(tria);
      dof_2.distribute_dofs(fe_2);
      dof_2.distribute_mg_dofs();

      std::vector<const DoFHandler<dim, dim> *> dh_ptrs{&dof_1, &dof_2};

      constexpr int n_q_points_1d = std::max(fe_degree_1, fe_degree_2) + 1;
      do_test<dim, fe_degree_1, fe_degree_2, n_q_points_1d, double>(dh_ptrs);
    }
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);
  mpi_initlog();
  deallog << std::setprecision(4);

  // 2 blocks
  {
    deallog.push("2d");
    test<2, 1, 1>();
    deallog << std::endl;
    test<2, 1, 3>();
    deallog << std::endl;
    test<2, 3, 1>();
    deallog.pop();
    deallog << std::endl;
    deallog.push("3d");
    test<3, 1, 1>();
    deallog << std::endl;
    test<3, 1, 2>();
    deallog << std::endl;
    test<3, 2, 1>();
    deallog.pop();
  }
}
