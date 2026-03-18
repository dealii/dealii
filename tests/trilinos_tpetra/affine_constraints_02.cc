// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2013 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

// check instantiations for some functions inside AffineConstraints
// like affine_constraints_01, but for block matrices and vectors

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/trilinos_tpetra_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_tpetra_block_vector.h>
#include <deal.II/lac/trilinos_tpetra_precondition.h>
#include <deal.II/lac/trilinos_tpetra_solver_direct.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"



template <int dim>
class RightHandSide : public Function<dim>
{
public:
  RightHandSide()
    : Function<dim>()
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const
  {
    double return_value = 0;
    for (unsigned int i = 0; i < dim; ++i)
      return_value += 2 * std::pow(p[i], 2);

    return return_value;
  }
};



template <int dim>
class BoundaryValues : public Function<dim>
{
public:
  BoundaryValues()
    : Function<dim>()
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const
  {
    return p.square();
  }
};



template <int dim, class MatrixType, class VectorType>
void
test()
{
  parallel::distributed::Triangulation<dim> triangulation(
    MPI_COMM_WORLD,
    typename Triangulation<dim>::MeshSmoothing(
      Triangulation<dim>::smoothing_on_refinement |
      Triangulation<dim>::smoothing_on_coarsening));

  FESystem<dim>   fe(FE_Q<dim>(1));
  DoFHandler<dim> dof_handler(triangulation);

  AffineConstraints<double> constraints;
  SparsityPattern           sparsity_pattern;

  GridGenerator::hyper_cube(triangulation, -1, 1);
  triangulation.refine_global(3);


  dof_handler.distribute_dofs(fe);

  constraints.clear();
  VectorTools::interpolate_boundary_values(dof_handler,
                                           0,
                                           BoundaryValues<dim>(),
                                           constraints);
  constraints.close();

  const IndexSet &locally_owned_dofs = dof_handler.locally_owned_dofs();
  const IndexSet  locally_relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(dof_handler);

  // compute the various partitionings between processors and blocks
  // of vectors and matrices
  const std::vector<types::global_dof_index> dofs_per_block = DoFTools::count_dofs_per_fe_block (dof_handler,
                                                                            {{0}});

  const std::vector<IndexSet> system_partitioning = locally_owned_dofs.split_by_block(dofs_per_block);
  const std::vector<IndexSet> system_relevant_partitioning = locally_relevant_dofs.split_by_block(dofs_per_block);

  BlockDynamicSparsityPattern dsp(system_relevant_partitioning);
  DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);
  
  SparsityTools::distribute_sparsity_pattern(dsp,
                                             locally_owned_dofs,
                                             MPI_COMM_WORLD,
                                             locally_relevant_dofs);

  MatrixType system_matrix;
  system_matrix.reinit(dsp);

  VectorType solution;
  solution.reinit(system_partitioning, system_relevant_partitioning, MPI_COMM_WORLD);

  VectorType system_rhs;
  system_rhs.reinit(system_partitioning, system_relevant_partitioning, MPI_COMM_WORLD);

  QGauss<dim> quadrature_formula(fe.degree + 1);

  const RightHandSide<dim> right_hand_side;

  FEValues<dim> fe_values(fe,
                          quadrature_formula,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const unsigned int n_q_points    = quadrature_formula.size();

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          cell_matrix = 0;
          cell_rhs    = 0;

          for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              {
                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                  cell_matrix(i, j) +=
                    (fe_values.shape_grad(i, q_point) *
                     fe_values.shape_grad(j, q_point) * fe_values.JxW(q_point));

                cell_rhs(i) +=
                  (fe_values.shape_value(i, q_point) *
                   right_hand_side.value(fe_values.quadrature_point(q_point)) *
                   fe_values.JxW(q_point));
              }

          cell->get_dof_indices(local_dof_indices);
          constraints.distribute_local_to_global(cell_matrix,
                                                 local_dof_indices,
                                                 system_matrix);

          constraints.distribute_local_to_global(cell_rhs,
                                                 local_dof_indices,
                                                 system_rhs);
        }
    }
  system_matrix.compress(VectorOperation::add);
  system_rhs.compress(VectorOperation::add);

  deallog << "Ok " << std::endl;
}


int
main(int argc, char **argv)
{
  initlog();

  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  test<
    2,
    LinearAlgebra::TpetraWrappers::BlockSparseMatrix<double, MemorySpace::Host>,
    LinearAlgebra::TpetraWrappers::BlockVector<double, MemorySpace::Host>>();
  test<
    3,
    LinearAlgebra::TpetraWrappers::BlockSparseMatrix<double, MemorySpace::Host>,
    LinearAlgebra::TpetraWrappers::BlockVector<double, MemorySpace::Host>>();

    // test<2,LinearAlgebra::TpetraWrappers::BlockSparseMatrix<double,
    // MemorySpace::Default>, LinearAlgebra::TpetraWrappers::BlockVector<double,
    // MemorySpace::Default>>();
    // test<3,LinearAlgebra::TpetraWrappers::BlockSparseMatrix<double,
    // MemorySpace::Default>, LinearAlgebra::TpetraWrappers::BlockVector<double,
    // MemorySpace::Default>>();
}
