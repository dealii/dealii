// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test matrix-free Laplace and Mass operators with MappingQEulerian
// on a refined mesh that triggered the bug below with 4 MPI
// ranks. This is due to a smaller set of ghost indices used by
// interpolate_to_mg.

// clang-format off
/*
An error occurred in line <1638> of file <../include/deal.II/lac/la_parallel_vector.h> in function
    Number dealii::LinearAlgebra::distributed::Vector<Number, MemorySpace>::operator()(dealii::LinearAlgebra::distributed::Vector<Number, MemorySpace>::size_type) const [with Number = double; MemorySpace = dealii::MemorySpace::Host; dealii::LinearAlgebra::distributed::Vector<Number, MemorySpace>::size_type = unsigned int]
The violated condition was:
    partitioner->in_local_range(global_index) || partitioner->ghost_indices().is_element(global_index)
*/
// clang-format on

#include <deal.II/base/function.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/thread_management.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/mapping_q_eulerian.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/identity_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_memory.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>

#include <deal.II/multigrid/mg_transfer_matrix_free.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <iostream>
#include <vector>

#include "../tests.h"


const std::vector<std::string> refined = {
  "0_0:",   "10_0:",  "11_0:",  "12_0:",  "13_0:",  "14_0:",  "15_0:",
  "1_0:",   "21_0:",  "23_0:",  "26_0:",  "27_0:",  "28_0:",  "29_0:",
  "2_0:",   "30_0:",  "31_0:",  "32_0:",  "33_0:",  "34_0:",  "34_1:0",
  "34_1:1", "34_1:2", "34_1:3", "35_0:",  "35_1:0", "35_1:2", "35_1:3",
  "36_0:",  "37_0:",  "38_0:",  "38_1:2", "39_0:",  "3_0:",   "40_0:",
  "40_1:0", "40_1:1", "40_1:2", "40_1:3", "41_0:",  "41_1:0", "41_1:1",
  "41_1:2", "41_1:3", "42_0:",  "42_1:0", "42_1:1", "42_1:2", "42_1:3",
  "43_0:",  "43_1:0", "43_1:1", "43_1:2", "43_1:3", "44_0:",  "44_1:0",
  "44_1:1", "44_1:2", "44_1:3", "45_0:",  "46_0:",  "46_1:0", "46_1:1",
  "46_1:2", "46_1:3", "47_0:",  "47_1:0", "47_1:2", "48_0:",  "49_0:",
  "4_0:",   "50_0:",  "51_0:",  "52_0:",  "53_0:",  "54_0:",  "55_0:",
  "56_0:",  "57_0:",  "58_0:",  "59_0:",  "60_0:",  "61_0:",  "62_0:",
  "63_0:",  "64_0:",  "65_0:",  "66_0:",  "67_0:",  "68_0:",  "69_0:",
  "6_0:",   "70_0:",  "71_0:",  "72_0:",  "73_0:",  "74_0:",  "75_0:",
  "76_0:",  "77_0:",  "78_0:",  "79_0:",  "80_0:",  "80_1:0", "80_1:1",
  "80_1:2", "80_1:3", "81_0:",  "81_1:0", "81_1:1", "81_1:2", "81_1:3",
  "82_0:",  "82_1:0", "82_1:1", "82_1:2", "82_1:3", "83_0:",  "83_1:0",
  "83_1:1", "83_1:2", "83_1:3", "84_0:",  "84_1:0", "84_1:1", "84_1:2",
  "84_1:3", "85_0:",  "85_1:0", "86_0:",  "86_1:0", "86_1:1", "86_1:2",
  "86_1:3", "87_0:",  "88_0:",  "89_0:",  "8_0:",   "92_0:",  "9_0:"};


// Displacement to be applied
template <int dim>
class Displacement : public Function<dim>
{
public:
  Displacement()
    : Function<dim>(dim)
  {}

  double
  value(const Point<dim> &p, const unsigned int component) const
  {
    return (0.1 + 2. * component) + (0.2 + 3. * component) * p[component];
  }

  template <typename NumberType>
  Tensor<1, dim, VectorizedArray<NumberType>>
  shift_value(const Point<dim, VectorizedArray<NumberType>> &p_vec) const
  {
    Tensor<1, dim, VectorizedArray<NumberType>> shift_vec;
    Point<dim>                                  p;
    for (unsigned int v = 0; v < VectorizedArray<NumberType>::size(); ++v)
      {
        for (unsigned int d = 0; d < dim; ++d)
          p[d] = p_vec[d][v];

        for (unsigned int d = 0; d < dim; ++d)
          shift_vec[d][v] = this->value(p, d);
      }

    return shift_vec;
  }
};


template <int dim,
          int fe_degree            = 2,
          int n_q_points           = fe_degree + 1,
          typename NumberType      = double,
          typename LevelNumberType = NumberType>
void
test(const unsigned int n_ref = 0)
{
  Displacement<dim>  displacement_function;
  const unsigned int euler_fe_degree = 2;

  MPI_Comm     mpi_communicator(MPI_COMM_WORLD);
  unsigned int myid    = Utilities::MPI::this_mpi_process(mpi_communicator);
  unsigned int numproc = Utilities::MPI::n_mpi_processes(mpi_communicator);

  parallel::distributed::Triangulation<dim> triangulation(
    mpi_communicator,
    Triangulation<dim>::limit_level_difference_at_vertices,
    parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy);

  GridGenerator::subdivided_hyper_cube(triangulation, 10, -5, 5);
  while (true)
    {
      deallog << "refine..." << std::endl;
      bool changes = false;

      for (const auto &cell : triangulation.active_cell_iterators())

        if (!cell->is_artificial() &&
            std::find(std::begin(refined),
                      std::end(refined),
                      cell->id().to_string()) != std::end(refined))
          {
            cell->set_refine_flag();
            changes = true;
          }

      if (Utilities::MPI::max(changes ? 1 : 0, MPI_COMM_WORLD) == 0)
        break;

      triangulation.execute_coarsening_and_refinement();
    }
  deallog << "Checksum: " << triangulation.get_checksum() << '\n'
          << "n_active_cells: " << triangulation.n_active_cells() << '\n'
          << "n_levels: " << triangulation.n_levels() << std::endl;

  FE_Q<dim>       fe(fe_degree);
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);
  dof_handler.distribute_mg_dofs();

  // quadrature for MatrixFree, not related to the degree in the
  // FE displacement field.
  QGauss<1> quadrature_formula(n_q_points);


  FESystem<dim>   fe_euler(FE_Q<dim>(euler_fe_degree), dim);
  DoFHandler<dim> dof_handler_euler(triangulation);
  dof_handler_euler.distribute_dofs(fe_euler);
  dof_handler_euler.distribute_mg_dofs();

  // IndexSets and constraints
  const IndexSet &locally_owned_dofs_euler =
    dof_handler_euler.locally_owned_dofs();
  const IndexSet locally_relevant_dofs_euler =
    DoFTools::extract_locally_relevant_dofs(dof_handler_euler);

  const IndexSet &locally_owned_dofs = dof_handler.locally_owned_dofs();
  const IndexSet  locally_relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(dof_handler);

  // constraints:
  AffineConstraints<double> constraints_euler;
  constraints_euler.reinit(locally_owned_dofs_euler,
                           locally_relevant_dofs_euler);
  DoFTools::make_hanging_node_constraints(dof_handler_euler, constraints_euler);
  constraints_euler.close();

  AffineConstraints<double> constraints;
  constraints.reinit(locally_owned_dofs, locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  constraints.close();

  // MG constraints:
  MGConstrainedDoFs mg_constrained_dofs_euler;
  mg_constrained_dofs_euler.initialize(dof_handler_euler);

  // Displacement vector
  LinearAlgebra::distributed::Vector<NumberType> displacement;
  displacement.reinit(locally_owned_dofs_euler,
                      locally_relevant_dofs_euler,
                      mpi_communicator);

  VectorTools::project<dim,
                       LinearAlgebra::distributed::Vector<NumberType>,
                       dim>(dof_handler_euler,
                            constraints_euler,
                            QGauss<dim>(n_q_points),
                            displacement_function,
                            displacement);
  displacement.update_ghost_values();

  MGTransferMatrixFree<dim, LevelNumberType> mg_transfer_euler(
    mg_constrained_dofs_euler);
  mg_transfer_euler.build(dof_handler_euler);

  // now the core of the test:
  const unsigned int max_level =
    dof_handler.get_triangulation().n_global_levels() - 1;
  const unsigned int min_level = 0;
  MGLevelObject<LinearAlgebra::distributed::Vector<LevelNumberType>>
    displacement_level(min_level, max_level);

  // Important! This preallocation of the displacement vectors fixes
  // the bug observed.
  for (unsigned int level = min_level; level <= max_level; ++level)
    {
      const IndexSet relevant_mg_dofs =
        DoFTools::extract_locally_relevant_level_dofs(dof_handler_euler, level);
      displacement_level[level].reinit(dof_handler_euler.locally_owned_mg_dofs(
                                         level),
                                       relevant_mg_dofs,
                                       mpi_communicator);
    }
  mg_transfer_euler.interpolate_to_mg(dof_handler_euler,
                                      displacement_level,
                                      displacement);

  // First, check fine-level only:
  {
    MappingQEulerian<dim, LinearAlgebra::distributed::Vector<NumberType>>
      euler_fine(euler_fe_degree, dof_handler_euler, displacement);



    MatrixFree<dim, NumberType>                          matrix_free_euler;
    typename MatrixFree<dim, NumberType>::AdditionalData data;
    data.tasks_parallel_scheme =
      MatrixFree<dim, NumberType>::AdditionalData::partition_color;
    data.mapping_update_flags = update_values | update_gradients |
                                update_JxW_values | update_quadrature_points;
    matrix_free_euler.reinit(
      euler_fine, dof_handler, constraints, quadrature_formula, data);

    MatrixFree<dim, NumberType> matrix_free;
    matrix_free.reinit(
      MappingQ1<dim>{}, dof_handler, constraints, quadrature_formula, data);


    // test fine-level mapping:
    {
      FEEvaluation<dim, fe_degree, n_q_points, 1, NumberType> fe_eval_euler(
        matrix_free_euler);
      FEEvaluation<dim, fe_degree, n_q_points, 1, NumberType> fe_eval(
        matrix_free);
      const unsigned int n_cells = matrix_free_euler.n_cell_batches();
      Assert(matrix_free_euler.n_cell_batches() == matrix_free.n_cell_batches(),
             ExcInternalError());
      const unsigned int nqp = fe_eval.n_q_points;
      for (unsigned int cell = 0; cell < n_cells; ++cell)
        {
          fe_eval_euler.reinit(cell);
          fe_eval.reinit(cell);
          for (unsigned int q = 0; q < nqp; ++q)
            {
              const auto &v1 = fe_eval_euler.quadrature_point(q);
              const auto &qp = fe_eval.quadrature_point(q);
              const auto  v2 = qp + displacement_function.shift_value(qp);
              VectorizedArray<NumberType> dist = v1.distance(v2);
              for (unsigned int v = 0; v < VectorizedArray<NumberType>::size();
                   ++v)
                AssertThrow(dist[v] < 1e-8,
                            ExcMessage("distance: " + std::to_string(dist[v])));
            }
        }
    }
  }

  // now go through all GMG levels:
  const std::set<types::boundary_id> dirichlet_boundary_ids = {0};
  MGConstrainedDoFs                  mg_constrained_dofs;
  mg_constrained_dofs.initialize(dof_handler);
  mg_constrained_dofs.make_zero_boundary_constraints(dof_handler,
                                                     dirichlet_boundary_ids);

  for (unsigned int level = min_level; level <= max_level; ++level)
    {
      typename MatrixFree<dim, LevelNumberType>::AdditionalData
        mg_additional_data;
      mg_additional_data.tasks_parallel_scheme =
        MatrixFree<dim, LevelNumberType>::AdditionalData::partition_color;
      mg_additional_data.mg_level = level;
      mg_additional_data.mapping_update_flags =
        update_values | update_gradients | update_JxW_values |
        update_quadrature_points;

      AffineConstraints<double> level_constraints;
      const IndexSet            relevant_dofs =
        DoFTools::extract_locally_relevant_level_dofs(dof_handler, level);
      level_constraints.reinit(dof_handler.locally_owned_mg_dofs(level),
                               relevant_dofs);
      level_constraints.add_lines(
        mg_constrained_dofs.get_boundary_indices(level));
      level_constraints.close();

      MatrixFree<dim, LevelNumberType> mg_level_euler;

      displacement_level[level].update_ghost_values();

      MappingQEulerian<dim, LinearAlgebra::distributed::Vector<LevelNumberType>>
        euler_level(euler_fe_degree,
                    dof_handler_euler,
                    displacement_level[level],
                    level);

      mg_level_euler.reinit(euler_level,
                            dof_handler,
                            level_constraints,
                            quadrature_formula,
                            mg_additional_data);

      MatrixFree<dim, LevelNumberType> mg_level;
      mg_level.reinit(MappingQ1<dim>{},
                      dof_handler,
                      level_constraints,
                      quadrature_formula,
                      mg_additional_data);

      // go through all cells and quadrature points:
      {
        FEEvaluation<dim, fe_degree, n_q_points, 1, NumberType> fe_eval_euler(
          mg_level_euler);
        FEEvaluation<dim, fe_degree, n_q_points, 1, NumberType> fe_eval(
          mg_level);
        const unsigned int n_cells = mg_level_euler.n_cell_batches();
        Assert(mg_level_euler.n_cell_batches() == mg_level.n_cell_batches(),
               ExcInternalError());
        const unsigned int nqp = fe_eval.n_q_points;
        for (unsigned int cell = 0; cell < n_cells; ++cell)
          {
            fe_eval_euler.reinit(cell);
            fe_eval.reinit(cell);
            for (unsigned int q = 0; q < nqp; ++q)
              {
                const auto &v1 = fe_eval_euler.quadrature_point(q);
                const auto &qp = fe_eval.quadrature_point(q);
                const auto  v2 = qp + displacement_function.shift_value(qp);
                VectorizedArray<NumberType> dist = v1.distance(v2);
                for (unsigned int v = 0;
                     v < VectorizedArray<NumberType>::size();
                     ++v)
                  AssertThrow(dist[v] < 1e-8,
                              ExcMessage(
                                "Level " + std::to_string(level) +
                                " distance: " + std::to_string(dist[v])));
              }
          }
      }
    }

  deallog << "Ok" << std::endl;
}



int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  test<2>();
}
