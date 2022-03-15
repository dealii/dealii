// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2020 by the deal.II authors
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


// test matrix-free Laplace and Mass operators with MappingQEulerian
// by comparing to the results obtained from a deformed mesh.
// As a displacement function use:  x exp (-2|x|)


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

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

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
    return p[component] * std::exp(-2. * std::abs(p[component]));
  }
};


template <int dim,
          int fe_degree       = 2,
          int n_q_points      = fe_degree + 1,
          typename NumberType = double>
void
test()
{
  MPI_Comm     mpi_communicator(MPI_COMM_WORLD);
  unsigned int myid    = Utilities::MPI::this_mpi_process(mpi_communicator);
  unsigned int numproc = Utilities::MPI::n_mpi_processes(mpi_communicator);

  parallel::distributed::Triangulation<dim> triangulation(
    mpi_communicator,
    Triangulation<dim>::limit_level_difference_at_vertices,
    parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy);

  GridGenerator::subdivided_hyper_cube(triangulation, 8, -5, 5);

  FE_Q<dim>       fe(fe_degree);
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);
  dof_handler.distribute_mg_dofs();

  const unsigned int euler_fe_degree = 1;
  FESystem<dim>      fe_euler(FE_Q<dim>(euler_fe_degree), dim);
  DoFHandler<dim>    dof_handler_euler(triangulation);
  dof_handler_euler.distribute_dofs(fe_euler);
  dof_handler_euler.distribute_mg_dofs();

  // IndexSets and constraints
  const IndexSet &locally_owned_dofs_euler =
    dof_handler_euler.locally_owned_dofs();
  IndexSet locally_relevant_dofs_euler;
  DoFTools::extract_locally_relevant_dofs(dof_handler_euler,
                                          locally_relevant_dofs_euler);

  const IndexSet &locally_owned_dofs = dof_handler.locally_owned_dofs();
  IndexSet        locally_relevant_dofs;
  DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);

  // constraints:
  AffineConstraints<double> constraints_euler;
  constraints_euler.reinit(locally_relevant_dofs_euler);
  DoFTools::make_hanging_node_constraints(dof_handler_euler, constraints_euler);
  constraints_euler.close();

  AffineConstraints<double> constraints;
  constraints.reinit(locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  constraints_euler.close();

  // Displacement vector
  LinearAlgebra::distributed::Vector<NumberType> displacement;
  displacement.reinit(locally_owned_dofs_euler,
                      locally_relevant_dofs_euler,
                      mpi_communicator);
  displacement = 0.;

  Displacement<dim> displacement_function;

  // first, move via mapping:
  VectorTools::interpolate(dof_handler_euler,
                           displacement_function,
                           displacement);
  displacement.compress(VectorOperation::insert);
  displacement.update_ghost_values();

  // The core : compute matrix-vector products with mass and laplace operators
  QGauss<1> quadrature_formula(n_q_points);

  MappingQEulerian<dim, LinearAlgebra::distributed::Vector<NumberType>>
    euler_mapping(euler_fe_degree, dof_handler_euler, displacement);

  std::shared_ptr<MatrixFree<dim, NumberType>> matrix_free(
    new MatrixFree<dim, NumberType>());
  typename MatrixFree<dim, NumberType>::AdditionalData data;
  data.tasks_parallel_scheme =
    MatrixFree<dim, NumberType>::AdditionalData::partition_color;
  data.mapping_update_flags =
    update_values | update_gradients | update_JxW_values;
  matrix_free->reinit(
    euler_mapping, dof_handler, constraints, quadrature_formula, data);

  LinearAlgebra::distributed::Vector<NumberType> src, dst1, dst2, dst3, dst4;
  matrix_free->initialize_dof_vector(src);
  matrix_free->initialize_dof_vector(dst1);
  matrix_free->initialize_dof_vector(dst2);
  matrix_free->initialize_dof_vector(dst3);
  matrix_free->initialize_dof_vector(dst4);

  for (unsigned int i = 0; i < src.local_size(); ++i)
    src.local_element(i) = random_value<NumberType>();

  MatrixFreeOperators::MassOperator<
    dim,
    fe_degree,
    n_q_points,
    1,
    LinearAlgebra::distributed::Vector<NumberType>>
    mf_mass;
  MatrixFreeOperators::LaplaceOperator<
    dim,
    fe_degree,
    n_q_points,
    1,
    LinearAlgebra::distributed::Vector<NumberType>>
    mf_laplace;
  mf_mass.initialize(matrix_free);
  mf_laplace.initialize(matrix_free);

  mf_mass.vmult(dst1, src);
  mf_laplace.vmult(dst3, src);

  // now move manually
  displacement = 0.;
  displacement.update_ghost_values();
  {
    typename DoFHandler<dim>::cell_iterator cell = dof_handler.begin_active(),
                                            endc = dof_handler.end();
    std::vector<bool> vertex_touched(triangulation.n_vertices(), false);
    for (cell = dof_handler.begin_active(); cell != endc; ++cell)
      for (const unsigned int vertex_no : GeometryInfo<dim>::vertex_indices())
        if (vertex_touched[cell->vertex_index(vertex_no)] == false)
          {
            Point<dim> &   v = cell->vertex(vertex_no);
            Tensor<1, dim> d;
            for (unsigned int i = 0; i < dim; ++i)
              d[i] = displacement_function.value(v, i);

            v += d;
            vertex_touched[cell->vertex_index(vertex_no)] = true;
          }
  }
  // minimize the data that is re-computed
  data.initialize_indices = false;
  matrix_free->reinit(
    euler_mapping, dof_handler, constraints, quadrature_formula, data);

  mf_mass.vmult(dst2, src);
  mf_laplace.vmult(dst4, src);

  deallog << "Mass operator: " << std::endl
          << "l2_norm:       " << std::abs(dst1.l2_norm() - dst2.l2_norm())
          << std::endl
          << "l1_norm:       " << std::abs(dst1.l1_norm() - dst2.l1_norm())
          << std::endl
          << "linfty_norm:   "
          << std::abs(dst1.linfty_norm() - dst2.linfty_norm()) << std::endl
          << "Laplace operator: " << std::endl
          << "l2_norm:       " << std::abs(dst3.l2_norm() - dst4.l2_norm())
          << std::endl
          << "l1_norm:       " << std::abs(dst3.l1_norm() - dst4.l1_norm())
          << std::endl
          << "linfty_norm:   "
          << std::abs(dst3.linfty_norm() - dst4.linfty_norm()) << std::endl;
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  mpi_initlog();

  deallog.push("2d");
  test<2, 2>();
  deallog.pop();

  deallog.push("3d");
  test<3, 2>();
  deallog.pop();
}
