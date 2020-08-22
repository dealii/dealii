// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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


// Test PolynomialsTet on quadrature points returned by QGaussTet.

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q_generic.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/la_sm_vector.h>

#include <deal.II/matrix_free/matrix_free.h>

#include "../tests.h"

using namespace dealii;

template <typename Number, int dim>
void
test(const int n_refinements, const int degree, const int group_size)
{
  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::subdivided_hyper_cube(tria, 2);

  DoFHandler<dim> dof_handler(tria);
  FE_Q<dim>       fe(degree);
  dof_handler.distribute_dofs(fe);

  AffineConstraints<Number> constraint;
  constraint.close();

  QGauss<1> quad(degree + 1);

  MappingQGeneric<dim> mapping(1);

  typename dealii::MatrixFree<dim, Number>::AdditionalData additional_data;
  additional_data.mapping_update_flags =
    update_gradients | update_JxW_values | update_quadrature_points;

  dealii::MatrixFree<dim, Number> matrix_free;
  matrix_free.reinit(mapping, dof_handler, constraint, quad, additional_data);

  const MPI_Comm comm = MPI_COMM_WORLD;

  MPI_Comm comm_sm = MPI_COMM_SELF;

  using VectorType = LinearAlgebra::SharedMPI::Vector<Number>;

  LinearAlgebra::SharedMPI::Vector<Number> dst;
  matrix_free.initialize_dof_vector(dst, comm_sm);

  const unsigned int local_size   = dst.local_size();
  unsigned int       local_offset = 0;
  unsigned int       n_ghost_elements =
    matrix_free.get_vector_partitioner(comm_sm)->n_ghost_indices();

  MPI_Exscan(
    &local_size, &local_offset, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);

  for (unsigned int i = 0; i < local_size; ++i)
    dst.local_element(i) = i + local_offset;

  dst.update_ghost_values();

  for (unsigned int i = 0; i < local_size + n_ghost_elements; ++i)
    deallog << static_cast<int>(dst.local_element(i)) << " ";
  deallog << std::endl;

  dst.set_ghost_state(false);
  dst.compress(VectorOperation::add);

  for (unsigned int i = 0; i < local_size; ++i)
    deallog << static_cast<int>(dst.local_element(i)) << " ";
  deallog << std::endl;
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  test<double, 2>(1, 1, 1);
}