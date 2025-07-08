// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test Portable::MatrixFree::initialize_dof_vector.

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/affine_constraints.h>

#include <deal.II/matrix_free/portable_matrix_free.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


template <typename Number, typename MemorySpaceType>
void
check(const LinearAlgebra::distributed::Vector<Number, MemorySpaceType> &vector,
      const Utilities::MPI::Partitioner &reference_partitioner)
{
  Assert(vector.get_partitioner()->locally_owned_range() ==
           reference_partitioner.locally_owned_range(),
         ExcInternalError());
  Assert(vector.get_partitioner()->ghost_indices() ==
           reference_partitioner.ghost_indices(),
         ExcInternalError());
}

template <int dim, int fe_degree, typename VectorType>
void
test()
{
  using Number = double;

  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(2);

  FE_Q<dim>       fe(fe_degree);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);

  const IndexSet &owned_set    = dof.locally_owned_dofs();
  const IndexSet  relevant_set = DoFTools::extract_locally_relevant_dofs(dof);

  deallog << "locally owned dofs :" << std::endl;
  owned_set.print(deallog.get_file_stream());

  deallog << "locally relevant dofs :" << std::endl;
  relevant_set.print(deallog.get_file_stream());

  AffineConstraints<double> constraints(owned_set, relevant_set);
  constraints.close();

  MappingQ<dim>                     mapping(fe_degree);
  Portable::MatrixFree<dim, Number> mf_data;
  const QGauss<1>                   quad(fe_degree + 1);
  typename Portable::MatrixFree<dim, Number>::AdditionalData additional_data;
  mf_data.reinit(mapping, dof, constraints, quad, additional_data);

  VectorType vector;
  mf_data.initialize_dof_vector(vector);

  Utilities::MPI::Partitioner reference_partitioner(owned_set,
                                                    relevant_set,
                                                    MPI_COMM_WORLD);
  check(vector, reference_partitioner);

  deallog << "OK" << std::endl;
}


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  deallog.push(Utilities::int_to_string(myid));

  MPILogInitAll mpi_inilog;

  test<2, 1, LinearAlgebra::distributed::Vector<double, MemorySpace::Host>>();

  test<2,
       1,
       LinearAlgebra::distributed::Vector<double, MemorySpace::Default>>();

  return 0;
}
