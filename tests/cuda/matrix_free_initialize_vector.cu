// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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



// Test CUDAWrappers::MatrixFree::initialize_dof_vector.

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/cuda_vector.h>

#include <deal.II/matrix_free/cuda_matrix_free.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


template <typename Number>
void
check(
  const LinearAlgebra::distributed::Vector<Number, MemorySpace::CUDA> &vector,
  const Utilities::MPI::Partitioner &reference_partitioner)
{
  Assert(vector.get_partitioner()->locally_owned_range() ==
           reference_partitioner.locally_owned_range(),
         ExcInternalError());
  Assert(vector.get_partitioner()->ghost_indices() ==
           reference_partitioner.ghost_indices(),
         ExcInternalError());
}

template <typename Number>
void
check(const LinearAlgebra::CUDAWrappers::Vector<Number> &vector,
      const Utilities::MPI::Partitioner &                reference_partitioner)
{
  AssertDimension(vector.size(), reference_partitioner.size());
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

  IndexSet owned_set = dof.locally_owned_dofs();
  IndexSet relevant_set;
  DoFTools::extract_locally_relevant_dofs(dof, relevant_set);

  deallog << "locally owned dofs :" << std::endl;
  owned_set.print(deallog.get_file_stream());

  deallog << "locally relevant dofs :" << std::endl;
  relevant_set.print(deallog.get_file_stream());

  AffineConstraints<double> constraints(relevant_set);
  constraints.close();

  MappingQ<dim>                         mapping(fe_degree);
  CUDAWrappers::MatrixFree<dim, Number> mf_data;
  const QGauss<1>                       quad(fe_degree + 1);
  typename CUDAWrappers::MatrixFree<dim, Number>::AdditionalData
    additional_data;
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

  init_cuda(true);
  MPILogInitAll mpi_inilog;

  test<2, 1, LinearAlgebra::distributed::Vector<double, MemorySpace::CUDA>>();
  if (Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) == 1)
    test<2, 1, LinearAlgebra::CUDAWrappers::Vector<double>>();
}
