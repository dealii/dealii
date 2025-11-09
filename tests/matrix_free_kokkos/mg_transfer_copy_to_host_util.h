// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/function_lib.h>

#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/multigrid/mg_transfer_global_coarsening.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"



template <typename MeshType, typename Number, typename MemorySpace>
void
initialize_dof_vector(
  LinearAlgebra::distributed::Vector<Number, MemorySpace> &vec,
  const MeshType                                          &dof_handler,
  const unsigned int                                       level)
{
  IndexSet locally_relevant_dofs;
  if (level == numbers::invalid_unsigned_int)
    locally_relevant_dofs =
      DoFTools::extract_locally_relevant_dofs(dof_handler);
  else
    locally_relevant_dofs =
      DoFTools::extract_locally_relevant_level_dofs(dof_handler, level);


  const parallel::TriangulationBase<MeshType::dimension> *dist_tria =
    dynamic_cast<const parallel::TriangulationBase<MeshType::dimension> *>(
      &(dof_handler.get_triangulation()));

  MPI_Comm comm =
    dist_tria != nullptr ? dist_tria->get_mpi_communicator() : MPI_COMM_SELF;

  vec.reinit(level == numbers::invalid_unsigned_int ?
               dof_handler.locally_owned_dofs() :
               dof_handler.locally_owned_mg_dofs(level),
             locally_relevant_dofs,
             comm);
}

template <typename Number>
void
print(const LinearAlgebra::distributed::Vector<Number, MemorySpace::Host> &vec)
{
  for (const auto &v : vec)
    deallog << v << " ";
  deallog << std::endl;
}

template <typename Number>
void
print_if_non_zero(const LinearAlgebra::distributed::Vector<Number> &vec,
                  const Number                                      tolerance)
{
  for (const auto &v : vec)
    if (std::abs(v) > tolerance)
      deallog << v << " ";
    else
      deallog << 0.0 << " ";
  deallog << std::endl;
}

template <typename Number>
void
copy_to_host(
  LinearAlgebra::distributed::Vector<Number, dealii::MemorySpace::Host> &dst,
  const LinearAlgebra::distributed::Vector<Number, dealii::MemorySpace::Default>
    &src)
{
  LinearAlgebra::ReadWriteVector<Number> rw_vector(
    src.get_partitioner()->locally_owned_range());
  rw_vector.import_elements(src, VectorOperation::insert);

  dst.reinit(src.get_partitioner());
  dst.import_elements(rw_vector, VectorOperation::insert);
}

template <typename Number>
void
copy_from_host(
  LinearAlgebra::distributed::Vector<Number, dealii::MemorySpace::Default> &dst,
  const LinearAlgebra::distributed::Vector<Number, dealii::MemorySpace::Host>
    &src)
{
  LinearAlgebra::ReadWriteVector<Number> rw_vector(
    src.get_partitioner()->locally_owned_range());
  rw_vector.import_elements(src, VectorOperation::insert);


  dst.reinit(src.get_partitioner());
  dst.import_elements(rw_vector, VectorOperation::insert);
}


template <int dim, typename Number, typename MemorySpace, typename MeshType>
void
test_copy_to_host_transfer_operator(
  const MGTwoLevelTransferCopyToHost<
    dim,
    LinearAlgebra::distributed::Vector<Number, MemorySpace>> &transfer,
  const MeshType                                             &dof_handler_fine,
  const MeshType    &dof_handler_coarse,
  const unsigned int mg_level_fine   = numbers::invalid_unsigned_int,
  const unsigned int mg_level_coarse = numbers::invalid_unsigned_int)
{
  AffineConstraints<Number> constraint_fine(
    dof_handler_fine.locally_owned_dofs(),
    DoFTools::extract_locally_relevant_dofs(dof_handler_fine));
  DoFTools::make_hanging_node_constraints(dof_handler_fine, constraint_fine);
  constraint_fine.close();

  // perform prolongation
  LinearAlgebra::distributed::Vector<Number, MemorySpace> src, dst;

  // Most operations are not safe on device vectors, we will have to copy the
  // result back to the host for most processing
  LinearAlgebra::distributed::Vector<Number, dealii::MemorySpace::Host>
    src_host, dst_host;

  initialize_dof_vector(dst, dof_handler_fine, mg_level_fine);
  initialize_dof_vector(src, dof_handler_coarse, mg_level_coarse);

  initialize_dof_vector(dst_host, dof_handler_fine, mg_level_fine);
  initialize_dof_vector(src_host, dof_handler_coarse, mg_level_coarse);

  // test prolongation
  {
    src = 0.0;
    src = 1.0;
    dst = 0.0;
    transfer.prolongate_and_add(dst, src);

    copy_to_host(dst_host, dst);
    copy_to_host(src_host, src);
    // transfer operator sets only non-constrained dofs -> update the rest
    // via constraint matrix
    constraint_fine.distribute(dst_host);

    // print norms
    if (true)
      {
        deallog << dst_host.l2_norm() << std::endl;
      }

    // print vectors
    if (true)
      {
        print(src_host);
        print(dst_host);
      }
  }

  // test restriction
  {
    dst = 1.0;
    src = 0.0;
    transfer.restrict_and_add(src, dst);

    // print norms
    // This operation is safe on the device
    if (true)
      {
        deallog << src.l2_norm() << std::endl;
      }

    copy_to_host(dst_host, dst);
    copy_to_host(src_host, src);
    // print vectors
    if (true)
      {
        print(dst_host);
        print(src_host);
      }
  }
}
