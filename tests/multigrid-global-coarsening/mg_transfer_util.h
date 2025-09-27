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

#ifndef dealii_multigrid_transfer_tests_h
#define dealii_multigrid_transfer_tests_h

#include <deal.II/base/function_lib.h>

#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/multigrid/mg_transfer_global_coarsening.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"



template <typename MeshType, typename Number>
void
initialize_dof_vector(LinearAlgebra::distributed::Vector<Number> &vec,
                      const MeshType                             &dof_handler,
                      const unsigned int                          level)
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
print(const LinearAlgebra::distributed::Vector<Number> &vec)
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


template <int dim, typename Number, typename MeshType>
void
test_transfer_operator(
  const MGTwoLevelTransfer<dim, LinearAlgebra::distributed::Vector<Number>>
                    &transfer,
  const MeshType    &dof_handler_fine,
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
  LinearAlgebra::distributed::Vector<Number> src, dst;

  initialize_dof_vector(dst, dof_handler_fine, mg_level_fine);
  initialize_dof_vector(src, dof_handler_coarse, mg_level_coarse);

  // test prolongation
  {
    src = 0.0;
    src = 1.0;
    dst = 0.0;
    transfer.prolongate_and_add(dst, src);

    // transfer operator sets only non-constrained dofs -> update the rest
    // via constraint matrix
    constraint_fine.distribute(dst);

    // print norms
    if (true)
      {
        deallog << dst.l2_norm() << std::endl;
      }

    // print vectors
    if (true)
      {
        print(src);
        print(dst);
      }

    // print full prolongation matrix
    if (Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) == 1 && false)
      {
        FullMatrix<Number> prolongation_matrix(dst.size(), src.size());
        for (unsigned int i = 0; i < src.size(); ++i)
          {
            src    = 0.0;
            src[i] = 1.0;
            dst    = 0.0;

            transfer.prolongate_and_add(dst, src);

            for (unsigned int j = 0; j < dst.size(); ++j)
              prolongation_matrix[j][i] = dst[j];
          }

        prolongation_matrix.print_formatted(
          deallog.get_file_stream(), 2, false, 5, "", 1, 1e-5);
      }
  }

  // test restriction
  {
    dst = 1.0;
    src = 0.0;
    transfer.restrict_and_add(src, dst);

    // print norms
    if (true)
      {
        deallog << src.l2_norm() << std::endl;
      }

    // print vectors
    if (true)
      {
        print(dst);
        print(src);
      }

    // print full restriction matrix
    if (Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) == 1 && false)
      {
        FullMatrix<Number> restriction_matrix(src.size(), dst.size());
        for (unsigned int i = 0; i < dst.size(); ++i)
          {
            dst    = 0.0;
            dst[i] = 1.0;
            src    = 0.0;

            transfer.restrict_and_add(src, dst);

            for (unsigned int j = 0; j < src.size(); ++j)
              restriction_matrix[j][i] = src[j];
          }

        restriction_matrix.print_formatted(
          deallog.get_file_stream(), 2, false, 5, "", 1, 1e-5);
      }
  }
}



template <int dim, typename Number, typename MeshType>
void
test_non_nested_transfer(
  const MGTwoLevelTransferBase<LinearAlgebra::distributed::Vector<Number>>
                              &transfer,
  const MeshType              &dof_handler_fine,
  const MeshType              &dof_handler_coarse,
  const Function<dim, Number> &function =
    Functions::ZeroFunction<dim, Number>(),
  const unsigned int mg_level_fine   = numbers::invalid_unsigned_int,
  const unsigned int mg_level_coarse = numbers::invalid_unsigned_int)
{
  AffineConstraints<Number> constraint_fine(
    dof_handler_fine.locally_owned_dofs(),
    DoFTools::extract_locally_relevant_dofs(dof_handler_fine));
  DoFTools::make_hanging_node_constraints(dof_handler_fine, constraint_fine);
  constraint_fine.close();

  // perform prolongation
  LinearAlgebra::distributed::Vector<Number> src, dst;

  initialize_dof_vector(dst, dof_handler_fine, mg_level_fine);
  initialize_dof_vector(src, dof_handler_coarse, mg_level_coarse);

  // test prolongation
  {
    src = 0.0;
    VectorTools::interpolate(dof_handler_coarse, function, src); // src = 1.0
    dst = 0.0;
    transfer.prolongate_and_add(dst, src);

    // transfer operator sets only non-constrained dofs -> update the rest
    // via constraint matrix
    constraint_fine.distribute(dst);

    // print norms
    if (true)
      {
        deallog << dst.l2_norm() << std::endl;
      }

    // print vectors
    if (true)
      {
        print(src);
        print(dst);
      }

    // print full prolongation matrix
    if (Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) == 1 && false)
      {
        FullMatrix<Number> prolongation_matrix(dst.size(), src.size());
        for (unsigned int i = 0; i < src.size(); ++i)
          {
            src    = 0.0;
            src[i] = 1.0;
            dst    = 0.0;

            transfer.prolongate_and_add(dst, src);

            for (unsigned int j = 0; j < dst.size(); ++j)
              prolongation_matrix[j][i] = dst[j];
          }

        prolongation_matrix.print_formatted(
          deallog.get_file_stream(), 2, false, 5, "", 1, 1e-5);
      }
  }

  // test restriction
  {
    VectorTools::interpolate(dof_handler_fine, function, dst); // dst = 1.0
    src = 0.0;
    transfer.restrict_and_add(src, dst);

    // print norms
    if (true)
      {
        deallog << src.l2_norm() << std::endl;
      }

    // print vectors
    if (true)
      {
        print(dst);
        print(src);
      }

    // print full restriction matrix
    if (Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) == 1 && false)
      {
        FullMatrix<Number> restriction_matrix(src.size(), dst.size());
        for (unsigned int i = 0; i < dst.size(); ++i)
          {
            dst    = 0.0;
            dst[i] = 1.0;
            src    = 0.0;

            transfer.restrict_and_add(src, dst);

            for (unsigned int j = 0; j < src.size(); ++j)
              restriction_matrix[j][i] = src[j];
          }

        restriction_matrix.print_formatted(
          deallog.get_file_stream(), 2, false, 5, "", 1, 1e-5);
      }
  }
}

#endif
