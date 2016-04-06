// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2015 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



// document problem in PETSc block system ("inserting nonzero")

#define USE_PETSC_LA

#include "../tests.h"
#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/base/index_set.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/base/function.h>

#include "gla.h"

template <class LA, int dim>
void test ()
{
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);

  if (myid==0)
    deallog << "numproc=" << numproc << std::endl;

  parallel::distributed::Triangulation<3> triangulation(MPI_COMM_WORLD);
  GridGenerator::hyper_cube (triangulation, -1.0, 1.0);
  triangulation.refine_global(2);

  FESystem<3> fe(FE_Q<3> (1), 2);
  DoFHandler<3> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);
  DoFRenumbering::block_wise(dof_handler);

  std::vector<types::global_dof_index> dofs_per_block (fe.n_blocks());
  std::vector<unsigned int> sub_blocks (fe.n_blocks(), 0);
  sub_blocks[1] = 1;
  DoFTools::count_dofs_per_block (dof_handler, dofs_per_block, sub_blocks);

  deallog << "size: " << dofs_per_block[0] << " " << dofs_per_block[1] << std::endl;

  std::vector<IndexSet> locally_relevant_partitioning;
  std::vector<IndexSet> locally_owned_partitioning;

  IndexSet locally_relevant_dofs;
  DoFTools::extract_locally_relevant_dofs (dof_handler, locally_relevant_dofs);
  locally_relevant_partitioning.push_back (locally_relevant_dofs.get_view (0, dofs_per_block[0]));
  locally_relevant_partitioning.push_back (locally_relevant_dofs.get_view (dofs_per_block[0], dofs_per_block[0] + dofs_per_block[1]));

  IndexSet locally_owned_dofs = dof_handler.locally_owned_dofs ();
  locally_owned_partitioning.push_back (locally_owned_dofs.get_view(0, dofs_per_block[0]));
  locally_owned_partitioning.push_back (locally_owned_dofs.get_view(dofs_per_block[0], dofs_per_block[0] + dofs_per_block[1]));

  deallog << "owned: ";
  locally_owned_dofs.print(deallog);
  deallog << "relevant: ";
  locally_relevant_dofs.print(deallog);


  ConstraintMatrix constraints(locally_relevant_dofs);
  constraints.close();

  BlockCompressedSimpleSparsityPattern bcsp (locally_relevant_partitioning);
  DoFTools::make_sparsity_pattern (dof_handler, bcsp, constraints, false);
  SparsityTools::distribute_sparsity_pattern (bcsp, dof_handler.locally_owned_dofs_per_processor (), MPI_COMM_WORLD, locally_relevant_dofs);

  typedef typename LA::MPI::BlockSparseMatrix::value_type number;
  typename LA::MPI::BlockSparseMatrix A;
  A.reinit(locally_owned_partitioning, bcsp, MPI_COMM_WORLD);

  QGauss<3> quadrature (3);
  FEValues<3> fe_values (fe, quadrature, update_values);

  std::vector<types::global_dof_index>  local_dof_indices (fe.dofs_per_cell);
  FullMatrix<number> local_matrix(fe.dofs_per_cell, fe.dofs_per_cell);

  for (DoFHandler<3>::active_cell_iterator cell = dof_handler.begin_active (); cell != dof_handler.end (); ++cell)
    if (cell->is_locally_owned ())
      {
        fe_values.reinit (cell);
        local_matrix = number();

        for (unsigned int q_point = 0; q_point < fe_values.n_quadrature_points; ++q_point)
          {
            for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
              {
                for (unsigned int j = 0; j < fe.dofs_per_cell; ++j)
                  local_matrix (i, j) += fe_values.shape_value(i, q_point) * fe_values.shape_value(j, q_point);
              }
          }

        cell->get_dof_indices (local_dof_indices);
        constraints.distribute_local_to_global (local_matrix, local_dof_indices, A);
      }

  A.compress(VectorOperation::add);
  A.print(deallog.get_file_stream());

  if (myid==0)
    deallog << "OK" << std::endl;
}

template <int dim>
void test_alt()
{
  typedef LA_Trilinos LA;
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);

  if (myid==0)
    deallog << "numproc=" << numproc << std::endl;

  parallel::distributed::Triangulation<3> triangulation(MPI_COMM_WORLD);
  GridGenerator::hyper_cube (triangulation, -1.0, 1.0);
  triangulation.refine_global(2);

  FESystem<3> fe(FE_Q<3> (1), 2);
  DoFHandler<3> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);
  DoFRenumbering::block_wise(dof_handler);

  std::vector<types::global_dof_index> dofs_per_block (fe.n_blocks());
  std::vector<unsigned int> sub_blocks (fe.n_blocks(), 0);
  sub_blocks[1] = 1;
  DoFTools::count_dofs_per_block (dof_handler, dofs_per_block, sub_blocks);

  deallog << "size: " << dofs_per_block[0] << " " << dofs_per_block[1] << std::endl;

  std::vector<IndexSet> locally_relevant_partitioning;
  std::vector<IndexSet> locally_owned_partitioning;

  IndexSet locally_relevant_dofs;
  DoFTools::extract_locally_relevant_dofs (dof_handler, locally_relevant_dofs);
  locally_relevant_partitioning.push_back (locally_relevant_dofs.get_view (0, dofs_per_block[0]));
  locally_relevant_partitioning.push_back (locally_relevant_dofs.get_view (dofs_per_block[0], dofs_per_block[0] + dofs_per_block[1]));

  IndexSet locally_owned_dofs = dof_handler.locally_owned_dofs ();
  locally_owned_partitioning.push_back (locally_owned_dofs.get_view(0, dofs_per_block[0]));
  locally_owned_partitioning.push_back (locally_owned_dofs.get_view(dofs_per_block[0], dofs_per_block[0] + dofs_per_block[1]));

  deallog << "owned: ";
  locally_owned_dofs.print(deallog);
  deallog << "relevant: ";
  locally_relevant_dofs.print(deallog);


  ConstraintMatrix constraints(locally_relevant_dofs);
  constraints.close();

  TrilinosWrappers::BlockSparsityPattern sp (locally_owned_partitioning,
                                             MPI_COMM_WORLD);
  DoFTools::make_sparsity_pattern (dof_handler, sp,
                                   constraints, false,
                                   Utilities::MPI::
                                   this_mpi_process(MPI_COMM_WORLD));
  sp.compress();
  typename LA::MPI::BlockSparseMatrix A;
  A.reinit(sp);

  QGauss<3> quadrature (3);
  FEValues<3> fe_values (fe, quadrature, update_values);

  std::vector<types::global_dof_index>  local_dof_indices (fe.dofs_per_cell);
  FullMatrix<double> local_matrix(fe.dofs_per_cell, fe.dofs_per_cell);

  for (DoFHandler<3>::active_cell_iterator cell = dof_handler.begin_active (); cell != dof_handler.end (); ++cell)
    if (cell->is_locally_owned ())
      {
        fe_values.reinit (cell);
        local_matrix = 0.0;

        for (unsigned int q_point = 0; q_point < fe_values.n_quadrature_points; ++q_point)
          {
            for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
              {
                for (unsigned int j = 0; j < fe.dofs_per_cell; ++j)
                  local_matrix (i, j) += fe_values.shape_value(i, q_point) * fe_values.shape_value(j, q_point);
              }
          }

        cell->get_dof_indices (local_dof_indices);
        constraints.distribute_local_to_global (local_matrix, local_dof_indices, A);
      }

  A.compress(VectorOperation::add);
  //A.print(deallog.get_file_stream());

  if (myid==0)
    deallog << "OK" << std::endl;
}


int main (int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, 1);
  MPILogInitAll log;

  {
    deallog.push("PETSc");
    test<LA_PETSc,3>();
    deallog.pop();
    deallog.push("Trilinos");
    test<LA_Trilinos,3>();
    deallog.pop();
    deallog.push("Trilinos_alt");
    test_alt<3>();
    deallog.pop();
  }

  // compile, don't run
  //if (myid==9999)
  //  test<LA_Dummy>();


}
