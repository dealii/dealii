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

// check BlockCSRMatrix::Tmmult() and Tr_Tmmult() in parallel using 1D example

#include <deal.II/base/logstream.h>

#include <deal.II/lac/block_csr_matrix.h>
#include <deal.II/lac/lapack_full_matrix.h>

#include <fstream>
#include <iostream>

#include "bcsr_helper.h"


using namespace dealii;

void
test()
{
  MPI_Comm           mpi_communicator(MPI_COMM_WORLD);
  const unsigned int myid =
    dealii::Utilities::MPI::this_mpi_process(mpi_communicator);
  const unsigned int n_proc =
    dealii::Utilities::MPI::n_mpi_processes(mpi_communicator);


  std::vector<IndexSet> local_support;
  std::vector<IndexSet> global_support;
  IndexSet              locally_owned_dofs;
  IndexSet              locally_relevant_dofs;
  IndexSet              column_partitioning;

  setup_1d_sparsity(global_support,
                    locally_owned_dofs,
                    locally_relevant_dofs,
                    column_partitioning,
                    mpi_communicator);

  for (auto g : global_support)
    local_support.push_back(g & locally_owned_dofs);

  const std::shared_ptr<dealii::Utilities::MPI::Partitioner> partitioner =
    std::make_shared<dealii::Utilities::MPI::Partitioner>(locally_owned_dofs,
                                                          locally_relevant_dofs,
                                                          mpi_communicator);

  // setup 2 row blocks for local row partitioning
  const std::vector<unsigned int> row_blocks_local = {
    {locally_owned_dofs.n_elements() / 2,
     locally_owned_dofs.n_elements() - locally_owned_dofs.n_elements() / 2}};

  deallog << "Column partitioning:" << std::endl;
  column_partitioning.print(deallog.get_file_stream());

  // partition 10 columns as follows
  std::vector<unsigned int> col_blocks_serial = {{3, 2, 2, 3}};
  Assert(local_support.size() == 10, ExcNotImplemented());
  const std::vector<unsigned int> col_blocks_local =
    n_proc == 1 ? col_blocks_serial :
                  get_local_col_blocks(column_partitioning, 2);

  std::vector<unsigned int> col_blocks;
  IndexSet                  owned_col_blocks;

  setup_column_blocks(owned_col_blocks,
                      col_blocks,
                      col_blocks_local,
                      mpi_communicator);

  deallog << "Column blocks:" << std::endl;
  for (const auto &b : col_blocks)
    deallog << " " << b;
  deallog << std::endl;

  deallog << "Owned column blocks:" << std::endl;
  owned_col_blocks.print(deallog.get_file_stream());

  // setup sparsity pattern
  DynamicSparsityPattern local_dsp_A;

  setup_bcsr_sparsity(row_blocks_local,
                      col_blocks,
                      local_support,
                      locally_owned_dofs,
                      local_dsp_A);

  // BlockIndices:
  std::shared_ptr<BlockIndices> rb =
    std::make_shared<BlockIndices>(row_blocks_local);
  std::shared_ptr<BlockIndices> cb = std::make_shared<BlockIndices>(col_blocks);
  std::shared_ptr<BlockIndices> cb_local =
    std::make_shared<BlockIndices>(col_blocks_local);

  DynamicSparsityPattern dsp_C_ghost(col_blocks.size(), col_blocks.size());

  dsp_C_ghost.compute_Tmmult_pattern(local_dsp_A, local_dsp_A);

  // print ghost dsp_C:
  {
    SparsityPattern sp;
    sp.copy_from(dsp_C_ghost);

    const std::string filename =
      "sparsity_C_ghost_" + std::to_string(myid) + ".svg";
    std::ofstream f(filename.c_str());
    sp.print_svg(f);
  }

  DynamicSparsityPattern dsp_C_global;
  gather_sparsity(dsp_C_global, dsp_C_ghost, mpi_communicator);

  // print local dsp_C_global:
  if (myid == 0)
    {
      SparsityPattern sp;
      sp.copy_from(dsp_C_global);

      const std::string filename = "sparsity_C_global.svg";
      std::ofstream     f(filename.c_str());
      sp.print_svg(f);
    }

  std::shared_ptr<dealii::Utilities::MPI::Partitioner> col_partitioner;
  DynamicSparsityPattern                               dsp_C;

  setup_col_partitioner_and_sparsity(col_partitioner,
                                     dsp_C,
                                     dsp_C_global,
                                     owned_col_blocks,
                                     local_dsp_A,
                                     mpi_communicator);

  // print local dsp_C:
  {
    SparsityPattern sp;
    sp.copy_from(dsp_C);

    const std::string filename =
      "sparsity_C_local_" + std::to_string(myid) + ".svg";
    std::ofstream f(filename.c_str());
    sp.print_svg(f);
  }


  // setup matrices
  BlockCSRMatrix<double> A, B, C;
  A.reinit(local_dsp_A, rb, cb, partitioner);
  B.reinit(local_dsp_A, rb, cb, partitioner);
  C.reinit(dsp_C, cb_local, cb, col_partitioner);

  init_bcsr(A, 0.01, 0.3, -0.12);
  init_bcsr(B, 0.02, 0.7, 0.22);

  // call Tmmult:
  A.Tmmult(C, B, false);

  const double trace = A.Tr_Tmmult(B);

  // now compare to full matrices
  const auto full_M = locally_owned_dofs.size();
  const auto full_N = column_partitioning.size();

  LAPACKFullMatrix<double> full_A(full_M, full_N), full_B(full_M, full_N),
    full_C(full_N, full_N), full_C_check(full_N, full_N);

  A.copy_to(full_A);
  B.copy_to(full_B);
  C.copy_to(full_C);

  full_A.Tmmult(full_C_check, full_B);

  const double trace_full = full_C_check.trace();

  full_C_check.add(-1, full_C);

  deallog << "diff norm: " << full_C_check.frobenius_norm() << std::endl;
  deallog << "diff trace: " << std::fabs(trace - trace_full) << std::endl;

  deallog << "Ok" << std::endl;
}

int
main(int argc, char **argv)
{
  dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  unsigned int myid = dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  const unsigned int n_procs =
    dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  std::string   deallogname = "output" + dealii::Utilities::int_to_string(myid);
  std::ofstream logfile(deallogname);
  dealii::deallog.attach(logfile, /*do not print job id*/ false);
  dealii::deallog.depth_console(0);

  test();

  logfile.close();

  MPI_Barrier(MPI_COMM_WORLD);

  if (myid == 0)
    for (unsigned int p = 0; p < n_procs; ++p)
      {
        std::string deallogname =
          "output" + dealii::Utilities::int_to_string(p);
        std::ifstream f(deallogname);
        std::string   line;
        while (std::getline(f, line))
          std::cout << p << ":" << line << std::endl;
      }

  return 0;
}
