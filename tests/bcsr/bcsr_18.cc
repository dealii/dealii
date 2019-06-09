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

// test operator-= for a matrix with ghosts.
// partitioning is same as in bcsr_16.

// row (NON-BLOCK) partitioning and ghost blocks are as follows:

// rows (100)
// rank :-- owned / relevant
// block sizes
//
// 1:-- [0,32] / [0,35]
// 16
// 17
// 2:-- [33,65] / [0,1],[30,68]
// 16
// 17
// 3:-- [66,99] / [0,1],[63,99]
// 17
// 17

// column blocks are
// {3, 2, 2, 3}

// global sparsity in blocks is:
//
// x x 0 0
// x x 0 0
// x x x 0
// 0 x x x
// 0 x x x
// 0 0 x x

#include <deal.II/base/logstream.h>
#include <deal.II/lac/lapack_full_matrix.h>

#include "bcsr_helper.h"
#include <deal.II/lac/block_csr_matrix.h>

#include <fstream>
#include <iostream>


using namespace dealii;

void test()
{
  MPI_Comm mpi_communicator(MPI_COMM_WORLD);
  const unsigned int myid =
    dealii::Utilities::MPI::this_mpi_process(mpi_communicator);

  std::vector<IndexSet> local_support;
  std::vector<IndexSet> global_support;
  IndexSet locally_owned_dofs;
  IndexSet locally_relevant_dofs;
  IndexSet column_partitioning;

  setup_1d_sparsity(global_support,
                    locally_owned_dofs,
                    locally_relevant_dofs,
                    column_partitioning,
                    mpi_communicator);

  // add 0 and 1 on all processors
  locally_relevant_dofs.add_index(0);
  locally_relevant_dofs.add_index(1);

  for (auto g : global_support)
    local_support.push_back(g & locally_owned_dofs);

  std::shared_ptr<dealii::Utilities::MPI::Partitioner> partitioner_row =
    std::make_shared<dealii::Utilities::MPI::Partitioner>(
      locally_owned_dofs, locally_relevant_dofs, mpi_communicator);

  // setup 2 row blocks for local partitioning
  const std::vector<unsigned int> row_blocks = {
    {locally_owned_dofs.n_elements() / 2,
     locally_owned_dofs.n_elements() - locally_owned_dofs.n_elements() / 2}};

  //
  // setup block partitioner
  //

  IndexSet locally_owned_blocks(6);
  locally_owned_blocks.add_index(myid*2);
  locally_owned_blocks.add_index(myid*2+1);
  IndexSet locally_relevant_blocks(locally_owned_blocks);
  if (myid==0)
    {
      locally_relevant_blocks.add_index(2);
    }
  else if (myid == 1)
    {
      locally_relevant_blocks.add_index(1);
      locally_relevant_blocks.add_index(4);
    }
  else if (myid == 2)
    {
      locally_relevant_blocks.add_index(3);
    }

  std::shared_ptr<dealii::Utilities::MPI::Partitioner> partitioner =
    std::make_shared<dealii::Utilities::MPI::Partitioner>(
      locally_owned_blocks, locally_relevant_blocks, mpi_communicator);

  std::shared_ptr<dealii::Utilities::MPI::Partitioner> partitioner_owned =
    std::make_shared<dealii::Utilities::MPI::Partitioner>(
      locally_owned_blocks, mpi_communicator);

  // partition 10 columns as follows
  Assert(local_support.size() == 10, ExcNotImplemented());
  const std::vector<unsigned int> col_blocks = {{3, 2, 2, 3}};

  // setup sparsity pattern
  DynamicSparsityPattern local_dsp_A, global_dsp_A;

  setup_bcsr_sparsity(
    row_blocks, col_blocks, local_support, locally_owned_dofs, local_dsp_A);

  IndexSet all_dofs(100);
  all_dofs.add_range(0, 100);

  std::vector<unsigned int> all_row_blocks;
  {
    const auto gathered_vectors =
      Utilities::MPI::all_gather(mpi_communicator, row_blocks);
    for (auto v : gathered_vectors)
      for (auto el : v)
        all_row_blocks.push_back(el);
  }

  Assert(std::accumulate(all_row_blocks.begin(), all_row_blocks.end(), 0) ==
           100,
         ExcInternalError());

  setup_bcsr_sparsity(
    all_row_blocks, col_blocks, global_support, all_dofs, global_dsp_A);

  if (myid == 0)
    {
      SparsityPattern global_sp_A;
      global_sp_A.copy_from(global_dsp_A);

      const std::string filename = "sparsity_A.svg";
      std::ofstream f(filename.c_str());
      global_sp_A.print_svg(f);
    }

  // BlockIndices:
  std::shared_ptr<BlockIndices> rb = std::make_shared<BlockIndices>(row_blocks);
  std::shared_ptr<BlockIndices> cb = std::make_shared<BlockIndices>(col_blocks);

  SparsityPattern local_sp_A;
  local_sp_A.copy_from(local_dsp_A);

  deallog << "Local sparsity pattern of A:" << std::endl;
  local_sp_A.print(deallog.get_file_stream());
  deallog << "Local row blocks of A:" << std::endl
          << rb->to_string() << std::endl;

  // setup matrices
  BlockCSRMatrix<double> A, A_owned;
  A.reinit(local_dsp_A, rb, cb, partitioner);
  A_owned.reinit(local_dsp_A, rb, cb, partitioner_owned);

  const auto &ghost_sp_A = A.get_sparsity_pattern();
  deallog << "Ghost sparsity pattern of A:" << std::endl;
  ghost_sp_A.print(deallog.get_file_stream());
  deallog << "Row blocks with ghosts of A:" << std::endl
          << A.get_row_blocks()->to_string() << std::endl
          << "has_ghost_elements: " << A.has_ghost_elements() << std::endl
          << "is_block_partitioned: " << A.is_block_partitioned() << std::endl;

  // now set elements to something
  const auto start = partitioner_row->local_range().first;

  auto set_elements = [&](BlockCSRMatrix<double> &A_,
                          const double row_mult = 1.,
                          const double col_mult = 1000.,
                          const double shift = 0) -> void {
    for (unsigned int r = 0; r < rb->size(); ++r)
      {
        const auto end = A_.end_local(r);
        const auto row_start = start + rb->block_start(r);
        const auto row_size = rb->block_size(r);
        for (auto it = A_.begin_local(r); it != end; ++it)
          {
            const auto c = it->column();
            const auto col_start = cb->block_start(c);
            const auto col_size = cb->block_size(c);

            for (unsigned int ii = 0; ii < row_size; ++ii)
              for (unsigned int jj = 0; jj < col_size; ++jj)
                *(it->data() + BlockCSRMatrix<double>::local_index(
                                 ii, jj, row_size, col_size)) =
                  (row_start + ii + 1) * row_mult +
                  (col_start + jj + 1) * col_mult + shift;
          }
      }
  };

  set_elements(A, 1, 1000, 0);
  set_elements(A_owned, 1, 999, -1);
  const double row_mult_expect = 0;
  const double col_mult_expect = 1;
  const double shift_expect = 1;

  A.update_ghost_values();

  const std::ios::fmtflags old_flags = deallog.get_file_stream().flags();
  deallog.get_file_stream().setf(std::ios::fixed, std::ios::floatfield);

  deallog << "has_ghost_elements: " << A.has_ghost_elements() << std::endl;
  deallog << "before operator-=" << std::endl;
  A.print(deallog.get_file_stream(), 6, 0);

  A-= A_owned;

  deallog <<"after operator-=" << std::endl;
  A.print(deallog.get_file_stream(), 6, 0);

  deallog.get_file_stream().flags(old_flags);

  const BlockCSRMatrix<double> &A_const = A;
  const auto &rb_ghost = A_const.get_row_blocks();
  for (unsigned int r = 0; r < rb_ghost->size(); ++r)
    {
      const auto end = A_const.end_local(r);
      const auto data = A_const.get_block_data(r);
      std::vector<types::global_dof_index> row_indices;
      for (const auto &range : data.second)
        for (unsigned int ind = range.first; ind < range.second; ++ind)
          row_indices.push_back(data.first + ind);

      const auto row_size = rb_ghost->block_size(r);
      Assert(row_indices.size() == row_size, ExcInternalError());
      for (auto it = A_const.begin_local(r); it != end; ++it)
        {
          const auto c = it->column();
          const auto col_start = cb->block_start(c);
          const auto col_size = cb->block_size(c);

          for (unsigned int ii = 0; ii < row_size; ++ii)
            for (unsigned int jj = 0; jj < col_size; ++jj)
              {
                const double val = *(
                  it->data() + BlockCSRMatrix<double>::local_index(
                                 ii, jj, row_size, col_size));
                const double expect =
                  (row_indices[ii] + 1) * row_mult_expect + (col_start + jj + 1) * col_mult_expect + shift_expect;
                AssertThrow(val == expect,
                            ExcMessage(std::to_string(val) +
                                       " != " + std::to_string(expect)));
              }
        }
    }

  deallog << "Ok" << std::endl;
}

int main(int argc, char **argv)
{
  dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  unsigned int myid = dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  const unsigned int n_procs =
    dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  std::string deallogname = "output" + dealii::Utilities::int_to_string(myid);
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
        std::string line;
        while (std::getline(f, line))
          std::cout << p << ":" << line << std::endl;
      }

  return 0;
}
