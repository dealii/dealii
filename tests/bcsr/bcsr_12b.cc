// test compress() for a matrix with MPI partitioner via rows.
// same layout as in bcsr_12.cc

// row partitioning and ghost blocks are as follows:

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

// therefore ghost blocks on each process are expected to be:
//
// 1:
// [33,35]            == 3
//
// 2:
// [0,1]              == 2
// [30,32]            == 3
// [66,68]            == 3
//
// 3:
// [0,1]              == 2
// [63,65]            == 3


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
#include <RFAStDFT/block_csr_matrix.h>

#include <fstream>
#include <iostream>

using namespace RealFAStDFT;
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

  std::shared_ptr<dealii::Utilities::MPI::Partitioner> partitioner =
    std::make_shared<dealii::Utilities::MPI::Partitioner>(
      locally_owned_dofs, locally_relevant_dofs, mpi_communicator);

  // setup 2 row blocks for local partitioning
  const std::vector<unsigned int> row_blocks = {
    {locally_owned_dofs.n_elements() / 2,
     locally_owned_dofs.n_elements() - locally_owned_dofs.n_elements() / 2}};

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
  BlockCSRMatrix<double> A;
  A.reinit(local_dsp_A, rb, cb, partitioner);

  const auto &ghost_sp_A = A.get_sparsity_pattern();
  deallog << "Ghost sparsity pattern of A:" << std::endl;
  ghost_sp_A.print(deallog.get_file_stream());
  deallog << "Row blocks with ghosts of A:" << std::endl
          << A.get_row_blocks()->to_string() << std::endl;
  deallog << "has_ghost_elements: " << A.has_ghost_elements() << std::endl;

  //===============================
  // here we differ from bcsr_12.cc
  //===============================

  // 1. set all elements including ghosts and run
  // compress(insert)
  const auto &rb_ghost = A.get_row_blocks();

  auto set_row = [&](const unsigned int r, const unsigned int shift = 0) {
    const auto end = A.end_local(r);
    const auto row_start =
      partitioner->local_to_global(rb_ghost->block_start(r));
    const auto row_size = rb_ghost->block_size(r);
    for (auto it = A.begin_local(r); it != end; ++it)
      {
        const auto c = it->column();
        const auto col_start = cb->block_start(c);
        const auto col_size = cb->block_size(c);

        for (unsigned int ii = 0; ii < row_size; ++ii)
          for (unsigned int jj = 0; jj < col_size; ++jj)
            *(it->data() + RealFAStDFT::BlockCSRMatrix<double>::local_index(
                             ii, jj, row_size, col_size)) =
              (row_start + ii + 1) + (col_start + jj + 1) * 1000 + shift;
      }
  };

  for (unsigned int r = 0; r < rb_ghost->size(); ++r)
    set_row(r);

  A.compress(VectorOperation::insert);

  // nothing should change so we can simply check all locally owned values:
  const auto start = partitioner->local_range().first;
  for (unsigned int r = 0; r < rb->size(); ++r)
    {
      const auto end = A.end_local(r);
      const auto row_start = start + rb->block_start(r);
      const auto row_size = rb->block_size(r);
      for (auto it = A.begin_local(r); it != end; ++it)
        {
          const auto c = it->column();
          const auto col_start = cb->block_start(c);
          const auto col_size = cb->block_size(c);

          for (unsigned int ii = 0; ii < row_size; ++ii)
            for (unsigned int jj = 0; jj < col_size; ++jj)
              AssertThrow(
                *(it->data() + RealFAStDFT::BlockCSRMatrix<double>::local_index(
                                 ii, jj, row_size, col_size)) ==
                  (row_start + ii + 1) + (col_start + jj + 1) * 1000,
                ExcInternalError());
        }
    }

  deallog.get_file_stream().setf(std::ios::fixed, std::ios::floatfield);

  deallog << "compress(insert):" << std::endl;
  A.print(deallog.get_file_stream(), 6, 0);

  // 2. reset matrix to zero and set only ghost elements to something
  // then call compress(add)
  A = 0.;

  // set owned the same on all cores
  for (unsigned int r = 0; r < rb->size(); ++r)
    set_row(r);

  // set ghost to different values on different cores
  const unsigned int p_mult = 100000;
  for (unsigned int r = rb->size(); r < rb_ghost->size(); ++r)
    set_row(r, (myid + 1) * p_mult);

  A.compress(VectorOperation::add);

  deallog << "compress(add):" << std::endl;
  A.print(deallog.get_file_stream(), 8, 0);

  // check values
  for (unsigned int r = 0; r < rb->size(); ++r)
    {
      const auto end = A.end_local(r);
      const auto row_start = start + rb->block_start(r);
      const auto row_size = rb->block_size(r);
      for (auto it = A.begin_local(r); it != end; ++it)
        {
          const auto c = it->column();
          const auto col_start = cb->block_start(c);
          const auto col_size = cb->block_size(c);

          for (unsigned int ii = 0; ii < row_size; ++ii)
            for (unsigned int jj = 0; jj < col_size; ++jj)
              {
                const double &val = *(
                  it->data() + RealFAStDFT::BlockCSRMatrix<double>::local_index(
                                 ii, jj, row_size, col_size));

                const double global_val =
                  (row_start + ii + 1) + (col_start + jj + 1) * 1000;
                // her we do a bit of hard-coding for the given partitioner and
                // 3 MPI processes:
                if (myid == 0)
                  {
                    if (r == 0 && ii < 2)
                      // summed from two processes (2nd and 3rd)
                      {
                        AssertThrow(val == (3 * global_val + 5 * p_mult),
                                    ExcInternalError());
                      }
                    else if (r == 1 && ii >= 14)
                      // summed from 2nd process
                      {
                        AssertThrow(val == (2 * global_val + 2 * p_mult),
                                    ExcInternalError());
                      }
                    else
                      // rest should be the same
                      {
                        AssertThrow(val == global_val, ExcInternalError());
                      }
                  }
                else if (myid == 1)
                  {
                    if (r == 0 && ii < 3)
                      // added from 1st proc
                      {
                        AssertThrow(val == (2 * global_val + 1 * p_mult),
                                    ExcInternalError());
                      }
                    else if (r == 1 && ii >= 14)
                      // added from 3rd proc
                      {
                        AssertThrow(val == (2 * global_val + 3 * p_mult),
                                    ExcInternalError());
                      }
                    else
                      // rest should be the same
                      {
                        AssertThrow(val == global_val, ExcInternalError());
                      }
                  }
                else
                  {
                    Assert(myid == 2, ExcNotImplemented());
                    if (r == 0 && ii < 3)
                      // added from 2st proc
                      {
                        AssertThrow(val == (2 * global_val + 2 * p_mult),
                                    ExcInternalError());
                      }
                    else
                      // rest should be zero
                      {
                        AssertThrow(val == global_val, ExcInternalError());
                      }
                  }
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
