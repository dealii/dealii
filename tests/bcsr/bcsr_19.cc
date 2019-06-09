// test frobenius_norm() for a matrix with MPI partitioner via rows.
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

  // BlockIndices:
  std::shared_ptr<BlockIndices> rb = std::make_shared<BlockIndices>(row_blocks);
  std::shared_ptr<BlockIndices> cb = std::make_shared<BlockIndices>(col_blocks);

  SparsityPattern local_sp_A;
  local_sp_A.copy_from(local_dsp_A);

  // setup matrices
  BlockCSRMatrix<double> A;
  A.reinit(local_dsp_A, rb, cb, partitioner);

  // now set elements to something
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
              *(it->data() + RealFAStDFT::BlockCSRMatrix<double>::local_index(
                               ii, jj, row_size, col_size)) =
                (row_start + ii + 1) * 0.25 + (col_start + jj + 1) * 0.37;
        }
    }

  // update ghost values to make sure frobenius_norm() on each process
  // does not take into account ghosts
  A.update_ghost_values();

  FullMatrix<double> A_full(100,10);
  A.copy_to(A_full);

  const auto bcsr_frob = A.frobenius_norm();
  const auto full_frob = A_full.frobenius_norm();
  const auto diff  = std::abs(bcsr_frob - full_frob);
  deallog << "BCSR: " << bcsr_frob << std::endl
          << "Full: " << full_frob << std::endl
          << "Diff: " << diff << std::endl;

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
