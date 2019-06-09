// test compress() for a matrix with MPI partitioner via blocks.
// same as bcsr_16.cc but tests compress()

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
          << A.get_row_blocks()->to_string() << std::endl
          << "has_ghost_elements: " << A.has_ghost_elements() << std::endl
          << "is_block_partitioned: " << A.is_block_partitioned() << std::endl;

  //===============================
  // here we differ from bcsr_16.cc
  //===============================

  const auto &rb_ghost = A.get_row_blocks();

  auto set_row = [&](const unsigned int r,
                     const unsigned int shift = 0,
                     const unsigned int mult = 1,
                     const bool test = false) {
    const auto end = A.end_local(r);
    const auto data = A.get_block_data(r);
    std::vector<types::global_dof_index> row_indices;
    for (const auto &range : data.second)
      for (unsigned int ind = range.first; ind < range.second; ++ind)
        row_indices.push_back(data.first + ind);

    const auto row_size = rb_ghost->block_size(r);
    Assert(row_indices.size() == row_size, ExcInternalError());
    for (auto it = A.begin_local(r); it != end; ++it)
      {
        const auto c = it->column();
        const auto col_start = cb->block_start(c);
        const auto col_size = cb->block_size(c);

        for (unsigned int ii = 0; ii < row_size; ++ii)
          for (unsigned int jj = 0; jj < col_size; ++jj)
            {
              double &val =
                *(it->data() + RealFAStDFT::BlockCSRMatrix<double>::local_index(
                                 ii, jj, row_size, col_size));
              const double set_val =
                (row_indices[ii] + 1) + (col_start + jj + 1) * 1000;
              const double expect = set_val * mult + shift;
              if (test)
                {
                  AssertThrow(val == expect,
                              ExcMessage(std::to_string(val) +
                                         " != " + std::to_string(expect)));
                }
              else
                {
                  val = expect;
                }
            }
      }
  };

  // 1. set all elements including ghosts and run
  // compress(insert)
  for (unsigned int r = 0; r < rb_ghost->size(); ++r)
    set_row(r);

  A.compress(VectorOperation::insert);

  deallog.get_file_stream().setf(std::ios::fixed, std::ios::floatfield);

  deallog << "compress(insert):" << std::endl;
  deallog << "has_ghost_elements: " << A.has_ghost_elements() << std::endl;
  A.print(deallog.get_file_stream(), 6, 0);

  // nothing should change so we can simply check all locally owned values:
  for (unsigned int r = 0; r < rb->size(); ++r)
    set_row(r, 0, 1, true);

  // 2. reset matrix to zero and set elements to something
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
  deallog << "has_ghost_elements: " << A.has_ghost_elements() << std::endl;
  A.print(deallog.get_file_stream(), 8, 0);

  // check values
  if (myid == 0)
    {
      set_row(0, 0, 1, true); // row 0 unchanged
      set_row(1, p_mult*2, 2, true); // row 1 added from proc 1
    }
  else if (myid == 1)
    {
      set_row(0, p_mult, 2, true); // row 0 added from proc 0
      set_row(1, 3*p_mult, 2, true); // row 0 added from proc 2
    }
  else
    {
      set_row(0, 2*p_mult, 2, true); // row 0 added from proc 1
      set_row(1, 0, 1, true); // row 1 unchanged
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
