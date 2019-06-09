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

// check BlockCSRMatrix::Tmmult() in parallel with MPI

#include <deal.II/base/logstream.h>
#include <deal.II/lac/lapack_full_matrix.h>

#include "bcsr_helper.h"
#include <deal.II/lac/block_csr_matrix.h>

#include <fstream>
#include <iostream>


using namespace dealii;

void test ()
{
  MPI_Comm mpi_communicator(MPI_COMM_WORLD);
  const unsigned int myid =
    dealii::Utilities::MPI::this_mpi_process(mpi_communicator);

  // number of blocks:
  const unsigned int M = 100;
  const unsigned int N = 50;
  const unsigned int O = 75;

  std::vector<unsigned int> M_blocks(M);
  std::vector<unsigned int> N_blocks(N, numbers::invalid_unsigned_int);
  std::vector<unsigned int> O_blocks(O, numbers::invalid_unsigned_int);

  // randomize block sizes
  const auto randomize_blocks = [](std::vector<unsigned int> &blocks) {
    const double s = 1.;
    for (auto &el : blocks)
      {
        const double v = Utilities::generate_normal_random_number(0, s);
        if (v > s)
          el = 11;
        else if (v > 0)
          el = 4;
        else if (v < -s)
          el = 1;
        else
          el = 6;
      }
  };

  randomize_blocks(M_blocks);
  if (myid == 0)
    {
      randomize_blocks(N_blocks);
      randomize_blocks(O_blocks);
    }

  int ierr = MPI_Bcast(
    N_blocks.data(), N_blocks.size(), MPI_UNSIGNED, 0, mpi_communicator);
  AssertThrowMPI(ierr);

  ierr = MPI_Bcast(
    O_blocks.data(), O_blocks.size(), MPI_UNSIGNED, 0, mpi_communicator);
  AssertThrowMPI(ierr);

  DynamicSparsityPattern dsp_A(M,N);
  DynamicSparsityPattern dsp_B(M,O);
  DynamicSparsityPattern dsp_C_ghost(N,O);

  const auto randomize_sp = [](DynamicSparsityPattern &sp) {
    for (unsigned int i = 0; i < sp.n_rows(); ++i)
      for (unsigned int j = 0; j < sp.n_cols(); ++j)
        if (Utilities::generate_normal_random_number(0, 0.2) > 0)
          {
            sp.add(i,j);
          }
  };

  randomize_sp(dsp_A);
  randomize_sp(dsp_B);
  dsp_C_ghost.compute_Tmmult_pattern(dsp_A, dsp_B);

  IndexSet owned_columns;
  std::vector<unsigned int> N_blocks_local;
  get_owned_columns(owned_columns,
                    N_blocks_local,
                    N_blocks,
                    mpi_communicator);

  // get local sparsity based on Tmmult from each MPI process
  DynamicSparsityPattern dsp_C;
  gather(dsp_C, dsp_C_ghost, owned_columns, mpi_communicator);

  std::shared_ptr<BlockIndices> Nb_local =
    std::make_shared<BlockIndices>(N_blocks_local);
  std::shared_ptr<BlockIndices> Nb =
    std::make_shared<BlockIndices>(N_blocks);
  std::shared_ptr<BlockIndices> Mb =
    std::make_shared<BlockIndices>(M_blocks);
  std::shared_ptr<BlockIndices> Ob =
    std::make_shared<BlockIndices>(O_blocks);

  const IndexSet owned_rows =
    get_owned_dofs(Mb->total_size(), mpi_communicator);
  auto bcsr_row_part = std::make_shared<dealii::Utilities::MPI::Partitioner>(
    owned_rows, mpi_communicator);

  // columns we need to know about on this process
  IndexSet ghost_columns = columns(dsp_A);
  ghost_columns.add_indices(owned_columns);

  auto bcsr_block_part = std::make_shared<dealii::Utilities::MPI::Partitioner>(
    owned_columns, ghost_columns, mpi_communicator);

  // setup matrices
  BlockCSRMatrix<double> A, B, C;
  A.reinit(dsp_A, Mb, Nb, bcsr_row_part);
  B.reinit(dsp_B, Mb, Ob, bcsr_row_part);
  C.reinit(dsp_C, Nb_local, Ob, bcsr_block_part);

  // randomize content of matrices
  const auto randomize_mat = [](BlockCSRMatrix<double> &mat) {
    for (unsigned int i = 0; i < mat.n_local_row_blocks(); ++i)
      {
        const auto M = mat.get_row_blocks()->block_size(i);
        for (auto it = mat.begin_local(i); it != mat.end_local(i); ++it)
          {
            const auto j = it->column();
            const auto N = mat.get_col_blocks()->block_size(j);
            unsigned int index = 0;
            for (unsigned int ii = 0; ii < M; ++ii)
              for (unsigned int jj = 0; jj < N; ++jj, ++index)
                *(it->data() + index) = Utilities::generate_normal_random_number(0, 0.2);
          }
      }
  };

  // both A and B can be random locally
  randomize_mat(A);
  randomize_mat(B);

  A.Tmmult(C,B,false);

  // now compare to full matrices
  const auto full_M = dealii::Utilities::MPI::sum(
    std::accumulate(M_blocks.begin(), M_blocks.end(), 0), mpi_communicator);

  const auto full_N = std::accumulate(N_blocks.begin(), N_blocks.end(), 0);
  const auto full_O = std::accumulate(O_blocks.begin(), O_blocks.end(), 0);

  LAPACKFullMatrix<double> full_A(full_M, full_N), full_B(full_M,full_O), full_C(full_N,full_O), full_C_check(full_N,full_O);

  A.copy_to(full_A);
  B.copy_to(full_B);
  C.copy_to(full_C);

  full_A.Tmmult(full_C_check, full_B);

  full_C_check.add(-1, full_C);

  const auto diff =
    Utilities::MPI::max(full_C_check.frobenius_norm(), mpi_communicator);
  if (myid == 0)
    deallog << "diff norm: " << diff << std::endl;
}


int main(int argc, char **argv)
{
  dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  const unsigned int myid =
    dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  std::ofstream logfile(myid == 0 ? "output"
                                  : "output_" + std::to_string(myid));
  dealii::deallog.attach(logfile,/*do not print job id*/false);
  dealii::deallog.depth_console(0);

  test ();
}
