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

// check BlockCSRMatrix::mmult() with square matrix B: C=A*B
// and MPI parallelization.
// otherwise similar to bcsr_03

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
  const unsigned int n_proc =
    dealii::Utilities::MPI::n_mpi_processes(mpi_communicator);

  // number of blocks:
  const unsigned int M = 33;
  const unsigned int N = 150;

  std::vector<unsigned int> M_blocks(M, numbers::invalid_unsigned_int);
  std::vector<unsigned int> N_blocks(N);

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

  randomize_blocks(N_blocks);
  if (myid == 0)
    randomize_blocks(M_blocks);

  const int ierr = MPI_Bcast(
    M_blocks.data(), M_blocks.size(), MPI_UNSIGNED, 0, mpi_communicator);
  AssertThrowMPI(ierr);

  DynamicSparsityPattern dsp_A(N,M);
  DynamicSparsityPattern dsp_B_global(M,M);
  DynamicSparsityPattern dsp_C(N,M);

  const auto randomize_sp = [](DynamicSparsityPattern &sp) {
    for (unsigned int i = 0; i < sp.n_rows(); ++i)
      for (unsigned int j = 0; j < sp.n_cols(); ++j)
        if (Utilities::generate_normal_random_number(0, 0.2) > 0)
          {
            sp.add(i,j);
          }
  };

  randomize_sp(dsp_A);
  if (myid == 0)
    randomize_sp(dsp_B_global);

  // broadcast global sparsity of B from root to others
  bcast_sp(dsp_B_global, mpi_communicator, 0);
  dsp_C.compute_mmult_pattern(dsp_A, dsp_B_global);

  // now do some partitioning of columns and take a view of the global sparsity
  // pattern
  IndexSet owned_columns;
  std::vector<unsigned int> M_blocks_local;
  get_owned_columns(owned_columns,
                    M_blocks_local,
                    M_blocks,
                    mpi_communicator);

  // columns we need to know about on this process
  IndexSet ghost_columns = columns(dsp_A);
  ghost_columns.add_indices(owned_columns);

  DynamicSparsityPattern dsp_B_local;
  get_view(dsp_B_local, dsp_B_global, owned_columns);

  std::shared_ptr<BlockIndices> Nb =
    std::make_shared<BlockIndices>(N_blocks);
  std::shared_ptr<BlockIndices> Mb =
    std::make_shared<BlockIndices>(M_blocks);
  std::shared_ptr<BlockIndices> Mb_local =
    std::make_shared<BlockIndices>(M_blocks_local);

  auto bcsr_block_part = std::make_shared<dealii::Utilities::MPI::Partitioner>(
    owned_columns, ghost_columns, mpi_communicator);

  const IndexSet owned_rows =
    get_owned_dofs(Nb->total_size(), mpi_communicator);
  auto bcsr_row_part = std::make_shared<dealii::Utilities::MPI::Partitioner>(
    owned_rows, mpi_communicator);

  // setup matrices
  BlockCSRMatrix<double> A, B, C;
  A.reinit(dsp_A, Nb, Mb, bcsr_row_part);
  B.reinit(dsp_B_local, Mb_local, Mb, bcsr_block_part);
  C.reinit(dsp_C, Nb, Mb, bcsr_row_part);

  // randomize content of matrices
  const auto randomize_mat = [](BlockCSRMatrix<double> &mat) {
    const auto & sp = mat.get_sparsity_pattern();
    for (unsigned int i = 0; i < sp.n_rows(); ++i)
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

  // A can be random on each MPI process
  randomize_mat(A);
  // B needs to be consistent
  init_bcsr(B, 0.01, 0.77, -22.0);

  A.mmult(C, B, false);

  // now compare to full matrices
  const auto full_M = std::accumulate(M_blocks.begin(), M_blocks.end(), 0);
  const auto full_N = dealii::Utilities::MPI::sum(
    std::accumulate(N_blocks.begin(), N_blocks.end(), 0), mpi_communicator);

  LAPACKFullMatrix<double> full_A(full_N, full_M), full_B(full_M,full_M), full_C(full_N,full_M), full_C_check(full_N,full_M);

  A.copy_to(full_A);
  B.copy_to(full_B);
  C.copy_to(full_C);

  full_A.mmult(full_C_check, full_B);

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
  dealii::deallog.attach(logfile, /*do not print job id*/ false);
  dealii::deallog.depth_console(0);

  test ();
}
