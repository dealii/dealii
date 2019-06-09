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

// check BlockCSRMatrix::Tr_Tmmult() when the result of Tmmult is a square matrix
// similar to bcsr_07b but simplified to small manually built sparsity patterns,
// that are easy to debug. In particular we construct A and B such that
// On 0-th row only the last element matches in sparsities
// On 1-st row there are no matches
// On 2-nd row one of the matrices is empty

#include <deal.II/base/logstream.h>
#include <deal.II/lac/lapack_full_matrix.h>

#include <deal.II/lac/block_csr_matrix.h>

#include <fstream>
#include <iostream>
#include <numeric>


using namespace dealii;

void test ()
{
  // number of blocks:
  const unsigned int M = 3;
  const unsigned int N = 5;

  std::vector<unsigned int> M_blocks = {2,3,5};
  std::vector<unsigned int> N_blocks = {2,3,4,5,6};

  const auto full_M = std::accumulate(M_blocks.begin(), M_blocks.end(), 0);
  const auto full_N = std::accumulate(N_blocks.begin(), N_blocks.end(), 0);

  DynamicSparsityPattern dsp_A(M,N);
  DynamicSparsityPattern dsp_B(M,N);
  DynamicSparsityPattern dsp_C(N,N);

  // A:
  //     0  1  2  3  4
  // 0      x     x  x
  // 1  x      x     x
  // 2  x
  dsp_A.add(0,1);
  dsp_A.add(0,3);
  dsp_A.add(0,4);
  dsp_A.add(1,0);
  dsp_A.add(1,2);
  dsp_A.add(1,4);
  dsp_A.add(2,0);

  // B:
  //     0  1  2  3  4
  // 0   x     x     x
  // 1      x     x
  // 2
  dsp_B.add(0,0);
  dsp_B.add(0,2);
  dsp_B.add(0,4);
  dsp_B.add(1,1);
  dsp_B.add(1,3);

  dsp_C.compute_Tmmult_pattern(dsp_A, dsp_B);

  std::shared_ptr<BlockIndices> Nb =
    std::make_shared<BlockIndices>(N_blocks);
  std::shared_ptr<BlockIndices> Mb =
    std::make_shared<BlockIndices>(M_blocks);

  auto bcsr_row_part =
    std::make_shared<dealii::Utilities::MPI::Partitioner>(Mb->total_size());

  auto bcsr_block_part =
    std::make_shared<dealii::Utilities::MPI::Partitioner>(Nb->size());

  // setup matrices
  BlockCSRMatrix<double> A, B, C;
  A.reinit(dsp_A, Mb, Nb, bcsr_row_part);
  B.reinit(dsp_B, Mb, Nb, bcsr_row_part);
  C.reinit(dsp_C, Nb, Nb, bcsr_block_part);

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

  randomize_mat(A);
  randomize_mat(B);

  A.Tmmult(C,B,false);

  deallog << "symmetric A: " << A.is_symmetric() << std::endl
          << "symmetric B: " << B.is_symmetric() << std::endl
          << "symmetric C: " << C.is_symmetric() << std::endl;

  const double trace = A.Tr_Tmmult(B);

  // now compare to full matrices
  LAPACKFullMatrix<double> full_A(full_M, full_N), full_B(full_M,full_N), full_C(full_N,full_N), full_C_check(full_N,full_N);

  A.copy_to(full_A);
  B.copy_to(full_B);
  C.copy_to(full_C);

  full_A.Tmmult(full_C_check, full_B);

  const double trace_full = full_C_check.trace();

  full_C_check.add(-1, full_C);

  deallog << "diff norm:  " << full_C_check.frobenius_norm() << std::endl;

  deallog << "diff trace: " << std::abs(trace_full - trace)/std::abs(std::max(trace_full, trace)) << std::endl;
}


int main(int argc, char **argv)
{
  dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  std::ofstream logfile("output");
  dealii::deallog.attach(logfile,/*do not print job id*/false);
  dealii::deallog.depth_console(0);

  test ();
}
