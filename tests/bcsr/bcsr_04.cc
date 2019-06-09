// check BlockCSRMatrix::Tmmult() for symmetric matrices

#include <deal.II/base/logstream.h>
#include <deal.II/lac/lapack_full_matrix.h>

#include <RFAStDFT/block_csr_matrix.h>

#include <fstream>
#include <iostream>
#include <numeric>

using namespace RealFAStDFT;
using namespace dealii;

void test ()
{
  // number of blocks:
  const unsigned int M = 100;
  const unsigned int N = 53;

  std::vector<unsigned int> M_blocks(M);
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

  randomize_blocks(M_blocks);
  randomize_blocks(N_blocks);

  const auto full_M = std::accumulate(M_blocks.begin(), M_blocks.end(), 0);
  const auto full_N = std::accumulate(N_blocks.begin(), N_blocks.end(), 0);

  DynamicSparsityPattern dsp_A(M,N);
  DynamicSparsityPattern dsp_C(N,N);

  const auto randomize_sp = [](DynamicSparsityPattern &sp) {
    for (unsigned int i = 0; i < sp.n_rows(); ++i)
      for (unsigned int j = 0; j < sp.n_cols(); ++j)
        if (Utilities::generate_normal_random_number(0, 0.2) > 0)
          {
            sp.add(i,j);
          }
  };

  randomize_sp(dsp_A);
  dsp_C.compute_Tmmult_pattern(dsp_A, dsp_A);

  std::shared_ptr<BlockIndices> Nb =
    std::make_shared<BlockIndices>(N_blocks);
  std::shared_ptr<BlockIndices> Mb =
    std::make_shared<BlockIndices>(M_blocks);

  auto bcsr_block_part =
    std::make_shared<dealii::Utilities::MPI::Partitioner>(Nb->size());

  auto bcsr_row_part =
    std::make_shared<dealii::Utilities::MPI::Partitioner>(Mb->total_size());

  // setup matrices
  BlockCSRMatrix<double> A, B, C;
  A.reinit(dsp_A, Mb, Nb, bcsr_row_part);
  B.reinit(dsp_A, Mb, Nb, bcsr_row_part);
  C.reinit(dsp_C, Nb, Nb, bcsr_block_part, true);

  // randomize content of matrices
  const auto randomize_mat = [](BlockCSRMatrix<double> &matA,
                                BlockCSRMatrix<double> &matB) {
    const auto & sp = matA.get_sparsity_pattern();
    for (unsigned int i = 0; i < sp.n_rows(); ++i)
      {
        const auto M = matA.get_row_blocks()->block_size(i);
        auto itA = matA.begin_local(i);
        auto itB = matB.begin_local(i);
        for (; itA != matA.end_local(i); ++itA, ++itB)
          {
            const auto j = itA->column();
            const auto N = matA.get_col_blocks()->block_size(j);
            unsigned int index = 0;
            for (unsigned int ii = 0; ii < M; ++ii)
              for (unsigned int jj = 0; jj < N; ++jj, ++index)
                {
                  const double v = Utilities::generate_normal_random_number(0, 0.2);
                  *(itA->data() + index) = v;
                  *(itB->data() + index) = v * 1.25;
                }
          }
      }
  };

  randomize_mat(A, B);

  A.Tmmult(C,B,false);

  deallog << "symmetric A: " << A.is_symmetric() << std::endl
          << "symmetric B: " << B.is_symmetric() << std::endl
          << "symmetric C: " << C.is_symmetric() << std::endl;

  // now compare to full matrices
  LAPACKFullMatrix<double> full_A(full_M, full_N), full_B(full_M,full_N), full_C(full_N,full_N), full_C_check(full_N,full_N);

  A.copy_to(full_A);
  B.copy_to(full_B);
  C.copy_to(full_C);

  full_A.Tmmult(full_C_check, full_B);

  full_C_check.add(-1, full_C);

  deallog << "diff norm: " << full_C_check.frobenius_norm() << std::endl;
}


int main(int argc, char **argv)
{
  dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  std::ofstream logfile("output");
  dealii::deallog.attach(logfile,/*do not print job id*/false);
  dealii::deallog.depth_console(0);

  test ();
}
