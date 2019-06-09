// test operator() access to actual stored numbers

#include <deal.II/base/logstream.h>
#include <deal.II/lac/lapack_full_matrix.h>

#include <RFAStDFT/block_csr_matrix.h>

#include <fstream>
#include <iostream>
#include <numeric>

using namespace RealFAStDFT;
using namespace dealii;

void test()
{
  // number of blocks:
  const std::vector<unsigned int> row_blocks = {{3, 2, 1}};
  const std::vector<unsigned int> col_blocks = {{2, 2, 1}};
  const unsigned int M = row_blocks.size();
  const unsigned int N = col_blocks.size();

  std::vector<dealii::types::global_dof_index> row_offset;
  std::vector<dealii::types::global_dof_index> col_offset;

  auto setup_offset = [](const std::vector<unsigned int> &blocks,
                         std::vector<dealii::types::global_dof_index> &offset) {
      offset.resize(blocks.size()+1, 0);
      std::partial_sum(blocks.begin(), blocks.end(), ++offset.begin());
    };

  setup_offset(row_blocks, row_offset);
  setup_offset(col_blocks, col_offset);

  deallog << "row offset:";
  for (auto el : row_offset)
    deallog << " " << el;
  deallog << std::endl;

  deallog << "col offset:";
  for (auto el : col_offset)
    deallog << " " << el;
  deallog << std::endl;

  DynamicSparsityPattern dsp(M, N);
  dsp.add(0, 0);
  dsp.add(0, 2);
  dsp.add(1, 1);
  dsp.add(2, 0);
  dsp.add(2, 2);

  std::shared_ptr<BlockIndices> rb =
    std::make_shared<BlockIndices>(row_blocks);
  std::shared_ptr<BlockIndices> cb =
    std::make_shared<BlockIndices>(col_blocks);

  auto bcsr_block_part =
    std::make_shared<dealii::Utilities::MPI::Partitioner>(rb->size());

  // setup matrices
  BlockCSRMatrix<double> A;
  A.reinit(dsp, rb, cb, bcsr_block_part);

  // setup
  {
    A(row_offset[0]+0,col_offset[0]+0) = 11;
    A(row_offset[0]+1,col_offset[0]+0) = 21;
    A(row_offset[0]+2,col_offset[0]+0) = 31;

    A(row_offset[0]+0,col_offset[0]+1) = 12;
    A(row_offset[0]+1,col_offset[0]+1) = 22;
    A(row_offset[0]+2,col_offset[0]+1) = 32;

    A(row_offset[0]+0,col_offset[2]+0) = 15;
    A(row_offset[0]+1,col_offset[2]+0) = 25;
    A(row_offset[0]+2,col_offset[2]+0) = 35;

    A(row_offset[1]+0,col_offset[1]+0) = 43;
    A(row_offset[1]+1,col_offset[1]+0) = 54;

    A(row_offset[1]+0,col_offset[1]+1) = 44;
    A(row_offset[1]+1,col_offset[1]+1) = 55;

    A(row_offset[2]+0,col_offset[0]+0) = 61;
    A(row_offset[2]+0,col_offset[0]+1) = 62;

    A(row_offset[2]+0,col_offset[2]+0) = 65;
  }

  deallog << "nonconst:" << std::endl;
  for (auto it = dsp.begin(); it != dsp.end(); ++it)
    {
      const auto &r = it->row();
      const auto &c = it->column();
      deallog << "block(" << r << "," << c << ")" << std::endl;
      for (unsigned int i = 0; i < row_blocks[r]; ++i)
        {
          for (unsigned int j = 0; j < col_blocks[c]; ++j)
            deallog << " " << A(row_offset[r] + i, col_offset[c] + j);

          deallog << std::endl;
        }

      deallog << std::endl;
    }

  deallog << "const:" << std::endl;
  const BlockCSRMatrix<double> &A_const = A;
  for (auto it = dsp.begin(); it != dsp.end(); ++it)
    {
      const auto &r = it->row();
      const auto &c = it->column();
      deallog << "block(" << r << "," << c << ")" << std::endl;
      for (unsigned int i = 0; i < row_blocks[r]; ++i)
        {
          for (unsigned int j = 0; j < col_blocks[c]; ++j)
            deallog << " " << A_const(row_offset[r] + i, col_offset[c] + j);

          deallog << std::endl;
        }

      deallog << std::endl;
    }

  deallog << "el():" << std::endl;
  const auto full_M =
    std::accumulate(row_blocks.begin(), row_blocks.end(), (unsigned int)0);
  const auto full_N =
    std::accumulate(col_blocks.begin(), col_blocks.end(), (unsigned int)0);
  for (unsigned int i = 0; i < full_M; ++i)
    {
      for (unsigned int j = 0; j < full_N; ++j)
        deallog << " " << A_const.el(i, j);

      deallog << std::endl;
    }

  deallog << "Ok" << std::endl;
}

int main(int argc, char **argv)
{
  dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  std::ofstream logfile("output");
  dealii::deallog.attach(logfile, /*do not print job id*/ false);
  dealii::deallog.depth_console(0);

  test();
}
