// test iterators in the presence of empty rows

//     12   34    5  678
//      2    2    1   3
//  3                       123
//  2        x        x     45
//  1                       6
//  2   x    x              78

#include <deal.II/base/logstream.h>
#include <deal.II/lac/lapack_full_matrix.h>

#include <RFAStDFT/block_csr_matrix.h>

#include <fstream>
#include <iostream>
#include <numeric>

using namespace RealFAStDFT;
using namespace dealii;

template <typename Matrix>
void print_const(const Matrix &matrix)
{
  for (typename Matrix::const_iterator i = matrix.begin(); i != matrix.end(); ++i)
    // use both i-> and (*i)
    deallog << i->row() << ' ' << i->column() << ' ' << *(i->data())
            << std::endl;
  deallog << std::endl;
}

template <typename Matrix>
void print(Matrix &matrix)
{
  for (typename Matrix::iterator i = matrix.begin(); i != matrix.end(); ++i)
    // use both i-> and (*i)
    deallog << i->row() << ' ' << i->column() << ' ' << *(i->data())
            << std::endl;
  deallog << std::endl;
}

template <typename Matrix>
void print_const_rows(const Matrix &matrix)
{
  const auto n_rows = matrix.get_sparsity_pattern().n_rows();
  for (unsigned int r = 0; r < n_rows; ++r)
    for (typename Matrix::const_iterator i = matrix.begin_local(r); i != matrix.end_local(r); ++i)
    // use both i-> and (*i)
    deallog << i->row() << ' ' << i->column() << ' ' << *(i->data())
            << std::endl;
  deallog << std::endl;
}


template <typename Matrix>
void print_rows(Matrix &matrix)
{
  const auto n_rows = matrix.get_sparsity_pattern().n_rows();
  for (unsigned int r = 0; r < n_rows; ++r)
    for (typename Matrix::iterator i = matrix.begin_local(r); i != matrix.end_local(r); ++i)
    // use both i-> and (*i)
    deallog << i->row() << ' ' << i->column() << ' ' << *(i->data())
            << std::endl;
  deallog << std::endl;
}

void test()
{
  // number of blocks:
  const std::vector<unsigned int> row_blocks = {{3, 2, 1, 2}};
  const std::vector<unsigned int> col_blocks = {{2, 2, 1, 3}};
  const unsigned int M = row_blocks.size();
  const unsigned int N = col_blocks.size();

  std::vector<dealii::types::global_dof_index> row_offset;
  std::vector<dealii::types::global_dof_index> col_offset;

  auto setup_offset = [](const std::vector<unsigned int> &blocks,
                         std::vector<dealii::types::global_dof_index> &offset) {
    offset.resize(blocks.size() + 1, 0);
    std::partial_sum(blocks.begin(), blocks.end(), ++offset.begin());
  };

  setup_offset(row_blocks, row_offset);
  setup_offset(col_blocks, col_offset);

  deallog << "row blocks:";
  for (auto el: row_blocks)
    deallog << " " << el;
  deallog << std::endl;

  deallog << "col blocks:";
  for (auto el: col_blocks)
    deallog << " " << el;
  deallog << std::endl;

  deallog << "row offset:";
  for (auto el : row_offset)
    deallog << " " << el;
  deallog << std::endl;

  deallog << "col offset:";
  for (auto el : col_offset)
    deallog << " " << el;
  deallog << std::endl;

  DynamicSparsityPattern dsp(M, N);
  dsp.add(1, 1);
  dsp.add(1, 3);
  dsp.add(3, 0);
  dsp.add(3, 1);

  std::shared_ptr<BlockIndices> rb =
    std::make_shared<BlockIndices>(row_blocks);
  std::shared_ptr<BlockIndices> cb =
    std::make_shared<BlockIndices>(col_blocks);

  auto bcsr_block_part =
    std::make_shared<dealii::Utilities::MPI::Partitioner>(rb->size());

  // setup matrices
  BlockCSRMatrix<double> A;
  const BlockCSRMatrix<double> &A_const = A;
  A.reinit(dsp, rb, cb, bcsr_block_part);

  // setup
  {
    A(row_offset[1]+0,col_offset[1]+0) = 43;
    A(row_offset[1]+1,col_offset[1]+0) = 53;

    A(row_offset[1]+0,col_offset[1]+1) = 44;
    A(row_offset[1]+1,col_offset[1]+1) = 54;

    A(row_offset[1]+0,col_offset[3]+0) = 46;
    A(row_offset[1]+0,col_offset[3]+1) = 47;
    A(row_offset[1]+0,col_offset[3]+2) = 48;

    A(row_offset[1]+1,col_offset[3]+0) = 56;
    A(row_offset[1]+1,col_offset[3]+1) = 57;
    A(row_offset[1]+1,col_offset[3]+2) = 58;

    A(row_offset[3]+0,col_offset[0]+0) = 71;
    A(row_offset[3]+0,col_offset[0]+1) = 72;

    A(row_offset[3]+1,col_offset[0]+0) = 81;
    A(row_offset[3]+1,col_offset[0]+1) = 82;

    A(row_offset[3]+0,col_offset[1]+0) = 73;
    A(row_offset[3]+0,col_offset[1]+1) = 74;

    A(row_offset[3]+1,col_offset[1]+0) = 83;
    A(row_offset[3]+1,col_offset[1]+1) = 84;
  }

  deallog << "m: " << A.m() << std::endl << "n: " << A.n() << std::endl;
  deallog << "initial:" << std::endl;
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

  // now test:
  // print first element in each block:
  deallog << "begin():" << std::endl;
  print(A);
  deallog << "begin() const:" << std::endl;
  print_const(A);
  deallog << "begin(r):" << std::endl;
  print_rows(A);
  deallog << "begin(r) const:" << std::endl;
  print_const_rows(A);

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
