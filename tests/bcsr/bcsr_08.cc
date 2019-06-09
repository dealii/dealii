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

// test read/write iterators for BCSR
// in particular
//
// ++it
// it++
// (*it)
// it->...
// ==
// !=
// <
// >
// it2=it+num

//     12   34    5  678
//      2    2    1   3
//  3   x         x         123
//  2        x        x     45
//  1   x         x   x     6
//  2   x    x              78

#include <deal.II/base/logstream.h>
#include <deal.II/lac/lapack_full_matrix.h>

#include <deal.II/lac/block_csr_matrix.h>

#include <fstream>
#include <iostream>
#include <numeric>


using namespace dealii;

template <typename Matrix>
void print_const(const Matrix &matrix)
{
  for (typename Matrix::const_iterator i = matrix.begin(); i != matrix.end(); ++i)
    // use both i-> and (*i)
    deallog << i->row() << ' ' << i->column() << ' ' << *((*i).data())
            << std::endl;
  deallog << std::endl;
}

template <typename Matrix>
void print(Matrix &matrix)
{
  for (typename Matrix::iterator i = matrix.begin(); i != matrix.end(); ++i)
    // use both i-> and (*i)
    deallog << i->row() << ' ' << i->column() << ' ' << *((*i).data())
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
  dsp.add(1, 3);
  dsp.add(2, 0);
  dsp.add(2, 2);
  dsp.add(2, 3);
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
    A(row_offset[1]+1,col_offset[1]+0) = 53;

    A(row_offset[1]+0,col_offset[1]+1) = 44;
    A(row_offset[1]+1,col_offset[1]+1) = 54;

    A(row_offset[1]+0,col_offset[3]+0) = 46;
    A(row_offset[1]+0,col_offset[3]+1) = 47;
    A(row_offset[1]+0,col_offset[3]+2) = 48;

    A(row_offset[1]+1,col_offset[3]+0) = 56;
    A(row_offset[1]+1,col_offset[3]+1) = 57;
    A(row_offset[1]+1,col_offset[3]+2) = 58;

    A(row_offset[2]+0,col_offset[0]+0) = 61;
    A(row_offset[2]+0,col_offset[0]+1) = 62;

    A(row_offset[2]+0,col_offset[2]+0) = 65;

    A(row_offset[2]+0,col_offset[3]+0) = 66;
    A(row_offset[2]+0,col_offset[3]+1) = 67;
    A(row_offset[2]+0,col_offset[3]+2) = 68;

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
  // ++it  --- print first element in each block:
  deallog << "first element in each block:" << std::endl;
  print(A);

  // Multiply by 10 each first element in the block
  for (auto i = A.begin(); i != A.end(); ++i)
    *(i->data()) *= 10;

  deallog << "multiply first element in each block by 10:" << std::endl;
  print_const(A);

  // Subtract 1 from each first element of each block in row 1
  for (auto i = A.begin(1); i != A.end(1); ++i)
    *(i->data()) -= 1.;

  //  Double each first element of each block in row 3 (last row)
  for (auto i = A.begin(3); i != A.end(3); ++i)
    *(i->data()) *= 2;

  deallog << "subtract 1 from blockrow 1 and multiply by 2 blockrow 3:" << std::endl;
  print_const_rows(A);

  // it++, ==, !=
  deallog << "Test it++ and ++it ..." << std::endl;
  {
    auto b = A_const.begin();
    auto b2 = b;
    auto b3 = b++;
    Assert(b3 == b2, ExcInternalError());
    Assert(b3 != b, ExcInternalError());
    auto b4 = ++b2;
    Assert(b4 == b, ExcInternalError());
  }

  deallog << "Test >, < and it+num ..." << std::endl;
  {
    const auto b = A_const.begin();
    auto b2 = b + 2;
    auto b3 = b;
    ++b3;
    ++b3;
    Assert(b2 > b, ExcInternalError());
    Assert(b < b2, ExcInternalError());
    Assert(b2 == b3, ExcInternalError());
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
