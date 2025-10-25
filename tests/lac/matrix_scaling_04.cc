// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2001 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test matrix scaling on block sparse linear systems. Check correctness of
// solution

#include <deal.II/base/numbers.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/matrix_scaling.h>
#include <deal.II/lac/sparse_direct.h>

#include "../tests.h"

using namespace dealii;

int
main(int argc, char **argv)
{
  initlog();

  const unsigned int                   n_blocks   = 5;
  const unsigned int                   block_size = 2;
  const unsigned int                   dim        = n_blocks * block_size;
  std::vector<types::global_dof_index> block_sizes(n_blocks, block_size);
  BlockIndices                         row_block_indices(block_sizes);
  BlockIndices                         col_block_indices(block_sizes);

  BlockDynamicSparsityPattern dsp(row_block_indices, col_block_indices);

  for (unsigned int block = 0; block < n_blocks; ++block)
    {
      for (unsigned int i = 0; i < block_size; ++i)
        dsp.block(block, block).add(i, i);
    }

  BlockSparsityPattern sp;
  sp.copy_from(dsp);
  BlockSparseMatrix<double> B(sp);
  BlockVector<double>       x(n_blocks, block_size), rhs(n_blocks, block_size);

  rhs = 1.0;

  for (unsigned int block = 0; block < n_blocks; ++block)
    {
      B.block(block, block).set(0, 0, 2 * block + 1.0);
      B.block(block, block).set(1, 1, 2 * block + 2.0);
    }

  MatrixScaling::AdditionalData control;
  MatrixScaling                 scaler(control);

  deallog << "Block Sparse Matrix: " << std::endl;

  std::ostringstream oss1, oss2;
  B.print(oss1);
  deallog << oss1.str();

  SparseDirectUMFPACK Binv;
  Binv.initialize(B);
  Binv.vmult(x, rhs);

  scaler.find_scaling_and_scale_linear_system(B, rhs);
  scaler.scale_system_solution(rhs);

  deallog << "Scaled Block Sparse Matrix: " << std::endl;
  B.print(oss2);
  deallog << oss2.str();

  const Vector<double> &row_scaling    = scaler.get_row_scaling();
  const Vector<double> &column_scaling = scaler.get_column_scaling();

  deallog << "Reciprocal of scaling vectors squared " << std::endl;

  for (unsigned int i = 0; i < dim; i++)
    deallog << 1.0 / (row_scaling[i] * row_scaling[i]) << " ";
  deallog << std::endl;
  for (unsigned int i = 0; i < dim; i++)
    deallog << 1.0 / (column_scaling[i] * column_scaling[i]) << " ";
  deallog << std::endl;

  deallog << "Solution of Ax=b: " << std::endl;
  for (unsigned int i = 0; i < dim; i++)
    deallog << x[i] << " ";
  deallog << std::endl;

  deallog << "Solution of Ax=b scaled: " << std::endl;
  for (unsigned int i = 0; i < dim; i++)
    deallog << rhs[i] << " ";
  deallog << std::endl;
}
