// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test the function projection_onto_orthogonal_tensors

#include <deal.II/base/tensor.h>

#include "../tests.h"

template <int dim, typename number>
void
print_type(const Tensor<2, dim, number> &)
{
  if (std::is_same_v<number, float>)
    deallog << " Tensor<2, " << dim << ", float>" << std::endl;
  else if (std::is_same_v<number, double>)
    deallog << " Tensor<2, " << dim << ", double>" << std::endl;
  else
    DEAL_II_NOT_IMPLEMENTED();
}

template <int dim, typename number>
void
test(const Tensor<2, dim, number> &A)
{
  print_type(A);
  deallog << "original tensor A = " << A << std::endl;
  deallog << "det A = " << determinant(A) << std::endl;
  deallog << "A A^T = " << (A * transpose(A)) << std::endl;
  const Tensor<2, dim, number> Q = project_onto_orthogonal_tensors(A);
  deallog << "projected to Q = " << Q << std::endl;
  deallog << "det Q = " << determinant(Q) << std::endl;
  deallog << "Q Q^T = " << (Q * transpose(Q)) << std::endl;
  deallog << std::endl;
}

int
main()
{
  initlog();
  // not orthogonal
  test(Tensor<2, 3, float>{{{1, 2, 3}, {2, 1, 4}, {3, 4, 1}}});
  test(Tensor<2, 3, double>{{{1, 2, 3}, {2, 1, 4}, {3, 4, 1}}});
  // already orthogonal
  test(Tensor<2, 3, float>{{{1, 0, 0}, {0, 0, 1}, {0, 1, 0}}});
  test(Tensor<2, 3, double>{{{1, 0, 0}, {0, 0, 1}, {0, 1, 0}}});
  // not orthogonal, but det = 1
  test(Tensor<2, 3, double>{{{1, 2, 0}, {0, 1, 0}, {0, 0, 1}}});
  // 2D not orthogonal
  test(Tensor<2, 2, float>{{{1, 2}, {3, 4}}});
  test(Tensor<2, 2, double>{{{1, 2}, {3, 4}}});
  // 2D already orthogonal
  test(Tensor<2, 2, float>{{{0, 1}, {1, 0}}});
  test(Tensor<2, 2, double>{{{0, 1}, {1, 0}}});
  // 1D not orthogonal
  test(Tensor<2, 1, float>{{{2.5}}});
  test(Tensor<2, 1, double>{{{2.5}}});
  // 1D already orthogonal
  test(Tensor<2, 1, float>{{{1.}}});
  test(Tensor<2, 1, double>{{{1.}}});

  deallog << "OK" << std::endl;
}
